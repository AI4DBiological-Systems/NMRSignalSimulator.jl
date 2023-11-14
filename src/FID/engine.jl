# based on core/engine.jl.

function fitfidproxies(
    As::Vector{HAM.SHType{T}},
    λ0::T,
    config::FIDSurrogateConfig{T};
    ) where T <: AbstractFloat

    SST = SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}

    N = length(As)
    Cs = Vector{MoleculeType{T,SST,FIDOperationRange{T}}}(undef, N)
    itp_samps = Vector{Vector{FIDInterpolationSamples{T}}}(undef, N)

    for n in eachindex(As)

        Cs[n], itp_samps[n] = fitfidproxy(SST, As[n], λ0, config)
    end

    shifts, phases, T2s, Δc_bars = setupSSvars(As, Cs)
    mixSS = FIDMixtureSpinSys(
        shifts,
        phases,
        T2s,
        Δc_bars,
        λ0,
    )

    return Cs, mixSS, itp_samps
end

function fitfidproxy(
    ::Type{SST},
    A::HAM.SHType{T},
    λ0::T,
    config::FIDSurrogateConfig,
    ) where {T,SST}

    #hz2ppmfunc = uu->(uu - A.ν_0ppm)*A.SW/A.fs
    #ppm2hzfunc = pp->(A.ν_0ppm + pp*A.fs/A.SW)

    # allocate
    N_coherence_vars_sys, N_resonance_groups_sys = getNcoherencesys(A.Δc_bar)
    ss_params = setupSSParamsparams(
        SST,
        N_coherence_vars_sys,
        N_resonance_groups_sys,
        λ0,
    )

    # prepare configuration parameters.
    # Δcs_max::Vector{T} = collect(
    #     config.Δcs_max_scalar for i in eachindex(A.N_spins_sys)
    # )
    # d_max::Vector{T} = ppm2hzfunc.(Δcs_max) .- ppm2hzfunc(zero(T))

    C = CLSurrogateSpinSysInputs(A.parts, A.αs, A.Ωs, zero(T), zero(T))
    
    itp_samps = getfiditpsamples(C, config)

    # in the current version, we don't have the gradient implemented for the FID surrogate, unlike setupclmoleculepartitionitp().
    qs = setupfidmoleculepartitionitp(
        itp_samps,
        ss_params.phase.cos_β,
        ss_params.phase.sin_β,
        C,
    )

    core = MoleculeType(
        qs,
        ss_params,
        FIDOperationRange(config.t_lb, config.Δt, config.t_ub),
        λ0,
    )

    return core, itp_samps
end

function setupfidmoleculepartitionitp(
    itp_samps::Vector{FIDInterpolationSamples{T}},
    cos_β::Vector{Vector{T}},
    sin_β::Vector{Vector{T}},
    C::CLSurrogateSpinSysInputs,
    ) where T <: AbstractFloat

    parts, αs = C.parts, C.αs

    N_sys = length(αs)

    qs = Vector{Vector{Function}}(undef, N_sys) # surrogates for lorentzian model.
    
    for i in eachindex(parts) # over elements in a spin group.

        N_groups = length(parts[i])
        qs[i] = Vector{Function}(undef, N_groups)
        
        samples = itp_samps[i].samples
        t_lb = itp_samps[i].t_lb
        Δt = itp_samps[i].Δt
        t_ub = itp_samps[i].t_ub

        A_t = step2LinRange(t_lb, Δt, t_ub)

        for k = 1:N_groups
            
            # get interpolation surrogates.
            A = samples[k]
            real_sitp, imag_sitp = setupfidpartitionitp(A, A_t)
            
            # construct simulator for this resonance group.
            qs[i][k] = (tt::T, rr::T)->evalqfid(
                real_sitp, imag_sitp, tt, rr, cos_β[i][k], sin_β[i][k],
            )::Complex{T}
        end
    end

    return qs
end

# doesn't include eps(-λ*t).
function evalqfid(real_sitp, imag_sitp, t::T, r::T, cos_β, sin_β)::Complex{T} where T <: AbstractFloat

    # Complex(a,b) approximates sum( α[l]*cis(Ω[l]*t) for l in eachindex(α) ) # see evalfidpartitionelement()

    a = real_sitp(t)
    b = imag_sitp(t)
    return Complex(a,b)*Complex(cos_β, sin_β)*cis(r*t)
end


function getfiditpsamples(
    C::CLSurrogateSpinSysInputs{T},
    config::FIDSurrogateConfig{T},
    )::Vector{FIDInterpolationSamples{T}} where T <: AbstractFloat

    parts, αs, Ωs = C.parts, C.αs, C.Ωs
    N_sys = length(αs)
    
    itp_samps = Vector{FIDInterpolationSamples{T}}(undef, N_sys)

    # loop over spin systems.
    for i in eachindex(parts)

        N_groups = length(parts[i])
        
        # set up default, which is for the singlet case.
        itp_samps[i] = FIDInterpolationSamples(T, length(N_groups))

        #A_r, A_λ = getfiditplocations(d_max[i], λ0, config)
        A_t = getfiditplocations(config)

        samples = Vector{Vector{Complex{T}}}(undef, N_groups)
        for k in eachindex(samples)

            inds = parts[i][k]
            α = αs[i][inds]
            Ω = Ωs[i][inds]
            samples[k] = computefidsamples(α, Ω, A_t)
        end

        # overwrite the default.
        itp_samps[i] = FIDInterpolationSamples(config, samples)

    end

    return itp_samps
end

function getfiditplocations(C::FIDSurrogateConfig{T}) where T

    t_lb, Δt, t_ub = C.t_lb, C.Δt, C.t_ub

    # LinRange is type stable.
    A_t = step2LinRange(t_lb, Δt, t_ub)

    return A_t
end

function computefidsamples(α::Vector{T}, Ω::Vector{T}, A_t)::Vector{Complex{T}} where T <: AbstractFloat

    # complex.
    # create complex-valued oracle and get values from it.
    f = tt->evalfidpartitionelement(tt, α, Ω)
    #A = [f(x1,x2) for x1 in A_r, x2 in A_ξ]
    A = collect( f(t) for t in A_t)

    return A
end

function setupfidpartitionitp(
    A::Vector{Complex{T}},
    A_t,
    ) where T <: AbstractFloat

    real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    #real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    real_sitp = Interpolations.scale(real_itp, A_t)
    real_setp = Interpolations.extrapolate(real_sitp, zero(T)) # zero outside interp range.

    imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    #imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    imag_sitp = Interpolations.scale(imag_itp, A_t)
    imag_setp = Interpolations.extrapolate(imag_sitp, zero(T)) # zero outside interp range.
    
    return real_setp, imag_setp
end

# # I am here. need to do a version of importmodel!() for FID's MoleculeType.
# function importmodel!(p::MixtureModelParameters)

#     updatespinsystems!(p.MSS, p.var_flat, p.systems_mapping)
#     resolvespinsystems!(p.MSS)