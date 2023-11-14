
"""
```
fitclproxies(
    As::Vector{HAM.SHType{T}},
    λ0::T,
    config::CLSurrogateConfig{T};
    names::Vector{String} = Vector{String}(undef, 0),
    ) where T <: AbstractFloat
```

Create surrogate for NMR spectrum, given the simulated resonance comopnent results from NMRHamiltonian.

Inputs:

- `As` -- a 1-D array of compound resonance simulations. See `NMRHamiltonian.simulate`.

- `λ0` -- the estimated T2* decay for the 0 ppm resonance component.

- `config` -- the configuration for building the surrogate. See `CLSurrogateConfig`.
"""
function fitclproxies(
    As::Vector{HAM.SHType{T}},
    λ0::T,
    config::CLSurrogateConfig{T};
    ) where T <: AbstractFloat
    
    ##type_SSParams = NMRSignalSimulator.getSpinSysParamsdatatype(NMRSignalSimulator.SharedShift{Float64})
    #type_SSParams = NMRSignalSimulator.getSpinSysParamsdatatype(NMRSignalSimulator.CoherenceShift{Float64})

    return fitclproxies(
        SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}},
        As,
        λ0,
        config,
    )
end

# """
# - `::Type{SpinSysParams{ST,PT,T2T}}`, set this to the surrogate model type you want to build. For now, only 
# ```
# NMRSignalSimulator.getSpinSysParamsdatatype(NMRSignalSimulator.CoherenceShift{T})
# ```
# is implemented, so pass that in as the input here. `T` is the concrete floating point number data type of your choice, e.g. `T = Float64`.

# """
function fitclproxies(
    ::Type{SpinSysParams{ST,PT,T2T}},
    As::Vector{HAM.SHType{T}},
    λ0::T,
    config::CLSurrogateConfig{T};
    )::Tuple{Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T}}},
    MixtureSpinSys{T,ST,PT,T2T},
    Vector{Vector{InterpolationSamples{T}}}} where {T,ST,PT,T2T}

    N = length(As)

    Bs = Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T}}}(undef, N)

    srs = Vector{Vector{Vector{Function}}}(undef, N)
    sis = Vector{Vector{Vector{Function}}}(undef, N)

    ∇srs! = Vector{Vector{Vector{Function}}}(undef, N)
    ∇sis! = Vector{Vector{Vector{Function}}}(undef, N)

    itp_samps = Vector{Vector{InterpolationSamples{T}}}(undef, N)

    for n in eachindex(As)
        A = As[n]

        # fit surrogate, save into `core`.
        Bs[n], srs[n], sis[n], ∇srs![n], ∇sis![n], itp_samps[n] = fitclproxy(
            SpinSysParams{ST,PT,T2T},
            A,
            λ0,
            config,
        )
    end

    # alternative specification for derivatives.
    shifts, phases, T2s, Δc_bars = setupSSvars(As, Bs)
    mixSS = MixtureSpinSys(
        srs,
        sis,
        ∇srs!,
        ∇sis!,
        shifts,
        phases,
        T2s,
        Δc_bars,
        λ0,
    )

    return Bs, mixSS, itp_samps
end


function fitclproxy(
    ::Type{SST},
    A::HAM.SHType{T},
    λ0::T,
    config::CLSurrogateConfig;
    ) where {T,SST}

    hz2ppmfunc = uu->(uu - A.ν_0ppm)*A.SW/A.fs
    ppm2hzfunc = pp->(A.ν_0ppm + pp*A.fs/A.SW)

    # allocate
    N_coherence_vars_sys, N_resonance_groups_sys = getNcoherencesys(A.Δc_bar)
    ss_params = setupSSParamsparams(
        SST,
        N_coherence_vars_sys,
        N_resonance_groups_sys,
        λ0,
    )

    # prepare configuration parameters.
    Δcs_max::Vector{T} = collect(
        config.Δcs_max_scalar for i in eachindex(A.N_spins_sys)
    )
    d_max::Vector{T} = ppm2hzfunc.(Δcs_max) .- ppm2hzfunc(zero(T))

    # threshold and partition the resonance components.
    Ωs_ppm = hz2ppmfunc.( combinevectors(A.Ωs) ./ twopi(T) )
    min_ppm = minimum(Ωs_ppm) - config.ppm_padding
    max_ppm = maximum(Ωs_ppm) + config.ppm_padding
    u_min = ppm2hzfunc(min_ppm)
    u_max = ppm2hzfunc(max_ppm)

    C = CLSurrogateSpinSysInputs(A.parts, A.αs, A.Ωs, u_min, u_max)
    
    itp_samps = getitpsamples(
        d_max,
        C,
        λ0,
        config,
    )

    qs, sr, si, ∇sr!, ∇si! = setupclmoleculepartitionitp(
        itp_samps,
        ss_params.phase.cos_β,
        ss_params.phase.sin_β,
        C,
    )

    core = MoleculeType(
        qs,
        ss_params,
        CLOperationRange(u_min, u_max, d_max, config),
        λ0,
    )

    return core, sr, si, ∇sr!, ∇si!, itp_samps
end



function getitpsamples(
    d_max::Vector{T},
    C::CLSurrogateSpinSysInputs{T},
    λ0::T,
    config::CLSurrogateConfig{T},
    )::Vector{InterpolationSamples{T}} where T <: AbstractFloat

    parts, αs, Ωs = C.parts, C.αs, C.Ωs
    N_sys = length(αs)
    
    itp_samps = Vector{InterpolationSamples{T}}(undef, N_sys)

    # loop over spin systems.
    for i in eachindex(parts)

        N_groups = length(parts[i])
        
        # set up default, which is for the singlet case.
        itp_samps[i] = InterpolationSamples(T, length(N_groups))

        r_min, r_max, A_r, A_λ = getitplocations(
            d_max[i],
            C.u_min,
            C.u_max,
            λ0,
            config,
        )

        samples = Vector{Matrix{Complex{T}}}(undef, N_groups)
        for k in eachindex(samples)
            inds = parts[i][k]
            α = αs[i][inds]
            Ω = Ωs[i][inds]
            samples[k] = computesamples(α, Ω, A_r, A_λ)
        end

        # overwrite the default.
        itp_samps[i] = InterpolationSamples(config, samples, r_min, r_max, λ0)

    end

    return itp_samps
end


function setupclmoleculepartitionitp(
    itp_samps::Vector{InterpolationSamples{T}},
    cos_β::Vector{Vector{T}},
    sin_β::Vector{Vector{T}},
    C::CLSurrogateSpinSysInputs,
    ) where T <: AbstractFloat

    parts, αs = C.parts, C.αs

    N_sys = length(αs)

    qs = Vector{Vector{Function}}(undef, N_sys) # surrogates for lorentzian model.
    
    sr = Vector{Vector{Function}}(undef, N_sys)
    si = Vector{Vector{Function}}(undef, N_sys)

    ∇sr! = Vector{Vector{Function}}(undef, N_sys)
    ∇si! = Vector{Vector{Function}}(undef, N_sys)
 
    for i in eachindex(parts) # over elements in a spin group.

        N_groups = length(parts[i])

        qs[i] = Vector{Function}(undef, N_groups)
        
        sr[i] = Vector{Function}(undef, N_groups)
        si[i] = Vector{Function}(undef, N_groups)

        ∇sr![i] = Vector{Function}(undef, N_groups)
        ∇si![i] = Vector{Function}(undef, N_groups)

        # parse.
        samples = itp_samps[i].samples
        r_min = itp_samps[i].r_min
        Δr = itp_samps[i].Δr
        r_max = itp_samps[i].r_max
        κ_λ_lb = itp_samps[i].κ_λ_lb
        Δκ_λ = itp_samps[i].Δκ_λ
        κ_λ_ub = itp_samps[i].κ_λ_ub
        λ0 = itp_samps[i].λ0

        # # StepRangeLen is type-unstable in Julia v 1.9.
        # A_r = r_min:Δr:r_max
        # A_λ = (κ_λ_lb:Δκ_λ:κ_λ_ub) .* λ0

        A_r = step2LinRange(r_min, Δr, r_max)
        A_λ = step2LinRange(κ_λ_lb*λ0, Δκ_λ*λ0, κ_λ_ub*λ0)

        # loop through resonance groups.
        for k = 1:N_groups

            # get interpolation surrogates.
            A = samples[k]
            real_sitp, imag_sitp = setupclpartitionitp(A, A_r, A_λ)
            
            # construct simulator for this resonance group.
            qs[i][k] = (rr::T, λλ::T)->evalq(real_sitp, imag_sitp, rr, λλ, cos_β[i][k], sin_β[i][k])::Complex{T}

            sr[i][k] = (rr::T, λλ::T)->real_sitp(rr,λλ)::T
            si[i][k] = (rr::T, λλ::T)->imag_sitp(rr,λλ)::T

            ∇sr![i][k] = (gg,rr,λλ)->Interpolations.gradient!(gg, real_sitp, rr, λλ)
            ∇si![i][k] = (gg,rr,λλ)->Interpolations.gradient!(gg, imag_sitp, rr, λλ)
        end

    
    end

    return qs, sr, si, ∇sr!, ∇si!
end


# function evalq(real_sitp, imag_sitp, r::T, λ::T, b::Vector{T}, c)::Complex{T} where T <: AbstractFloat

#     return (real_sitp(r,λ)+im*imag_sitp(r,λ))*cis(dot(b, c))
# end

function evalq(real_sitp, imag_sitp, r::T, λ::T, c, d)::Complex{T} where T <: AbstractFloat

    a = real_sitp(r,λ)
    b = imag_sitp(r,λ)
    return Complex(a*c-b*d, a*d+b*c)
end

function computesamples(α::Vector{T}, Ω::Vector{T}, A_r, A_λ)::Matrix{Complex{T}} where T <: AbstractFloat

    # complex.
    # f = (rr,λλ)->evalclpart(rr, α, Ω, λλ*λ0)
    # A = [f(x1,x2) for x1 in A_r, x2 in A_ξ]
    f = (rr,λλ)->evalclpart(rr, α, Ω, λλ)
    A = [f(x1,x2) for x1 in A_r, x2 in A_λ]

    return A
end

function setupclpartitionitp(
    A::Matrix{Complex{T}},
    A_r,
    A_λ,
    ) where T <: AbstractFloat

    real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    #real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    real_sitp = Interpolations.scale(real_itp, A_r, A_λ)
    real_setp = Interpolations.extrapolate(real_sitp, zero(T)) # zero outside interp range.

    imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    #imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    imag_sitp = Interpolations.scale(imag_itp, A_r, A_λ)
    imag_setp = Interpolations.extrapolate(imag_sitp, zero(T)) # zero outside interp range.
    
    return real_setp, imag_setp
end

function getitplocations(
    d_max::T,
    u_min::T,
    u_max::T,
    λ0::T,
    C::CLSurrogateConfig{T},
    ) where T
    
    κ_λ_lb, κ_λ_ub, Δr, Δκ_λ = C.κ_λ_lb, C.κ_λ_ub, C.Δr, C.Δκ_λ
    
    r_min = twopi(T)*(u_min - d_max)
    r_max = twopi(T)*(u_max + d_max)

    # # set up sample grid range.

    # # StepRangeLen is type-unstable.
    # A_r = r_min:Δr:r_max # large range.
    # A_ξ = κ_λ_lb:Δκ_λ:κ_λ_ub
    # A_λ = A_ξ .* λ0

    # LinRange is type stable.
    A_r = step2LinRange(r_min, Δr, r_max) # large range.
    A_λ = step2LinRange(κ_λ_lb*λ0, Δκ_λ*λ0, κ_λ_ub*λ0)

    return r_min, r_max, A_r, A_λ
end

function step2LinRange(x1::T, s::T, x2::T)::LinRange{T,Int} where T
    N_steps = floor(Int, (x2+s-x1)/s)
    y1 = x1
    y2 = convert(T, x1 + s*(N_steps-1))
    return LinRange(y1, y2, N_steps)
end
# x1, s, x2 = (80692.28f0, 1.0f0, 87574.445f0)
# r = x1:s:x2
# N_steps = floor(Int, (x2+s-x1)/s)
# y1 = x1
# y2 = convert(T, x1 + s*(N_steps-1))
# q = LinRange(y1,y2,N_steps)
# @assert norm(collect(r) - collect(q)) < 1e-8@show