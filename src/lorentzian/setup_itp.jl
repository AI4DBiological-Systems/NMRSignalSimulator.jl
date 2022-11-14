
function fitclproxies(
    ::Type{SpinSysParams{ST,PT,T2T}},
    As::Vector{SHType{T}},
    λ0::T;
    names::Vector{String} = Vector{String}(undef, 0),
    config_path::String = "",
    Δcs_max_scalar_default = 0.2,
    κ_λ_lb_default = 0.5,
    κ_λ_ub_default = 2.5,
    u_min = Inf,
    u_max = Inf,
    Δr_default = 1.0,
    Δκ_λ_default = 0.05,
    default_ppm_padding = 0.5,
    )::Tuple{Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T}}},
    MixtureSpinSys{T,ST,PT,T2T},
    MixtureSinglets{T}} where {T,ST,PT,T2T}

    config_dict = Dict()
    if ispath(config_path)

        # TODO add error-handling if name is not found in the dictionary or filename does not exist.
        #config_dict = JSON.parsefile(config_path)
        config_dict = JSON3.read(read(config_path))
    end

    N = length(As)

    Bs = Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T}}}(undef, N)

    srs = Vector{Vector{Vector{Function}}}(undef, N)
    sis = Vector{Vector{Vector{Function}}}(undef, N)

    ∇srs! = Vector{Vector{Vector{Function}}}(undef, N)
    ∇sis! = Vector{Vector{Vector{Function}}}(undef, N)

    for n in eachindex(As)
        A = As[n]

        # fit surrogate, save into `core`.
        Bs[n], srs[n], sis[n], ∇srs![n], ∇sis![n] = fitclproxy(
            SpinSysParams{ST,PT,T2T},
            A,
            λ0,
            config_dict;
            molecule_name = names[n],
            Δcs_max_scalar_default = Δcs_max_scalar_default,
            κ_λ_lb_default = κ_λ_lb_default,
            κ_λ_ub_default = κ_λ_ub_default,
            u_min = u_min,
            u_max = u_max,
            Δr_default = Δr_default,
            Δκ_λ_default = Δκ_λ_default,
            default_ppm_padding = default_ppm_padding,
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
    )

    MS = MixtureSinglets(
        collect( Bs[n].κs_λ_singlets for n in eachindex(Bs) ),
        collect( Bs[n].β_singlets for n in eachindex(Bs) ),
        collect( Bs[n].d_singlets for n in eachindex(Bs) ),
        collect( As[n].αs_singlets for n in eachindex(As) ),
        collect( As[n].Ωs_singlets for n in eachindex(As) ),
        λ0,
    )
    return Bs, mixSS, MS
end

function fitclproxy(
    ::Type{SST},
    A::SHType{T},
    λ0::T,
    config_dict;
    molecule_name::String = "",
    Δcs_max_scalar_default = 0.2,
    κ_λ_lb_default = 0.5,
    κ_λ_ub_default = 2.5,
    u_min = Inf,
    u_max = Inf,
    Δr_default = 1.0,
    Δκ_λ_default = 0.05,
    default_ppm_padding = 0.5,
    ) where {T,SST}

    hz2ppmfunc = uu->(uu - A.ν_0ppm)*A.SW/A.fs
    ppm2hzfunc = pp->(A.ν_0ppm + pp*A.fs/A.SW)

    # allocate `core` data structure.
    N_singlets = length(A.αs_singlets)
    κs_λ_singlets = ones(T, N_singlets)
    λ_singlets = λ0 .* κs_λ_singlets
    κs_β_singlets = zeros(T, N_singlets)
    d_singlets = zeros(T, N_singlets)

    #N_β_vars_sys = A.N_spins_sys # no equivalence used.
    N_coherence_vars_sys::Vector{Int} = collect( length(A.Δc_bar[i][begin]) for i in eachindex(A.Δc_bar) )
    N_resonance_groups_sys::Vector{Int} = collect( length(A.Δc_bar[i]) for i in eachindex(A.Δc_bar) )

    # proxy placeholder.
    #qs = Vector{Vector{Function}}(undef, length(N_spins_sys))

    SSParams_obj = setupSSParamsparams(
        SST,
        N_coherence_vars_sys,
        N_resonance_groups_sys,
        #λ0,
    )

    # prepare configuration parameters.
    κ_λ_lb = κ_λ_lb_default
    κ_λ_ub = κ_λ_ub_default
    Δκ_λ = Δκ_λ_default
    Δcs_max::Vector{T} = collect( Δcs_max_scalar_default for i in eachindex(A.N_spins_sys))
    Δr = Δr_default

    d_max::Vector{T} = ppm2hzfunc.(Δcs_max) .- ppm2hzfunc(zero(T))

    if !isempty(config_dict) && !isempty(molecule_name)
        dict = config_dict[Symbol(molecule_name)] # TODO graceful error-handle.

        κ_λ_lb = dict[Symbol("κ_λ lower bound")]
        κ_λ_ub = dict[Symbol("κ_λ upper bound")]
        Δκ_λ = dict[Symbol("κ_λ surrogate sampling step size")]
        if length(dict[Symbol("Δcs_max")]) == length(A.N_spins_sys)
            Δcs_max = dict[Symbol("Δcs_max")]
        # else
        #     println("Warning: problem an entry's Δcs_max value in config file. Using default scalar value for Δcs_max")
        end
        Δr = dict[Symbol("radians surrogate sampling step size")]

        d_max = collect( ppm2hzfunc(Δcs_max[i])-ppm2hzfunc(0.0) for i in eachindex(Δcs_max) )
    end

    # threshold and partition the resonance components.
    if !isfinite(u_min) || !isfinite(u_max)
        Ωs_ppm = hz2ppmfunc.( combinevectors(A.Ωs) ./ (2*π) )
        tmp = hz2ppmfunc.( A.Ωs_singlets ./ (2*π) )
        push!(Ωs_ppm, tmp...)

        min_ppm = minimum(Ωs_ppm) - default_ppm_padding
        max_ppm = maximum(Ωs_ppm) + default_ppm_padding
        u_min = ppm2hzfunc(min_ppm)
        u_max = ppm2hzfunc(max_ppm)
    end

    qs, sr, si, ∇sr!, ∇si! = setupclmoleculepartitionitp(
        d_max,
        #SSParams_obj.κs_β,
        SSParams_obj.phase.cos_β,
        SSParams_obj.phase.sin_β,
        A.Δc_bar,
        A.part_inds_molecule,
        A.αs, A.Ωs,
        λ0, 
        u_min, u_max;
        κ_λ_lb = κ_λ_lb,
        κ_λ_ub = κ_λ_ub,
        Δr = Δr,
        Δκ_λ = Δκ_λ,
    )

    core = MoleculeType(
        qs,
        SSParams_obj,
        κs_λ_singlets,
        λ_singlets,
        κs_β_singlets,
        d_singlets,
        Δcs_max,
        λ0
    )

    return core, sr, si, ∇sr!, ∇si!
end

function setupclmoleculepartitionitp(
    d_max::Vector{T},
    #κs_β::Vector{Vector{T}},
    cos_β::Vector{Vector{T}},
    sin_β::Vector{Vector{T}},
    Δc_bar::Vector{Vector{Vector{T}}},
    part_inds_molecule::Vector{Vector{Vector{Int}}},
    αs::Vector{Vector{T}},
    Ωs::Vector{Vector{T}},
    λ0::T,
    u_min::T,
    u_max::T;
    κ_λ_lb = 0.5,
    κ_λ_ub = 2.5,
    Δr = 1.0,
    Δκ_λ = 0.05,
    ) where T <: AbstractFloat

    N_sys = length(αs)

    qs = Vector{Vector{Function}}(undef, N_sys) # surrogates for lorentzian model.
    
    sr = Vector{Vector{Function}}(undef, N_sys)
    si = Vector{Vector{Function}}(undef, N_sys)

    ∇sr! = Vector{Vector{Function}}(undef, N_sys)
    ∇si! = Vector{Vector{Function}}(undef, N_sys)
    
    #gs = Vector{Vector{Function}}(undef, length(αs)) # surrogates for FID model.

    for i in eachindex(part_inds_molecule) # over elements in a spin group.

        N_partition_elements = length(part_inds_molecule[i])

        qs[i] = Vector{Function}(undef, N_partition_elements)
        
        sr[i] = Vector{Function}(undef, N_partition_elements)
        si[i] = Vector{Function}(undef, N_partition_elements)

        ∇sr![i] = Vector{Function}(undef, N_partition_elements)
        ∇si![i] = Vector{Function}(undef, N_partition_elements)

        #gs[i] = Vector{Function}(undef, N_partition_elements)

        for k = 1:N_partition_elements
            #println("i,k", (i,k))

            inds = part_inds_molecule[i][k]
            α = αs[i][inds]
            Ω = Ωs[i][inds]

            real_sitp, imag_sitp = setupclpartitionitp(α, Ω,
            d_max[i], λ0, u_min, u_max; κ_λ_lb = κ_λ_lb, κ_λ_ub = κ_λ_ub,
            Δr = Δr, Δκ_λ = Δκ_λ)

            #qs[i][k] = (rr, ξξ)->evalq(real_sitp, imag_sitp, rr, ξξ, κs_β[i], Δc_bar[i][k])
            qs[i][k] = (rr::T, ξξ::T)->evalq(real_sitp, imag_sitp, rr, ξξ, cos_β[i][k], sin_β[i][k])::Complex{T}

            sr[i][k] = (rr::T, ξξ::T)->real_sitp(rr,ξξ)::T
            si[i][k] = (rr::T, ξξ::T)->imag_sitp(rr,ξξ)::T

            ∇sr![i][k] = (gg,rr,ξξ)->Interpolations.gradient!(gg, real_sitp, rr, ξξ)
            ∇si![i][k] = (gg,rr,ξξ)->Interpolations.gradient!(gg, imag_sitp, rr, ξξ)
        end
    end

    return qs, sr, si, ∇sr!, ∇si!
end

# function evalq(real_sitp, imag_sitp, r::T, ξ::T, b::Vector{T}, c)::Complex{T} where T <: AbstractFloat

#     return (real_sitp(r,ξ)+im*imag_sitp(r,ξ))*cis(dot(b, c))
# end

function evalq(real_sitp, imag_sitp, r::T, ξ::T, c, d)::Complex{T} where T <: AbstractFloat

    a = real_sitp(r,ξ)
    b = imag_sitp(r,ξ)
    return Complex(a*c-b*d, a*d+b*c)
end

function setupclpartitionitp(
    α::Vector{T},
    Ω::Vector{T},
    d_max::T,
    λ0::T,
    u_min::T,
    u_max::T;
    κ_λ_lb = 0.5,
    κ_λ_ub = 2.5,
    Δr = 1.0,
    Δκ_λ = 0.05,
    ) where T <: AbstractFloat

    # A_x1 = 1:.1:10
    # A_x2 = 1:.5:20
    # f(x1, x2) = log(x1+x2)
    # A = [f(x1,x2) for x1 in A_x1, x2 in A_x2]
    # itp = Interpolations.interpolate(A, BSpline(Cubic(Line(OnGrid()))))
    # sitp = Interpolations.scale(itp, A_x1, A_x2)
    # sitp(5., 10.) # exactly log(5 + 10)
    # sitp(5.6, 7.1) # approximately log(5.6 + 7.1)

    # set up bounds.

    # r_min = 2*π*u_min - d_max
    # r_max = 2*π*u_max + d_max
    r_min = 2*π*(u_min - d_max)
    r_max = 2*π*(u_max + d_max)

    # set up samples.
    A_r = r_min:Δr:r_max # large range.
    #A_r = r_min:20.0:r_max # large range.
    A_ξ = κ_λ_lb:Δκ_λ:κ_λ_ub
    # println("length(A_r) = ", length(A_r))
    # println("length(A_ξ) = ", length(A_ξ))

    # complex.
    f = (rr,ξξ)->evalclpart(rr, α, Ω, ξξ*λ0)
    #f = (rr,ξξ)->evalclpart(rr, α, Ω, λ[ind])
    A = [f(x1,x2) for x1 in A_r, x2 in A_ξ]

    real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    #real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    real_sitp = Interpolations.scale(real_itp, A_r, A_ξ)
    real_setp = Interpolations.extrapolate(real_sitp, 0.0) # zero outside interp range.

    imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    #imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    imag_sitp = Interpolations.scale(imag_itp, A_r, A_ξ)
    imag_setp = Interpolations.extrapolate(imag_sitp, 0.0) # zero outside interp range.


    # @show α
    # @show Ω
    # @show A_r
    # @show A_ξ
    # @assert 1==2

    #return real_sitp, imag_sitp
    return real_setp, imag_setp
end