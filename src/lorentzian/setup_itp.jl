
function fitclproxies(As::Vector{SHType{T}},
    dummy_SSFID::SST,
    λ0::T;
    names::Vector{String} = Vector{String}(undef, 0),
    config_path::String = "",
    Δcs_max_scalar_default = 0.2,
    κ_λ_lb_default = 0.5,
    κ_λ_ub_default = 2.5,
    u_min = Inf,
    u_max = Inf,
    Δr_default = 1.0,
    Δκ_λ_default = 0.05) where {T,SST}

    config_dict = Dict()
    if ispath(config_path)

        # TODO add error-handling if name is not found in the dictionary or filename does not exist.
        #config_dict = JSON.parsefile(config_path)
        config_dict = JSON3.read(read(config_path))
    end

    cores = Vector{CompoundType{T,SST}}(undef, length(As))

    for n = 1:length(As)
        A = As[n]

        # fit surrogate, save into `core`.
        cores[n] = fitclproxy(dummy_SSFID, A, λ0, config_dict;
        compound_name = names[n],
        Δcs_max_scalar_default = Δcs_max_scalar_default,
        κ_λ_lb_default = κ_λ_lb_default,
        κ_λ_ub_default = κ_λ_ub_default,
        u_min = u_min,
        u_max = u_max,
        Δr_default = Δr_default,
        Δκ_λ_default = Δκ_λ_default)
    end

    return cores
end


function fitclproxy(dummy_SSFID::SST,
    A::SHType{T},
    λ0::T,
    config_dict;
    compound_name::String = "",
    Δcs_max_scalar_default = 0.2,
    κ_λ_lb_default = 0.5,
    κ_λ_ub_default = 2.5,
    u_min = Inf,
    u_max = Inf,
    Δr_default = 1.0,
    Δκ_λ_default = 0.05)::CompoundType{T,SST} where {T,SST}

    hz2ppmfunc = uu->(uu - A.ν_0ppm)*A.SW/A.fs
    ppm2hzfunc = pp->(A.ν_0ppm + pp*A.fs/A.SW)

    # allocate `core` data structure.
    N_singlets = length(A.αs_singlets)
    κs_λ_singlets = ones(T, N_singlets)
    κs_β_singlets = zeros(T, N_singlets)
    d_singlets = zeros(T, N_singlets)

    #N_β_vars_sys = A.N_spins_sys # no equivalence used.
    N_β_vars_sys::Vector{Int} = collect( length(A.Δc_bar[i][1]) for i = 1:length(A.Δc_bar) )

    # proxy placeholder.
    #qs = Vector{Vector{Function}}(undef, length(N_spins_sys))

    SSFID_obj = setupSSFIDparams(dummy_SSFID, A.part_inds_compound, N_β_vars_sys)

    # prepare configuration parameters.
    κ_λ_lb = κ_λ_lb_default
    κ_λ_ub = κ_λ_ub_default
    Δκ_λ = Δκ_λ_default
    Δcs_max::Vector{T} = collect( Δcs_max_scalar_default for i = 1:length(A.N_spins_sys))
    Δr = Δr_default

    d_max::Vector{T} = ppm2hzfunc.(Δcs_max) .- ppm2hzfunc(zero(T))

    if !isempty(config_dict) && !isempty(compound_name)
        dict = config_dict[Symbol(compound_name)] # TODO graceful error-handle.

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

        min_ppm = minimum(Ωs_ppm) - 0.5
        max_ppm = maximum(Ωs_ppm) + 0.5
        u_min = ppm2hzfunc(min_ppm)
        u_max = ppm2hzfunc(max_ppm)
    end

    qs = setupclcompoundpartitionitp(d_max,
        SSFID_obj.κs_β,
        A.Δc_bar,
        A.part_inds_compound,
        A.αs, A.Ωs,
        λ0, u_min, u_max;
        κ_λ_lb = κ_λ_lb,
        κ_λ_ub = κ_λ_ub,
        Δr = Δr,
        Δκ_λ = Δκ_λ)

    core = CompoundType(qs, SSFID_obj, κs_λ_singlets, κs_β_singlets, d_singlets,
        Δcs_max, λ0)

    return core
end

function setupclpartitionitp(α::Vector{T}, Ω::Vector{T}, d_max::T, λ0::T,
    u_min::T, u_max::T;
    κ_λ_lb = 0.5,
    κ_λ_ub = 2.5,
    Δr = 1.0,
    Δκ_λ = 0.05) where T <: Real

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
    f = (rr,ξξ)->evalclpartitionelement(rr, α, Ω, ξξ*λ0)
    A = [f(x1,x2) for x1 in A_r, x2 in A_ξ]

    #real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    real_itp = Interpolations.interpolate(real.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    real_sitp = Interpolations.scale(real_itp, A_r, A_ξ)
    real_setp = Interpolations.extrapolate(real_sitp, 0.0) # zero outside interp range.

    #imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid()))))
    imag_itp = Interpolations.interpolate(imag.(A), Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))))
    imag_sitp = Interpolations.scale(imag_itp, A_r, A_ξ)
    imag_setp = Interpolations.extrapolate(imag_sitp, 0.0) # zero outside interp range.

    #return real_sitp, imag_sitp
    return real_setp, imag_setp
end

function setupclcompoundpartitionitp(d_max::Vector{T},
    κs_β::Vector{Vector{T}},
    #Δc_m_compound::Vector{Vector{Vector{T}}},
    Δc_bar::Vector{Vector{Vector{T}}},
    part_inds_compound::Vector{Vector{Vector{Int}}},
    αs::Vector{Vector{T}}, Ωs::Vector{Vector{T}},
    λ0::T, u_min::T, u_max::T;
    κ_λ_lb = 0.5,
    κ_λ_ub = 2.5,
    Δr = 1.0,
    Δκ_λ = 0.05) where T

    qs = Vector{Vector{Function}}(undef, length(αs)) # surrogates for lorentzian model.
    #gs = Vector{Vector{Function}}(undef, length(αs)) # surrogates for FID model.

    for i = 1:length(part_inds_compound) # over elements in a spin group.

        N_partition_elements = length(part_inds_compound[i])
        qs[i] = Vector{Function}(undef, N_partition_elements)
        #gs[i] = Vector{Function}(undef, N_partition_elements)

        for k = 1:N_partition_elements
            #println("i,k", (i,k))

            inds = part_inds_compound[i][k]
            α = αs[i][inds]
            Ω = Ωs[i][inds]

            real_sitp, imag_sitp = setupclpartitionitp(α, Ω,
            d_max[i], λ0, u_min, u_max; κ_λ_lb = κ_λ_lb, κ_λ_ub = κ_λ_ub,
            Δr = Δr, Δκ_λ = Δκ_λ)

            qs[i][k] = (rr, ξξ)->evalq(real_sitp, imag_sitp, rr, ξξ, κs_β[i], Δc_bar[i][k])
        end
    end

    return qs
end

function evalq(real_sitp, imag_sitp, r::T, ξ::T, b::Vector{T}, c)::Complex{T} where T <: Real

    return (real_sitp(r,ξ)+im*imag_sitp(r,ξ))*cis(dot(b, c))
end
