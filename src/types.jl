
# reduce clutter in code.
#SHType = NMRHamiltonian.SHType
const SHType{T} = NMRHamiltonian.SHType{T} where T

### model.

# This is without the compensation amplitude parameter, κ_α, denoted κs_α in code.
struct CompoundType{T,SST} # parameters for surrogate model.

    # non-singlet spin systems.
    qs::Vector{Vector{Function}} # spin group, partition element index.

    # κs_λ::Vector{T} # multipliers for each spin group.
    # κs_β::Vector{Vector{T}} # coefficients for each (spin group, partition element).
    # d::Vector{T}
    ss_params::SST

    # singlets.
    κs_λ_singlets::Vector{T}
    β_singlets::Vector{T}
    d_singlets::Vector{T}

    # misc.
    Δcs_max::Vector{T}
    λ0::T
end

# This is with the compensation amplitude parameter, κ_α, denoted κs_α in code.
mutable struct καCompoundType{T,SST}
    κs_α::Vector{Vector{T}} # spin group, partition element index.
    κs_α_singlets::Vector{T}
    core::CompoundType{T,SST}
end

function καCompoundType(core::CompoundType{T,SST}) where {T,SST}

    N_spins = length(core.qs)
    κs_α = Vector{Vector{T}}(undef, N_spins)
    for i = 1:N_spins
        κs_α[i] = ones(T, length(core.qs[i]))
    end

    κs_α_singlets = ones(T, length(core.d_singlets))

    return καCompoundType(κs_α, κs_α_singlets, core)
end

### different parameterizations of the spin system FID parameters.
struct SpinSysParamsType2{T}
    κs_λ::Vector{T} # a multiplier for each (spin group.
    κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
    d::Vector{Vector{T}} # a multiplier for each (spin group, partition element).
    κs_d::Vector{Vector{T}} # same size as κs_β. # intermediate buffer for d.
end

function SpinSysParamsType2(x::T) where T
    return SpinSysParamsType2(Vector{T}(undef,0), Vector{Vector{T}}(undef, 0), Vector{Vector{T}}(undef, 0), Vector{Vector{T}}(undef, 0))
end

struct SpinSysParamsType1{T}
    κs_λ::Vector{T} # a common multiplier for each spin group. length: number of spin groups.
    κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
    d::Vector{T} # a common multiplier for each spin group. length: number of spin groups.
end

function SpinSysParamsType1(x::T) where T
    return SpinSysParamsType1(Vector{T}(undef,0), Vector{Vector{T}}(undef, 0), Vector{T}(undef, 0))
end

function constructorSSFID(x::SpinSysParamsType1{T}, y...)::SpinSysParamsType1{T} where T
    return SpinSysParamsType1(y...)
end

function constructorSSFID(x::SpinSysParamsType2{T}, y...)::SpinSysParamsType2{T} where T
    return SpinSysParamsType2(y...)
end

########### more elaborate constructors.

function setupSSFIDparams(dummy_SSFID::SpinSysParamsType1{T}, part_inds_compound, N_β_vars_sys)::SpinSysParamsType1{T} where T
    L = length(part_inds_compound)

    κs_λ = ones(T, L)
    κs_β = collect( zeros(T, N_β_vars_sys[i]) for i = 1:length(N_β_vars_sys))
    #d = rand(length(αs))
    d = zeros(T, L)

    return constructorSSFID(dummy_SSFID, κs_λ, κs_β, d)
end

function setupSSFIDparams(dummy_SSFID::SpinSysParamsType2{T}, part_inds_compound::Vector{Vector{Vector{Int}}}, N_β_vars_sys)::SpinSysParamsType2{T} where T

    N_sys = length(part_inds_compound)
    κs_λ = ones(T, N_sys)
    d = Vector{Vector{T}}(undef, N_sys)

    for i = 1:length(d)
        N_partition_elements = length(part_inds_compound[i])

        d[i] = zeros(T, N_partition_elements)
    end

    κs_β = collect( zeros(T, N_β_vars_sys[i]) for i = 1:length(N_β_vars_sys))
    κs_d = collect( zeros(T, N_β_vars_sys[i]) for i = 1:length(N_β_vars_sys))

    return constructorSSFID(dummy_SSFID, κs_λ, κs_β, d, κs_d)
end
