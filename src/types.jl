
# reduce clutter in code.
#SHType = NMRHamiltonian.SHType
const SHType{T} = NMRHamiltonian.SHType{T} where T

### model.

# This is without the compensation amplitude parameter, κ_α, denoted κs_α in code.
struct MoleculeType{T,SST} # parameters for surrogate model.

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
mutable struct καMoleculeType{T,SST}
    κs_α::Vector{Vector{T}} # spin group, partition element index.
    κs_α_singlets::Vector{T}
    core::MoleculeType{T,SST}
end

function καMoleculeType(core::MoleculeType{T,SST}) where {T,SST}

    N_spins = length(core.qs)
    κs_α = Vector{Vector{T}}(undef, N_spins)
    for i = 1:N_spins
        κs_α[i] = ones(T, length(core.qs[i]))
    end

    κs_α_singlets = ones(T, length(core.d_singlets))

    return καMoleculeType(κs_α, κs_α_singlets, core)
end

### different parameterizations of the spin system FID parameters.
abstract type SpinSysParams{T<:AbstractFloat} end


struct CoherenceShift{T} <: SpinSysParams{T}
    κs_λ::Vector{T} # a multiplier for each (spin group.
    κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
    d::Vector{Vector{T}} # a multiplier for each (spin group, partition element).
    κs_d::Vector{Vector{T}} # same size as κs_β. # intermediate buffer for d.
end

function CoherenceShift(x::T) where T
    return CoherenceShift(Vector{T}(undef,0), Vector{Vector{T}}(undef, 0), Vector{Vector{T}}(undef, 0), Vector{Vector{T}}(undef, 0))
end

struct SharedShift{T} <: SpinSysParams{T}
    κs_λ::Vector{T} # a common multiplier for each spin group. length: number of spin groups.
    κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
    d::Vector{T} # a common multiplier for each spin group. length: number of spin groups.
end

function SharedShift(x::T) where T
    return SharedShift(Vector{T}(undef,0), Vector{Vector{T}}(undef, 0), Vector{T}(undef, 0))
end

function constructorSSParams(x::SharedShift{T}, y...)::SharedShift{T} where T
    return SharedShift(y...)
end

function constructorSSParams(x::CoherenceShift{T}, y...)::CoherenceShift{T} where T
    return CoherenceShift(y...)
end

########### more elaborate constructors.

function setupSSParamsparams(dummy_SSParams::SharedShift{T}, part_inds_molecule, N_β_vars_sys)::SharedShift{T} where T
    L = length(part_inds_molecule)

    κs_λ = ones(T, L)
    κs_β = collect( zeros(T, N_β_vars_sys[i]) for i in eachindex(N_β_vars_sys))
    #d = rand(length(αs))
    d = zeros(T, L)

    return constructorSSParams(dummy_SSParams, κs_λ, κs_β, d)
end

function setupSSParamsparams(dummy_SSParams::CoherenceShift{T}, part_inds_molecule::Vector{Vector{Vector{Int}}}, N_β_vars_sys)::CoherenceShift{T} where T

    N_sys = length(part_inds_molecule)
    κs_λ = ones(T, N_sys)
    d = Vector{Vector{T}}(undef, N_sys)

    for i in eachindex(d)
        N_partition_elements = length(part_inds_molecule[i])

        d[i] = zeros(T, N_partition_elements)
    end

    κs_β = collect( zeros(T, N_β_vars_sys[i]) for i in eachindex(N_β_vars_sys))
    κs_d = collect( zeros(T, N_β_vars_sys[i]) for i in eachindex(N_β_vars_sys))

    return constructorSSParams(dummy_SSParams, κs_λ, κs_β, d, κs_d)
end
