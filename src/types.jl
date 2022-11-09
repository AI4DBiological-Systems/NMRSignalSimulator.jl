
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
    κs_λ_singlets::Vector{T} # multiplier wrt λ0.
    λ_singlets::Vector{T} # actual T2.
    β_singlets::Vector{T}
    d_singlets::Vector{T}

    # misc.
    Δcs_max::Vector{T}
    λ0::T
end


### different parameterizations of the spin system FID parameters.
abstract type ShiftParms{T<:AbstractFloat} end
abstract type PhaseParms{T<:AbstractFloat} end
abstract type T2Parms{T<:AbstractFloat} end


# struct CoherenceShift{T} <: SpinSysParams{T}
#     κs_λ::Vector{T} # a multiplier for each (spin group.
#     κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
#     d::Vector{Vector{T}} # a multiplier for each (spin group, partition element).
#     κs_d::Vector{Vector{T}} # same size as κs_β. # intermediate buffer for d.
# end



# struct SharedShift{T} <: SpinSysParams{T}
#     κs_λ::Vector{T} # a common multiplier for each spin group. length: number of spin groups.
#     κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
#     d::Vector{T} # a common multiplier for each spin group. length: number of spin groups.
# end




struct SharedT2{T} <: T2Parms{T}
    κs_λ::Vector{T} # multiplier wrt some λ0. length: number of spin groups.
    #λ::Vector{T} # actual decay. length: number of spin groups.
end

function SharedT2(
    ::Type{T},
    N_sys::Int,
    #λ0::T
    )::SharedT2{T} where T <: AbstractFloat

    κs_λ = ones(T, N_sys)
    #λ = λ0 .* κs_λ

    #return SharedT2(κs_λ, λ)
    return SharedT2(κs_λ)
end

struct SharedShift{T} <: ShiftParms{T}
    d::Vector{T} # length: number of spin groups.
end

function SharedShift(
    ::Type{T},
    N_sys::Int,
    )::SharedShift{T} where T <: AbstractFloat

    d = zeros(T, N_sys)

    return SharedShift(d)
end

struct CoherenceShift{T} <: ShiftParms{T}
    κs_d::Vector{Vector{T}} # first index for spin systems, second for coherence dimension.
    d::Vector{Vector{T}} # first index for spin systems, second for resonance groups.
end

function CoherenceShift(
    ::Type{T},
    N_coherence_vars_sys::Vector{Int},
    N_resonance_groups_sys::Vector{Int},
    )::CoherenceShift{T} where T <: AbstractFloat

    κs_d = collect( zeros(T, N_coherence_vars_sys[i]) for i in eachindex(N_coherence_vars_sys))
    d = collect( zeros(T, N_resonance_groups_sys[i]) for i in eachindex(N_resonance_groups_sys))

    return CoherenceShift(κs_d, d)
end


struct CoherencePhase{T} <: PhaseParms{T}
    κs_β::Vector{Vector{T}} # first index for spin systems, second for coherence dimension.
    cos_β::Vector{Vector{T}} # first index for spin systems, second for resonance groups.
    sin_β::Vector{Vector{T}}
end

function CoherencePhase(
    ::Type{T},
    N_coherence_vars_sys::Vector{Int},
    N_resonance_groups_sys::Vector{Int},
    )::CoherencePhase{T} where T <: AbstractFloat

    κs_β = collect( zeros(T, N_coherence_vars_sys[i]) for i in eachindex(N_coherence_vars_sys))
    cos_β = collect( ones(T, N_resonance_groups_sys[i]) for i in eachindex(N_resonance_groups_sys))
    sin_β = collect( zeros(T, N_resonance_groups_sys[i]) for i in eachindex(N_resonance_groups_sys))

    return CoherencePhase(κs_β, cos_β, sin_β)
end



struct SpinSysParams{ST,PT,DT}
    shift::ST
    phase::PT
    T2::DT
end

function setupSSParamsparams(
    ::Type{SpinSysParams{SharedShift{T}, CoherencePhase{T}, SharedT2{T}}},
    N_coherence_vars_sys::Vector{Int},
    N_resonance_groups_sys::Vector{Int},
    )::SpinSysParams{SharedShift{T}, CoherencePhase{T}, SharedT2{T}} where T <: AbstractFloat

    N_sys = length(N_coherence_vars_sys)

    return SpinSysParams(
        SharedShift(T, N_sys),
        CoherencePhase(T, N_coherence_vars_sys, N_resonance_groups_sys),
        SharedT2(T, N_sys),
    )
end

function setupSSParamsparams(
    ::Type{SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}},
    N_coherence_vars_sys::Vector{Int},
    N_resonance_groups_sys::Vector{Int},
    )::SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}} where T <: AbstractFloat

    N_sys = length(N_coherence_vars_sys)

    return SpinSysParams(
        CoherenceShift(T, N_coherence_vars_sys, N_resonance_groups_sys),
        CoherencePhase(T, N_coherence_vars_sys, N_resonance_groups_sys),
        SharedT2(T, N_sys),
    )
end

function getSpinSysParamsdatatype(
    ::Type{CoherenceShift{T}}
    )::DataType where T <: AbstractFloat

    return SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}
end

function getSpinSysParamsdatatype(
    ::Type{SharedShift{T}}
    )::DataType where T <: AbstractFloat

    return SpinSysParams{SharedShift{T}, CoherencePhase{T}, SharedT2{T}}
end