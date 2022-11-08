
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


struct SpinSysParams{ST,PT,DT}
    shift::ST
    phase::PT
    T2::DT
end

struct SharedT2{T} <: T2Parms{T}
    κs_λ::Vector{T} # multiplier wrt some λ0. length: number of spin groups.
    λ::Vector{T} # actual decay. length: number of spin groups.
end

function SharedT2(::Type{T}) where T
    return SharedT2(
        Vector{T}(undef, 0),
        Vector{T}(undef, 0),
    )
end

struct SharedShift{T} <: ShiftParms{T}
    d::Vector{T} # length: number of spin groups.
end

function SharedShift(::Type{T}) where T
    return CoherenceShift(Vector{T}(undef, 0))
end

struct CoherenceShift{T} <: ShiftParms{T}
    d::Vector{Vector{T}} # first index for spin systems, second for resonance groups.
    κs_d::Vector{Vector{T}} # first index for spin systems, second for coherence dimension.
end

function CoherenceShift(::Type{T}) where T
    return CoherenceShift(
        Vector{Vector{T}}(undef, 0),
        Vector{Vector{T}}(undef, 0)
        )
end

struct CoherencePhase{T} <: PhaseParms{T}
    β::Vector{Vector{T}} # first index for spin systems, second for resonance groups.
    κs_β::Vector{Vector{T}} # first index for spin systems, second for coherence dimension.
end

function CoherencePhase(::Type{T}) where T
    return CoherencePhase(
        Vector{Vector{T}}(undef, 0),
        Vector{Vector{T}}(undef, 0)
        )
end