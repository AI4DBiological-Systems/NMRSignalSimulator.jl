

### model.

# # free-induction decay (FID).
@kwdef struct FIDSurrogateConfig{T}
    Δt::T = convert(T, 1e-5)
    t_lb::T = zero(T)
    t_ub::T = convert(T, 3.0)
end

abstract type OperationRange end

struct FIDOperationRange{T} <: OperationRange
    t_lb::T
    t_ub::T
    Δt::T
end

struct FIDInterpolationSamples{T <: AbstractFloat}

    samples::Vector{Vector{Complex{T}}} # [resonance group index][sample index]

    t_lb::T
    Δt::T
    t_ub::T

    #λ0::T
end

# default to invalid values. At deserialization, check value for λ0. If is negative, then we know it is the singlet case; no interpolation surrogates are used.
function FIDInterpolationSamples(::Type{T}, N_groups::Int) where T
    return FIDInterpolationSamples(
        collect( ones(Complex{T}, 1) for _ = 1:N_groups ),
        (-ones(T, 3))...
    )
end

function FIDInterpolationSamples(
    C::FIDSurrogateConfig{T},
    s::Vector{Vector{Complex{T}}},
    )::FIDInterpolationSamples{T} where T <: AbstractFloat
    #
    return FIDInterpolationSamples(s, C.t_lb, C.Δt, C.t_ub)
end


# # complex Lorentzian (CL).
"""
    struct CLSurrogateConfig{T}
        Δr::T = convert(T, 1.0)
        Δκ_λ::T = convert(T, 0.05)
        Δcs_max_scalar::T = convert(T, 0.2)
        κ_λ_lb::T = convert(T, 0.5)
        κ_λ_ub::T = convert(T, 2.5)

        use_compound_freqs::Bool = false
        default_ppm_padding::T = convert(T , 0.5)
    end

- `κ_λ_lb`, `κ_λ_ub` -- the default lower and upper bounds, respectively, for the κ_λ input of the surrogate.

- `Δcs_max_scalar`, `default_ppm_padding` -- specifies the surrogate's default frequency operating range, [u_min - delta_u, u_max + delta_u], where delta_u is Δcs_max_scalar expressed in Hz.
u_min, u_max are in Hz, Δcs_max_scalar is in ppm. `u_min` is set to the smallest resonance frequency from `As`, and `u_max` is set to the largest resonance frequency from `As`.
`Δcs_max_scalar` is in units of ppm. It is the border that is added to u_min and u_max (once they are converted to ppm) to get the final frequency interval for which the surrogate of a spin system is fitted to.
The extrapolation of the surrogate for frequency queries outisde this interval is set to return zero, so the surrogate is the zero signal outside this interval for that spin system.

- `use_compound_freqs`, if set to `true`, the surrogate for each compound is valid only over +/- `Δcs_max_scalar` (ppm units) from the comopund's minimum and maximum resonance frequencies.
If set to `false`, the surrogate for all compounds is valid between the minimum and maximum resonance frequencies of the entire mixture.

- `Δr`, `Δκ_λ` -- the sampling increment for the frequency input r and T2 multiplier input κ_λ for generating samples to fit the surrogate. Smaller means the surrogate is more accurate, but slower to construct the surrogate.
"""
@kwdef struct CLSurrogateConfig{T}
    Δr::T = convert(T, 1.0) # the samples used to build the surrogate is taken every `Δr` radian on the frequency axis. Decrease for improved accuracy at the expense of computation resources.
    Δκ_λ::T = convert(T, 0.05) # the samples used to build thes urrogate for κ_λ are taken at this sampling spacing. Decrease for improved accuracy at the expense of computation resources.
    Δcs_max_scalar::T = convert(T, 0.2) # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.
    κ_λ_lb::T = convert(T, 0.5) # interpolation lower limit for κ_λ.
    κ_λ_ub::T = convert(T, 2.5) # interpolation upper limit for κ_λ.

    use_compound_freqs::Bool = false # if true, do not fit surrogate over the ensemble of frequency interval of the entire mixture, only individual compounds.
    # if true, save memory and faster to fit.

    ppm_padding::T = convert(T , 0.5)
end

# the ingredients for making a surrogate spin system.
struct CLSurrogateSpinSysInputs{T <: AbstractFloat}
    parts::Vector{Vector{Vector{Int}}}
    αs::Vector{Vector{T}}
    Ωs::Vector{Vector{T}}

    u_min::T
    u_max::T
end


"""
```
struct CLOperationRange{T}
    u_min::T # in Hz.
    u_max::T # in Hz.
    
    # d_max = ppm2hzfunc.(Δcs_max) .- ppm2hzfunc(zero(T))
    d_max::Vector{T} # in Hz. ζ is in rad.

    κ_λ_lb::T
    κ_λ_ub::T
    Δr::T
    Δκ_λ::T
end
```
Operation range for the spline surrogate q(r,λ) of a resonance group.
Bounds on `r` for groups in the `i`-th spin system:
```
r_min = 2*π*(u_min - d_max[i])
r_max = 2*π*(u_max + d_max[i])
```
Lower and upper bounds on `λ` are `κ_λ_lb` and `κ_λ_ub`, respectively.
Able to serialize/deserialize.
"""
struct CLOperationRange{T} <: OperationRange
    u_min::T # in Hz.
    u_max::T # in Hz.
    
    # d_max = ppm2hzfunc.(Δcs_max) .- ppm2hzfunc(zero(T))
    d_max::Vector{T} # in Hz. ζ is in rad.

    κ_λ_lb::T
    κ_λ_ub::T
    Δr::T
    Δκ_λ::T
end


function CLOperationRange(u_min::T, u_max::T, d_max::Vector{T}, C::CLSurrogateConfig{T})::CLOperationRange{T} where T <: AbstractFloat
    return CLOperationRange(u_min, u_max, d_max, C.κ_λ_lb, C.κ_λ_ub, C.Δr, C.Δκ_λ)
end

# for one spin system.
# disregard is λ0 is a negative value.
struct InterpolationSamples{T <: AbstractFloat}

    samples::Vector{Matrix{Complex{T}}} # [resonance group index][sample index]

    #A_r = r_min:Δr:r_max
    #A_λ = (κ_λ_lb:Δκ_λ:κ_λ_ub) .* λ0
    r_min::T
    Δr::T
    r_max::T

    κ_λ_lb::T
    Δκ_λ::T
    κ_λ_ub::T
    λ0::T
end

# default to invalid values. At deserialization, check value for λ0. If is negative, then we know it is the singlet case; no interpolation surrogates are used.
function InterpolationSamples(::Type{T}, N_groups::Int) where T
    return InterpolationSamples(
        collect( ones(Complex{T}, 1,1) for _ = 1:N_groups ),
        (-ones(T, 7))...
    )
end

function InterpolationSamples(
    C::CLSurrogateConfig{T},
    s::Vector{Matrix{Complex{T}}},
    r_min::T,
    r_max::T,
    λ0::T,
    )::InterpolationSamples{T} where T <: AbstractFloat
    #
    return InterpolationSamples(s, r_min, C.Δr, r_max, C.κ_λ_lb, C.Δκ_λ, C.κ_λ_ub, λ0)
end

# This is without the compensation amplitude parameter, κ_α, denoted κs_α in code.
struct MoleculeType{T,SST, OT <: OperationRange} # parameters for surrogate model.

    # simulator for each resonance group. [i][k] is i-th spin system, k-th group.
    qs::Vector{Vector{Function}} # spin group, partition element index. range is complex numbers.

    # parameters.
    ss_params::SST # include singlets. SST <: SpinSysParams
    
    # operation range of spline surrogate.
    op_range::OT

    # base T2. formula: λ = ξ*λ0.
    λ0::T
end


### different parameterizations of the spin system FID parameters.

abstract type MoleculeParams end

# abstract type ShiftParms{T} <: MoleculeParms{T} end
# abstract type PhaseParms{T} <: MoleculeParms{T} end
# abstract type T2Parms{T} <: MoleculeParms{T} end

# struct CoherenceShift{T} <: SpinSysParams{T}
#     κs_λ::Vector{T} # a multiplier for each (spin group.
#     κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
#     d::Vector{Vector{T}} # a multiplier for each (spin group, partition element).
#     κs_ζ::Vector{Vector{T}} # same size as κs_β. # intermediate buffer for d.
# end

# struct SharedParams <: T2Parms{T}
#     var::Vector{T} # multiplier wrt some λ0. length: number of spin groups.

# end
# function SharedParams(
#     ::Type{T},
#     N_sys::Int;
#     default_value::T = one(T),
#     )::SharedParams where T <: AbstractFloat

#     var = ones(T, N_sys)
#     fill!(var, default_value)

#     return SharedParams(var)
# end
abstract type SharedParams <: MoleculeParams end
abstract type CoherenceParams <: MoleculeParams end

struct SharedT2{T} <: SharedParams
    var::Vector{T} # multiplier wrt some λ0. length: number of spin groups.
    λ::Vector{T} # actual T2. length: number of spin groups.
end

function SharedT2(
    #::Type{T},
    λ0::T,
    N_sys::Int,
    )::SharedT2{T} where T <: AbstractFloat

    κs_λ = ones(T, N_sys)
    #λ = ones(T, N_sys)
    λ = λ0 .* κs_λ

    return SharedT2(κs_λ, λ)
    #return SharedT2(κs_λ)
end

struct CoherenceT2{T} <: CoherenceParams
    var::Vector{Vector{T}} # [spin system][coherence dimension], multiplier.
    λ::Vector{Vector{T}} # [spin system][resonance group]. actual T2.
end

function CoherenceT2(
    λ0::T,
    N_coherence_vars_sys::Vector{Int},
    N_resonance_groups_sys::Vector{Int},
    )::CoherenceT2{T} where T <: AbstractFloat

    κs_λ = collect( ones(T, N_coherence_vars_sys[i]) for i in eachindex(N_coherence_vars_sys))
    λ = collect( λ0 .* ones(T, N_resonance_groups_sys[i]) for i in eachindex(N_resonance_groups_sys))

    return CoherenceT2(κs_λ, λ)
end

struct SharedShift{T} <: SharedParams
    var::Vector{T} # length: number of spin groups.
    #d::Vector{T} # length: number of spin groups.
    #d::Vector{Vector{T}}
end

function SharedShift(
    ::Type{T},
    N_sys::Int,
    )::SharedShift{T} where T <: AbstractFloat

    ζ = zeros(T, N_sys)

    return SharedShift(ζ)
end

struct CoherenceShift{T} <: CoherenceParams
    var::Vector{Vector{T}} # [spin system][coherence dimension]
    #κs_ζ::Vector{Vector{T}} # first index for spin systems, second for coherence dimension.
    ζ::Vector{Vector{T}} # [spin system][resonance group]
end

function CoherenceShift(
    ::Type{T},
    N_coherence_vars_sys::Vector{Int},
    N_resonance_groups_sys::Vector{Int},
    )::CoherenceShift{T} where T <: AbstractFloat

    κs_ζ = collect( zeros(T, N_coherence_vars_sys[i]) for i in eachindex(N_coherence_vars_sys))
    ζ = collect( zeros(T, N_resonance_groups_sys[i]) for i in eachindex(N_resonance_groups_sys))

    return CoherenceShift(κs_ζ, ζ)
end

function CoherenceShift(Δc_bar::Vector{Vector{Vector{T}}})::CoherenceShift{T} where T

    N_coherence_vars_sys, N_resonance_groups_sys = getNcoherencesys(Δc_bar)

    return CoherenceShift(
        T,
        N_coherence_vars_sys,
        N_resonance_groups_sys,
    )
end

struct CoherencePhase{T} <: CoherenceParams
    var::Vector{Vector{T}} # first index for spin systems, second for coherence dimension.
    β::Vector{Vector{T}} # first index for spin systems, second for coherence dimension.
    cos_β::Vector{Vector{T}} # first index for spin systems, second for resonance groups.
    sin_β::Vector{Vector{T}}
end

function CoherencePhase(
    ::Type{T},
    N_coherence_vars_sys::Vector{Int},
    N_resonance_groups_sys::Vector{Int},
    )::CoherencePhase{T} where T <: AbstractFloat

    κs_β = collect( zeros(T, N_coherence_vars_sys[i]) for i in eachindex(N_coherence_vars_sys))
    β = collect( zeros(T, N_resonance_groups_sys[i]) for i in eachindex(N_resonance_groups_sys))
    cos_β = collect( ones(T, N_resonance_groups_sys[i]) for i in eachindex(N_resonance_groups_sys))
    sin_β = collect( zeros(T, N_resonance_groups_sys[i]) for i in eachindex(N_resonance_groups_sys))

    return CoherencePhase(κs_β, β, cos_β, sin_β)
end

# include singlets.
struct SpinSysParams{ST,PT,DT}
    shift::ST
    phase::PT
    T2::DT
end

function setupSSParamsparams(
    ::Type{SpinSysParams{SharedShift{T}, CoherencePhase{T}, SharedT2{T}}},
    N_coherence_vars_sys::Vector{Int},
    N_resonance_groups_sys::Vector{Int},
    λ0,
    )::SpinSysParams{SharedShift{T}, CoherencePhase{T}, SharedT2{T}} where T <: AbstractFloat

    N_sys = length(N_coherence_vars_sys)

    return SpinSysParams(
        SharedShift(T, N_sys),
        CoherencePhase(T, N_coherence_vars_sys, N_resonance_groups_sys),
        SharedT2(λ0, N_sys),
    )
end

function setupSSParamsparams(
    ::Type{SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}},
    N_coherence_vars_sys::Vector{Int},
    N_resonance_groups_sys::Vector{Int},
    λ0,
    )::SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}} where T <: AbstractFloat

    N_sys = length(N_coherence_vars_sys)

    return SpinSysParams(
        CoherenceShift(T, N_coherence_vars_sys, N_resonance_groups_sys),
        CoherencePhase(T, N_coherence_vars_sys, N_resonance_groups_sys),
        SharedT2(λ0, N_sys),
    )
end

function setupSSParamsparams(
    ::Type{SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, CoherenceT2{T}}},
    N_coherence_vars_sys::Vector{Int},
    N_resonance_groups_sys::Vector{Int},
    λ0,
    )::SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, CoherenceT2{T}} where T <: AbstractFloat

    N_sys = length(N_coherence_vars_sys)

    return SpinSysParams(
        CoherenceShift(T, N_coherence_vars_sys, N_resonance_groups_sys),
        CoherencePhase(T, N_coherence_vars_sys, N_resonance_groups_sys),
        CoherenceT2(λ0, N_coherence_vars_sys, N_resonance_groups_sys),
    )
end

function getNcoherencesys(Δc_bar::Vector{Vector{Vector{T}}})::Tuple{Vector{Int},Vector{Int}} where T

    N_coherence_vars_sys = collect( length(Δc_bar[i][begin]) for i in eachindex(Δc_bar) )
    N_resonance_groups_sys = collect( length(Δc_bar[i]) for i in eachindex(Δc_bar) )

    return N_coherence_vars_sys, N_resonance_groups_sys
end

"""
```
function getSpinSysParamsdatatype(
    ::Type{CoherenceShift{T}}
    )::DataType where T <: AbstractFloat
```

Convinence constructor for SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}
"""
function getSpinSysParamsdatatype(
    ::Type{CoherenceShift{T}}
    )::DataType where T <: AbstractFloat

    return SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}
end

"""
```
function getSpinSysParamsdatatype(
    ::Type{SharedShift{T}}
    )::DataType where T <: AbstractFloat
```

Convinence constructor for SpinSysParams{SharedShift{T}, CoherencePhase{T}, SharedT2{T}}
"""
function getSpinSysParamsdatatype(
    ::Type{SharedShift{T}}
    )::DataType where T <: AbstractFloat

    return SpinSysParams{SharedShift{T}, CoherencePhase{T}, SharedT2{T}}
end


### fetch intermediate parameters.

function fetchphase(
    phases::Vector{CoherencePhase{T}},
    n::Int,
    i::Int,
    k::Int,
    )::Tuple{T,T} where T

    return phases[n].cos_β[i][k], phases[n].sin_β[i][k]
end

function fetchshift(
    shifts::Vector{CoherenceShift{T}},
    n::Int,
    i::Int,
    k::Int,
    )::T where T

    return shifts[n].ζ[i][k]
end

function fetchshift(
    shifts::Vector{SharedShift{T}},
    n::Int,
    i::Int,
    args...
    )::T where T

    return shifts[n].var[i]
end

function fetchT2(
    T2s::Vector{SharedT2{T}},
    n::Int,
    i::Int,
    args...
    )::T where T

    return T2s[n].λ[i]
end

abstract type SystemsTrait end
struct Coherence <: SystemsTrait end
struct Shared <: SystemsTrait end

# # no Default behavior.
# SystemsTrait(::Type) = Nothing

# gradient buffers for mixtures.
#SystemsTrait(::Type{<:Vector{T}}) where T <: AbstractFloat = Shared()
#SystemsTrait(::Type{<:Vector{Vector{T}}}) where T <: AbstractFloat = Coherence()


# parameters for mixtures.
SystemsTrait(::Y) where Y <: MoleculeParams = Shared() # default for super type MoleculeParams.
SystemsTrait(::Y) where Y <: SharedParams = Shared()
SystemsTrait(::Y) where Y <: CoherenceParams = Coherence()



# parameters for mixtures. won't work since Vector{abstract type}.
# SystemsTrait(::Type{<:Vector{Y}}) where Y <: MoleculeParams = Shared() # default for super type MoleculeParams.
# SystemsTrait(::Type{<:Vector{Y}}) where Y <: SharedParams = Shared()
# SystemsTrait(::Type{<:Vector{Y}}) where Y <: CoherenceParams = Coherence()

######## for construction a table for visualization and analysis.
abstract type TableConstructionTrait end
struct ProcessShifts <: TableConstructionTrait end
