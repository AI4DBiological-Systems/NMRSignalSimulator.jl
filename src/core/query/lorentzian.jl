


#### top-level

function evalclmolecule(u_rad, A::HAM.SHType{T}, B::MoleculeType{T,SST})::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out_sys = evalclspinsystem(u_rad, A.αs, A.Ωs, B.ss_params,
    B.λ0, A.Δc_bar, A.parts)

    # out_singlets = evalclsinglets(u_rad, B.ζ_singlets, A.αs_singlets, A.Ωs_singlets,
    # B.β_singlets, B.λ0, B.κs_λ_singlets)

    return out_sys #+ out_singlets
end

"""
```
function evalclmixture(
    u_rad,
    As::Vector{HAM.SHType{T}},
    Bs::Vector{MoleculeType{T,SST}};
    w::Vector{T} = ones(T, length(As)),
)::Complex{T} where {T <: Real, SST}
```

Evaluates the preliminary model at radial frequency (in radians) `u_rad`.
Inputs:

- `Bs` is the surrogate model.
- `As` is the resonance component data structure.
- `w` is relative concentration of each compound entry.
"""
function evalclmixture(
    u_rad,
    As::Vector{HAM.SHType{T}},
    Bs::Vector{MoleculeType{T,SST}};
    w::Vector{T} = ones(T, length(As)),
    )::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out = zero(Complex{T})
    for n in eachindex(As)

        out += w[n]*evalclmolecule(u_rad, As[n], Bs[n])
    end

    return out
end



#### complex Lorentzian spin system.

# """
# Evaluates a part from the partition of a spin group, i.e., a resonance group of a spin group.
# r := 2*π*u-d, d is chem shift of this element in radians.
# Does not evaluate the complex phase β.
# """
function evalclpart(r,
    α::Vector{T}, Ω::Vector{T}, λ::T)::Complex{T} where T <: Real

    out = sum( α[l]/(λ+im*(r-Ω[l])) for l in eachindex(α) )

    return out
end

function evalclspinsystem(
    u_rad,
    αs::Vector{Vector{T}},
    Ωs::Vector{Vector{T}},
    x::SpinSysParams{SharedShift{T}, CoherencePhase{T}, SharedT2{T}},
    λ0::T,
    c,
    parts)::Complex{T} where T <: Real


    #u_rad = 2*π*u

    out = zero(Complex{T})
    for i in eachindex(αs)
        r = u_rad - x.shift.var[i]

        λ = x.T2.λ[i]
        for k in eachindex(parts[i])
            inds = parts[i][k]

            out += evalclpart(r, αs[i][inds],
                Ωs[i][inds], λ)*cis(dot(x.phase.var[i], c[i][k]))
        end
    end

    return out
end

function evalclspinsystem(u_rad,
    αs::Vector{Vector{T}}, Ωs::Vector{Vector{T}},
    x::SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}},
    λ0::T,
    c, parts)::Complex{T} where T <: Real

    #u_rad = 2*π*u

    out = zero(Complex{T})
    for i in eachindex(αs)

        λ = x.T2.λ[i]
        for k in eachindex(parts[i])
            
            r = u_rad - x.shift.ζ[i][k]
            inds = parts[i][k]

            out += evalclpart(r, αs[i][inds],
                Ωs[i][inds], λ)*cis(dot(x.phase.var[i], c[i][k]))
        end
    end

    return out
end


#### complex Lorentzian singlets.
# function evalclsinglets(
#     u_rad::T,
#     d::Vector{T},
#     αs_singlets::Vector{T},
#     Ωs_singlets,
#     βs_singlets,
#     λ0::T,
#     λ_multipliers::Vector{T}
#     ) where T <: Real

#     #u_rad = 2*π*u

#     out = zero(Complex{T})
#     for i in eachindex(αs_singlets)
#         ζ = u_rad - d[i]

#         λ = λ0*λ_multipliers[i]
#         Ω = Ωs_singlets[i]
#         out += αs_singlets[i]*cis(βs_singlets[i])/(λ+im*(ζ-Ω))
#     end
#     return out
# end

# arbitrary phase.
function evalcomplexlorentzian(
    r::T,
    λ::T,
    α::T,
    Ω::T,
    cos_β::T,
    sin_β::T,
    )::Complex{T} where T <: AbstractFloat

    #u_rad = 2*π*u

    # common quantities:
    τ = r-Ω
    denominator = λ^2 + τ^2

    # numerators and denoimators.
    real_n = α*( sin_β*τ + λ*cos_β )
    imag_n = α*( -cos_β*τ + λ*sin_β )

    return Complex(real_n/denominator, imag_n/denominator)
end

# zero phase version.
function evalcomplexlorentzian(
    r::T,
    λ::T,
    α::T,
    Ω::T,
    )::Complex{T} where T <: AbstractFloat

    return evalcomplexlorentzian(r, λ, α, Ω, one(T), zero(T))
end

## arbitrary phase and zero phase.
function evalabsorptionlorentzian(
    r::T,
    λ::T,
    α::T,
    Ω::T,
    cos_β::T,
    sin_β::T,
    )::T where T

    # common quantities:
    τ = r-Ω
    denominator = λ^2 + τ^2

    # numerators and denoimators.
    real_n = α*( sin_β*τ + λ*cos_β )
    #imag_n = α*( -cos_β*τ + λ*sin_β )

    return real_n/denominator
end

function evalabsorptionlorentzian(
    r::T,
    λ::T,
    α::T,
    Ω::T,
    )::T where T <: AbstractFloat

    return evalabsorptionlorentzian(r, λ, α, Ω, one(T), zero(T))
end

function evaldispersionlorentzian(
    r::T,
    λ::T,
    α::T,
    Ω::T,
    cos_β::T,
    sin_β::T,
    )::T where T <: AbstractFloat

    # common quantities:
    τ = r-Ω
    denominator = λ^2 + τ^2

    # numerators and denoimators.
    #real_n = α*( sin_β*τ + λ*cos_β )
    imag_n = α*( -cos_β*τ + λ*sin_β )

    return imag_n/denominator
end

function evaldispersionlorentzian(
    r::T,
    λ::T,
    α::T,
    Ω::T,
    )::T where T <: AbstractFloat

    return evaldispersionlorentzian(r, λ, α, Ω, one(T), zero(T))
end

function evalabsorptionlorentzianderivatives!(
    out::Vector{T},
    r::T,
    λ::T,
    α::T,
    Ω::T,
    c::T, #cos_β::T,
    s::T, #sin_β::T,
    ) where T

    # common quantities:
    #τ = Ω-r
    U = -r^2 - Ω^2 + λ^2 + 2*r*Ω
    Z =  r^2 + Ω^2 + λ^2 - 2*r*Ω
    factor = α/Z^2

    # numerators
    num_dq_dr =  U*s + 2*c*λ*( -r + Ω )
    num_dq_dλ = -U*c + 2*s*λ*( -r + Ω )

    out[begin] = num_dq_dr*factor
    out[begin+1] = num_dq_dλ*factor

    return nothing
end

function evalabsorptionlorentzianderivatives!(
    out::Vector{T},
    r::T,
    λ::T,
    α::T,
    Ω::T,
    ) where T

    return evalabsorptionlorentzianderivatives!(out, r, λ, α, Ω, one(T), zero(T))
end

function evaldispersionlorentzianderivatives!(
    out::Vector{T},
    r::T,
    λ::T,
    α::T,
    Ω::T,
    c::T, #cos_β::T,
    s::T, #sin_β::T,
    ) where T <: AbstractFloat

    # common quantities:
    #τ = Ω-r
    U = -r^2 - Ω^2 + λ^2 + 2*r*Ω
    Z =  r^2 + Ω^2 + λ^2 - 2*r*Ω
    factor = α/Z^2

    # numerators
    num_dq_dr = -U*c + 2*s*λ*( -r + Ω )
    num_dq_dλ = -U*s - 2*c*λ*( -r + Ω )

    out[begin] = num_dq_dr*factor
    out[begin+1] = num_dq_dλ*factor
    
    return nothing
end

function evaldispersionlorentzianderivatives!(
    out::Vector{T},
    r::T,
    λ::T,
    α::T,
    Ω::T,
    ) where T <: AbstractFloat

    return evaldispersionlorentzianderivatives!(out, r, λ, α, Ω, one(T), zero(T))
end

## without phase (this is actually used as input to evalq)
