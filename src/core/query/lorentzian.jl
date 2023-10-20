


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
