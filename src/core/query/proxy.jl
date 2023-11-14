################### spin system.

function evalclproxysys(
    qs::Vector{Vector{Function}},
    u_rad::T,
    x::SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}},
    )::Complex{T} where T

    ζ = x.shift.ζ
    λ = x.T2.λ

    @assert length(ζ) == length(qs)

    out = zero(Complex{T})
    for i in eachindex(qs)

        for k in eachindex(qs[i])
            r = u_rad - ζ[i][k]

            out += qs[i][k](r, λ[i])
        end
    end

    return out
end

###################### front end.

"""
```
function evalclproxymixture(
    u_rad,
    As::Vector{HAM.SHType{T}},
    Bs::Vector{MoleculeType{T,SST,CLOperationRange{T}}};
    w::Vector{T} = ones(T, length(As)),
)::Complex{T} where {T <: Real,SST}
```

Evaluates the surrogate model at radial frequency (in radians) `u_rad`.
Inputs:

- `Bs` is the surrogate model.
- `As` is the resonance component data structure.
- `w` is relative concentration of each compound entry.
"""
function evalclproxymixture(
    u_rad,
    Bs::Vector{MoleculeType{T,SST,CLOperationRange{T}}};
    w::Vector{T} = ones(T, length(Bs)),
    )::Complex{T} where {T <: Real,SST}

    out = zero(Complex{T})

    for n in eachindex(Bs)
        out += w[n]*evalclproxymolecule(u_rad, Bs[n])
    end

    return out
end

# with proxy.
function evalclproxymolecule(u_rad, B::MoleculeType{T,SST,CLOperationRange{T}})::Complex{T} where {T <: Real, SST}

    return evalclproxysys(B.qs, u_rad, B.ss_params)
end


function findsimplecoherences(c2::Vector{T};
    atol = 1e-6)::Vector{Int} where T

    c2a = collect( abs.(c2[i]) for i in eachindex(c2))
    sum_c2a = collect(sum(c2a[i]) for i in eachindex(c2a))
    flags = isapprox.(sum_c2a, 1.0, atol = atol)

    return collect(1:length(c2))[flags]
end
