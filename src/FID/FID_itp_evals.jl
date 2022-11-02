
function evalFIDproxymixture(t, As::Vector{SHType{T}},
    Bs::Vector{MoleculeType{T,SST}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where {T <: Real,SST}

    out = zero(Complex{T})

    for n in eachindex(As)
        out += w[n]*evalFIDproxymolecule(t, As[n], Bs[n])
    end

    return out
end

function evalFIDproxymolecule(t, A::SHType{T}, B::MoleculeType{T,SST})::Complex{T} where {T <: Real, SST}

    out_sys = evalFIDproxysys(B.qs, t, B.ss_params, B.λ0)

    out_singlets = evalFIDsinglets(t, B.d_singlets, A.αs_singlets, A.Ωs_singlets,
    B.β_singlets, B.λ0, B.κs_λ_singlets)

    return out_sys + out_singlets
end


function evalFIDproxysys(qs::Vector{Vector{Function}},
    t::T, x::SharedShift{T}, λ0)::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    for i in eachindex(qs)
        r = d[i]

        sys_sum = zero(Complex{T})
        for k in eachindex(qs[i])

            sys_sum += qs[i][k](r, t)
        end
        out += sys_sum*exp(-x.κs_λ[i]*λ0*t)
    end

    return out
end

function evalFIDproxysys(qs::Vector{Vector{Function}},
    t::T, x::CoherenceShift{T}, λ0)::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    for i in eachindex(qs)
        sys_sum = zero(Complex{T})

        for k in eachindex(qs[i])
            r = d[i][k]

            sys_sum += qs[i][k](r, t)
        end
        out += sys_sum*exp(-x.κs_λ[i]*λ0*t)
    end

    return out
end
