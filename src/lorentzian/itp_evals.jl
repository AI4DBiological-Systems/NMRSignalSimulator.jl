################### spin system.

function evalclproxysys(
    qs::Vector{Vector{Function}},
    u_rad::T,
    x::SpinSysParams{SharedShift{T}, CoherencePhase{T}, SharedT2{T}},
    )::Complex{T} where T

    d = x.shift.var
    κs_λ = x.T2.var
    #κs_β = x.κs_β

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    for i in eachindex(qs)
        r = u_rad - d[i]

        for k in eachindex(qs[i])

            out += qs[i][k](r, κs_λ[i])
        end
    end

    return out
end


function evalclproxysys(
    qs::Vector{Vector{Function}},
    u_rad::T,
    x::SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}},
    )::Complex{T} where T

    d = x.shift.d
    κs_λ = x.T2.var
    #κs_β = x.κs_β

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    #u_rad = 2*π*u
    for i in eachindex(qs)

        for k in eachindex(qs[i])
            r = u_rad - d[i][k]

            out += qs[i][k](r, κs_λ[i])
        end
    end

    return out
end

###################### front end.

# without κ_α compensation.
function evalclproxymixture(u_rad, As::Vector{SHType{T}},
    Bs::Vector{MoleculeType{T,SST}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where {T <: Real,SST}

    #u_rad = 2*π*u

    out = zero(Complex{T})

    for n in eachindex(As)
        out += w[n]*evalclproxymolecule(u_rad, As[n], Bs[n])
    end

    return out
end

# with proxy.
function evalclproxymolecule(u_rad, A::SHType{T}, B::MoleculeType{T,SST})::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out_sys = evalclproxysys(B.qs, u_rad, B.ss_params)

    out_singlets = evalclsinglets(u_rad, B.d_singlets, A.αs_singlets, A.Ωs_singlets,
    B.β_singlets, B.λ0, B.κs_λ_singlets)

    return out_sys + out_singlets
end


function findsimplecoherences(c2::Vector{T};
    atol = 1e-6)::Vector{Int} where T

    c2a = collect( abs.(c2[i]) for i in eachindex(c2))
    sum_c2a = collect(sum(c2a[i]) for i in eachindex(c2a))
    flags = isapprox.(sum_c2a, 1.0, atol = atol)

    return collect(1:length(c2))[flags]
end
