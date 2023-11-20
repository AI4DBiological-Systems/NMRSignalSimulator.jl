
function evalfidmolecule(t, A::HAM.SHType{T}, B::MoleculeType{T,SST,FIDOperationRange{T}})::Complex{T} where {T <: Real, SST}

    out_sys = evalfidspinsystem(t, A.αs, A.Ωs, B.ss_params, A.Δc_bar, A.parts)

    return out_sys
end

function evalfidmixture(
    t,
    As::Vector{HAM.SHType{T}},
    Bs::Vector{MoleculeType{T,SST, FIDOperationRange{T}}};
    w::Vector{T} = ones(T, length(As)),
    )::Complex{T} where {T <: Real, SST}

    out = zero(Complex{T})
    for n in eachindex(As)

        out += w[n]*evalfidmolecule(t, As[n], Bs[n])
    end

    return out
end




#### FID spin system.

function evalfidpartitionelement(t, α::Vector{T}, Ω::Vector{T})::Complex{T} where T <: Real

    out = sum( α[l]*cis(Ω[l]*t) for l in eachindex(α) )
    return out
end

function evalfidspinsystem(
    t,
    αs::Vector{Vector{T}},
    Ωs::Vector{Vector{T}},
    x::SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}},
    c::Vector{Vector{Vector{T}}}, # A.Δc_bar.
    parts::Vector{Vector{Vector{Int64}}}, # A.parts
    )::Complex{T} where T <: Real

    rt = zero(T) # pre-allocate.
    inner_sum = zero(Complex{T}) # pre-allocate.

    out = zero(Complex{T})
    for i in eachindex(αs)

        λ = x.T2.λ[i]
        inner_sum = zero(Complex{T})
        for k in eachindex(parts[i])
            inds = parts[i][k]

            rt = x.shift.ζ[i][k]*t
            inner_sum += evalfidpartitionelement(
                t,
                αs[i][inds],
                Ωs[i][inds],
            )*cis(
                dot(
                    x.phase.var[i],
                    c[i][k],
                ) + rt,
            )
        end
        out += (exp(-λ*t)*inner_sum)
    end

    return out
end


function evalfidspinsystem(
    t,
    αs::Vector{Vector{T}},
    Ωs::Vector{Vector{T}},
    x::SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, CoherenceT2{T}},
    c::Vector{Vector{Vector{T}}}, # A.Δc_bar.
    parts::Vector{Vector{Vector{Int64}}}, # A.parts
    )::Complex{T} where T <: Real

    λ = x.T2.λ

    rt = zero(T) # pre-allocate.
    inner_sum = zero(Complex{T}) # pre-allocate.

    out = zero(Complex{T})
    for i in eachindex(αs)
        
        inner_sum = zero(Complex{T})
        for k in eachindex(parts[i])
            inds = parts[i][k]

            rt = x.shift.ζ[i][k]*t
            inner_sum += evalfidpartitionelement(
                t,
                αs[i][inds],
                Ωs[i][inds],
            )*cis(
                dot(
                    x.phase.var[i],
                    c[i][k],
                ) + rt,
            )*exp(-λ[i][k]*t)
        end
        out += inner_sum
    end

    return out
end

########## proxy


function evalfidproxymixture(
    t,
    Bs::Vector{MoleculeType{T,SST,FIDOperationRange{T}}};
    w::Vector{T} = ones(T, length(Bs)),
    )::Complex{T} where {T <: Real,SST}

    out = zero(Complex{T})

    for n in eachindex(Bs)
        out += w[n]*evalfidproxymolecule(t, Bs[n])#*exp(-Bs[n].λ0^2*t)
    end

    return out
end


function evalfidproxymolecule(t, B::MoleculeType{T,SST,FIDOperationRange{T}})::Complex{T} where {T <: Real, SST}

    return evalfidproxysys(B.qs, t, B.ss_params)
end

function evalfidproxysys(
    qs::Vector{Vector{Function}},
    t::T,
    x::SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}},
    )::Complex{T} where T

    ζ = x.shift.ζ
    λ = x.T2.λ

    @assert length(ζ) == length(qs)

    out = zero(Complex{T})
    for i in eachindex(qs)
        sys_sum = zero(Complex{T})

        for k in eachindex(qs[i])
            r = ζ[i][k]

            sys_sum += qs[i][k](t, r)
        end
        out += sys_sum*exp(-λ[i]*t)
    end

    return out
end

function evalfidproxysys(
    qs::Vector{Vector{Function}},
    t::T,
    x::SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, CoherenceT2{T}},
    )::Complex{T} where T

    ζ = x.shift.ζ
    λ = x.T2.λ

    @assert length(ζ) == length(qs)

    out = zero(Complex{T})
    for i in eachindex(qs)
        sys_sum = zero(Complex{T})

        for k in eachindex(qs[i])
            r = ζ[i][k]
            

            sys_sum += qs[i][k](t, r)*exp(-λ[i][k]*t)
        end
        out += sys_sum
    end

    return out
end