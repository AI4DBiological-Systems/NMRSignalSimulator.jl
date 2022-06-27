################### spin system.

function evalclproxysys(qs::Vector{Vector{Function}},
    u_rad::T, x::SpinSysParamsType1{T})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ
    #κs_β = x.κs_β

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    for i = 1:length(qs)
        r = u_rad - d[i]

        for k = 1:length(qs[i])

            out += qs[i][k](r, κs_λ[i])
        end
    end

    ## slower possibly due to r = u_rad - d[i] being evaluated every time qs is called.
    #out = sum( sum(qs[i][k](u_rad - d[i], κs_λ[i]) for k = 1:length(qs[i])) for i = 1:length(qs) )

    return out
end

function evalclproxysys(qs::Vector{Vector{Function}},
    u_rad::T, x::SpinSysParamsType1{T}, κs_α::Vector{Vector{T}})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ
    #κs_β = x.κs_β

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    for i = 1:length(qs)
        r = u_rad - d[i]

        for k = 1:length(qs[i])

            out += κs_α[i][k]*qs[i][k](r, κs_λ[i])
        end
    end

    return out
end

function evalclproxysys(qs::Vector{Vector{Function}},
    u_rad::T, x::SpinSysParamsType2{T})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ
    #κs_β = x.κs_β

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    #u_rad = 2*π*u
    for i = 1:length(qs)

        for k = 1:length(qs[i])
            r = u_rad - d[i][k]

            out += qs[i][k](r, κs_λ[i])
        end
    end

    return out
end

function evalκitpproxysys(κ_α::Vector{Vector{T}}, qs::Vector{Vector{Function}},
    u_rad::T, x::SpinSysParamsType1{T})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ
    #κs_β = x.κs_β

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    #u_rad = 2*π*u
    for i = 1:length(qs)
        r = u_rad - d[i]

        for k = 1:length(qs[i])

            out += κ_α[i][k]*qs[i][k](r, κs_λ[i])
        end
    end

    return out
end

function evalκitpproxysys(κ_α::Vector{Vector{T}}, qs::Vector{Vector{Function}},
    u_rad::T, x::SpinSysParamsType2{T})::Complex{T} where T

    d = x.d
    κs_λ = x.κs_λ

    @assert length(d) == length(qs)

    out = zero(Complex{T})

    #u_rad = 2*π*u
    for i = 1:length(qs)

        for k = 1:length(qs[i])
            r = u_rad - d[i][k]

            out += κ_α[i][k]*qs[i][k](r, κs_λ[i])
        end
    end

    return out
end


###################### front end.

# without κ_α compensation.
function evalclproxymixture(u_rad, As::Vector{SHType{T}},
    Bs::Vector{CompoundType{T,SST}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where {T <: Real,SST}

    #u_rad = 2*π*u

    out = zero(Complex{T})

    for n = 1:length(As)
        out += w[n]*evalclproxycompound(u_rad, As[n], Bs[n])
    end

    return out
end

# with proxy.
function evalclproxycompound(u_rad, A::SHType{T}, B::CompoundType{T,SST})::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out_sys = evalclproxysys(B.qs, u_rad, B.ss_params)

    out_singlets = evalclsinglets(u_rad, B.d_singlets, A.αs_singlets, A.Ωs_singlets,
    B.β_singlets, B.λ0, B.κs_λ_singlets)

    return out_sys + out_singlets
end


# with κ_α compensation.
function evalclproxymixture(u_rad, As::Vector{SHType{T}},
    Es::Vector{καCompoundType{T,SST}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where {T <: Real,SST}

    #u_rad = 2*π*u

    out = zero(Complex{T})

    for n = 1:length(Es)
        out += w[n]*evalclproxycompound(u_rad, As[n], Es[n])
    end

    return out
end

# with κ-proxy.
function evalclproxycompound(u_rad, A::SHType{T}, E::καCompoundType{T,SST})::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out_sys = evalclproxysys(E.core.qs, u_rad, E.core.ss_params, E.κs_α)

    out_singlets = evalclsinglets(u_rad, E.core.d_singlets, A.αs_singlets, A.Ωs_singlets,
    E.core.β_singlets, E.core.λ0, E.core.κs_λ_singlets, E.κs_α_singlets)

    return out_sys + out_singlets
end

function findsimplecoherences(c2::Vector{T};
    atol = 1e-6)::Vector{Int} where T

    c2a = collect( abs.(c2[i]) for i = 1:length(c2))
    sum_c2a = collect(sum(c2a[i]) for i = 1:length(c2a))
    flags = isapprox.(sum_c2a, 1.0, atol = atol)

    return collect(1:length(c2))[flags]
end
