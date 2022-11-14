


#### top-level

function evalclmolecule(u_rad, A::SHType{T}, B::MoleculeType{T,SST})::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out_sys = evalclspinsystem(u_rad, A.αs, A.Ωs, B.ss_params,
    B.λ0, A.Δc_bar, A.part_inds_molecule)

    out_singlets = evalclsinglets(u_rad, B.d_singlets, A.αs_singlets, A.Ωs_singlets,
    B.β_singlets, B.λ0, B.κs_λ_singlets)

    return out_sys + out_singlets
end

function evalclmixture(u_rad, As::Vector{SHType{T}}, Bs::Vector{MoleculeType{T,SST}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where {T <: Real, SST}

    #u_rad = 2*π*u

    out = zero(Complex{T})
    for n in eachindex(As)

        out += w[n]*evalclmolecule(u_rad, As[n], Bs[n])
    end

    return out
end



#### complex Lorentzian spin system.

"""
Evaluates a part from the partition of a spin group, i.e., a resonance group of a spin group.
r := 2*π*u-d, d is chem shift of this element in radians.
Does not evaluate the complex phase β.
"""
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
    part_inds_molecule)::Complex{T} where T <: Real


    #u_rad = 2*π*u

    out = zero(Complex{T})
    for i in eachindex(αs)
        r = u_rad - x.shift.var[i]

        λ = x.T2.var[i]*λ0
        for k in eachindex(part_inds_molecule[i])
            inds = part_inds_molecule[i][k]

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
    c, part_inds_molecule)::Complex{T} where T <: Real

    #u_rad = 2*π*u

    out = zero(Complex{T})
    for i in eachindex(αs)

        λ = x.T2.var[i]*λ0
        for k in eachindex(part_inds_molecule[i])
            
            r = u_rad - x.shift.d[i][k]
            inds = part_inds_molecule[i][k]

            out += evalclpart(r, αs[i][inds],
                Ωs[i][inds], λ)*cis(dot(x.phase.var[i], c[i][k]))
        end
    end

    return out
end


#### complex Lorentzian singlets.
function evalclsinglets(u_rad::T, d::Vector{T}, αs_singlets::Vector{T}, Ωs_singlets,
    βs_singlets, λ0::T, λ_multipliers::Vector{T}) where T <: Real

    #u_rad = 2*π*u

    out = zero(Complex{T})
    for i in eachindex(αs_singlets)
        τ = u_rad - d[i]

        λ = λ0*λ_multipliers[i]
        Ω = Ωs_singlets[i]
        out += αs_singlets[i]*cis(βs_singlets[i])/(λ+im*(τ-Ω))
    end
    return out
end

#= function evalclsinglets(u_rad::T, d::Vector{T},
    αs_singlets::Vector{T}, Ωs_singlets,
    βs_singlets, λ0::T, λ_multipliers::Vector{T},
    κ_α_singlets::Vector{T}) where T <: Real

    #u_rad = 2*π*u

    out = zero(Complex{T})
    for i in eachindex(αs_singlets)
        τ = u_rad - d[i]

        λ = λ0*λ_multipliers[i]
        Ω = Ωs_singlets[i]
        out += κ_α_singlets[i]*αs_singlets[i]*cis(βs_singlets[i])/(λ+im*(τ-Ω))
    end
    return out
end =#
