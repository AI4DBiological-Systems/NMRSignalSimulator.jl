
function evalFIDmolecule(t, A::SHType{T}, B::MoleculeType{T,SST})::Complex{T} where {T <: Real, SST}

    out_sys = evalFIDspinsystem(t, A.αs, A.Ωs, B.ss_params,
    B.λ0, A.Δc_bar, A.part_inds_molecule)

    out_singlets = evalFIDsinglets(t, B.d_singlets, A.αs_singlets, A.Ωs_singlets,
    B.β_singlets, B.λ0, B.κs_λ_singlets)

    return out_sys + out_singlets
end

function evalFIDmixture(t, As::Vector{SHType{T}}, Bs::Vector{MoleculeType{T,SST}};
    w::Vector{T} = ones(T, length(As)))::Complex{T} where {T <: Real, SST}

    out = zero(Complex{T})
    for n = 1:length(As)

        out += w[n]*evalFIDmolecule(t, As[n], Bs[n])
    end

    return out
end




#### FID spin system.

"""
Evaluates an element from the partition of a spin group.
r := 2*π*u-d, d is chem shift of this element in radians.
Does not evaluate the complex phase β.
"""
function evalFIDpartitionelement(t,
    α::Vector{T}, Ω::Vector{T}, r::T)::Complex{T} where T <: Real

    out = sum( α[l]*cis((Ω[l]-r)*t) for l = 1:length(α) )

    return out
end

function evalFIDpartitionelement(t,
    α::Vector{T}, Ω::Vector{T})::Complex{T} where T <: Real

    out = sum( α[l]*cis(Ω[l]*t) for l = 1:length(α) )

    return out
end


function evalFIDspinsystem(t,
    αs::Vector{Vector{T}}, Ωs::Vector{Vector{T}},
    x::SpinSysParamsType1{T}, λ0::T,
    c, part_inds_molecule)::Complex{T} where T <: Real

    out = zero(Complex{T})
    for i = 1:length(αs)

        rt = x.d[i]*t

        sys_sum = zero(Complex{T})

        for k = 1:length(part_inds_molecule[i])
            inds = part_inds_molecule[i][k]

            sys_sum += evalFIDpartitionelement(t, αs[i][inds],
                Ωs[i][inds])*cis(dot(x.κs_β[i], c[i][k])+rt)

        end

        out += sys_sum*exp(-x.κs_λ[i]*λ0*t)
    end

    return out
end

function evalFIDspinsystem(t,
    αs::Vector{Vector{T}}, Ωs::Vector{Vector{T}},
    x::SpinSysParamsType2{T}, λ0::T,
    c, part_inds_molecule)::Complex{T} where T <: Real

    rt = zero(T) # pre-allocate.
    inner_sum = zero(Complex{T}) # pre-allocate.

    out = zero(Complex{T})
    for i = 1:length(αs)

        λ = x.κs_λ[i]*λ0
        inner_sum = zero(Complex{T})
        for k = 1:length(part_inds_molecule[i])
            inds = part_inds_molecule[i][k]

            rt = x.d[i][k]*t
            #out += exp(-λ*t)*evalFIDpartitionelement(t, αs[i][inds],
            inner_sum += evalFIDpartitionelement(t, αs[i][inds],
                Ωs[i][inds])*cis(dot(x.κs_β[i], c[i][k])+rt)
        end
        out += (exp(-λ*t)*inner_sum)
    end

    return out
end


#### FID singlets.
function evalFIDsinglets(t::T, d::Vector{T}, αs_singlets::Vector{T}, Ωs_singlets,
    βs_singlets, λ0::T, λ_multipliers::Vector{T}) where T <: Real

    out = zero(Complex{T})
    for i = 1:length(αs_singlets)

        λ = λ0*λ_multipliers[i]
        Ω = Ωs_singlets[i] + d[i]
        out += αs_singlets[i]*cis(Ω*t+βs_singlets[i])*exp(-λ*t)
    end
    return out
end
