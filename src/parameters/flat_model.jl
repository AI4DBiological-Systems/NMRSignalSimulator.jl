################# move to model_types later. for forward model.

struct FlatMixtureModel{T}

    ## non-singlet spin systems.
    ds::Vector{T}
    cos_βs::Vector{T}
    sin_βs::Vector{T}
    ξs::Vector{T} # multiplier.
    
    # constants.
    srs::Vector{Function}
    sis::Vector{Function}
    ∇srs!::Vector{Function}
    ∇sis!::Vector{Function}

    ## singlets.
    singlet_λs::Vector{T} # actual T2s, not multiplier.
    singlet_βs::Vector{T}
    singlet_ds::Vector{T}

    # constaints.
    singlet_αs::Vector{T}
    singlet_Ωs::Vector{T}

    λ0::T
end

function FlatMixtureModel(
    As::Vector{SHType{T}},
    MSS::MixtureSpinSys{T,ST,PT,T2T},
    MS::MixtureSinglets{T},
    ) where {T,ST,PT,T2T}
    
    ## non-singlet spin systems.
    N_groups = getNgroups(As)

    # variables.
    cos_βs = ones(T, N_groups) .* -Inf
    sin_βs = ones(T, N_groups) .* -Inf
    ξs = ones(T, N_groups) .* -Inf
    ds = ones(T, N_groups) .* -Inf

    # constants.
    # Vector{Vector{Vector{Function}}} to Vector{Function}.
    srs = collect( Base.Iterators.flatten( Base.Iterators.flatten(MSS.srs) ) )
    sis = collect( Base.Iterators.flatten( Base.Iterators.flatten(MSS.sis) ) )
    ∇srs! = collect( Base.Iterators.flatten( Base.Iterators.flatten(MSS.∇srs!) ) )
    ∇srs! = collect( Base.Iterators.flatten( Base.Iterators.flatten(MSS.∇srs!) ) )

    ## singlets.
    N_singlets = getNsinglets(MS)
    
    # variables. update later.
    singlet_λs = ones(T, N_singlets) .* -Inf
    singlet_βs = ones(T, N_singlets) .* -Inf
    singlet_ds = ones(T, N_singlets) .* -Inf

    # constants. update now.
    singlet_αs = ones(T, N_singlets) .* -Inf
    singlet_Ωs = ones(T, N_singlets) .* -Inf
    resolvesingletparameters!(singlet_αs, MS.αs)
    resolvesingletparameters!(singlet_Ωs, MS.Ωs)

    FM = FlatMixtureModel(
        ds,
        ξs,
        cos_βs,
        sin_βs,
        srs,
        sis,
        ∇srs!,
        ∇srs!,
        singlet_λs,
        singlet_βs,
        singlet_ds,
        singlet_αs,
        singlet_Ωs,
        MS.λ0
    )

    # update variables.
    resolveflatten!(FM, MSS)
    resolveflatten!(FM, MS)

    return FM
end

function getNsinglets(
    MS::MixtureSinglets{T},
    ) where T

    return sum( length(MS.αs[n]) for n in eachindex(MS.αs) )
end

# mutates FM.
function resolveflatten!(
    FM::FlatMixtureModel{T},
    MSS::MixtureSpinSys{T,ST,PT,T2T},
    ) where {T,ST,PT,T2T}

    resolveparameters!(FM.cos_βs, FM.sin_βs, MSS.phases, MSS.Δc_bars)
    resolveparameters!(FM.ds, MSS.shifts, MSS.Δc_bars)
    resolveparameters!(FM.ξs, MSS.T2s, MSS.Δc_bars)

    return nothing
end

function resolveflatten!(
    FM::FlatMixtureModel{T},
    MS::MixtureSinglets{T},
    ) where T

    resolvesingletparameters!(FM.singlet_βs, MS.βs)
    resolvesingletparameters!(FM.singlet_ds, MS.ds)
    resolvesingletparameters!(FM.singlet_λs, MS.ξs, FM.λ0)

    return nothing
end


function evalforwardmodel(
    u_rad,
    srs,
    sis,
    ds::Vector{T},
    ξs,
    cos_βs,
    sin_βs,
    )::Complex{T} where T

    return sum( 
        evalq(
            
            srs[l],
            sis[l],
            u_rad - ds[l],
            ξs[l],
            cos_βs[l],
            sin_βs[l],
        )
    for l in eachindex(srs))
end

function evalforwardmodel(
    u_rad,
    FM::FlatMixtureModel{T},
    MS,
    )::Complex{T} where T
    
    # eval.
    # out_sys = zero(Complex{T})
    # for l in eachindex(FM.srs)

    #     out_sys += evalq(
            
    #         FM.srs[l],
    #         FM.sis[l],
    #         u_rad - FM.ds[l],
    #         FM.ξs[l],
    #         FM.cos_βs[l],
    #         FM.sin_βs[l],
    #     )
    # end
    out_sys = evalforwardmodel(
        u_rad,
        FM.srs,
        FM.sis,
        FM.ds,
        FM.ξs,
        FM.cos_βs,
        FM.sin_βs,
    )
    out_singlets = evalclsinglets(
        u_rad,
        FM.singlet_ds,
        FM.singlet_αs,
        FM.singlet_Ωs,
        FM.singlet_βs,
        FM.singlet_λs,
    )

    return out_sys + out_singlets
end

function evalclsinglets(
    u_rad::T, 
    ds::Vector{T}, 
    αs::Vector{T}, 
    Ωs,
    βs,
    λs::Vector{T}
    )::Complex{T} where T <: Real

    out = zero(Complex{T})
    for i in eachindex(αs)
        τ = u_rad - ds[i]

        out += αs[i]*cis(βs[i])/(λs[i]+im*(τ-Ωs[i]))
    end
    return out
end



##### flat version. For every resonance group.

function resolveparameters!(
    cos_βs::Vector{T},
    sin_βs::Vector{T},
    phases::Vector{CoherencePhase{T}},
    Δc_bars::Vector{Vector{Vector{Vector{T}}}},
    ) where T

    #@assert length(cos_β) == length(sin_β) == getNgroups(As)

    offset = 0
    for n in eachindex(phases)
        x = phases[n]

        #@assert length(x_var) == length(x.cos_β) == length(x.sin_β) # same number of non-singlet spin systems.

        for i in eachindex(x.cos_β)
            for k in eachindex(x.cos_β[i])

                β = dot(Δc_bars[n][i][k], x.var[i])

                cos_βs[begin+offset] = cos(β)
                sin_βs[begin+offset] = sin(β)
                offset += 1
            end
        end
    end

    return offset
end

function resolveparameters!(
    ds::Vector{T},
    shifts::Vector{CoherenceShift{T}},
    Δc_bars::Vector{Vector{Vector{Vector{T}}}},
    ) where T

    offset = 0
    for n in eachindex(shifts)
        x = shifts[n]

        #@assert length(x_var) == length(x.d) # same number of non-singlet spin systems.

        for i in eachindex(x.d)    

            for k in eachindex(x.d[i])

                ds[begin+offset] = dot(Δc_bars[n][i][k], x.var[i])
                offset += 1
            end
        end
    end

    return offset
end

function resolveparameters!(
    ds::Vector{T},
    shifts::Vector{SharedShift{T}},
    Δc_bars::Vector{Vector{Vector{Vector{T}}}},
    ) where T

    offset = 0
    for n in eachindex(shifts)
        x = shifts[n]

        #@assert length(x_var) == length(x.d) # same number of non-singlet spin systems.

        for i in eachindex(x.var)    

            for k in eachindex(Δc_bars[n][i])

                ds[begin+offset] = x.var[i]
                offset += 1
            end
        end
    end

    return offset
end

function resolveparameters!(
    ξs::Vector{T},
    T2s::Vector{SharedT2{T}},
    Δc_bars::Vector{Vector{Vector{Vector{T}}}},
    ) where T

    offset = 0
    for n in eachindex(T2s)
        x = T2s[n]

        #@assert length(x_var) == length(x.d) # same number of non-singlet spin systems.

        for i in eachindex(x.var)    

            for k in eachindex(Δc_bars[n][i])

                #λs[begin+offset] = x.var[i]*λ0
                ξs[begin+offset] = x.var[i]
                offset += 1
            end
        end
    end

    return offset
end

function resolvesingletparameters!(
    λs::Vector{T},
    nested_ξs::Vector{Vector{T}},
    λ0::T,
    ) where T

    offset = 0
    for n in eachindex(nested_ξs)
        for i in eachindex(nested_ξs[n])    

            λs[begin+offset] = nested_ξs[n][i]*λ0
            offset += 1
        end
    end

    return offset
end

function resolvesingletparameters!(
    ps::Vector{T},
    nested_ps::Vector{Vector{T}},
    ) where T

    offset = 0
    for n in eachindex(nested_ps)
        for i in eachindex(nested_ps[n])    

            ps[begin+offset] = nested_ps[n][i]
            offset += 1
        end
    end

    return offset
end