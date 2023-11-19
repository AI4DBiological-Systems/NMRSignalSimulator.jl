###
# - updatesparameters!() transfer contents from a flat storage to a parameter data structure.
# exportparameters!() transfer the reverse direction in comparison to updateparameters!().
# - before evaluating the model with updated parameter data structures, do resolvespinsystems!() to compute the intermediate quantities used by the model evaluation code.
# No intermediates are used for evaluating singlets, so a resolve procedure isn't needed.


###### update and export routines for parameters.

function updateparameters!(
    θs::Vector{ST}, # mutates.
    p,
    θs_mapping::MoleculeParamsMapping,
    ) where ST<:SharedParams

    for n in eachindex(θs)

        for i in eachindex(θs[n].var)

            st_ind = θs_mapping.st[n][i]
            
            θs[n].var[i] = p[st_ind]
            # fin_ind = θs_mapping.fin[n][i]
            # @assert fin_ind - st_ind +1 == length(θs[n].var[i])

            # for j in eachindex(θs[n].var[i])
            #     θs[n].var[i][j] = p[st_ind+j-1]
            # end
        end
    end

    return nothing
end

# reverse of updateparameters!()
function exportparameters!(
    p, # mutates.
    θs::Vector{ST},
    θs_mapping::MoleculeParamsMapping,
    ) where ST<:SharedParams

    for n in eachindex(θs)

        for i in eachindex(θs[n].var)

            st_ind = θs_mapping.st[n][i]
            
            p[st_ind] = θs[n].var[i]
        end
    end

    return nothing
end

function updateparameters!(
    θs::Vector{CT}, # mutates.
    p,
    θs_mapping::MoleculeParamsMapping,
    ) where CT<:CoherenceParams

    for n in eachindex(θs)
        x = θs[n]

        for i in eachindex(x.var)

            st_ind = θs_mapping.st[n][i]
            fin_ind = θs_mapping.fin[n][i]

            x.var[i][:] = p[begin+st_ind-1:begin+fin_ind-1]
        end
    end

    return nothing
end

# reverse of updateparameters!()
function exportparameters!(
    p, # mutates.
    θs::Vector{CT},
    θs_mapping::MoleculeParamsMapping,
    ) where CT<:CoherenceParams

    for n in eachindex(θs)
        x = θs[n]

        for i in eachindex(x.var)

            st_ind = θs_mapping.st[n][i]
            fin_ind = θs_mapping.fin[n][i]

            p[begin+st_ind-1:begin+fin_ind-1] = x.var[i]#[:]
        end
    end

    return nothing
end



# update MSS' var fields with contents in p.
function updatespinsystems!(
    MSS::MixtureSpinSys, # mutates.
    p,
    mapping::ParamsMapping,
    )

    updateparameters!(MSS.T2s, p, mapping.T2)
    updateparameters!(MSS.shifts, p, mapping.shift)
    updateparameters!(MSS.phases, p, mapping.phase)

    return nothing
end

# update p with contents of MSS' var fields.
function exportspinsystems!(
    p, # mutates.
    MSS::MixtureSpinSys,
    mapping::ParamsMapping,
    )

    exportparameters!(p, MSS.T2s, mapping.T2)
    exportparameters!(p, MSS.shifts, mapping.shift)
    exportparameters!(p, MSS.phases, mapping.phase)

    return nothing
end

##### resolve non=singlet spin systems.

# update MSS variables using the contents in MSS' var fields.
function resolvespinsystems!(MSS::MixtureSpinSys)

    resolveparameters!(MSS.T2s, MSS.λ0)
    resolveparameters!(MSS.shifts, MSS.Δc_bars)
    resolveparameters!(MSS.phases, MSS.Δc_bars)

    return nothing
end

function resolveparameters!(
    phases::Vector{CoherencePhase{T}},
    Δc_bars::Vector{Vector{Vector{Vector{T}}}},
    ) where T <: AbstractFloat

    for n in eachindex(phases)
        x = phases[n]

        @assert length(x.cos_β) == length(x.sin_β) # same number of non-singlet spin systems.

        for i in eachindex(x.cos_β)
            for k in eachindex(x.cos_β[i])

                β = dot(Δc_bars[n][i][k], x.var[i])

                x.β[i][k] = β
                x.cos_β[i][k] = cos(β)
                x.sin_β[i][k] = sin(β)
            end
        end
        #resolveparameters!(phases[n], Δc_bars[n])
    end

    return nothing
end

function resolveparameters!(
    shifts::Vector{CoherenceShift{T}},
    Δc_bars::Vector{Vector{Vector{Vector{T}}}},
    ) where T <: AbstractFloat

    for n in eachindex(shifts)
        x = shifts[n]

        for i in eachindex(x.ζ)    

            for k in eachindex(x.ζ[i])

                x.ζ[i][k] = dot(Δc_bars[n][i][k], x.var[i])
            end
        end
    end

    return nothing
end

function resolveparameters!(::Vector{SharedShift{T}}, args...) where T <: AbstractFloat
    return nothing
end

function resolveparameters!(T2s::Vector{SharedT2{T}}, λ0::T, args...) where T <: AbstractFloat
    
    for n in eachindex(T2s)
        @assert length(T2s[n].var) == length(T2s[n].λ)
        
        for i in eachindex(T2s[n].var)
            T2s[n].λ[i] = T2s[n].var[i]*λ0
        end
    end

    return nothing
end
