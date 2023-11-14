###
# - updatesparameters!() transfer contents from a flat storage to a parameter data structure.
# exportparameters!() transfer the reverse direction in comparison to updateparameters!().
# - before evaluating the model with updated parameter data structures, do resolvespinsystems!() to compute the intermediate quantities used by the model evaluation code.
# No intermediates are used for evaluating singlets, so a resolve procedure isn't needed.


###### update and export routines for parameters.

function updateparameters!(
    θs::Vector{ST},
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
# don't loop over θs_mapping.st or θs_mapping.fin because they might have undefined references.
# instead, loop over θs, which has empty collections instead of undefined references.
# TODO: change θs_mapping.st and θs_mapping.fin such that they don't have this issue later.
function exportparameters!(
    p,
    θs::Vector{ST},
    θs_mapping::MoleculeParamsMapping,
    ) where ST<:SharedParams

    for n in eachindex(θs)

        for i in eachindex(θs[n].var)

            st_ind = θs_mapping.st[n][i]
            
            p[st_ind] = θs[n].var[i]
            # fin_ind = θs_mapping.fin[n][i]
            # @show (fin_ind - st_ind +1)
            # @show length(θs[n].var[i])
            # @assert (fin_ind - st_ind +1) == #length(θs[n].var[i])

            # for j in eachindex(θs[n].var[i])
            #     p[st_ind+j-1] = θs[n].var[i][j]
            # end
        end
    end

    return nothing
end

function updateparameters!(
    θs::Vector{CT},
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
    p,
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
    MSS::CLMixtureSpinSys{T,SST},
    p,
    mapping::ParamsMapping,
    ) where {T,SST}

    updateparameters!(MSS.T2s, p, mapping.T2)
    updateparameters!(MSS.shifts, p, mapping.shift)
    updateparameters!(MSS.phases, p, mapping.phase)

    return nothing
end

# update p with contents of MSS' var fields.
function exportspinsystems!(
    p,
    MSS::CLMixtureSpinSys{T,SST},
    mapping::ParamsMapping,
    ) where {T,SST}

    exportparameters!(p, MSS.T2s, mapping.T2)
    exportparameters!(p, MSS.shifts, mapping.shift)
    exportparameters!(p, MSS.phases, mapping.phase)

    return nothing
end

##### resolve non=singlet spin systems.

# update MSS variables using the contents in MSS' var fields.
function resolvespinsystems!(
    MSS::CLMixtureSpinSys{T,SST},
    ) where {T,SST}

    resolveparameters!(MSS.T2s, MSS.λ0)
    resolveparameters!(MSS.shifts, MSS.Δc_bars)
    resolveparameters!(MSS.phases, MSS.Δc_bars)

    return nothing
end

function resolveparameters!(
    phases::Vector{CoherencePhase{T}},
    Δc_bars::Vector{Vector{Vector{Vector{T}}}},
    ) where T

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
    ) where T

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

function resolveparameters!(shifts::Vector{SharedShift{T}}, args...) where T
    return nothing
end

function resolveparameters!(T2s::Vector{SharedT2{T}}, λ0::T, args...) where T
    
    for n in eachindex(T2s)
        @assert length(T2s[n].var) == length(T2s[n].λ)
        
        for i in eachindex(T2s[n].var)
            T2s[n].λ[i] = T2s[n].var[i]*λ0
        end
    end

    return nothing
end
