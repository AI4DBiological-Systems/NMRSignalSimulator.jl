
"""
```
updateparameters!(A::CoherencePhase{T})
```

Updates the other fields of A using the contents of `A.var`.
"""
function updateparameters!(A::CoherencePhase{T}, Δc_bar::Vector{Vector{Vector{T}}}) where T <: AbstractFloat
    
    @assert length(Δc_bar) == length(A.var)

    β = NaN

    for i in eachindex(A.var)
        @assert length(Δc_bar[i]) == length(A.cos_β[i]) == length(A.sin_β[i])

        for k in eachindex(A.var[i])
            β = dot(A.var[i], Δc_bar[i][k])
            A.cos_β[i][k] = cos(β)
            A.sin_β[i][k] = sin(β)
        end
    end

    return nothing
end

"""
```
updateparameters!(A::CoherencePhase{T})
```

Updates the other fields of A using the contents of `A.var`.
"""
function updateparameters!(A::CoherenceShift{T}, Δc_bar::Vector{Vector{Vector{T}}}) where T <: AbstractFloat
    
    @assert length(Δc_bar) == length(A.var)

    β = NaN

    for i in eachindex(A.var)
        @assert length(Δc_bar[i]) == length(A.d[i])

        for k in eachindex(A.var[i])
            A.d[i][k] = dot(A.var[i], Δc_bar[i][k])
        end
    end

    return nothing
end



function updatemodel!(
    Bs,
    #Bs::Vector{MoleculeType{T,CoherenceShift{T}}},
    #mapping::CoherenceShiftMapping,
    x::Vector{T},
    ) where T <: Real

    #N_d, N_κs_d, N_κs_β = getnumvars(Bs)

    # These are 1-indexing arrays.
    κs_d_st = mapping.κs_d_st
    κs_d_fin = mapping.κs_d_fin
    κs_β_st = mapping.κs_β_st
    κs_β_fin = mapping.κs_β_fin

    # shift κs_d
    l = 0
    for n in eachindex(Bs)
        
        for i in eachindex(Bs[n].ss_params.phase.var)
            l += 1
            #st_ind +=1
            #fin_ind = st_ind + length(Bs[n].ss_params.phase.var[i]) - 1

            #Bs[n].ss_params.shift.var[i][:] = x[begin+ranges_κs_d[n][i].st-1:begin+ranges_κs_d[n][i].fin-1]

            Bs[n].ss_params.shift.var[i][:] = x[begin+κs_d_st[l]-1:begin+κs_d_fin[l]-1]

            # update shift.
            for k in eachindex(Bs[n].ss_params.shift.var[i])
                Bs[n].ss_params.shift.var[i][k] = dot(As[n].Δc_bar[i][k], Bs[n].ss_params.shift.var[i])
            end

        end

    end

    # phase κs_β
    l = 0
    for n in eachindex(Bs)
        
        for i in eachindex(Bs[n].ss_params.phase.var)
            l += 1
            #st_ind +=1
            #fin_ind = st_ind + length(Bs[n].ss_params.phase.var[i]) - 1

            #Bs[n].ss_params.shift.var[i][:] = x[begin+ranges_κs_d[n][i].st-1:begin+ranges_κs_d[n][i].fin-1]

            Bs[n].ss_params.phase.var[i][:] = x[begin+κs_β_st[l]-1:begin+κs_β_fin[l]-1]

            #ranges_κs_β[n][i] = st_ind:fin_ind
            # update shift.
            for k in eachindex(Bs[n].ss_params.β[i])
                Bs[n].ss_params.β[i][k] = dot(As[n].Δc_bar[i][k], Bs[n].ss_params.phase.var[i])
            end
        end
    end

    return nothing
end