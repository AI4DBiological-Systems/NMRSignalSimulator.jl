# Builds nested shifts, no copying.
function extractshifts(Bs::Vector{MoleculeType{T,CoherenceShift{T}}}) where T

    out_d = Vector{Vector{Vector{T}}}(undef, length(Bs))
    out_β = Vector{Vector{Vector{T}}}(undef, length(Bs))
    
    for n in eachindex(Bs)
        out_d[n] = Vector{Vector{T}}(undef, length(Bs[n].ss_params.d))
        out_β[n] = Vector{Vector{T}}(undef, length(Bs[n].ss_params.κs_β))

        for i in eachindex(Bs[n].ss_params.d)
            
            out_d[n][i] = Bs[n].ss_params.d[i]
            out_β[n][i] = Bs[n].ss_params.κs_β[i]

            # out_d[n][i] = Vector{T}(undef, length(Bs[n].ss_params.d))
            # for k in eachindex(Bs[n].ss_params.d[i])
            #     out_d[n][i][k] = Bs[n].ss_params.d[i]
            # end
        end
    end

    return out_d, out_β
end

function getnumvars(Bs::Vector{MoleculeType{T,CoherenceShift{T}}}) where T

    N_d::Int = 0
    N_κs_d::Int = 0
    N_κs_β::Int = 0

    for n in eachindex(Bs) # molecule index.
        for i in eachindex(Bs[n].ss_params.d) # spin system index.
            for k in eachindex(Bs[n].ss_params.d[i]) # resonance group index.

                N_d += length(Bs[n].ss_params.d[i][k])
            end

            for j in eachindex(Bs[n].ss_params.κs_β[i]) # effective chemical shift index.

                N_κs_d += length(Bs[n].ss_params.κs_d[i][j])
                N_κs_β += length(Bs[n].ss_params.κs_β[i][j])
            end
        end
    end

    return N_d, N_κs_d, N_κs_β
end

# for a mixture.
struct CoherenceShiftMapping
    #κs_d_st::Vector{Vector{Int}}
    #κs_d_fin::Vector{Vector{Int}}
    κs_d_st::Vector{Int}
    κs_d_fin::Vector{Int}
    
    #κs_β::Vector{Vector{T}}
    κs_β_st::Vector{Int}
    κs_β_fin::Vector{Int}
end

function buildshiftmapping(
    Bs::Vector{MoleculeType{T,SharedShift{T}}};
    offset_ind::Int = 0,
    ) where T

    ranges_κs_d_st = Vector{Int}(undef, length(Bs))
    ranges_κs_d_fin = Vector{Int}(undef, length(Bs))

    st_ind = 0 + offset_ind
    for n in eachindex(Bs)
        
        st_ind +=1
        fin_ind = st_ind + length(Bs[n].ss_params.κs_d) - 1

        ranges_κs_d_st[n] = st_ind
        ranges_κs_d_fin[n] = fin_ind
    end

    return ranges_κs_d_st, ranges_κs_d_fin
end

function buildshiftmapping(
    Bs::Vector{MoleculeType{T,CoherenceShift{T}}};
    offset_ind::Int = 0,
    ) where T

    return buildβrmapping(Bs; offset_ind = offset_ind)
end

function buildβrmapping(
    Bs::Vector{MoleculeType{T,CoherenceShift{T}}};
    offset_ind::Int = 0,
    ) where T

    ranges_κs_β_st = Vector{Vector{Int}}(undef, length(Bs))
    ranges_κs_β_fin = Vector{Vector{Int}}(undef, length(Bs))

    st_ind = 0 + offset_ind
    for n in eachindex(Bs)
        ranges_κs_β_st[n] = Vector{Int}(undef, length(Bs[n].ss_params.κs_β))
        ranges_κs_β_fin[n] = Vector{Int}(undef, length(Bs[n].ss_params.κs_β))

        for i in eachindex(Bs[n].ss_params.κs_β)
            st_ind +=1
            fin_ind = st_ind + length(Bs[n].ss_params.κs_β[i]) - 1

            ranges_κs_β_st[n][i] = st_ind
            ranges_κs_β_fin[n][i] = fin_ind
        end
    end

    return ranges_κs_β_st, ranges_κs_β_fin
end


function updatemodel!(
    Bs,
    #Bs::Vector{MoleculeType{T,CoherenceShift{T}}},
    mapping::CoherenceShiftMapping,
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
        
        for i in eachindex(Bs[n].ss_params.κs_β)
            l += 1
            #st_ind +=1
            #fin_ind = st_ind + length(Bs[n].ss_params.κs_β[i]) - 1

            #Bs[n].ss_params.κs_d[i][:] = x[begin+ranges_κs_d[n][i].st-1:begin+ranges_κs_d[n][i].fin-1]

            Bs[n].ss_params.κs_d[i][:] = x[begin+κs_d_st[l]-1:begin+κs_d_fin[l]-1]

            # update shift.
            for k in eachindex(Bs[n].ss_params.d[i])
                Bs[n].ss_params.d[i][k] = dot(As[n].Δc_bar[i][k], Bs[n].ss_params.κs_d[i])
            end

        end

    end

    # phase κs_β
    l = 0
    for n in eachindex(Bs)
        
        for i in eachindex(Bs[n].ss_params.κs_β)
            l += 1
            #st_ind +=1
            #fin_ind = st_ind + length(Bs[n].ss_params.κs_β[i]) - 1

            #Bs[n].ss_params.κs_d[i][:] = x[begin+ranges_κs_d[n][i].st-1:begin+ranges_κs_d[n][i].fin-1]

            Bs[n].ss_params.κs_β[i][:] = x[begin+κs_β_st[l]-1:begin+κs_β_fin[l]-1]

            #ranges_κs_β[n][i] = st_ind:fin_ind
            # update shift.
            for k in eachindex(Bs[n].ss_params.β[i])
                Bs[n].ss_params.β[i][k] = dot(As[n].Δc_bar[i][k], Bs[n].ss_params.κs_β[i])
            end
        end
    end

    return nothing
end


#= function extractshiftranges(
    Bs::Vector{MoleculeType{T,CoherenceShift{T}}},
    offset_ind::Int,
    ) where T
    
    dtype = UnitRange{Int}

    ranges_d = Vector{Vector{dtype}}(undef, length(Bs))
    ranges_κs_d = Vector{Vector{dtype}}(undef, length(Bs))
    ranges_κs_β = Vector{Vector{dtype}}(undef, length(Bs))
    
    st_ind = 0 + offset_ind
    for n in eachindex(Bs)
        ranges_d[n] = Vector{dtype}(undef, length(Bs[n].ss_params.d))

        for i in eachindex(Bs[n].ss_params.d)
            st_ind +=1
            fin_ind = st_ind + length(Bs[n].ss_params.d[i]) - 1

            ranges_d[n][i] = st_ind:fin_ind
        end
    end

    st_ind = 0
    for n in eachindex(Bs)
        ranges_κs_d[n] = Vector{dtype}(undef, length(Bs[n].ss_params.κs_d))
        ranges_κs_β[n] = Vector{dtype}(undef, length(Bs[n].ss_params.κs_β))

        for i in eachindex(Bs[n].ss_params.κs_β)
            st_ind +=1
            fin_ind = st_ind + length(Bs[n].ss_params.κs_β[i]) - 1

            ranges_κs_d[n][i] = st_ind:fin_ind
            ranges_κs_β[n][i] = st_ind:fin_ind
        end
    end

    return buffer_d, buffer_κs_d, buffer_κs_β, ranges_d, ranges_κs_d, ranges_κs_β
end

function extractshiftviews(Bs::Vector{MoleculeType{T,CoherenceShift{T}}}) where T
    
    dtype = SubArray{T,1, Vector{T}, Tuple{UnitRange{Int}}, true}

    N_d, N_κs_d, N_κs_β = getnumvars(Bs)

    buffer_d = Vector{T}(undef, N_d)
    buffer_κs_d = Vector{T}(undef, N_κs_d)
    buffer_κs_β = Vector{T}(undef, N_κs_β)

    views_d = Vector{Vector{dtype}}(undef, length(Bs))
    views_κs_d = Vector{Vector{dtype}}(undef, length(Bs))
    views_κs_β = Vector{Vector{dtype}}(undef, length(Bs))
    
    st_ind = 0
    for n in eachindex(Bs)
        views_d[n] = Vector{dtype}(undef, length(Bs[n].ss_params.d))

        for i in eachindex(Bs[n].ss_params.d)
            st_ind +=1
            fin_ind = st_ind + length(Bs[n].ss_params.d[i]) - 1

            views_d[n][i] = @view buffer_d[begin+st_ind-1:begin+fin_ind-1]
        end
    end

    st_ind = 0
    for n in eachindex(Bs)
        views_κs_d[n] = Vector{dtype}(undef, length(Bs[n].ss_params.κs_d))
        views_κs_β[n] = Vector{dtype}(undef, length(Bs[n].ss_params.κs_β))

        for i in eachindex(Bs[n].ss_params.κs_β)
            st_ind +=1
            fin_ind = st_ind + length(Bs[n].ss_params.κs_β[i]) - 1

            views_κs_d[n][i] = @view buffer_κs_d[begin+st_ind-1:begin+fin_ind-1]
            views_κs_β[n][i] = @view buffer_κs_β[begin+st_ind-1:begin+fin_ind-1]
        end
    end

    return buffer_d, buffer_κs_d, buffer_κs_β, views_d, views_κs_d, views_κs_β
end =#