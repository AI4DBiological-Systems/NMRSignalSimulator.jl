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