# Builds nested shifts, no copying.
function extractshifts(
    Bs::Vector{MoleculeType{T, SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}}},
    ) where T

    out_d = Vector{Vector{Vector{T}}}(undef, length(Bs))
    out_β = Vector{Vector{Vector{T}}}(undef, length(Bs))
    
    for n in eachindex(Bs)
        out_d[n] = Vector{Vector{T}}(undef, length(Bs[n].ss_params.shift.var))
        out_β[n] = Vector{Vector{T}}(undef, length(Bs[n].ss_params.phase.var))

        for i in eachindex(Bs[n].ss_params.shift.var)
            
            out_d[n][i] = Bs[n].ss_params.shift.var[i]
            out_β[n][i] = Bs[n].ss_params.phase.var[i]

            # out_d[n][i] = Vector{T}(undef, length(Bs[n].ss_params.shift.var))
            # for k in eachindex(Bs[n].ss_params.shift.var[i])
            #     out_d[n][i][k] = Bs[n].ss_params.shift.var[i]
            # end
        end
    end

    return out_d, out_β
end


function setupvars(
    Bs::Vector{MoleculeType{T, SpinSysParams{ST, PT, T2T}}},
    ) where {T, ST, PT, T2T}

    N = length(Bs)
    shifts = Vector{ST}(undef, N)
    phases = Vector{PT}(undef, N)
    T2s = Vector{T2T}(undef, N)

    for n in eachindex(Bs)
        shifts[n] = Bs[n].ss_params.shift
        phases[n] = Bs[n].ss_params.phase
        T2s[n] = Bs[n].ss_params.T2
    end

    return shifts, phases, T2s
end

# I am here.
function getNvars(
    Vs::Vector{VT},
    ) where VT <: MoleculeParams

    N_vars::Int = 0

    for n in eachindex(Vs) # molecule index.
        for i in eachindex(Vs[n].var) # spin system index.

            for j in eachindex(Vs[n].var[i]) # effective chemical shift index.

                N_vars += length(Vs[n].var[i][j])
            end
        end
    end

    return N_vars
end

#abstract type MoleculeParamsMapping end

struct MoleculeParamsMapping
    st::Vector{Vector{Int}}
    fin::Vector{Vector{Int}}
end

# struct CoherenceMapping <: MoleculeParamsMapping
#     st::Vector{Vector{Int}}
#     fin::Vector{Vector{Int}}
# end

# struct SharedMapping <: MoleculeParamsMapping
#     st::Vector{Int}
#     fin::Vector{Int}
# end

# function getmapping(
#     ::Type{VT},
#     st::Vector{Vector{Int}},
#     fin::Vector{Vector{Int}}
#     )::CoherenceMapping where VT <: CoherenceParams

#     return CoherenceMapping(st, fin)
# end

# function getmapping(
#     ::Type{VT},
#     st::Vector{Int},
#     fin::Vector{Int}
#     )::SharedMapping where VT <: SharedParams

#     return SharedMapping(st, fin)
# end

# for a mixture.
# struct ParamsMapping{ST<:MoleculeParamsMapping,PT<:MoleculeParamsMapping,T2T<:MoleculeParamsMapping}
#     #κs_d_st::Vector{Vector{Int}}
#     #κs_d_fin::Vector{Vector{Int}}
#     shift::ST
#     phase::PT
#     T2::T2T
# end
struct ParamsMapping
    #κs_d_st::Vector{Vector{Int}}
    #κs_d_fin::Vector{Vector{Int}}
    shift::MoleculeParamsMapping
    phase::MoleculeParamsMapping
    T2::MoleculeParamsMapping
end

function getParamsMapping(
    shifts::Vector{ST},
    phases::Vector{PT},
    T2s::Vector{T2T},
    ) where {ST,PT,T2T}

    N_κs_d = getNvars(shifts)
    N_κs_β = getNvars(phases)
    N_κs_λ = getNvars(T2s)
    
    κs_d_st, κs_d_fin, fin_ind = buildmapping(shifts; offset_ind = 0)
    #@assert κs_d_fin[end] -κs_d_st[begin] +1 == N_κs_d

    κs_β_st, κs_β_fin, fin_ind = buildmapping(phases; offset_ind = N_κs_d)
    #@assert κs_β_fin[end] -κs_β_st[begin] +1 == N_κs_β

    κs_λ_st, κs_λ_fin, fin_ind = buildmapping(T2s; offset_ind = N_κs_d+N_κs_β)
    #@assert κs_λ_fin[end] -κs_λ_st[begin] +1 == N_κs_λ

    return ParamsMapping(
        MoleculeParamsMapping(κs_d_st, κs_d_fin),
        MoleculeParamsMapping(κs_β_st, κs_β_fin),
        MoleculeParamsMapping(κs_λ_st, κs_λ_fin),
    )
end

# function buildmapping(
#     Vs::Vector{VT};
#     offset_ind::Int = 0,
#     ) where VT <: SharedParams

#     sts = Vector{Int}(undef, length(Vs))
#     fins = Vector{Int}(undef, length(Vs))

#     fill!(sts, -1)
#     fill!(fins, -1)

#     #st_ind = 0 + offset_ind
#     fin_ind = 0 + offset_ind
#     for n in eachindex(Vs)

#         sts[n], fin_ind = buildsharedmapping(Vs[n].var, fin_ind)
#         fins[n] = fin_ind
#     end

#     return sts, fins, fin_ind
# end

function buildvarmapping(κs_λ::Vector{T}, fin_ind::Int) where T

    st_ind = fin_ind + 1
    fin_ind = st_ind + length(κs_λ) - 1

    return [st_ind;], [fin_ind;], fin_ind
end


# identical to buildphasemapping() for CoherencePhase, but this dispatches on CohereneShift.
function buildmapping(
    Vs::Vector{VT};
    offset_ind::Int = 0,
    ) where VT <: MoleculeParams #CoherenceParams

    κs_d_st = Vector{Vector{Int}}(undef, length(Vs))
    κs_d_fin = Vector{Vector{Int}}(undef, length(Vs))

    #st_ind = 0 + offset_ind
    fin_ind = 0 + offset_ind
    for n in eachindex(Vs)

        if !isempty(Vs[n].var)
                
            κs_d_st[n], κs_d_fin[n], fin_ind = buildvarmapping(
                Vs[n].var,
                fin_ind,
            )
        end
    end

    return κs_d_st, κs_d_fin, fin_ind
end

function buildvarmapping(κs_β::Vector{Vector{T}}, fin_ind::Int) where T

    κs_β_st = Vector{Int}(undef, length(κs_β))
    κs_β_fin = Vector{Int}(undef, length(κs_β))

    for i in eachindex(κs_β)
        st_ind = fin_ind + 1
        fin_ind = st_ind + length(κs_β[i]) - 1

        κs_β_st[i] = st_ind
        κs_β_fin[i] = fin_ind
    end

    return κs_β_st, κs_β_fin, fin_ind
end

#=
function buildβrmapping(
    phases::Vector{CoherencePhase{T}};
    offset_ind::Int = 0,
    ) where T

    N = length(phases)

    ranges_κs_β_st = Vector{Vector{Int}}(undef, N)
    ranges_κs_β_fin = Vector{Vector{Int}}(undef, N)

    st_ind = 0 + offset_ind
    fin_ind = st_ind
    for n in eachindex(Bs)

        ranges_κs_β_st[n], ranges_κs_β_fin[n], fin_ind = buildcoherencemapping(
            phases[n].var,
            fin_ind,
        )
    end

    return ranges_κs_β_st, ranges_κs_β_fin
end



 function extractshiftranges(
    Bs::Vector{MoleculeType{T,CoherenceShift{T}}},
    offset_ind::Int,
    ) where T
    
    dtype = UnitRange{Int}

    ranges_d = Vector{Vector{dtype}}(undef, length(Bs))
    ranges_κs_d = Vector{Vector{dtype}}(undef, length(Bs))
    ranges_κs_β = Vector{Vector{dtype}}(undef, length(Bs))
    
    st_ind = 0 + offset_ind
    for n in eachindex(Bs)
        ranges_d[n] = Vector{dtype}(undef, length(Bs[n].ss_params.shift.var))

        for i in eachindex(Bs[n].ss_params.shift.var)
            st_ind +=1
            fin_ind = st_ind + length(Bs[n].ss_params.shift.var[i]) - 1

            ranges_d[n][i] = st_ind:fin_ind
        end
    end

    st_ind = 0
    for n in eachindex(Bs)
        ranges_κs_d[n] = Vector{dtype}(undef, length(Bs[n].ss_params.shift.var))
        ranges_κs_β[n] = Vector{dtype}(undef, length(Bs[n].ss_params.phase.var))

        for i in eachindex(Bs[n].ss_params.phase.var)
            st_ind +=1
            fin_ind = st_ind + length(Bs[n].ss_params.phase.var[i]) - 1

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
        views_d[n] = Vector{dtype}(undef, length(Bs[n].ss_params.shift.var))

        for i in eachindex(Bs[n].ss_params.shift.var)
            st_ind +=1
            fin_ind = st_ind + length(Bs[n].ss_params.shift.var[i]) - 1

            views_d[n][i] = @view buffer_d[begin+st_ind-1:begin+fin_ind-1]
        end
    end

    st_ind = 0
    for n in eachindex(Bs)
        views_κs_d[n] = Vector{dtype}(undef, length(Bs[n].ss_params.shift.var))
        views_κs_β[n] = Vector{dtype}(undef, length(Bs[n].ss_params.phase.var))

        for i in eachindex(Bs[n].ss_params.phase.var)
            st_ind +=1
            fin_ind = st_ind + length(Bs[n].ss_params.phase.var[i]) - 1

            views_κs_d[n][i] = @view buffer_κs_d[begin+st_ind-1:begin+fin_ind-1]
            views_κs_β[n][i] = @view buffer_κs_β[begin+st_ind-1:begin+fin_ind-1]
        end
    end

    return buffer_d, buffer_κs_d, buffer_κs_β, views_d, views_κs_d, views_κs_β
end =#