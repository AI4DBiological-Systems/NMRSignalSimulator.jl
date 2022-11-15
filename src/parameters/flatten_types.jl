

#abstract type MoleculeParamsMapping end

struct MoleculeParamsMapping
    st::Vector{Vector{Int}}
    fin::Vector{Vector{Int}}
end

struct ParamsMapping
    #κs_d_st::Vector{Vector{Int}}
    #κs_d_fin::Vector{Vector{Int}}
    shift::MoleculeParamsMapping
    phase::MoleculeParamsMapping
    T2::MoleculeParamsMapping
end

function getsingletsParamsMapping(
    MS::MixtureSinglets{T};
    offset_ind::Int = 0,
    ) where T

    ξs, βs, ds = MS.ξs, MS.βs, MS.ds

    N = length(ξs)
    @assert length(βs) == length(ds)

    #st_ind = 0 + offset_ind
    fin_ind = 0 + offset_ind
    
    # shift.
    ds_st = Vector{Vector{Int}}(undef, N)
    ds_fin = Vector{Vector{Int}}(undef, N)
    for n in eachindex(ds)

        if !isempty(ds[n])
            ds_st[n], ds_fin[n], fin_ind = buildvarmapping(ds[n], fin_ind)
        end
    end

    # phase.
    βs_st = Vector{Vector{Int}}(undef, N)
    βs_fin = Vector{Vector{Int}}(undef, N)
    for n in eachindex(βs)

        if !isempty(βs[n])
            βs_st[n], βs_fin[n], fin_ind = buildvarmapping(βs[n], fin_ind)
        end
    end

    # T2.
    ξs_st = Vector{Vector{Int}}(undef, N)
    ξs_fin = Vector{Vector{Int}}(undef, N)
    for n in eachindex(ξs)

        if !isempty(ξs[n])
            ξs_st[n], ξs_fin[n], fin_ind = buildvarmapping(ξs[n], fin_ind)
        end
    end

    return ParamsMapping(
        MoleculeParamsMapping(ds_st, ds_fin),
        MoleculeParamsMapping(βs_st, βs_fin),
        MoleculeParamsMapping(ξs_st, ξs_fin),
    )
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


############# utility functions.
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

function extractshifts(
    Bs::Vector{MoleculeType{T, SpinSysParams{SharedShift{T}, CoherencePhase{T}, SharedT2{T}}}},
    ) where T

    out_d = Vector{Vector{T}}(undef, length(Bs))
    out_β = Vector{Vector{Vector{T}}}(undef, length(Bs))
    
    for n in eachindex(Bs)
        out_d[n] = Bs[n].ss_params.shift.var
        out_β[n] = Vector{Vector{T}}(undef, length(Bs[n].ss_params.phase.var))

        for i in eachindex(Bs[n].ss_params.phase.var)
            
            out_β[n][i] = Bs[n].ss_params.phase.var[i]

            # out_d[n][i] = Vector{T}(undef, length(Bs[n].ss_params.shift.var))
            # for k in eachindex(Bs[n].ss_params.shift.var[i])
            #     out_d[n][i][k] = Bs[n].ss_params.shift.var[i]
            # end
        end
    end

    return out_d, out_β
end

# creates referenecs to the nested array objects in Bs, As.
function setupSSvars(
    As::Vector{SHType{T}},
    Bs::Vector{MoleculeType{T, SpinSysParams{ST, PT, T2T}}},
    ) where {T, ST, PT, T2T}

    N = length(Bs)
    @assert length(As) == N

    shifts = Vector{ST}(undef, N)
    phases = Vector{PT}(undef, N)
    T2s = Vector{T2T}(undef, N)
    Δc_bars = Vector{Vector{Vector{Vector{T}}}}(undef, N)

    for n in eachindex(Bs)
        shifts[n] = Bs[n].ss_params.shift
        phases[n] = Bs[n].ss_params.phase
        T2s[n] = Bs[n].ss_params.T2

        Δc_bars[n] = As[n].Δc_bar
    end

    return shifts, phases, T2s, Δc_bars
end

# for non-singlet spin systems.
function getNvars(
    Vs::Vector{VT},
    )::Int where VT <: MoleculeParams

    N_vars = 0

    for n in eachindex(Vs) # molecule index.
        for i in eachindex(Vs[n].var) # spin system index.

            for j in eachindex(Vs[n].var[i]) # effective chemical shift index.

                N_vars += length(Vs[n].var[i][j])
            end
        end
    end

    return N_vars
end

# number of resonance groups. For non-singlet spin systems.
function getNgroups(
    As::Vector{SHType{T}},
    )::Int where T

    N = 0

    for n in eachindex(As) # molecule index.
        for i in eachindex(As[n].Δc_bar) # spin system index.

            for k in eachindex(As[n].Δc_bar[i]) # effective chemical shift index.

                N += 1
            end
        end
    end

    return N
end