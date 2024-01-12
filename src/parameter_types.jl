
abstract type MixtureSpinSys end

#abstract type MoleculeParamsMapping end

struct MixtureSinglets{T}
    
    # variables.
    ξs::Vector{Vector{T}} # multiplier wrt λ0.
    βs::Vector{Vector{T}}
    ζs::Vector{Vector{T}}

    # constants.
    αs::Vector{Vector{T}}
    Ωs::Vector{Vector{T}}

    # misc.
    #λ_singlets::Vector{Vector{T}} # actual T2.
    λ0::T
end

# each vector here have length equal to the number of compound entries.
struct CLMixtureSpinSys{T, ST,PT,T2T} <: MixtureSpinSys
    
    # The range of the functions here are all the real number line. This is in contrast to qs, which maps to the complex numbers.
    # Therefore, we need to do cis(phase angle) after using the functions here. see constructdesignmatrix!(), BatchEvalBuffer{T}. in b_models.jl.
    # spline surrogate portion of qs, and their derivatives. For computing qs' derivatives.
    srs::Vector{Vector{Vector{Function}}}
    sis::Vector{Vector{Vector{Function}}}

    ∇srs!::Vector{Vector{Vector{Function}}}
    ∇sis!::Vector{Vector{Vector{Function}}}

    # variables.
    shifts::Vector{ST}
    phases::Vector{PT}
    T2s::Vector{T2T}

    # misc.
    Δc_bars::Vector{Vector{Vector{Vector{T}}}}
    λ0::T
end

# We don't need the explicit surrogate functions for the FID model, since we don't implement the derivatives for the time-domain fit.
# For now, this data structure is used for updating the model parameters via importmodel!().
struct FIDMixtureSpinSys{T, ST,PT,T2T} <: MixtureSpinSys

    # variables.
    shifts::Vector{ST}
    phases::Vector{PT}
    T2s::Vector{T2T}

    # misc.
    Δc_bars::Vector{Vector{Vector{Vector{T}}}}
    λ0::T
end

"""
```
struct MoleculeParamsMapping
    st::Vector{Vector{Int}}
    fin::Vector{Vector{Int}}
end
```
"""
struct MoleculeParamsMapping
    st::Vector{Vector{Int}}
    fin::Vector{Vector{Int}}
end

function getvarrange(A::MoleculeParamsMapping)
    return A.st[begin][begin]:A.fin[end][end]
end

"""
```
struct ParamsMapping
    shift::MoleculeParamsMapping
    phase::MoleculeParamsMapping
    T2::MoleculeParamsMapping
end
```
"""
struct ParamsMapping
    shift::MoleculeParamsMapping
    phase::MoleculeParamsMapping
    T2::MoleculeParamsMapping
end

"""
```
function getshiftrange(A::ParamsMapping)
```

Returns the range for the shift parameters (ζ).
"""
function getshiftrange(A::ParamsMapping)
    return getvarrange(A.shift)
end

"""
```
function getphaserange(A::ParamsMapping)
```

Returns the range for the phase parameters (κ_β).
"""
function getphaserange(A::ParamsMapping)
    return getvarrange(A.phase)
end

"""
```
function getT2range(A::ParamsMapping)
```

Returns the range for the T2 parameters (κ_λ).
"""
function getT2range(A::ParamsMapping)
    return getvarrange(A.T2)
end

# alternative to MoleculeType, for mixtures.
struct MixtureModelParameters{T,ST <: MixtureSpinSys}

    # parameters used for computing shift (ζ), phase (β), T2 (λ).
    var_flat::Vector{T} # flat storage.
    
    MSS::ST # spin system parameters.
    systems_mapping::ParamsMapping # for addressing the flat storage.

    w::Vector{T}
end

"""
    function MixtureModelParameters(
        MSS::MT;
        w = ones(T, getNentries(MS)),
    )::MixtureModelParameters{T,MT} where {T <: AbstractFloat, MT <: MixtureSpinSys}

Convinence constructor for `MixtureModelParameters()`. Does not create a copy of of the inputs.
"""
function MixtureModelParameters(
    MSS::MT,
    w::Vector{T},
    )::MixtureModelParameters{T,MT} where {T <: AbstractFloat, MT <: MixtureSpinSys}

    @assert length(w) == getNentries(MSS)

    # set up.
    shifts, phases, T2s = MSS.shifts, MSS.phases, MSS.T2s

    N_κs_ζ = getNvars(shifts)
    N_κs_β = getNvars(phases)
    N_κs_λ = getNvars(T2s)
    #N_SS_vars = N_κs_ζ + N_κs_β + N_κs_λ

    #N_singlet_vars = getNvars(MS)
    N_singlet_vars = 0

    # spin system mapping.
    systems_mapping = ParamsMapping(shifts, phases, T2s)

    # singlet mapping.
    #singlets_mapping = getsingletsParamsMapping(MS, offset_ind = N_SS_vars)

    # flat storage for the variable parameters.
    var_flat = Vector{T}(undef, N_κs_ζ+N_κs_β+N_κs_λ+N_singlet_vars)

    # weights.
    @assert length(w) == getNentries(MSS)

    # allocate.
    out = MixtureModelParameters(var_flat, MSS, systems_mapping, w)

    # initialize var_flat with the parameters from MSS and MS.
    exportmodel!(out)

    return out
end

function ParamsMapping(MSS::MixtureSpinSys)
    return ParamsMapping(MSS.shifts, MSS.phases, MSS.T2s)
end

"""
```
function ParamsMapping(
    shifts::Vector{ST},
    phases::Vector{PT},
    T2s::Vector{T2T},
) where {ST,PT,T2T}
```

Convinence constructor for the type `ParamsMapping`.
"""
function ParamsMapping(
    shifts::Vector{ST},
    phases::Vector{PT},
    T2s::Vector{T2T},
    ) where {ST <: MoleculeParams, PT <: MoleculeParams, T2T <: MoleculeParams}

    N_κs_ζ = getNvars(shifts)
    N_κs_β = getNvars(phases)
    #N_κs_λ = getNvars(T2s)
    
    κs_ζ_st, κs_ζ_fin, fin_ind = buildmapping(shifts; offset_ind = 0)
    #@assert κs_ζ_fin[end] -κs_ζ_st[begin] +1 == N_κs_ζ

    κs_β_st, κs_β_fin, fin_ind = buildmapping(phases; offset_ind = N_κs_ζ)
    #@assert κs_β_fin[end] -κs_β_st[begin] +1 == N_κs_β

    κs_λ_st, κs_λ_fin, fin_ind = buildmapping(T2s; offset_ind = N_κs_ζ+N_κs_β)
    #@assert κs_λ_fin[end] -κs_λ_st[begin] +1 == N_κs_λ

    return ParamsMapping(
        MoleculeParamsMapping(κs_ζ_st, κs_ζ_fin),
        MoleculeParamsMapping(κs_β_st, κs_β_fin),
        MoleculeParamsMapping(κs_λ_st, κs_λ_fin),
    )
end

function buildvarmapping(
    ::Shared,
    κs_λ::Vector{T},
    prev_fin_ind::Int,
    ) where T

    st_ind = prev_fin_ind + 1
    fin_ind = st_ind + length(κs_λ) - 1
    # return [st_ind;], [fin_ind;], fin_ind

    sts = collect(st_ind:fin_ind)
    fins = sts

    return sts, fins, fin_ind
end


# identical to buildphasemapping() for CoherencePhase, but this dispatches on CohereneShift.
function buildmapping(
    Vs::Vector{VT};
    offset_ind::Int = 0,
    ) where VT <: MoleculeParams #CoherenceParams

    κs_ζ_st = Vector{Vector{Int}}(undef, length(Vs))
    κs_ζ_fin = Vector{Vector{Int}}(undef, length(Vs))

    #st_ind = 0 + offset_ind
    fin_ind = 0 + offset_ind
    for n in eachindex(Vs)

        if !isempty(Vs[n].var)
                
            κs_ζ_st[n], κs_ζ_fin[n], fin_ind = buildvarmapping(
                SystemsTrait(Vs[n]),
                Vs[n].var,
                fin_ind,
            )
        end
    end

    return κs_ζ_st, κs_ζ_fin, fin_ind
end

function buildvarmapping(
    ::Coherence,
    κs_β::Vector{Vector{T}},
    fin_ind::Int,
    ) where T

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
    Bs::Vector{MoleculeType{T, SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}, OT}},
    ) where {T, OT <: OperationRange}

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
    Bs::Vector{MoleculeType{T, SpinSysParams{SharedShift{T}, CoherencePhase{T}, SharedT2{T}}, OT}},
    ) where {T, OT <: OperationRange}

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
    As::Vector{HAM.SHType{T}},
    Bs::Vector{MoleculeType{T, SpinSysParams{ST, PT, T2T}, OT}},
    ) where {T, ST, PT, T2T, OT <: OperationRange}

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

function getNsinglets(αs::Vector{Vector{T}}) where T

    return sum( length(αs[n]) for n in eachindex(αs) )
end

# function getNsinglets(MS::MixtureSinglets{T}) where T
#     return getNsinglets(MS.αs)
# end

# function getNvars(MS::MixtureSinglets{T})::Int where T

#     return getNsinglets(MS)*3 # three types of variables: ζs, βs, ξs.
# end

function getNvars(model_params::MixtureModelParameters)::Int
    
    return getNvars(model_params.MSS)
end

function getNvars(MSS::MixtureSpinSys)::Int

    return getNvars(MSS.shifts) + getNvars(MSS.phases) + getNvars(MSS.T2s)
end

"""
```
function getNentries(MSS::MixtureSpinSys)::Int
```

Returns the number of molecule entries associated with `MSS`
"""
function getNentries(MSS::MixtureSpinSys)::Int

    return length(MSS.Δc_bars)
end

"""
```
function getNentries(model_params::MixtureModelParameters)::Int
```

Return the number of molecule entries associated with `model_params`
"""
function getNentries(model_params::MixtureModelParameters)::Int

    return getNentries(model_params.MSS)
end

# number of resonance groups. for all spin systems.
function getNgroups(
    As::Vector{HAM.SHType{T}},
    )::Int where T

    N = 0

    for n in eachindex(As) # molecule index.
        for i in eachindex(As[n].Δc_bar) # spin system index.

            for _ in eachindex(As[n].Δc_bar[i]) # effective chemical shift index.

                N += 1
            end
        end
    end

    return N
end
function getNsystems(
    As::Vector{HAM.SHType{T}},
    )::Int where T

    N = 0

    for n in eachindex(As) # molecule index.
        for i in eachindex(As[n].Δc_bar) # spin system index.

            N += 1
        end
    end

    return N
end