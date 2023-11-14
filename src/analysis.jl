
############# model conversion: SharedShift to CoherenceShift. Implement the reverse direction later.

function modelconversion(
    Xs::Vector{ST},
    Δc_bars::Vector{Vector{Vector{Vector{T}}}}, # MSS.Δc_bars
    )::Vector{CoherenceShift{T}} where {T <: AbstractFloat, ST <: SharedParams}
    @assert length(Xs) == length(Δc_bars)

    Ys = collect( CoherenceShift(Δc_bars[n]) for n in eachindex(Δc_bars) )
    modelconversion!(Ys, Xs)

    return Ys
end

function modelconversion!(
    Ys::Vector{CT},
    Xs::Vector{ST},
    ) where {CT<:CoherenceParams, ST<:SharedParams}

    @assert length(Ys) == length(Xs)

    for n in eachindex(Ys)
        for i in eachindex(Ys[n].var)
            fill!(Ys[n].var[i], Xs[n].var[i])
        end
    end

    return nothing
end

# approach:
# shared to coherence params.
# coherence params to Phys.
# Phys to flat vectors.
function convertmodel(
    ssp::Vector{SpinSysParams{SharedShift{T},CT,T2T}},
    Δc_bars::Vector{Vector{Vector{Vector{T}}}}, # MSS.Δc_bars
    )::Vector{SpinSysParams{CoherenceShift{T},CT,T2T}} where {T,CT,T2T}

    shifts = collect( ssp[n].shift for n in eachindex(ssp) )
    shifts_c = modelconversion(shifts, Δc_bars)

    return collect( SpinSysParams(shifts_c[n], ssp[n].phase, ssp[n].T2) for n in eachindex(ssp) )
end

function convertmodel(
    ssp::Vector{SpinSysParams{SharedShift{T},CT,T2T}},
    As::Vector{HAM.SHType{T}},
    )::Vector{SpinSysParams{CoherenceShift{T},CT,T2T}} where {T,CT,T2T}

    Δc_bars = collect( As[n].Δc_bar for n in eachindex(As) )
    return convertmodel(ssp, Δc_bars)
end

####### preparation routines for analysis, visualize, and storage tasks.
"""
```
function createtablecolumns(Phys::Vector{NMRHamiltonian.PhysicalParamsType{T}})::Tuple{Vector{Vector{Int}},Vector{Tuple{Int,Int}},Vector{T}} where T
```

Extract information from `Phys` to produce three 1-D arrays that is useful for constructing tables for visualization purposes.

Outputs:

- `nucleus_ID_set` -- nucleus ID numbers, as recorded from `Phys[n].H_IDs`, i.e., the original not the re-labelled IDs.

- `location_set` -- pairs of the format `(n,i)`, `n` is the compound entry index, `i` is the spin system index.

- `cs_set` -- chemical shift numbers for each nucleus, in units ppm.
"""
function createtablecolumns(Phys::Vector{NMRHamiltonian.PhysicalParamsType{T}})::Tuple{Vector{Vector{Int}},Vector{Tuple{Int,Int}},Vector{T}} where T

    nucleus_ID_set = Vector{Vector{Int}}(undef, 0)
    location_set = Vector{Tuple{Int,Int}}(undef, 0) # contains 1-indexing info of molecule entries and spin systems.
    cs_set = Vector{T}(undef, 0)

    for n in eachindex(Phys)

        for i in eachindex(Phys[n].cs_sys)
            
            @assert length(Phys[n].cs_sys[i]) == length(Phys[n].H_inds_sys[i])

            if isempty(Phys[n].ME[i])

                # case: every chemical shift is only associated with a nucleus, and should have a separate table row entry.
                for j in eachindex(Phys[n].cs_sys[i])
                    
                    ind = Phys[n].H_inds_sys[i][j]
                    nucleus_ID = Phys[n].H_IDs[ind]

                    cs = Phys[n].cs_sys[i][j]

                    location = (n,i)

                    push!(nucleus_ID_set, [nucleus_ID;])
                    push!(location_set, location)
                    push!(cs_set, cs)
                end
            else

                # case: every chemical shift is associated with the multiple nucleus in the magnetic equivalence group, and should have a separate table row entry.
                for l in eachindex(Phys[n].ME[i])
                    
                    inds = Phys[n].ME[i][l]
                    nucleus_IDs = Phys[n].H_IDs[inds]

                    
                    tmp = Phys[n].cs_sys[i][inds]
                    cs = tmp[begin]
                    @assert norm(tmp - ones(length(tmp)) .* cs ) < 1e-12 # sanity-check: tmp should be a vector of the same value. 

                    location = (n,i)

                    push!(nucleus_ID_set, nucleus_IDs)
                    push!(location_set, location)
                    push!(cs_set, cs)
                end
            end
        end

        for i in eachindex(Phys[n].cs_singlets)
                
            inds = Phys[n].H_inds_singlets[i]
            nucleus_IDs = Phys[n].H_IDs[inds]

            cs = Phys[n].cs_singlets[i]

            spin_sys_ind = i + length(Phys[n].cs_sys) # singlet spin system numbering starts after the non-singlet spin systems.

            location = (n,spin_sys_ind)

            push!(nucleus_ID_set, nucleus_IDs)
            push!(location_set, location)
            push!(cs_set, cs)
        end
    end

    return nucleus_ID_set, location_set, cs_set
end



############ table construction.


# hz2ppmfunc = uu->hz2ppm(uu, experiment)
function hz2ppm(u::T, fs, SW, ν_0ppm)::T where T
    return (u - ν_0ppm)*SW/fs
end

# ppm2hzfunc = pp->ppm2hz(pp, experiment)
function ppm2hz(p::T, fs, SW, ν_0ppm)::T where T
    return (ν_0ppm + p*fs/SW)
end

# ppm2radfunc = pp->ppm2hzfunc(pp)*2*π
function ppm2rad(p::T, fs, SW, ν_0ppm)::T where T
    return ppm2hz(p, fs, SW, ν_0ppm)*2*π
end

# rad2ppmfunc = ww->hz2ppmfunc(ww/(2*π))
function rad2ppm(ω::T, fs, SW, ν_0ppm)::T where T
    return hz2ppm(ω/(2*π), fs, SW, ν_0ppm)
end


function Δppm2Δrad(Δp::T, fs, SW, ν_0ppm)::T where T
    return ppm2rad(Δp, fs, SW, ν_0ppm) - ppm2rad(zero(T), fs, SW, ν_0ppm)
end

# inverse of Δppm2Δrad().
function Δrad2Δppm(Δω::T, fs, SW, ν_0ppm)::T where T
    return rad2ppm(Δω + ppm2rad(zero(T), fs, SW, ν_0ppm), fs, SW, ν_0ppm)
end


"""
```
exportphysicalparams(
    MSS::CLMixtureSpinSys,
    Phys::Vector{HAM.PhysicalParamsType{T}}, # template.
    fs::T,
    SW::T,
    ν_0ppm::T,
) where T
```

Outputs: out_shift, out_phase, out_T2.

The fields `cs_sys` and `cs_singlets` for `out_shift`, `out_phase`, and `out_T2` does not contain chemical shift. It contains the shift difference (in ppm), phase variable κ_β (in radians), and the T2 variable ξ (a multiplier, dimensionless), respectively.
"""
function exportphysicalparams(
    MSS::CLMixtureSpinSys,
    Phys::Vector{HAM.PhysicalParamsType{T}}, # template.
    fs::T,
    SW::T,
    ν_0ppm::T,
    ) where T

    shifts, phases, T2s = MSS.shifts, MSS.phases, MSS.T2s

    # shift..
    out_shift = deepcopy(Phys)
    
    for n in eachindex(Phys)
        cs_n = HAM.readbasechemshifts(Phys[n])

        @assert length(cs_n) == length(shifts[n].var)
        for i in eachindex(shifts[n].var)

            @assert length(cs_n[i]) == length(shifts[n].var[i])
            for l in eachindex(shifts[n].var[i])

                cs_n[i][l] = Δrad2Δppm(
                    shifts[n].var[i][l],
                    fs,
                    SW,
                    ν_0ppm,
                )
            end
        end

        HAM.writebasechemshifts!(out_shift[n], cs_n)
    end

    # phase..
    out_phase = deepcopy(Phys)

    for n in eachindex(Phys)
        cs_n = HAM.readbasechemshifts(Phys[n])

        @assert length(cs_n) == length(phases[n].var)
        for i in eachindex(phases[n].var)

            @assert length(cs_n[i]) == length(phases[n].var[i])
            for l in eachindex(phases[n].var[i])

                cs_n[i][l] = phases[n].var[i][l] # no transofrmation required. in Radians.
            end
        end

        HAM.writebasechemshifts!(out_phase[n], cs_n)
    end

    # phase..
    out_phase = deepcopy(Phys)

    for n in eachindex(Phys)
        cs_n = HAM.readbasechemshifts(Phys[n])

        @assert length(cs_n) == length(phases[n].var)
        for i in eachindex(phases[n].var)

            @assert length(cs_n[i]) == length(phases[n].var[i])
            for l in eachindex(phases[n].var[i])

                cs_n[i][l] = phases[n].var[i][l] # no transofrmation required. in Radians.
            end
        end

        HAM.writebasechemshifts!(out_phase[n], cs_n)
    end
    
    # T2..
    out_T2 = deepcopy(Phys)

    for n in eachindex(Phys)
        cs_n = HAM.readbasechemshifts(Phys[n])

        @assert length(cs_n) == length(T2s[n].var)
        for i in eachindex(T2s[n].var)

            fill!(cs_n[i], T2s[n].var[i]) # no transofrmation required. in Radians.
        end

        HAM.writebasechemshifts!(out_T2[n], cs_n)
    end
    
    return out_shift, out_phase, out_T2
end

