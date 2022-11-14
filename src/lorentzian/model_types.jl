

struct MixtureModel{T,SST}

    var_flat::Vector{T}

    # shift_flat::Vector{T}
    # phase_flat::Vector{T}
    # T2_flat::Vector{T}
    
    mapping::ParamsMapping

    MSS::MixtureSpinSys{T,SST} # mixture parameters, non-singlet spin systems.
    MS::MixtureSinglets{T} # mixture parameters, singlets.
end

function MixtureModel(MSS::MixtureSpinSys{T,SST}, MS) where {T,SST}

    #
    shifts, phases, T2s = MSS.shifts, MSS.phases, MSS.T2s
    mapping = getParamsMapping(shifts, phases, T2s)
    
    N_κs_d = getNvars(shifts)
    N_κs_β = getNvars(phases)
    N_κs_λ = getNvars(T2s)
    var_flat = zeros(T, N_κs_d+N_κs_β+N_κs_λ)

    return MixtureModel(var_flat, mapping, MSS, MS)
end

# function evalforwardmodel2(
#     u_rad,
#     #κ_d,
#     #κ_β,
#     #κ_λ,
#     p::MixtureModel{T,SST},
#     MSS,
#     ) where {T,SST}
    
#     # setup.
#     mapping, MSS = p.MSS, p.mapping
    
#     # eval.
#     out = zero(Complex{T})

#     for n in eachindex(phases)
#         out += evalforwardmodel(
#             u_rad,
#             MSS.srs[n],
#             MSS.sis[n],
#             MSS.shifts[n],
#             MSS.phases[n],
#             MSS.T2s[n],
#         )
#     end

#     return out
# end

# no error checking on the lengths of fieldnames in MSS.
function evalforwardmodel(
    u_rad,
    MSS::MixtureSpinSys{T,SST},
    #As,
    #Bs,
    MS
    )::Complex{T} where {T,SST}
    
    # eval.
    out_sys = zero(Complex{T})

    for n in eachindex(MSS.shifts)
        out_sys += evalforwardmodel(
            u_rad,
            MSS.srs[n],
            MSS.sis[n],
            MSS.shifts[n],
            MSS.phases[n],
            MSS.T2s[n],
        )
    end

    out_singlets = zero(Complex{T})
    # for n in eachindex(Bs)
    #     A = As[n]
    #     B = Bs[n]
    #     out_singlets += NMRSignalSimulator.evalclsinglets(
    #         u_rad,
    #         B.d_singlets,
    #         A.αs_singlets,
    #         A.Ωs_singlets,
    #         B.β_singlets,
    #         B.λ0,
    #         B.κs_λ_singlets
    #     )
    # end
    for n in eachindex(MS.ds)

        out_singlets += NMRSignalSimulator.evalclsinglets(
            u_rad,
            MS.ds[n],
            MS.αs[n],
            MS.Ωs[n],
            MS.βs[n],
            MS.λ0,
            MS.ξs[n]
        )
    end

    return out_sys + out_singlets
end

function evalforwardmodel(
    u_rad::T,
    #ξ::T,
    sr,
    si,
    shift::CoherenceShift{T}, # contains d.
    phase::CoherencePhase{T}, # contains β.
    T2::SharedT2{T}, # contains ξ.
    ) where T <: AbstractFloat

    out = zero(Complex{T})

    # the number of spin systems should be the same.
    @assert length(phase.cos_β) == length(phase.sin_β) == length(shift.d) == length(T2.var)

    for i in eachindex(phase.cos_β)

        # the number of resonance groups should be the same.
        @assert length(phase.cos_β[i]) == length(phase.sin_β[i]) == length(shift.d[i])

        ξ = T2.var[i]

        for k in eachindex(phase.cos_β[i])
            r = u_rad - shift.d[i][k]

            out += evalq(
                sr[i][k],
                si[i][k],
                r,
                ξ,
                phase.cos_β[i][k],
                phase.sin_β[i][k],
            )
        
        end
    end

    return out
end