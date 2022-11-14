

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

    var_flat = Vector{T}(undef, N_κs_d+N_κs_β+N_κs_λ)


    return MixtureModel(var_flat, mapping, MSS, MS)
end


function updatemodel!(
    p::Vector{T},
    MSS,
    MS,
    mapping,
    ) where T <: AbstractFloat

    for n in eachindex(MSS.phases)

        updateparameters!(MSS.T2s, p, mapping.T2)
        updateparameters!(MSS.shifts, p, mapping.shift)
        updateparameters!(MSS.phases, p, mapping.phase)

        #updateparameters!(MSS.phases, p, mapping.phase)
    end

    return nothing
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
    MS::MixtureSinglets{T},
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
    for n in eachindex(MS.ds)

        # TODO: make evalclsinglets utilize cos_β_singlets to avoid complex multiply for speed.
        # only save ~10% execution time though.
        out_singlets += evalclsinglets(
            u_rad,
            MS.ds[n],
            MS.αs[n],
            MS.Ωs[n],
            MS.βs[n],
            MS.λ0,
            MS.ξs[n],
        )
    end

    return out_sys + out_singlets
end

function evalforwardmodel(
    u_rad::T,
    sr,
    si,
    shift::SharedShift{T}, # contains d.
    phase::CoherencePhase{T}, # contains β.
    T2::SharedT2{T}, # contains ξ.
    ) where T <: AbstractFloat

    out = zero(Complex{T})

    # # tmp storage.
    r = zero(T)
    ξ = one(T)

    for i in eachindex(phase.cos_β)
        
        ξ = T2.var[i]
        r = u_rad - shift.var[i]

        for k in eachindex(phase.cos_β[i])

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
    # out_r = zero(T)
    # out_i = zero(T)

    # # tmp storage.
    # real_part = zero(T)
    # imag_part = zero(T)
    #out_buf = zeros(T,2)
    r = zero(T)

    # the number of spin systems should be the same.
    #@assert length(phase.cos_β) == length(phase.sin_β) == length(shift.d) == length(T2.var)

    for i in eachindex(phase.cos_β)

        # the number of resonance groups should be the same.
        #@assert length(phase.cos_β[i]) == length(phase.sin_β[i]) == length(shift.d[i])

        ξ = T2.var[i]

        for k in eachindex(phase.cos_β[i])
            r = u_rad - shift.d[i][k]

            # 300 ms.
            out += evalq(
                sr[i][k],
                si[i][k],
                r,
                ξ,
                phase.cos_β[i][k],
                phase.sin_β[i][k],
            )

            # # slower. 400 ms.
            # real_part, imag_part = evalq2(
            #     sr[i][k](r,ξ),
            #     si[i][k](r,ξ),
            #     phase.cos_β[i][k],
            #     phase.sin_β[i][k],
            # )
            # out_r += real_part
            # out_i += imag_part
            
            # ## very slow. 700 ms.
            # real_part, imag_part = evalq3(
            #     sr[i][k],
            #     si[i][k],
            #     r,
            #     ξ,
            #     phase.cos_β[i][k],
            #     phase.sin_β[i][k],
            # )
            # out_r += real_part
            # out_i += imag_part
            
            # # slower 420 ms
            # evalq2!(
            #     out_buf,
            #     sr[i][k](r,ξ),
            #     si[i][k](r,ξ),
            #     phase.cos_β[i][k],
            #     phase.sin_β[i][k],
            # )
            
            # # second fastest. 309ms
            # evalq3!(
            #     out_buf,
            #     sr[i][k],
            #     si[i][k],
            #     r,
            #     ξ,
            #     phase.cos_β[i][k],
            #     phase.sin_β[i][k],
            # )

            ## very slow. 1 sec.
            # out_r += evalqr(
            #     sr[i][k](r,ξ),
            #     si[i][k](r,ξ),
            #     phase.cos_β[i][k],
            #     phase.sin_β[i][k],
            # )
            # out_i += evalqi(
            #     sr[i][k](r,ξ),
            #     si[i][k](r,ξ),
            #     phase.cos_β[i][k],
            #     phase.sin_β[i][k],
            # )

        end
    end
    #out = Complex(out_r, out_i)
    #out = Complex(out_buf[begin], out_buf[begin+1])

    return out
end

function evalq2(fr_eval, fi_eval, cos_β, sin_β::T)::Tuple{T,T} where T
    return fr_eval*cos_β - fi_eval*sin_β, fr_eval*sin_β + fi_eval*cos_β
end

function evalq2!(out::Vector{T}, fr_eval::T, fi_eval::T, cos_β::T, sin_β::T) where T <: AbstractFloat
    
    out[begin] += fr_eval*cos_β - fi_eval*sin_β
    out[begin+1] += fr_eval*sin_β + fi_eval*cos_β

    return nothing
end

function evalq3(fr, fi, r::T, ξ::T, c::T, d::T)::Tuple{T,T} where T <: AbstractFloat

    a = fr(r,ξ)
    b = fi(r,ξ)
    return a*c-b*d, a*d+b*c
end

function evalq3!(out::Vector{T}, fr, fi, r::T, ξ::T, c::T, d::T) where T <: AbstractFloat

    a = fr(r,ξ)
    b = fi(r,ξ)
    out[begin] += a*c-b*d
    out[begin+1] += a*d+b*c
    
    return nothing
end

function evalqr(fr_eval, fi_eval, cos_β, sin_β::T)::T where T
    return fr_eval*cos_β - fi_eval*sin_β
end

function evalqi(fr_eval, fi_eval, cos_β, sin_β::T)::T where T
    return fr_eval*sin_β + fi_eval*cos_β
end