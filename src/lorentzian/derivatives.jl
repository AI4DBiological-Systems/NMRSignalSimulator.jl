# derivatives for itp.

# coherenceshift, coherencephase, sharedT2.
# ∂fR_∂x in notes.
function derivative!(
    buf_r::Vector{T},
    buf_i::Vector{T},
    x::T,
    κ_λ::T,
    λ0,
    #Bs::Vector{MoleculeType{T, SpinSysParams{CoherenceShift{T}, CoherenceParams{T}, SharedT2{T}}}},
    ∂fr!,
    ∂fi!,
    fr, # used for derivatives of β.
    fi,
    #β,
    # shift::CoherenceShift{T},
    phases::Vector{CoherencePhase{T}},
    # T2::SharedT2{T},
    n::Int,
    i::Int,
    k::Int,
    ) where T

    @assert length(buf_r) == 2 == length(buf_i)
    
    ∂fr!(buf_r, x, κ_λ)
    ∂fi!(buf_i, x, κ_λ)

    sin_β = phases[n].sin_β[i][k] #sin(β)
    cos_β = phases[n].cos_β[i][k] #cos(β)

    ∂qr_∂x = buf_r[begin]*cos_β - buf_i[begin]*sin_β
    ∂qr_∂κ_λ = buf_r[begin+1]*λ0*cos_β - buf_i[begin+1]*λ0*sin_β
    ∂qr_∂β = -fr(x, κ_λ)*sin_β -fi(x, κ_λ)*cos_β

    #∂qr_∂κ_β = 

    return ∂qr_∂x, ∂qr_∂κ_λ, ∂qr_∂β
    
    # for i in eachindex(qs)
    #     for k in eachindex(qs[i])

    #         ∂fr[i][k](x, κ_λ)
    #     end

    # end

    #return nothing
end

# I am here.
function forwardeval(
    MSS,
    u_rad,
    d,
    ξ,
    β
    )::Complex{T} where T <: AbstractFloat

    

    for n in eachindex(MSS.sr)
        for i in eachindex(MSS.sr[i])
            for k in eachindex(MSS.sr[i])
                
                sr, si = MSS.sr[n][i][k], MSS.sis[n][i][k]
                
                β = dot(MSS.Δc_bars[n][i][k], MSS.phases[n][i][k].κs_β)

                MSS.phases[n][i][k].cos_β = cos(β)
                MSS.phases[n][i][k].sin_β = sin(β)

                out += evalq(sr, si, r, ξ)

            end
        end
    end

    return nothing
end