

function combinevectors(x::Vector{Vector{T}})::Vector{T} where T

    if isempty(x)
        return Vector{T}(undef, 0)
    end

    N = sum(length(x[i]) for i in eachindex(x))

    y = Vector{T}(undef,N)

    st_ind = 0
    fin_ind = 0
    for i in eachindex(x)
        st_ind = fin_ind + 1
        fin_ind = st_ind + length(x[i]) - 1

        y[st_ind:fin_ind] = x[i]
    end

    return y
end



###### forward model and derivatives.

# This should be the same as evalclproxymixture()
function forwardmodelmolecule(
    u_rad::T,
    As,
    Bs,
    ) where T <: AbstractFloat

    out = zero(Complex{T})
    for n in eachindex(Bs)
        A = As[n]
        B = Bs[n]

        # d = x.d
        # κs_λ = x.κs_λ

        # for i in eachindex(Bs[n].qs)
            
        #     r = u_rad - d[i]
            
        #     for k in eachindex(Bs[n].qs[i])
        #         out += qs[i][k](r, κs_λ[i])
        #     end
        # end

        out_sys = NMRSignalSimulator.evalclproxysys(B.qs, u_rad, B.ss_params)

        out_singlets = NMRSignalSimulator.evalclsinglets(u_rad, B.d_singlets, A.αs_singlets, A.Ωs_singlets,
        B.β_singlets, B.λ0, B.κs_λ_singlets)

        out += out_sys + out_singlets
    end

    return out
end