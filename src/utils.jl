
"""
```
function Δcs2ζ(Δcs::T, ppm2hzfunc)::T where T
```

Convert ppm to radial frequency.
'''
# # test.
a = Δcs2ζ(0.1, ppm2hzfunc)
@show ζ2Δcs(a, ν_0ppm, hz2ppmfunc)
'''
"""
function Δcs2ζ(Δcs::T, ppm2hzfunc)::T where T

    return (ppm2hzfunc(Δcs)-ppm2hzfunc(zero(T)))*2*π
end

"""
```
function ζ2Δcs(ζ::T, ν_0ppm::T, hz2ppmfunc)::T where T
```

Convert radial frequency to ppm.
'''
# # test.
a = Δcs2ζ(0.1, ppm2hzfunc)
@show ζ2Δcs(a, ν_0ppm, hz2ppmfunc)
'''
"""
function ζ2Δcs(ζ::T, ν_0ppm::T, hz2ppmfunc)::T where T
    return hz2ppmfunc(ζ/(2*π) + ν_0ppm)
end


"""
```
function fetchbounds(
    p::MixtureModelParameters,
    Bs::Vector{MoleculeType{T, SST}};
    shift_proportion = 0.9, # between 0 to 1. Control the returned shift bounds as a proportion of the maximum allowed shift bounds used when the surrogates were created.
    phase_lb = convert(T, -π),
    phase_ub = convert(T, π),
)::Tuple{Vector{T},Vector{T}} where {T,SST}
```

Returns a `Vector` for lower bound and `Vector` for upper bound for each parameter variable. The order (first elements to last elements) of the returned vectors are: Shift (ζ), phase (κ_β), then T2 (κ_λ) parameters.

Input:

- `p` -- model parameters. Used to determine the size of the output.

- `Bs` -- Surrogate model. This function uses the `OperationRange` field from each element of `Bs`.

Optional inputs:

- `shift_proportion` -- between 0 to 1. Control the returned shift bounds as a proportion of the maximum allowed shift bounds used when the surrogates were created.

- `phase_lb` -- the fill value for the lower bound of phase parameter (κ_β).

- `phase_ub` -- the fill value for the lower bound of phase parameter (κ_β).

"""
function fetchbounds(
    p::MixtureModelParameters,
    Bs::Vector{MoleculeType{T, SST}};
    shift_proportion = 0.9,
    phase_lb = convert(T, -π),
    phase_ub = convert(T, π),
    )::Tuple{Vector{T},Vector{T}} where {T,SST}

    @assert zero(T) < shift_proportion < one(T)

    mapping = p.systems_mapping

    lbs = Vector{T}(undef, length(p.var_flat))
    ubs = Vector{T}(undef, length(p.var_flat))

    for n in eachindex(mapping.shift.st)

        d_max = Bs[n].op_range.d_max
        κ_λ_lb = Bs[n].op_range.κ_λ_lb
        κ_λ_ub = Bs[n].op_range.κ_λ_ub

        # shift.
        for i in eachindex(mapping.shift.st[n])

            #r_lb = 2*π*(u_min - d_max[i])
            #r_ub = 2*π*(u_max + d_max[i])
            ζ_max = d_max[i]*2*π*shift_proportion

            for l = mapping.shift.st[n][i]:mapping.shift.fin[n][i]
                lbs[l] = -ζ_max #r_lb
                ubs[l] = ζ_max #r_ub
            end
        end

        # phase.
        for i in eachindex(mapping.phase.st[n])

            for l = mapping.phase.st[n][i]:mapping.phase.fin[n][i]
                lbs[l] = phase_lb
                ubs[l] = phase_ub
            end
        end

        # T2 multiplier.
        for i in eachindex(mapping.T2.st[n])

            for l = mapping.T2.st[n][i]:mapping.T2.fin[n][i]
                lbs[l] = κ_λ_lb
                ubs[l] = κ_λ_ub
            end
        end
    end

    return lbs, ubs
end

function initializeparameter(αs::Vector{Vector{T}}, default_val::T) where T <: AbstractFloat
    out = similar(αs)
    
    for n in eachindex(αs)
        out[n] = Vector{T}(undef, length(αs[n]))
        fill!(out[n], default_val)
    end

    return out
end

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

# """
#     gettimerange(N::Int, fs::T) where T

# Returns the time stamps for a sequence, starting at time 0. Returns zero(T):Ts:(N-1)*Ts, Ts = 1/fs.
# """
function gettimerange(N::Int, fs::T) where T
    Ts::T = 1/fs

    return zero(T):Ts:(N-1)*Ts
end

# """
#     getDFTfreqrange(N::Int, fs::T) where T

# Returns the frequency stamps for a DFT sequence, computed by fft().
#     Starting at frequency 0 Hz. Returns LinRange(0, fs-fs/N, N).
# """
function getDFTfreqrange(N::Int, fs::T)::LinRange{T} where T
    a = zero(T)
    b = fs-fs/N

    return LinRange(a, b, N)
end



##### tests.
# mutates shifts, phases, T2s.

function testupdate!(MSS, MS::MixtureSinglets{T}) where T <: AbstractFloat

    # generate random input p.
    N_κs_ζ = getNvars(MSS.shifts)
    N_κs_β = getNvars(MSS.phases)
    N_κs_λ = getNvars(MSS.T2s)
    N_singlets = getNsinglets(MS)

    N_SS_vars = N_κs_ζ+N_κs_β+N_κs_λ
    N_singlet_vars = 3*N_singlets

    p = randn(N_SS_vars + N_singlet_vars)

    # spin system update.
    mapping = ParamsMapping(MSS.shifts, MSS.phases, MSS.T2s)
    updatespinsystems!(MSS, p, mapping)

    # singlets.
    s_mapping = getsingletsParamsMapping(MS, offset_ind = N_SS_vars)
    updatesinglets!(MS, p, s_mapping)
    
    # verify.
    
    discrepancy1, fin_ind = verifyupdatedvars(p, MSS.shifts; st_ind = 1)
    discrepancy2, fin_ind = verifyupdatedvars(p, MSS.phases; st_ind = fin_ind + 1)
    discrepancy3, fin_ind = verifyupdatedvars(p, MSS.T2s; st_ind = fin_ind + 1)
    
    discrepancy4, fin_ind = verifyupdatedvars(p, MS; st_ind = fin_ind+1)
    
    return discrepancy1 + discrepancy2 + discrepancy3 + discrepancy4, p
    #return discrepancy1, discrepancy2, discrepancy3, discrepancy4, p
end



# based on the hardcoded ordering of ParamsMapping(): shift, phase, T2.
# function verifyupdatedvars(
#     p,
#     shifts::Vector{CoherenceShift{T}};
#     st_ind = 1,
#     ) where T

#     @assert length(shifts) == length(phases) == length(T2s)

#     discrepancy = zero(T)
#     l = st_ind

#     for n in eachindex(shifts)
#         for i in eachindex(shifts[n].var)
#             for j in eachindex(shifts[n].var[i])
#                 l += 1
#                 discrepancy += abs(shifts[n].var[i][j] - p[begin+l-1])
#             end
#         end
#     end

#     return discrepancy, l
# end

function verifyupdatedvars(
    p,
    phases::Vector{PT};
    st_ind = 1,
    ) where PT <: CoherenceParams

    discrepancy = 0.0
    l = st_ind

    for n in eachindex(phases)
        for i in eachindex(phases[n].var)
            for j in eachindex(phases[n].var[i])
                
                discrepancy += abs(phases[n].var[i][j] - p[begin+l-1])
                l += 1
            end
        end
    end

    return discrepancy, l-1
end

function verifyupdatedvars(
    p,
    T2s::Vector{ST};
    st_ind = 1,
    ) where ST <: SharedParams

    discrepancy = 0.0
    l = st_ind 

    for n in eachindex(T2s)
        for i in eachindex(T2s[n].var)
            for j in eachindex(T2s[n].var[i])
                
                discrepancy += abs(T2s[n].var[i][j] - p[begin+l-1])
                l += 1
            end
        end
    end

    return discrepancy, l-1
end

function verifyupdatedvars(
    p,
    MS::MixtureSinglets{T};
    st_ind = 1,
    ) where T

    discrepancy = zero(T)
    l = st_ind

    for n in eachindex(MS.ζs)
        for i in eachindex(MS.ζs[n])
            
            discrepancy += abs(MS.ζs[n][i] - p[begin + l-1])
            l += 1
        end
    end

    # for n in eachindex(MS.βs)
    #     for i in eachindex(MS.βs[n])
            
            
    #         discrepancy += abs(MS.βs[n][i] - p[begin+l-1])
    #         l += 1
    #     end
    # end

    # for n in eachindex(MS.ξs)
    #     for i in eachindex(MS.ξs[n])
            
            
    #         discrepancy += abs(MS.ξs[n][i] - p[begin+l-1])
    #         l += 1
    #     end
    # end

    return discrepancy, l-1
end





############### tsting.

function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where T <: Real

    return (x-a)*(d-c)/(b-a)+c
end

function generateparameters(lbs::Vector{T}, ubs::Vector{T})::Vector{T} where T
    
    @assert length(lbs) == length(ubs)

    return collect( convertcompactdomain(rand(T), zero(T), one(T), lbs[i], ubs[i]) for i in eachindex(lbs) )
end

