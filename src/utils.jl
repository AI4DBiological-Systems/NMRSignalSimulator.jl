

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

"""
    gettimerange(N::Int, fs::T) where T

Returns the time stamps for a sequence, starting at time 0. Returns zero(T):Ts:(N-1)*Ts, Ts = 1/fs.
"""
function gettimerange(N::Int, fs::T) where T
    Ts::T = 1/fs

    return zero(T):Ts:(N-1)*Ts
end

"""
    getDFTfreqrange(N::Int, fs::T) where T

Returns the frequency stamps for a DFT sequence, computed by fft().
    Starting at frequency 0 Hz. Returns LinRange(0, fs-fs/N, N).
"""
function getDFTfreqrange(N::Int, fs::T)::LinRange{T} where T
    a = zero(T)
    b = fs-fs/N

    return LinRange(a, b, N)
end



##### tests.
# mutates shifts, phases, T2s.

function testupdate!(MSS, MS::MixtureSinglets{T}) where T <: AbstractFloat

    # generate random input p.
    N_κs_d = getNvars(MSS.shifts)
    N_κs_β = getNvars(MSS.phases)
    N_κs_λ = getNvars(MSS.T2s)
    N_singlets = getNsinglets(MS)

    N_SS_vars = N_κs_d+N_κs_β+N_κs_λ
    N_singlet_vars = 3*N_singlets

    p = randn(N_SS_vars + N_singlet_vars)

    # spin system update.
    mapping = getParamsMapping(MSS.shifts, MSS.phases, MSS.T2s)
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



# based on the hardcoded ordering of getParamsMapping(): shift, phase, T2.
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

    for n in eachindex(MS.ds)
        for i in eachindex(MS.ds[n])
            
            discrepancy += abs(MS.ds[n][i] - p[begin + l-1])
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