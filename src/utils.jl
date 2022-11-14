

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
function testupdateparameters1!(shifts, phases, T2s)

    mapping = getParamsMapping(shifts, phases, T2s)

    N_κs_d = getNvars(shifts)
    N_κs_β = getNvars(phases)
    N_κs_λ = getNvars(T2s)

    p = randn(N_κs_d+N_κs_β+N_κs_λ)

    updateparameters!(T2s, p, mapping.T2)
    updateparameters!(shifts, p, mapping.shift)
    updateparameters!(phases, p, mapping.phase)

    return p, mapping
end

# based on the hardcoded ordering of getParamsMapping(): shift, phase, T2.
function testupdateparameters2(
    p,
    shifts::Vector{CoherenceShift{T}},
    phases,
    T2s
    ) where T

    @assert length(shifts) == length(phases) == length(T2s)

    discrepancy = zero(T)

    l = 0
    for n in eachindex(shifts)
        for i in eachindex(shifts[n].var)
            for j in eachindex(shifts[n].var[i])
                l += 1
                discrepancy += abs(shifts[n].var[i][j] - p[l])
            end
        end
    end

    for n in eachindex(phases)
        for i in eachindex(phases[n].var)
            for j in eachindex(phases[n].var[i])
                l += 1
                discrepancy += abs(phases[n].var[i][j] - p[l])
            end
        end
    end

    for n in eachindex(T2s)
        for i in eachindex(T2s[n].var)
            for j in eachindex(T2s[n].var[i])
                l += 1
                discrepancy += abs(T2s[n].var[i][j] - p[l])
            end
        end
    end

    return discrepancy
end