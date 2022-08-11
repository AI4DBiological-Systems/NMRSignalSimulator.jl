

function combinevectors(x::Vector{Vector{T}})::Vector{T} where T

    if isempty(x)
        return Vector{T}(undef, 0)
    end

    N = sum(length(x[i]) for i = 1:length(x))

    y = Vector{T}(undef,N)

    st_ind = 0
    fin_ind = 0
    for i = 1:length(x)
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
