# show the approximation between the free-induction decay singlet and complex Lorentzian singlet.

include("../src/NMRSignalSimulator.jl")
import .NMRSignalSimulator

using LinearAlgebra
using FFTW
import PyPlot

import Random
Random.seed!(25)


include("helpers/utils.jl")

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

fs = 9615.38461538462
SW = 16.0196918511501
ν_0ppm = 6753.577042707225
λ0 = 3.4

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

t2 = NMRSignalSimulator.gettimerange(16384, fs)

Ω = ν_0ppm * 2*π
cl =  uu->(1/(λ0+im*(2*π*uu-Ω)))
fid = tt->exp(-λ0*tt)*exp(im*Ω*tt)

f_t = fid.(t2) ./fs
DFT_f = fft(f_t)

U_DFT = NMRSignalSimulator.getDFTfreqrange(length(t2), fs)
P_DFT = hz2ppmfunc.(U_DFT)

C_U = cl.(U_DFT)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(U_DFT, real.(DFT_f), label = "DFT f", linewidth = "2")
PyPlot.plot(U_DFT, real.(C_U), label = "cL C", "--", linewidth = "2")

PyPlot.legend()
PyPlot.xlabel("Hz")
PyPlot.ylabel("")
PyPlot.title("real part DFT-FID vs complex lorentzian, ν at $(ν_0ppm) Hz")
