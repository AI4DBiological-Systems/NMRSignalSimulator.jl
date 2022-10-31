import NMRHamiltonian

using DataDeps
import Tar

using LinearAlgebra
import PyPlot
#import JSON3

include("./helpers/data.jl")
include("./helpers/SH.jl")

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

### user inputs.

#molecule_entries = ["L-Methionine"; "L-Phenylalanine"; "DSS"; "Ethanol"; "L-Isoleucine"]
molecule_entries = ["alpha-D-Glucose"; "beta-D-Glucose"; "DSS"; "D2O"]

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665

# # machine values for the BMRB 500 MHz glucose experiment.
# ν_0ppm = 6752.490995937095
# SW = 16.0196917451925
# fs = 9615.38461538462

### end inputs.

As, Rs = runSH(molecule_entries)

#@assert 1==2

### visualize a target molecule and spin group, given a T-2 decay parameter.
λ0 = 3.4 # T-2 decay parameter.
molecule_select = 1
spin_system_select = 1
ppm_offset = 0.2 # for display padding.
N_viz = 50000

a = As[molecule_select].αs[spin_system_select]
F = As[molecule_select].Ωs[spin_system_select]

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

ΩS_ppm = hz2ppmfunc.( F ./ (2*π) )
ΩS_ppm_sorted = sort(ΩS_ppm)

u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - ppm_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + ppm_offset)

P_min = hz2ppmfunc(u_min)
P_max = hz2ppmfunc(u_max)

P = LinRange(P_min, P_max, N_viz)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

# absorption Lorentzian.

function evalzerophasecl1Darray(u_rad, αs::Vector{T}, Ωs::Vector{T}, λ::T)::Complex{T} where T <: AbstractFloat

    out = zero(Complex{T})
    for l = 1:length(αs)
        out += evalzerophaseclpartitionelement(u_rad, αs[l], Ωs[l], λ)
    end

    return out
end

function evalzerophaseclpartitionelement(r,
    α::T, Ω::T, λ::T)::Complex{T} where T <: AbstractFloat

    out = α/(λ+im*(r-Ω))

    return out
end

q = uu->evalzerophasecl1Darray(uu, a, F, λ0)
q_U = q.(U_rad)

# plot.
PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, q_U)
PyPlot.gca().invert_xaxis()

PyPlot.legend()
PyPlot.title("spectrum of $(molecule_entries[molecule_select]), spin system $(spin_system_select)")
