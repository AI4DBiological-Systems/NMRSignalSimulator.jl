
import NMRHamiltonian

# include("../src/NMRSignalSimulator.jl")
# import .NMRSignalSimulator
# ## remove later.
# import Interpolations, OffsetArrays
# ##

import NMRSignalSimulator

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
println("finished As")
println()

#@assert 1==2

λ0 = 3.4

Δr_default = 1.0 # the samples used to build the surrogate is taken every `Δr` radian on the frequency axis. Decrease for improved accuracy at the expense of computation resources.
Δκ_λ_default = 0.05 # the samples used to build thes urrogate for κ_λ are taken at this sampling spacing. Decrease for improved accuracy at the expense of computation resources.
Δcs_max_scalar_default = 0.2 # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.
κ_λ_lb_default = 0.5 # interpolation lower limit for κ_λ.
κ_λ_ub_default = 2.5 # interpolation upper limit for κ_λ.

surrogate_config_path = "/home/roy/Documents/repo/NMRData/input/select_molecules_surrogate_configs.json"


#dummy_SSFID = NMRSignalSimulator.SpinSysParamsType1(0.0)
dummy_SSFID = NMRSignalSimulator.SpinSysParamsType2(0.0)
# u_min = ppm2hzfunc(-0.5)
# u_max = ppm2hzfunc(4.0)

Bs = NMRSignalSimulator.fitclproxies(As, dummy_SSFID, λ0;
    names = molecule_entries,
    config_path = surrogate_config_path,
    Δcs_max_scalar_default = Δcs_max_scalar_default,
    κ_λ_lb_default = κ_λ_lb_default,
    κ_λ_ub_default = κ_λ_ub_default,
    # u_min = u_min,
    # u_max = u_max,
    Δr_default = Δr_default,
    Δκ_λ_default = Δκ_λ_default)

#
### plot.

# purposely distort the spectra by assigning random values to model parameters.
B = Bs[1]
if typeof(dummy_SSFID) <: NMRSignalSimulator.SpinSysParamsType1
    B.ss_params.d[:] = rand(length(B.ss_params.d))
elseif typeof(dummy_SSFID) <: NMRSignalSimulator.SpinSysParamsType2
    B.ss_params.d[:] = collect( rand(length(B.ss_params.d[i])) .* (2*π) for i = 1:length(B.ss_params.d) )
end
B.ss_params.κs_λ[:] = rand(length(B.ss_params.κs_λ)) .+ 1
B.ss_params.κs_β[:] = collect( rand(length(B.ss_params.κs_β[i])) .* (2*π) for i = 1:length(B.ss_params.κs_β) )


f = uu->NMRSignalSimulator.evalclmixture(uu, As, Bs)

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

# test params.
#ΩS_ppm = collect( hz2ppmfunc.( NMRSignalSimulator.combinevectors(A.Ωs) ./ (2*π) ) for A in As )
#ΩS_ppm_flat = NMRSignalSimulator.combinevectors(ΩS_ppm)


#P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
P = LinRange(-0.1, 6.0, 80000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

## parameters that affect qs.
# A.d, A.κs_λ, A.κs_β
# A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets
q = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs)

#Es = collect( NMRSignalSimulator.καMoleculeType(Bs[i]) for i = 1:length(Bs) )
#q = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Es)

f_U = f.(U_rad)
q_U = q.(U_rad)


discrepancy = abs.(f_U-q_U)
max_val, ind = findmax(discrepancy)
println("relative discrepancy = ", norm(discrepancy)/norm(f_U))
println("max discrepancy: ", max_val)
println()

# ## remove areas with low signal from plotting to reduce plot size.
# reduction_factor = 100
# threshold_factor =  α_relative_lower_threshold/10
# inds, keep_inds, inds_for_reducing = NMRSignalSimulator.prunelowsignalentries(q_U, threshold_factor, reduction_factor)
#
# q_U_display = q_U[inds]
# f_U_display = f_U[inds]
# P_display = P[inds]

q_U_display = q_U
f_U_display = f_U
P_display = P

## visualize.
PyPlot.figure(fig_num)
fig_num += 1

# PyPlot.plot(P, real.(f_U), label = "f")
# PyPlot.plot(P, real.(q_U), label = "q")
PyPlot.plot(P_display, real.(f_U_display), label = "f")
PyPlot.plot(P_display, real.(q_U_display), label = "q")
PyPlot.plot(P_display, real.(q_U_display), "x")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("f vs q")


## visualize. zoom in.

inds = findall(xx->(2.5<xx<3.9), P_display)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_display[inds], real.(f_U_display[inds]), label = "f")
PyPlot.plot(P_display[inds], real.(q_U_display[inds]), label = "q")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("f vs q")


# using BenchmarkTools
#
# m = 1
# A = As[1];
#
# println("qs[i][k], gs eval:")
# r0 = 2*π*U[m] - A.ss_params.d[1]
# @btime A.qs[1][1](r0, 1.0)
# @btime A.gs[1][1](r0, 1.0)
#
# println("q eval.")
# @btime q.(U_rad[m]);
#
# println("q_U eval")
# @btime q.(U_rad);
# println()
