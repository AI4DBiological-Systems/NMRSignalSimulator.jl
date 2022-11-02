
# Under construction. paused.

import NMRHamiltonian

include("../src/NMRSignalSimulator.jl")
import .NMRSignalSimulator

using LinearAlgebra
using FFTW
import PyPlot
import JSON

#import Clustering
import Statistics

import Random
Random.seed!(25)


include("helpers/utils.jl")

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

### user inputs.

tol_coherence = 1e-2 # resonances are pairs of eigenvalues of the Hamiltonian that have quantum numbers that differ by -1. This is the tolerance away from -1 that is allowed.
α_relative_lower_threshold = 0.05 # resonances with relative amplitude less than this factor compared to the maximum resonance in the spin group will be removed. Set to 0.0 to see every single resonance component.
Δc_partition_radius = 0.3 # determines how many resonances get grouped together. Larger number means less groups and thus more resonances per group.
λ0 = 3.4

Δr_default = 1.0 # the samples used to build the surrogate is taken every `Δr` radian on the frequency axis. Decrease for improved accuracy at the expense of computation resources.
Δκ_λ_default = 0.05 # the samples used to build thes urrogate for κ_λ are taken at this sampling spacing. Decrease for improved accuracy at the expense of computation resources.
Δcs_max_scalar_default = 0.2 # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.
κ_λ_lb_default = 0.5 # interpolation lower limit for κ_λ.
κ_λ_ub_default = 2.5 # interpolation upper limit for κ_λ.

SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_molecules_SH_configs.json"
surrogate_config_path = "/home/roy/Documents/repo/NMRData/input/surrogate_configs/select_molecules_SH_configs.json"

molecule_entries = ["L-Serine"; "L-Phenylalanine"; "Ethanol"; "L-Isoleucine"; "DSS"; "D2O"]
#molecule_entries = ["D-(+)-Glucose"; "DSS"]
#molecule_entries = ["L-Serine";]
#molecule_entries = ["L-Serine"; "D2O";]
#molecule_entries = ["D2O";]

# # machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
# fs = 14005.602240896402
# SW = 20.0041938620844
# ν_0ppm = 10656.011933076665

# # machine values for the BMRB 500 MHz glucose experiment.
# ν_0ppm = 6752.490995937095
# SW = 16.0196917451925
# fs = 9615.38461538462

# machine values for the NRC-2022 mixture experiment. 600 MHz.
ν_0ppm = 6753.577042707225
SW = 16.0196918511501
fs = 9615.38461538462

#t = LinRange(0, 1.693536, 16384)
t = NMRSignalSimulator.gettimerange(16384, fs)
#t = NMRSignalSimulator.gettimerange(16384*5, fs) # for smaller λ.

# path to the json file that provides the mapping from a molecule name to its spin system info file name.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_molecule_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/molecule_mapping/select_molecules.json")

### end inputs.

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

println("Timing: getphysicalparameters")
@time Phys = NMRHamiltonian.getphysicalparameters(molecule_entries,
    H_params_path,
    dict_molecule_to_filename;
    unique_cs_atol = 1e-6)

#
println("Timing: setupmixtureSH()")
@time As = NMRHamiltonian.setupmixtureSH(molecule_entries,
    fs, SW, ν_0ppm,
    Phys;
    config_path = SH_config_path,
    tol_coherence = tol_coherence,
    α_relative_lower_threshold = α_relative_lower_threshold,
    Δc_partition_radius = Δc_partition_radius)

#dummy_SSParams = NMRSignalSimulator.SharedShift(0.0)
dummy_SSParams = NMRSignalSimulator.CoherenceShift(0.0)
# u_min = ppm2hzfunc(-0.5)
# u_max = ppm2hzfunc(4.0)
#u_min = ppm2hzfunc(3.5)
#u_max = ppm2hzfunc(4.2)

println("fitclproxies():")
@time Bs_cl = NMRSignalSimulator.fitclproxies(As, dummy_SSParams, λ0;
    names = molecule_entries,
    config_path = surrogate_config_path,
    Δcs_max_scalar_default = Δcs_max_scalar_default,
    κ_λ_lb_default = κ_λ_lb_default,
    κ_λ_ub_default = κ_λ_ub_default,
    #u_min = u_min,
    #u_max = u_max,
    Δr_default = Δr_default,
    Δκ_λ_default = Δκ_λ_default)

#
#t_test = t[1:500] # 150 sec for glucose + DSS. v 1.
#t_test = t[1:8000] # v2.
t_test = t
delta_t = 1/fs

println("fitFIDproxies():")
@time Bs = NMRSignalSimulator.fitFIDproxies(As, dummy_SSParams, λ0;
    names = molecule_entries,
    config_path = surrogate_config_path,
    Δcs_max_scalar_default = 0.2,
    #t_lb_default = t[1],
    #t_ub_default = t[end],
    t_lb_default = t_test[1],
    t_ub_default = t_test[end],
    #u_min = u_min,
    #u_max = u_max,
    Δr_default = 1.0,
    #Δr_default = 0.1,
    Δt_default = delta_t*0.5)


### plot.

# # purposely distort the spectra by assigning random values to model parameters.

####  manual assignment of FID parameters.
# ## type 1.
# for n_select in eachindex(Bs)
#     tmp_d = 100.0 .* ones(length(Bs[n_select].ss_params.d))
#     Bs[n_select].ss_params.d[:] = tmp_d
#     Bs_cl[n_select].ss_params.d[:] = tmp_d
# end

# ## type 2.
for n_select in eachindex(Bs)
    for i_select in eachindex(Bs[n_select].ss_params.d)
        tmp_d = 100.0 .* ones(length(Bs[n_select].ss_params.d[i_select]))
        Bs[n_select].ss_params.d[i_select] = tmp_d
        Bs_cl[n_select].ss_params.d[i_select] = tmp_d
    end
end

## common to both types.
for n_select in eachindex(Bs)
    tmp_λ = rand(length(Bs[n_select].ss_params.κs_λ)) .+ 1
    Bs[n_select].ss_params.κs_λ[:] = tmp_λ
    Bs_cl[n_select].ss_params.κs_λ[:] = tmp_λ

    tmp_β = collect( rand(length(Bs[n_select].ss_params.κs_β[i])) .* (2*π) for i in eachindex(Bs[n_select].ss_params.κs_β) )
    Bs[n_select].ss_params.κs_β[:] = tmp_β
    Bs_cl[n_select].ss_params.κs_β[:] = tmp_β
end

# manual move the D2O singlet, assumed to be in the last position.
for n_select in eachindex(Bs)
    tmp_d_singlets = 20 .* rand(length(Bs[n_select].d_singlets))
    Bs[n_select].d_singlets[:] = tmp_d_singlets
    Bs_cl[n_select].d_singlets[:] = tmp_d_singlets

    # if λ is too small relative to the length of t, DFT-FID and CL models won't agree.
    # λ is large enough if the FID signal decays close to zero before end of t.
    tmp_λ_singlet = rand(length(Bs[n_select].κs_λ_singlets)) .+ 1
    #tmp_λ_singlet[1] = 0.2657313746474563 # small λ.
    Bs[n_select].κs_λ_singlets[:] = tmp_λ_singlet
    Bs_cl[n_select].κs_λ_singlets[:] = tmp_λ_singlet

    tmp_β_singlet = rand(length(Bs[n_select].β_singlets))
    Bs[n_select].β_singlets[:] = tmp_β_singlet
    Bs_cl[n_select].β_singlets[:] = tmp_β_singlet
end

#### end manual assignment of FID parameters.

w = rand(length(Bs))

B = Bs[1]
f = tt->NMRSignalSimulator.evalFIDmixture(tt, As, Bs; w = w)
f_t = f.(t_test)

q = tt->NMRSignalSimulator.evalFIDproxymixture(tt, As, Bs; w = w)
q_t = q.(t_test)
#q_t=f_t

discrepancy = norm(q_t-f_t)
println("discrepancy = ", discrepancy)
println()

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(t_test, real.(f_t), label = "f", linewidth = "2")
PyPlot.plot(t_test, real.(q_t), label = "q", linewidth = "2", "--")
PyPlot.plot(t_test, real.(f_t), "x")

PyPlot.legend()
PyPlot.xlabel("sec")
PyPlot.ylabel("intensity")
PyPlot.title("real part")


##### compare against complex Lorentzian.





U_DFT = NMRSignalSimulator.getDFTfreqrange(length(q_t), fs)
P_DFT = hz2ppmfunc.(U_DFT)

Q = uu->NMRSignalSimulator.evalclproxymixture(uu, As, Bs_cl; w = w)
Q_U = Q.(U_DFT .* (2*π))

DFT_q = fft(q_t) ./ fs

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_DFT, real.(DFT_q), label = "DFT q", linewidth = "2")
PyPlot.plot(P_DFT, real.(Q_U), label = "cL Q", "--", linewidth = "2")

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("DFT-FID vs. cL, 0 ppm at $(ν_0ppm) Hz")


import NMRDataSetup
offset_Hz = ν_0ppm - (ppm2hzfunc(0.3)-ppm2hzfunc(0.0))
_, U_q, U_inds_q = NMRDataSetup.getwraparoundDFTfreqs(length(q_t), fs, offset_Hz)

DFT_q_wrap = (fft(q_t) ./ fs)[U_inds_q]
P_DFT_wrap = hz2ppmfunc.(U_q)
Q_U_wrap = Q.(U_q .* (2*π))

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P_DFT_wrap, real.(DFT_q_wrap), label = "DFT q", linewidth = "2")
PyPlot.plot(P_DFT_wrap, real.(Q_U_wrap), label = "cL Q", "--", linewidth = "2")
PyPlot.gca().invert_xaxis()

PyPlot.legend()
PyPlot.xlabel("ppm")
PyPlot.ylabel("real")
PyPlot.title("wrapped, DFT-FID vs. cL, 0 ppm at $(ν_0ppm) Hz")


function getΩSppm(As::Vector{NMRHamiltonian.SHType{T}}, hz2ppmfunc) where T

    ΩS_ppm = Vector{Vector{T}}(undef, length(As))

    for (n,A) in enumerate(As)

        ΩS_ppm[n] = hz2ppmfunc.( combinevectors(A.Ωs) ./ (2*π) )

        tmp = hz2ppmfunc.( A.Ωs_singlets ./ (2*π) )
        push!(ΩS_ppm[n], tmp...)
    end

    return ΩS_ppm
end

ΩS_ppm = getΩSppm(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(NMRSignalSimulator.combinevectors(ΩS_ppm))


tmp = copy(Q_U_wrap)
tmp[U_inds_q] = Q_U_wrap
invDFT_Q = ifft(tmp .* fs) #.* fs #.* exp(im*2*π*ν_0ppm)

PyPlot.figure(fig_num)
fig_num += 1


PyPlot.plot(t_test, real.(q_t), label = "q", linewidth = "2")
PyPlot.plot(t_test, real.(invDFT_Q), label = "invDFT Q", "--", linewidth = "2")

PyPlot.legend()
PyPlot.xlabel("sec")
PyPlot.ylabel("real")
PyPlot.title("invDFT cL (wrapped) vs. FID")

println("q_t vs. invDFT_Q descrepancy: ", norm(real.(q_t)-real.(invDFT_Q)))
println()

### abs will be the same even if we don't wrap the frequency, since it amounts to a *cis(constant) for all time domain values.
PyPlot.figure(fig_num)
fig_num += 1


PyPlot.plot(t_test, abs.(q_t), label = "q", linewidth = "2")
PyPlot.plot(t_test, abs.(invDFT_Q), label = "invDFT Q", "--", linewidth = "2")

PyPlot.legend()
PyPlot.xlabel("sec")
PyPlot.ylabel("abs")
PyPlot.title("invDFT cL (wrapped) vs. FID")
