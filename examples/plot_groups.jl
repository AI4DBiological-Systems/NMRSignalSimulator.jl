
import NMRHamiltonian

include("../src/NMRSignalSimulator.jl")
import .NMRSignalSimulator

using LinearAlgebra
using FFTW
import PyPlot
import JSON

#import Clustering
import Statistics

import PlotlyJS
using Plots; plotly()

include("../examples/helpers/plot_partition.jl")
include("../examples/helpers/utils.jl")


import Random
Random.seed!(25)

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

### user inputs.

save_folder = "/home/roy/MEGAsync/outputs/NMR/groups"

tol_coherence = 1e-2 # resonances are pairs of eigenvalues of the Hamiltonian that have quantum numbers that differ by -1. This is the tolerance away from -1 that is allowed.
#α_relative_threshold = 0.05 # resonances with relative amplitude less than this factor compared to the maximum resonance in the spin group will be removed. Set to 0.0 to see every single resonance component.
α_relative_threshold = 0.0
Δc_partition_radius = 0.3 # determines how many resonances get grouped together. Larger number means less groups and thus more resonances per group.
λ0 = 3.4

Δr_default = 1.0 # the samples used to build the surrogate is taken every `Δr` radian on the frequency axis. Decrease for improved accuracy at the expense of computation resources.
Δκ_λ_default = 0.05 # the samples used to build thes urrogate for κ_λ are taken at this sampling spacing. Decrease for improved accuracy at the expense of computation resources.
Δcs_max_scalar_default = 0.2 # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.
κ_λ_lb_default = 0.5 # interpolation lower limit for κ_λ.
κ_λ_ub_default = 2.5 # interpolation upper limit for κ_λ.

SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs_reduce.json"
surrogate_config_path = "/home/roy/Documents/repo/NMRData/input/surrogate_configs/select_compounds_SH_configs.json"

#molecule_names = ["L-Valine"; ]
#molecule_names = ["L-Isoleucine"; ]
#molecule_names = ["L-Leucine"; ]
#molecule_names = ["alpha-D-Glucose"; ]
#molecule_names = ["beta-D-Glucose"; ]
#molecule_names = ["L-Phenylalanine"; ]
#molecule_names = ["L-Glutamine"; ]
#molecule_names = ["Ethanol"; ]
#molecule_names = ["L-Serine"; ]
#molecule_names = ["DSS"; ]
#molecule_names = ["ATP"; ] # idea: coherence is the compensation for intensity.
molecule_names = ["L(-)-Glutathione, oxidized"; ]
# molecule_names = ["L(-)-Glutathione, reduced"; ]
# molecule_names = ["beta-Alanine"; ]
# molecule_names = ["L-Alanine"; ]

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665

# # # machine values taken from the BMRB 700 MHz 20 mM ATP experiment.
# fs = 14005.602240896402
# SW = 20.0041938620844
# ν_0ppm = 10654.75043163717


# # machine values for the BMRB 500 MHz glucose experiment.
# ν_0ppm = 6752.490995937095
# SW = 16.0196917451925
# fs = 9615.38461538462

# path to the json file that provides the mapping from a compound name to its spin system info file name.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

### end inputs.

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

println("Timing: mag equivalence")
@time MEs = NMRHamiltonian.getmageqinfomixture(molecule_names,
    H_params_path,
    dict_compound_to_filename;
    unique_cs_atol = 1e-6)

#
println("Timing: setupmixtureSH()")
@time mixture_params = NMRHamiltonian.setupmixtureSH(molecule_names,
    H_params_path, dict_compound_to_filename, fs, SW,
    ν_0ppm;
    MEs = MEs,
    config_path = SH_config_path,
    # tol_coherence = tol_coherence,
    # α_relative_threshold = α_relative_threshold,
    # Δc_partition_radius = Δc_partition_radius,
    prune_combo_Δc_bar_flag = true)
As = mixture_params

##prunecombocoherences!(As[1], α_relative_threshold, tol_coherence, Δc_partition_radius)
#prunecombocoherencesbar!(As[1], α_relative_threshold, tol_coherence, Δc_partition_radius)

dummy_SSFID = NMRSignalSimulator.SpinSysParamsType1(0.0)

## frequency locations. For plotting.
ΩS_ppm = getPsnospininfo(mixture_params, hz2ppmfunc)
ΩS_ppm_sorted = sort(combinevectors(ΩS_ppm))

u_offset = 0.2 # in units ppm.
u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)

println("fitproxies():")
@time Bs = NMRSignalSimulator.fitproxies(As, dummy_SSFID, λ0;
    names = molecule_names,
    config_path = surrogate_config_path,
    Δcs_max_scalar_default = Δcs_max_scalar_default,
    κ_λ_lb_default = κ_λ_lb_default,
    κ_λ_ub_default = κ_λ_ub_default,
    u_min = u_min,
    u_max = u_max,
    Δr_default = Δr_default,
    Δκ_λ_default = Δκ_λ_default)

#
### plot.




# This is the frequency range that we shall work with.
P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

## parameters that affect qs.
# A.d, A.κs_λ, A.κs_β
# A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets
q = uu->NMRSignalSimulator.evalclproxymixture(uu, mixture_params, Bs)

# create the functions for each resonance group.
B = Bs[1]
A = As[1]
qs = collect( collect( ωω->B.qs[i][k](ωω-B.ss_params.d[i], B.ss_params.κs_λ[i]) for k = 1:length(B.qs[i]) ) for i = 1:length(B.qs) )
q_singlets = ωω->NMRSignalSimulator.evalclsinglets(ωω, B.d_singlets, A.αs_singlets, A.Ωs_singlets, B.β_singlets, B.λ0, B.κs_λ_singlets)

# create the function for the entire compound.
q = uu->NMRSignalSimulator.evalclproxymixture(uu, As[1:1], Bs[1:1])

# evaluate at the plotting positions.
q_U = q.(U_rad)

qs_U = collect( collect( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
q_singlets_U = q_singlets.(U_rad)


#### sanity check.
q_check_U = q_singlets_U
if !isempty(qs) # some compounds only have singlets.
    q_check_U += sum( sum( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
end

discrepancy = norm(q_check_U- q_U)
println("sanity check. This should be numerically zero: ", discrepancy)



#### plot.


# reduce the plotting positions for low signal regions. Otherwise the plot store size will be too large, and the time to load the plot will be long.
display_reduction_factor = 100
display_threshold_factor =  α_relative_threshold/10

inds, _ = prunelowsignalentries(q_U, display_threshold_factor, display_reduction_factor)
P_display = P[inds]
U_display = U[inds]

canvas_size = (1000, 400)


plots_save_path = joinpath(save_folder, "$(molecule_names[1])_groups_real.html")
title_string = "Resonance groups, real part: $(molecule_names[1])"
plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P_display, U_display, q, qs, q_singlets, real, P[1]; canvas_size = canvas_size)
Plots.savefig(plot_obj, plots_save_path)
display(plot_obj)

include("radius_strategies.jl")
