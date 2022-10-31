
# loop version of plot_groups.jl
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

#save_folder = "/home/roy/MEGAsync/outputs/NMR/groups/700MHz"
#SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_molecules_SH_configs_reduce.json"

save_folder = "/home/roy/MEGAsync/outputs/NMR/groups/700MHz_low_intensity_threshold"
SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_molecules_SH_configs_low_intensity_threshold.json"

isdir(save_folder) || mkpath(save_folder)

λ0 = 3.4
unique_cs_atol = 1e-6
prune_Δc_option = 5
u_offset = 0.2

surrogate_config_path = "/home/roy/Documents/repo/NMRData/input/surrogate_configs/select_molecules_SH_configs.json"

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665

# path to the json file that provides the mapping from a molecule name to its spin system info file name.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_molecule_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/molecule_mapping/select_molecules.json")

### end inputs.

function plotgroupsbulk(molecule_entries::Vector{String},
        dict_molecule_to_filename,
        H_params_path, SH_config_path,
        fs, SW, ν_0ppm;
        unique_cs_atol = 1e-6,
        prune_Δc_option = 5,
        u_offset = 0.2,
        display_reduction_factor = 100,
        display_threshold_factor =  0.01/10,
        canvas_size = (1000, 400) )

    for n = 1:length(molecule_entries)
        plotgroups(molecule_entries[n],
            dict_molecule_to_filename,
            H_params_path, SH_config_path,
            fs, SW, ν_0ppm;
            unique_cs_atol = unique_cs_atol,
            prune_Δc_option = prune_Δc_option,
            u_offset = u_offset,
            display_reduction_factor = display_reduction_factor,
            display_threshold_factor = display_threshold_factor,
            canvas_size = canvas_size)
    end
end


function plotgroups(name::String,
    dict_molecule_to_filename,
    H_params_path, SH_config_path,
    fs, SW, ν_0ppm;
    unique_cs_atol = 1e-6,
    prune_Δc_option = 5,
    u_offset = 0.2,
    display_reduction_factor = 100,
    display_threshold_factor =  0.01/10,
    canvas_size = (1000, 400) )

    println("Working on $(name)")
    hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
    ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

    molecule_entries = [name;]
    Phys = NMRHamiltonian.getphysicalparameters(molecule_entries,
        H_params_path,
        dict_molecule_to_filename;
        unique_cs_atol = unique_cs_atol)

    mixture_params = NMRHamiltonian.setupmixtureSH(molecule_entries,
        fs, SW, ν_0ppm,
        Phys;
        config_path = SH_config_path,
        prune_Δc_option = prune_Δc_option)
    As = mixture_params

    dummy_SSFID = NMRSignalSimulator.SpinSysParamsType1(0.0)

    ## frequency locations. For plotting.
    ΩS_ppm = getPsnospininfo(mixture_params, hz2ppmfunc)
    ΩS_ppm_sorted = sort(combinevectors(ΩS_ppm))

    u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
    u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)

    Bs = NMRSignalSimulator.fitclproxies(As, dummy_SSFID, λ0;
        names = molecule_entries,
        config_path = surrogate_config_path,
        u_min = u_min,
        u_max = u_max)

    # This is the frequency range that we shall work with.
    P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
    U = ppm2hzfunc.(P)
    U_rad = U .* (2*π)

    q = uu->NMRSignalSimulator.evalclproxymixture(uu, mixture_params, Bs)

    # create the functions for each resonance group.
    B = Bs[1]
    A = As[1]
    qs = collect( collect( ωω->B.qs[i][k](ωω-B.ss_params.d[i], B.ss_params.κs_λ[i]) for k = 1:length(B.qs[i]) ) for i = 1:length(B.qs) )
    q_singlets = ωω->NMRSignalSimulator.evalclsinglets(ωω, B.d_singlets, A.αs_singlets, A.Ωs_singlets, B.β_singlets, B.λ0, B.κs_λ_singlets)

    # create the function for the entire molecule.
    q = uu->NMRSignalSimulator.evalclproxymixture(uu, As[1:1], Bs[1:1])

    # evaluate at the plotting positions.
    q_U = q.(U_rad)

    qs_U = collect( collect( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
    q_singlets_U = q_singlets.(U_rad)


    #### sanity check.
    q_check_U = q_singlets_U
    if !isempty(qs) # some molecules only have singlets.
        q_check_U += sum( sum( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
    end

    discrepancy = norm(q_check_U- q_U)
    println("sanity check. This should be numerically zero: ", discrepancy)

    # reduce the plotting positions for low signal regions. Otherwise the plot store size will be too large, and the time to load the plot will be long.

    inds, _ = prunelowsignalentries(q_U, display_threshold_factor, display_reduction_factor)
    P_display = P[inds]
    U_display = U[inds]

    save_molecule_name = replace("$(molecule_entries[1])", ","=>"-", " "=>"-")
    plots_save_path = joinpath(save_folder, "$(save_molecule_name)_groups_real.html")
    title_string = "Resonance groups, real part: $(molecule_entries[1])"
    plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P_display, U_display, q, qs, q_singlets, real, P[1]; canvas_size = canvas_size)
    Plots.savefig(plot_obj, plots_save_path)

    save_molecule_name = replace("$(molecule_entries[1])", ","=>"-", " "=>"-")
    println("name = ", save_molecule_name)

    if length(A.N_spins_sys) > 0
        println("Number of non-singlet spin systems: ", length(A.N_spins_sys))
        println("Resonance group sizes for each system: ", collect(collect(length(qs[i]) for i = 1:length(qs) ) ))
        println("Number of spins for each system: ", A.N_spins_sys)
        println("Number of resonance components for each system: ", length.(A.αs))
        println("Number of resonance components in each resonance group (inner index), for each system (outer index): ", collect( length.(A.part_inds_molecule[i]) for i = 1:length(A.part_inds_molecule)))
        println("summed intensity of each system: ", sum.(A.αs))
    end

    if length(A.αs_singlets) > 0
        println("intensity of singlets: ", A.αs_singlets)
    end
    println()

end

# molecule_entries = ["beta-Alanine";
# "L-Alanine"; ]
molecule_entries = collect( key for (key,val) in dict_molecule_to_filename)

println("Timing: plotgroupsbulk()")
@time plotgroupsbulk(molecule_entries,
        dict_molecule_to_filename,
        H_params_path, SH_config_path,
        fs, SW, ν_0ppm;
        unique_cs_atol = 1e-6,
        prune_Δc_option = 5,
        u_offset = 0.2,
        display_reduction_factor = 100,
        display_threshold_factor =  0.01/10,
        canvas_size = (1000, 400) )
