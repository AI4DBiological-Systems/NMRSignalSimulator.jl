
import NMRHamiltonian
import NMRSignalSimulator
import NMRSpecifyRegions

using DataDeps
import Tar

using LinearAlgebra
import MakiePlots

using Parameters

include("./helpers/data.jl")
include("./helpers/utils.jl")
include("./helpers/config.jl")
include("./helpers/SH.jl")
include("./helpers/plot.jl")

# # User inputs.

root_data_path = getdatapath() # coupling values data repository root path

H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.

molecule_mapping_root_path = joinpath(root_data_path, "molecule_name_mapping")
molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "select_molecules.json")
#molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "GISSMO_names.json")


#=
Some example choices:
```
molecule_entries = ["D2O";]
molecule_entries = ["Glycerol";]
molecule_entries = ["Glycine";]
molecule_entries = ["DSS";]
molecule_entries = ["L-Glutathione oxidized";]
molecule_entries = ["L-Glutathione reduced";]
molecule_entries = ["beta-D-Glucose";]
molecule_entries = ["alpha-D-Glucose";]
molecule_entries = ["L-Glutamine";]
molecule_entries = ["Gamma-Aminobutyric acid";]
molecule_entries = ["Agmatine";]
molecule_entries = ["L-Phenylalanine";]
molecule_entries = ["L-Leucine";]
molecule_entries = ["L-Valine";]
molecule_entries = ["L-Isoleucine";]
molecule_entries = ["L-Histidine";]
```
See joinpath(getdatapath(), "molecule_name_mapping/select_molecules.json") for some more valid molecule entries.
=#
molecule_entries = ["L-Glutamine";]

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665

#=
Another example: achine values for the BMRB 500 MHz glucose experiment.
```
ν_0ppm = 6752.490995937095
SW = 16.0196917451925
fs = 9615.38461538462
```
=#

# I am here. make this run, then go to model.

# These are in units ppm.
u_offset = 0.2 # affects the frequency range to plot.
Δcs_padding = 0.02 # affects the interval creation.
min_window_cs = 0.06 # minimum interval length.

# Show each resonance group and singlets, but also show their sum.
show_sum = true

# The spin Hamiltonian simulation does not simulate T2 decay. The user needs to specify it.
λ0 = 4.0

# # Plot settings
# ## Specify physical dimension.
# Suppose we have a width constraint of 8.3 cm.
max_width_inches = 8.3 / 2.54

# Choose the width we want to use, and do a sanity check.
width_inches = 8 / 2.54
@assert width_inches < max_width_inches

# The aspect ratio of the entire figure.
aspect_fig = 0.54

size_inches = (width_inches, aspect_fig*width_inches)

# ## Other settings
# This is a function to apply to the complex-valued function output before displaying. The options are real, imag, abs for real, imaginary, magnitude spectrum, respectively. These options are functions from the Julia standard library.
postfunc = real

# Change this to the desired storage path. The plotting function will create the path if it doesn't exist.
save_folder_path = "./plots/resonance_groups" # makes path if doesn't exist.


# Assemble the configuration specified in the user settings section at the beginning of this code. The keyword struct from parameters.jl is used so that unspecified field names have default values. See the struct `ResonancePlotConfigType` for the field names and their default values in ./helpers/config.jl, or type `?ResonancePlotConfigType` in the Julia REPL.
# Add more explicit field name initialization here as desired.
config = ResonancePlotConfigType(
    save_folder_path = save_folder_path,
    size_inches = size_inches,
    postfunc = real,
)

# Currently, we're using a feature of NMRHamiltonian where we do an iterative search over γ Stopping condition for partitioning algorithm
max_partition_size_offset = 9 #2

# # Set up
# Get conversion functions.
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

# Get spin Hamiltonian simulation. This could be slow.
As, Rs = runSH(
    H_params_path,
    molecule_mapping_file_path,
    fs,
    SW,
    ν_0ppm,
    molecule_entries,
    max_partition_size_offset;
    #search_θ = true,
    #θ_default = 0.0,
    starting_manual_knn = 60,
    γ_base = 0.1,
    #γ_rate = 1.05,
    max_iters_γ = 100,
    #min_dynamic_range = 0.95,
    cc_gap_tol = 1e-8,
    cc_max_iters = 300,
    assignment_zero_tol = 1e-3,
)

# Fit surrogate model. Since we are not fitting against data, a `SharedShift` surrogate would suffice.
type_SSParams = NMRSignalSimulator.getSpinSysParamsdatatype(NMRSignalSimulator.SharedShift{Float64})

# Get the frequency range over which the surrogate acts on.
ΩS_ppm = getPsnospininfo(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(combinevectors(ΩS_ppm))

u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)

# Construct surrogate.
Bs = NMRSignalSimulator.fitclproxies(
    type_SSParams,
    As,
    λ0;
    names = molecule_entries,
    u_min = u_min,
    u_max = u_max
)

# This is the frequency range that we shall work with.
P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000) # ppm.
U = ppm2hzfunc.(P) # Hz.
U_rad = U .* (2*π) # radians.

# ## Simulated spectrum for the first molecule
# The following returns the simulated spectrum for the entire spectrum (i.e. summed groups and singlets) `q`, a nested array of all resonance group spectrums for all non-singlet spin systems in the molecule `qs`, and the spectrum of all singlets. Each spectrum is a function that takes in a frequency in radians, and outputs a complex number.
q, qs, q_singlets = setupresonanceplotfuncs(As, Bs)

# ### Optional: sanity checks.

# evaluate the summed spectrum at the plotting positions.
q_U = q.(U_rad)

# the singlets spectrum and the sum of the resonance groups' spectrum.
q_singlets_U = q_singlets.(U_rad)

aggregated_U = q_singlets_U
if !isempty(qs) # some molecules only have singlets.
    aggregated_U += sum( sum( qs[i][k].(U_rad) for k in eachindex(qs[i]) ) for i in eachindex(qs) )
end

# they should be the same.
discrepancy = norm(aggregated_U- q_U)
@assert discrepancy < 1e-14

#=
For new Julia users: Although agreegated_U was first assigned to the array q_singlets_U, it was reassigned in the `+=` line. This would not mutate (i.e. change) the contents of q_singlets_U. If we wanted to change the contents of q_singlets_U indirectly via working on aggregated_U, we'd do:
```
aggregated_U[:] += sum( sum( qs[i][k].(U_rad) for k in eachindex(qs[i]) ) for 
```
which tells Julia that we want aggregated_U to remain a reference/pointer-like variable to the array that is *currently* referenced by q_singlets_U, and assign the right-hand side to that array. We didn't use `[:]`, so Julia creates a new array on the right-hand side, and re-assigns aggregated_U to be another the reference/pointer-like variable that points to this new array, which does not modify the array that q_singlets_U references.

This is a large course of confusion for those who have issues viewing this as combining the dynamic assignment behavior of MATLAB with the reference to an array idea from C. The combination is actually very well-defined and becomes second-nature with practice. More details [here](https://m3g.github.io/JuliaNotes.jl/stable/assignment/).
=#

# ## get frequency intervals
# Find the minimum and maximum frequency, and the intervals within the range where the spectrum. An interval boundary point is created when there are no resonance components within `Δcs_padding` ppm units.

ΩS0 = getΩS(As)
ΩS0_ppm = getPs(ΩS0, hz2ppmfunc)

Δsys_cs = initializeΔsyscs(As, Δcs_padding)
exp_info = NMRSpecifyRegions.setupexperimentresults(molecule_entries, ΩS0_ppm, Δsys_cs; min_dist = Δcs_padding)
plot_lbs, plot_ubs = setupplotintervals(exp_info; min_window_cs = min_window_cs)

intervals = collect( (plot_lbs[i], plot_ubs[i]) for i in eachindex(plot_lbs) )

# # Save plot
# Plot resonance groups of all spin systems on a figure for the first molecule.

# title.
spectrometer_freq = round(fs/SW, digits = 2)
"$(molecule_entries[1]) at $(spectrometer_freq) MHz"

# save name. Does not check if there is already a file with this name to avoid overwrite!
formated_name = replace(
    molecule_entries[1],
    " "=> "",
    ":" => "-",
    ","=>"-",
    "'"=>"comma"
)
save_name = "$(formated_name)_$(round(Int,spectrometer_freq))_MHz.svg"



# `Qs` is the list of functions to plot.
# Line style options are `:solid`, `:dash`, `:dashdot`.

#g = pp->q(ppm2hzfunc(pp)*(2*π))
Qs, legend_labels, line_styles, plot_colours = setupresonanceplot(
    q,
    qs,
    q_singlets,
    ppm2hzfunc;
    show_sum = show_sum,
    show_singlets = !isempty(first(As).αs_singlets), # set false if there is no singlets.
    line_style_q = :dash, # line style for sum
    line_style_qs = :solid, # line style for each resonance group.
    line_style_q_singlets = :solid, # line style for all singlets.
    postfunc = config.postfunc,
)

# same width for all lines.
line_widths = ones(length(Qs)) .* config.line_width


fig = MakiePlots.plotmultiinterval1D(
    Qs,
    intervals;
    legend_placement = :bottom,
    legend_labels = legend_labels,
    save_name = save_name,
    color_list = plot_colours,
    line_style_list = line_styles,
    line_width_list = line_widths,
    title = "$(molecule_entries[1]), $(spectrometer_freq) MHz",

    title_font_size = config.title_font_size,
    width_padding_proportion = config.width_padding_proportion,
    reverse_x_axis = config.reverse_x_axis,
    grid_visible = config.grid_visible,

    x_label = config.x_label,
    x_tick_decimal_places = config.x_tick_decimal_places,
    x_label_font_size = config.x_label_font_size,
    x_tick_label_size = config.x_tick_label_size,

    y_label = config.y_label,
    y_tick_decimal_places = config.y_tick_decimal_places,
    y_label_font_size = config.y_label_font_size,
    y_tick_label_size = config.y_tick_label_size,

    use_Wilkinson_for_ticks = config.use_Wilkinson_for_ticks,

    legend_font_size = config.legend_font_size,
    legend_orientation = config.legend_orientation,
    legend_line_width = config.legend_line_width,

    save_folder_path = config.save_folder_path,
    size_inches = config.size_inches,
    pt_per_inch = config.pt_per_inch,
    font_size = config.misc_font_size, # default font size for text that isn't specified.
    )


