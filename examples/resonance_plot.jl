
# run a.jl first.

include("./helpers/data.jl")
include("./helpers/utils.jl")
include("./helpers/config.jl")
include("./helpers/SH.jl")


### user inputs.

#molecule_entries = ["L-Methionine"; "L-Phenylalanine"; "DSS"; "Ethanol"; "L-Isoleucine"]
#molecule_entries = ["alpha-D-Glucose"; "beta-D-Glucose"; "DSS"; "D2O"]

molecule_entries = ["beta-D-Glucose";]
molecule_entries = ["L-Glutathione oxidized";]
molecule_entries = ["L-Glutathione reduced";]

molecule_entries = ["alpha-D-Glucose";]

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665

# # machine values for the BMRB 500 MHz glucose experiment.
# ν_0ppm = 6752.490995937095
# SW = 16.0196917451925
# fs = 9615.38461538462


u_offset = 0.2 #in units ppm.
Δcs_padding = 0.02 #in units ppm.
min_window_cs = 0.06 #in units ppm.
show_sum = true

λ0 = 4.0


# # Specify physical dimension

save_folder_path = "./plots/resonance_groups" # makes path if doesn't exist.

# Suppose we have a width constraint of 8.3 cm.
max_width_inches = 8.3 / 2.54

# Choose the width we want to use, and do a sanity check.
width_inches = 8 / 2.54
@assert width_inches < max_width_inches

# The aspect ratio of the entire figure.
aspect_fig = 0.54

size_inches = (width_inches, aspect_fig*width_inches)

# options are real, imag, abs for real, imaginary, magnitude spectrum, respectively.
postfunc = real

# assemble.
# see the struct `ResonancePlotConfigType` for more fields in ./helpers/config.jl
config = ResonancePlotConfigType(
    save_folder_path = save_folder_path,
    size_inches = size_inches,
    postfunc = real,
)
### end inputs.

### set up.
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)


# slow.
As, Rs = runSH(molecule_entries)

### end set up.



dummy_SSFID = NMRSignalSimulator.SpinSysParamsType1(0.0)

## frequency locations. For plotting.
ΩS_ppm = getPsnospininfo(As, hz2ppmfunc)
ΩS_ppm_sorted = sort(combinevectors(ΩS_ppm))


u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)

Bs = NMRSignalSimulator.fitclproxies(As, dummy_SSFID, λ0;
names = molecule_entries,
#config_path = surrogate_config_path,
u_min = u_min,
u_max = u_max)

# This is the frequency range that we shall work with.
P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

## parameters that affect qs.
# A.d, A.κs_λ, A.κs_β
# A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets


# create the functions for each resonance group.
q, qs, q_singlets = setupresonanceplotfuncs(As, Bs)

# create the function for the entire molecule.
#q = uu->NMRSignalSimulator.evalclproxymixture(uu, As[1:1], Bs[1:1])

# evaluate at the plotting positions.
q_U = q.(U_rad)

qs_U = collect( collect( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
q_singlets_U = q_singlets.(U_rad)

### TODO make qs, and use ppm as input.

#### sanity check.
q_check_U = q_singlets_U
if !isempty(qs) # some molecules only have singlets.
    q_check_U += sum( sum( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
end

discrepancy = norm(q_check_U- q_U)
@assert discrepancy < 1e-14

## get intervals.
ΩS0 = getΩS(As)
ΩS0_ppm = getPs(ΩS0, hz2ppmfunc)

Δsys_cs = initializeΔsyscs(As, Δcs_padding)
exp_info = NMRSpecifyRegions.setupexperimentresults(molecule_entries, ΩS0_ppm, Δsys_cs; min_dist = Δcs_padding)
plot_lbs, plot_ubs = setupplotintervals(exp_info; min_window_cs = min_window_cs)

intervals = collect( (plot_lbs[i], plot_ubs[i]) for i in eachindex(plot_lbs) )

######  plot.
#spin_sys_select = 1

# save name. Does not check if there is already a file with this name to avoid overwrite!
formated_name = replace(
    molecule_entries[1],
    " "=> "",
    ":" => "-",
    ","=>"-",
    "'"=>"comma"
)
save_name = "$(formated_name).svg"


# title.
spectrometer_freq = round(fs/SW, digits = 2)
"$(molecule_entries[1]) at $(spectrometer_freq) MHz"


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


MakiePlots.plotmultiinterval1D(
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


# @assert 44==3


# NMR1Dbrokenfrequencyaxisplot(P, ys, legend_labels, plot_lbs, plot_ubs;
#     fig_size = fig_size,
#     fig_title_string = "$(molecule_entries[1]) at $(spectrometer_freq) MHz",
#     plot_title_string = "$(molecule_entries[1])\n$(spectrometer_freq) MHz",
#     subplot_spacing = subplot_spacing,
#     line_styles = line_styles,
#     legend_space_per_row = legend_space_per_row,
#     N_cols_legend = N_cols_legend,
#     legend_font = legend_font,
#     x_tick_font_size = x_tick_font_size,
#     line_width = line_width)