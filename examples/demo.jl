# # Load dependencies
using LinearAlgebra
using FFTW
using Printf

import Random
Random.seed!(25)

import PythonPlot as PLT #hide
import PlotlyLight as PLY
import DisplayAs

import PublicationDatasets as DS
import NMRDataSetup as DSU

import Statistics

import NMRSignalSimulator
SIG = NMRSignalSimulator
NMRHamiltonian = NMRSignalSimulator.NMRHamiltonian
JSON3 = NMRHamiltonian.JSON3
HAM = NMRHamiltonian

PLT.close("all")
fig_num = 1

#T = Float32 #hide
T = Float64;

# # User inputs for simulation
# Let's use a preset from a BMRB experiment at 700 MHz.
fs, SW, ν_0ppm = HAM.getpresetspectrometer(T, "700");

# Specify a 1/T2 inverse time constant.
λ0 = convert(T, 3.4);

# put the database coupling values into dictionary structures. You can supply your own coupling values; see the documentation website for NMRSignalSimulator.jl and NMRHamiltonian.jl.

root_data_path = DS.getdatapath(DS.NMR2023()) # coupling values data repository root path

H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.

molecule_mapping_root_path = joinpath(
    root_data_path,
    "molecule_name_mapping",
)
molecule_mapping_file_path = joinpath(
    molecule_mapping_root_path,
    "select_molecules.json",
);

# list of target compounds for simulation and surrogate construction.
molecule_entries = [
    "L-Glutamine";
    #"Gamma-Aminobutyric acid";
    #"L-Histidine";
];

# make up a set of relative concentration.
w_oracle = rand(T, length(molecule_entries))
w_oracle[end] = convert(T, 16.0); # make solvent very large.

# # Spin Hamiltonian simulation and resonance group computation
# generate the spin Hamiltonian simulation and cluster to get resonance groups. See NMRHamiltonian.jl documentation for further details for the following code.
config = HAM.SHConfig{T}(
    coherence_tol = convert(T, 0.01),
    relative_α_threshold = convert(T, 0.001),
    max_deviation_from_mean = convert(T, 0.05),
    acceptance_factor = convert(T, 0.99),
    total_α_threshold = zero(T),
)
unique_cs_digits = 6

Phys, As, MSPs = HAM.loadandsimulate(
    fs, SW, ν_0ppm,
    molecule_entries,
    H_params_path,
    molecule_mapping_file_path,
    config;
    unique_cs_digits = unique_cs_digits,
);

# # Frequency-domain surrogate construction

proxy_config = SIG.CLSurrogateConfig{T}(
    Δr = convert(T, 1.0), # radial frequency resolution: smaller means slower to build surrogate, but more accurate.
    Δκ_λ = convert(T, 0.05), # T2 multiplier resolution. smaller means slower to build surrogate, but more accurate.
    Δcs_max_scalar = convert(T, 0.2), # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.
    κ_λ_lb = convert(T, 0.5), # lower limit for κ_λ for which the surrogate is made from.
    κ_λ_ub = convert(T, 2.5), # upper limit for κ_λ for which the surrogate is made from.
    ppm_padding = convert(T , 0.5),
)

Bs, MSS, itp_samps = SIG.fitclproxies(As, λ0, proxy_config);
# Bs and MSS are linked. Modification to one of its fields will affect the other.

# # Visualize

# These are the conversion formulae for ppm frequency and Hz frequency.
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

# find a good display bound that contains all the compounds in our mixture. In this case, there is only L-Glutamine.
ΩS_ppm = collect( hz2ppmfunc.( SIG.combinevectors(A.Ωs) ./ SIG.twopi(T) ) for A in As )
ΩS_ppm_flat = SIG.combinevectors(ΩS_ppm)
P_max = maximum(ΩS_ppm_flat) + convert(T, 0.5)
P_min = minimum(ΩS_ppm_flat) - convert(T, 0.5)

P = LinRange(P_min, P_max, 80000)
U = ppm2hzfunc.(P)
U_rad = U .* SIG.twopi(T);

## parameters that affect qs.
q = uu->SIG.evalclproxymixture(uu, Bs; w = w_oracle)
q_U = q.(U_rad);

# Let's visualize the resonance groups of the only compound in our simulated mixture.
B = Bs[begin];


# Evaluate the plotting radial frequencies in `U_rad` for each resonance group surrogate.
qs = collect(
    collect(
        ww->SIG.evalclrg(ww, i, k, B)
        for k = 1:length(B.qs[i])
    )
    for i = 1:length(B.qs)
)
qs_U = collect( collect( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i in eachindex(qs) );

# The first index of `qs_U` is over the number of spin systems. Since alpha-D-Glucose has only one spin system, we set this index to 1. We only show the real part of the spectrum for this demo.
rgs_data = collect( real.(qs_U[begin][k]) for k in eachindex(qs_U[begin]))
sum_rg = sum(rgs_data);

# ## Visualize the 2.05 ppm to 2.2 ppm region of L-glutamine.
inds = findall(xx->(2.05<xx<2.2), P)
x_plot = P[inds]

y_plots = Vector{Vector{T}}(undef, length(rgs_data)+1)
y_plots[begin] = sum_rg[inds]
y_plots[begin+1:end] = collect( rgs_data[k][inds] for k in eachindex(rgs_data) )

line_styles = ["dash"; repeat(["solid";], length(rgs_data))]

name_strings = ["Summed signal";] # labels for each plotted line.
append!(name_strings, collect("Group $i" for i = 1:length(rgs_data)))

spec_freq_rounded = @sprintf("%.2f", fs/SW)
title_string = "Resonance groups of L-Glutamine at $spec_freq_rounded"
y_label = "real part of spectrum"
x_label = "ppm"
x_autorange = "reversed" # The convention for ppm in NMR literature is to reverse the axis.


plot_data_set = collect(
    PLY.Config(
        x = x_plot,
        y = y_plots[m],
        name = name_strings[m],
        type = "line-chart",
        line = PLY.Config(dash = line_styles[m]),
    ) for m in eachindex(y_plots)
)
plot_data = vcat(plot_data_set...)

layout = PLY.Config(
    #type = "line-chart",
    title = title_string,
    height = 400,
    width = 800,
    margin = 10,
    font = PLY.Config(
        family = "Lato",
        size = 16,
    ),
    xaxis = PLY.Config(
        title = x_label,
        autorange = x_autorange,
        titlefont = PLY.Config(
            color = "black",
            size = 12,
        )
    ),
    yaxis = PLY.Config(
        title = y_label,        
        titlefont = PLY.Config(
            color = "black",
            size = 12,
        )
    ),
)
ph = PLY.Plot(plot_data, layout);
#display(ph) # shows up in browser. #hide

import Cobweb #hide
Cobweb.save(
    Cobweb.Page(ph),
    #"../docs/src/generated/plot_1.html",
    "plot_1.html", # Literate.jl operates at the destination's folder. this'll show up in the generated/ folder.
);#hide
# You can view the rendered plot [here](plot_1.html)


# ## Visualize the 2.4 ppm to 2.5 ppm region of L-glutamine.
inds = findall(xx->(2.4<xx<2.5), P)
x_plot = P[inds]

y_plots = Vector{Vector{T}}(undef, length(rgs_data)+1)
y_plots[begin] = sum_rg[inds]
y_plots[begin+1:end] = collect( rgs_data[k][inds] for k in eachindex(rgs_data) )

line_styles = ["dash"; repeat(["solid";], length(rgs_data))]

name_strings = ["Summed signal";] # labels for each plotted line.
append!(name_strings, collect("Group $i" for i = 1:length(rgs_data)))

spec_freq_rounded = @sprintf("%.2f", fs/SW)
title_string = "Resonance groups of L-Glutamine at $spec_freq_rounded"
y_label = "real part of spectrum"
x_label = "ppm"
x_autorange = "reversed" # The convention for ppm in NMR literature is to reverse the axis.


plot_data_set = collect(
    PLY.Config(
        x = x_plot,
        y = y_plots[m],
        name = name_strings[m],
        type = "line-chart",
        line = PLY.Config(dash = line_styles[m]),
    ) for m in eachindex(y_plots)
)
plot_data = vcat(plot_data_set...)

layout = PLY.Config(
    #type = "line-chart",
    title = title_string,
    height = 400,
    width = 800,
    margin = 10,
    font = PLY.Config(
        family = "Lato",
        size = 16,
    ),
    xaxis = PLY.Config(
        title = x_label,
        autorange = x_autorange,
        titlefont = PLY.Config(
            color = "black",
            size = 12,
        )
    ),
    yaxis = PLY.Config(
        title = y_label,        
        titlefont = PLY.Config(
            color = "black",
            size = 12,
        )
    ),
)
ph = PLY.Plot(plot_data, layout);
display(ph) # shows up in browser. #hide

import Cobweb #hide
Cobweb.save(
    Cobweb.Page(ph),
    #"../docs/src/generated/plot_1.html",
    "plot_2.html", # Literate.jl operates at the destination's folder. this'll show up in the generated/ folder.
); #hide
# You can view the rendered plot [here](plot_2.html)



# ## Visualize the 3.72 ppm to 3.79 ppm region of L-glutamine.
inds = findall(xx->(3.72<xx<3.79), P)
x_plot = P[inds]

y_plots = Vector{Vector{T}}(undef, length(rgs_data)+1)
y_plots[begin] = sum_rg[inds]
y_plots[begin+1:end] = collect( rgs_data[k][inds] for k in eachindex(rgs_data) )

line_styles = ["dash"; repeat(["solid";], length(rgs_data))]

name_strings = ["Summed signal";] # labels for each plotted line.
append!(name_strings, collect("Group $i" for i = 1:length(rgs_data)))

spec_freq_rounded = @sprintf("%.2f", fs/SW)
title_string = "Resonance groups of L-Glutamine at $spec_freq_rounded"
y_label = "real part of spectrum"
x_label = "ppm"
x_autorange = "reversed" # The convention for ppm in NMR literature is to reverse the axis.


plot_data_set = collect(
    PLY.Config(
        x = x_plot,
        y = y_plots[m],
        name = name_strings[m],
        type = "line-chart",
        line = PLY.Config(dash = line_styles[m]),
    ) for m in eachindex(y_plots)
)
plot_data = vcat(plot_data_set...)

layout = PLY.Config(
    #type = "line-chart",
    title = title_string,
    height = 400,
    width = 800,
    margin = 10,
    font = PLY.Config(
        family = "Lato",
        size = 16,
    ),
    xaxis = PLY.Config(
        title = x_label,
        autorange = x_autorange,
        titlefont = PLY.Config(
            color = "black",
            size = 12,
        )
    ),
    yaxis = PLY.Config(
        title = y_label,        
        titlefont = PLY.Config(
            color = "black",
            size = 12,
        )
    ),
)
ph = PLY.Plot(plot_data, layout);
#display(ph) # shows up in browser. #hide

import Cobweb #hide
Cobweb.save(
    Cobweb.Page(ph),
    #"../docs/src/generated/plot_1.html",
    "plot_3.html", # Literate.jl operates at the destination's folder. this'll show up in the generated/ folder.
); #hide
# You can view the rendered plot [here](plot_3.html)

