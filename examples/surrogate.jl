# # Load dependencies
using LinearAlgebra
using FFTW
using Printf

import Random
Random.seed!(25)

import PythonPlot as PLT # hide
import PlotlyLight as PLY
#import PlotlyDocumenter
import Cobweb
plot_counter = 0

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
    #"alpha-D-Glucose";
    #"beta-D-Glucose";
    "Ethanol";
    "Singlet - 0 ppm";
    "Glycine";
    "L-Serine";
    "L-Histidine";
    "DSS"; # this simulates the 0 ppm frequency reference compound.
    "Singlet - 4.9 ppm"; # this simulates D2O, the solvent.
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


# # Modify the phase of a spin system. Fill it with random numbers.

# choose this spin system:
n_select = 3
sys_select = 1
println("We're modifying $(molecule_entries[n_select]), spin system $sys_select.")

# create model parameters.
model_params = SIG.MixtureModelParameters(MSS, copy(w_oracle));

# create random phase parameters.
x_oracle = copy(model_params.var_flat)

st = model_params.systems_mapping.phase.st[n_select][sys_select]
fin = model_params.systems_mapping.phase.fin[n_select][sys_select]
x_oracle[st:fin] = randn(T, fin-st+1);

# make the solvent heavy-tailed.
x_oracle[end] = 2.4;

# flip phase for solvent.
st = model_params.systems_mapping.phase.st[end][1]
fin = model_params.systems_mapping.phase.fin[end][1]
x_oracle[st:fin] = ones(T, 1) .* -π/2;

# load into model parameters data structure, then update model. model_params is linked to the surrogate model Bs and MSS.
#model_params.var_flat[:] = x_oracle
SIG.importmodel!(model_params, x_oracle);



# # Visualize

f = uu->SIG.evalclmixture(uu, As, Bs; w = w_oracle)

# These are the conversion formulae for ppm frequency and Hz frequency.
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

# let's see how much discrepancy there is between the target mixture of Lorentzian model and the surrogate.
ΩS_ppm = collect( hz2ppmfunc.( SIG.combinevectors(A.Ωs) ./ SIG.twopi(T) ) for A in As )
ΩS_ppm_flat = SIG.combinevectors(ΩS_ppm)
P_max = maximum(ΩS_ppm_flat) + convert(T, 0.5)
P_min = minimum(ΩS_ppm_flat) - convert(T, 0.5)

P = LinRange(P_min, P_max, 80000)
U = ppm2hzfunc.(P)
U_rad = U .* SIG.twopi(T);

q = uu->SIG.evalclproxymixture(uu, Bs; w = w_oracle)

f_U = f.(U_rad)
q_U = q.(U_rad);
S_U = copy(q_U) #hide

discrepancies = abs.(f_U-q_U)
max_val, ind = findmax(discrepancies)
println("Spline approximation for surrogate model:")
println("relative l-2 discrepancy = ", norm(f_U-q_U)/norm(f_U))
println("max l-2 discrepancy: ", max_val, ", at ppm $(P[ind]).")
println();


# # Visualize.
# First, visualize the surrogate `q` against its fit target, `f`.
q_U_display = q_U
f_U_display = f_U
P_display = P

PLT.figure(fig_num) #hide
fig_num += 1 #hide
PLT.plot(P_display, real.(f_U_display), label = "f") #hide
PLT.plot(P_display, real.(q_U_display), "--", label = "q") #hide
PLT.legend() #hide
PLT.xlabel("ppm") #hide
PLT.ylabel("real") #hide
PLT.gca().invert_xaxis() #hide
PLT.title("complex Lorentzian model (f) vs surrogate (q)") #hide
PLT.gcf() #hide

# test serialization. #hide

using Serialization #hide

blah = SIG.exportmixtureproxy(Bs) #hide
serialize("tmp", blah) #hide
ss_params_set2, op_range_set2, λ02 = deserialize("tmp") #hide
Bs2, MSS2 = SIG.recoverclproxies(itp_samps, ss_params_set2, op_range_set2, As, λ02) #hide

serialize("tmp", Bs) #hide
Bs2 = deserialize("tmp") #hide
q2 = uu->SIG.evalclproxymixture(uu, Bs2; w = w_oracle) #hide
q2_U = q.(U_rad) #hide

@assert norm(q_U - q2_U) < 1e-12 #hide
rm("tmp") #hide

nothing
