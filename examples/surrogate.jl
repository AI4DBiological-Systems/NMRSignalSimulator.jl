include("a.jl")


import Random
Random.seed!(25)

include("./helpers/data.jl")

PythonPlot.close("all")
fig_num = 1

T = Float32
#T = Float64

### user inputs.
fs, SW, ν_0ppm = HAM.getpresetspectrometer(T, "700")

## pull the sample coupling values into dictionary structures.

root_data_path = getdatapath() # coupling values data repository root path

H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.

molecule_mapping_root_path = joinpath(
    root_data_path,
    "molecule_name_mapping",
)
molecule_mapping_file_path = joinpath(
    molecule_mapping_root_path,
    "select_molecules.json",
)

config = HAM.SHConfig{T}(
    coherence_tol = convert(T, 0.01),
    relative_α_threshold = convert(T, 0.005),
    tol_radius_1D = convert(T, 0.1), # strictly between 0 and 1. The lower, the better the approximation, but would a larger partition (i.e. more resonance groups).
    nuc_factor = convert(T, 1.5),
)
unique_cs_atol = convert(T, 1e-6)

molecule_entries = [
    "Singlet - 0 ppm";
    "L-Serine";
    "alpha-D-Glucose";
    "beta-D-Glucose";
    "Ethanol";
    #"L-Methionine";     
    "L-Phenylalanine";
    #"L-Glutathione reduced";
    #"L-Glutathione oxidized";       
    #"Uridine";
    "L-Glutamine";
    #"L-Valine";
    "DSS";
]
### end user input.


w_oracle = ones(T, length(molecule_entries))
Phys, As, MSPs = HAM.loadandsimulate(
    T,
    "700",
    molecule_entries,
    H_params_path,
    molecule_mapping_file_path;
    config = config,
    unique_cs_atol = unique_cs_atol
)


###
proxy_config = SIG.CLSurrogateConfig{T}(
    λ0 = convert(T, 3.4),
    Δr = convert(T, 1.0), # radial frequency resolution: smaller means slower to build surrogate, but more accurate.
    Δκ_λ = convert(T, 0.05), # T2 multiplier resolution. smaller means slower to build surrogate, but more accurate.
    Δcs_max_scalar = convert(T, 0.2), # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.
    κ_λ_lb = convert(T, 0.5), # lower limit for κ_λ for which the surrogate is made from.
    κ_λ_ub = convert(T, 2.5), # upper limit for κ_λ for which the surrogate is made from.
    ppm_padding = convert(T , 0.5),
)

Bs, MSS, itp_samps = SIG.fitclproxies(As, proxy_config)
# Bs and MSS are linked. Modification to one of its fields will affect the other.



# # Modify the phase of a spin system. Fill it with random numbers.

# choose this spin system:
n_select = 1
sys_select = 1

# create model parameters.
model_params = SIG.MixtureModelParameters(MSS; w = copy(w_oracle))

# create random phase parameters.
x_oracle = copy(model_params.var_flat)

st = model_params.systems_mapping.phase.st[n_select][sys_select]
fin = model_params.systems_mapping.phase.fin[n_select][sys_select]
x_oracle[st:fin] = randn(T, fin-st+1)

# load into model parameters data structure, then update model. model_params is linked to the surrogate model Bs and MSS.
model_params.var_flat[:] = x_oracle
SIG.importmodel!(model_params)
###  end modification of phase.


### plot.

f = uu->SIG.evalclmixture(uu, As, Bs; w = w_oracle)

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

# test params.
ΩS_ppm = collect( hz2ppmfunc.( SIG.combinevectors(A.Ωs) ./ SIG.twopi(T) ) for A in As )
ΩS_ppm_flat = SIG.combinevectors(ΩS_ppm)
P_max = maximum(ΩS_ppm_flat) + convert(T, 0.5)
P_min = minimum(ΩS_ppm_flat) - convert(T, 0.5)

P = LinRange(P_min, P_max, 80000)
#P = LinRange(-0.2, 5.5, 80000)
U = ppm2hzfunc.(P)
U_rad = U .* SIG.twopi(T)

## parameters that affect qs.
q = uu->SIG.evalclproxymixture(uu, As, Bs; w = w_oracle)

f_U = f.(U_rad)
q_U = q.(U_rad)

# for downstream examples.
S_U = copy(q_U)

discrepancy = norm(f_U-q_U)
max_val, ind = findmax(discrepancy)
println("relative l-2 discrepancy = ", norm(discrepancy)/norm(f_U))
println("max l-2 discrepancy: ", abs(f_U[ind]-q_U[ind])/abs(f_U[ind]), ", at $(P[ind]) ppm.")
println()

q_U_display = q_U
f_U_display = f_U
P_display = P

## visualize.
PythonPlot.figure(fig_num)
fig_num += 1

PythonPlot.plot(P_display, real.(f_U_display), label = "f")
PythonPlot.plot(P_display, real.(q_U_display), label = "q")
PythonPlot.plot(P_display, real.(q_U_display), "x")

PythonPlot.legend()
PythonPlot.xlabel("ppm")
PythonPlot.ylabel("real")
PythonPlot.gca().invert_xaxis()
PythonPlot.title("complex Lorentzian model (f) vs surrogate (q)")


nothing
