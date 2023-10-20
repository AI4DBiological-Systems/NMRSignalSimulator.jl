var documenterSearchIndex = {"docs":
[{"location":"api/#Function-References","page":"Public API","title":"Function References","text":"","category":"section"},{"location":"api/#Surrogate","page":"Public API","title":"Surrogate","text":"","category":"section"},{"location":"api/","page":"Public API","title":"Public API","text":"CLSurrogateConfig","category":"page"},{"location":"api/#NMRSignalSimulator.CLSurrogateConfig","page":"Public API","title":"NMRSignalSimulator.CLSurrogateConfig","text":"struct CLSurrogateConfig{T}     λ0::T = convert(T, 3.4)     Δr::T = convert(T, 1.0)     Δκλ::T = convert(T, 0.05)     Δcsmaxscalar::T = convert(T, 0.2)     κλlb::T = convert(T, 0.5)     κλ_ub::T = convert(T, 2.5)\n\ndefault_ppm_padding::T = convert(T , 0.5)\n\nend\n\ndefault_ppm_padding – \nκ_λ_lb, κ_λ_ub – the default lower and upper bounds, respectively, for the κ_λ input of the surrogate.\nΔcs_max_scalar, default_ppm_padding – specifies the surrogate's default frequency operating range, [umin - deltau, umax + deltau], where deltau is Δcsmax_scalar expressed in Hz.\n\numin, umax are in Hz, Δcsmaxscalar is in ppm. u_min is set to the smallest resonance frequency from As, and u_max is set to the largest resonance frequency from As. Δcs_max_scalar is in units of ppm. It is the border that is added to umin and umax (once they are converted to ppm) to get the final frequency interval for which the surrogate of a spin system is fitted to. The extrapolation of the surrogate for frequency queries outisde this interval is set to return zero, so the surrogate is the zero signal outside this interval for that spin system.\n\nΔr, Δκ_λ – the sampling increment for the frequency input r and T2 multiplier input κ_λ for generating samples to fit the surrogate. Smaller means the surrogate is more accurate, but slower to construct the surrogate.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"Public API","title":"Public API","text":"fitclproxies","category":"page"},{"location":"api/#NMRSignalSimulator.fitclproxies","page":"Public API","title":"NMRSignalSimulator.fitclproxies","text":"fitclproxies(\n    As::Vector{HAM.SHType{T}},\n    config::CLSurrogateConfig{T};\n    names::Vector{String} = Vector{String}(undef, 0),\n    ) where T <: AbstractFloat\n\nCreate surrogate for NMR spectrum, given the simulated resonance comopnent results from NMRHamiltonian.\n\nInputs:\n\nAs – a 1-D array of compound resonance simulations. See NMRHamiltonian.simulate.\nconfig – the configuration for building the surrogate. See CLSurrogateConfig.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"evalclmixture","category":"page"},{"location":"api/#NMRSignalSimulator.evalclmixture","page":"Public API","title":"NMRSignalSimulator.evalclmixture","text":"function evalclmixture(\n    u_rad,\n    As::Vector{HAM.SHType{T}},\n    Bs::Vector{MoleculeType{T,SST}};\n    w::Vector{T} = ones(T, length(As)),\n)::Complex{T} where {T <: Real, SST}\n\nEvaluates the preliminary model at radial frequency (in radians) u_rad. Inputs:\n\nBs is the surrogate model.\nAs is the resonance component data structure.\nw is relative concentration of each compound entry.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"evalclproxymixture","category":"page"},{"location":"api/#NMRSignalSimulator.evalclproxymixture","page":"Public API","title":"NMRSignalSimulator.evalclproxymixture","text":"function evalclproxymixture(\n    u_rad,\n    As::Vector{HAM.SHType{T}},\n    Bs::Vector{MoleculeType{T,SST}};\n    w::Vector{T} = ones(T, length(As)),\n)::Complex{T} where {T <: Real,SST}\n\nEvaluates the surrogate model at radial frequency (in radians) u_rad. Inputs:\n\nBs is the surrogate model.\nAs is the resonance component data structure.\nw is relative concentration of each compound entry.\n\n\n\n\n\n","category":"function"},{"location":"api/#Constructors","page":"Public API","title":"Constructors","text":"","category":"section"},{"location":"api/","page":"Public API","title":"Public API","text":"getSpinSysParamsdatatype","category":"page"},{"location":"api/#NMRSignalSimulator.getSpinSysParamsdatatype","page":"Public API","title":"NMRSignalSimulator.getSpinSysParamsdatatype","text":"function getSpinSysParamsdatatype(\n    ::Type{CoherenceShift{T}}\n    )::DataType where T <: AbstractFloat\n\nConvinence constructor for SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}\n\n\n\n\n\nfunction getSpinSysParamsdatatype(\n    ::Type{SharedShift{T}}\n    )::DataType where T <: AbstractFloat\n\nConvinence constructor for SpinSysParams{SharedShift{T}, CoherencePhase{T}, SharedT2{T}}\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"MixtureModelParameters","category":"page"},{"location":"api/#NMRSignalSimulator.MixtureModelParameters","page":"Public API","title":"NMRSignalSimulator.MixtureModelParameters","text":"function MixtureModelParameters(\n    MSS::MixtureSpinSys{T,ST,PT,T2T};\n    w = ones(T, getNentries(MS)),\n) where {T,ST,PT,T2T}\n\nConvinence constructor for MixtureModelParameters(). Does not create a copy of of the inputs.\n\n\n\n\n\n","category":"type"},{"location":"api/#Model-updates","page":"Public API","title":"Model updates","text":"","category":"section"},{"location":"api/","page":"Public API","title":"Public API","text":"importmodel!","category":"page"},{"location":"api/#NMRSignalSimulator.importmodel!","page":"Public API","title":"NMRSignalSimulator.importmodel!","text":"function importmodel!(p::MixtureModelParameters)::Nothing\n\nUpdates p.MSS with the contents of p.var_flat.\n\n\n\n\n\nfunction importmodel!(p::MixtureModelParameters, x, inds::SubsetVarsIndices)::Nothing\n\nOnly update the entries as specified by inds.\n\n\n\n\n\nfunction importmodel!(p::MixtureModelParameters, x, ::AllVars)::Nothing\n\nCalls importmodel!(p, x).\n\n\n\n\n\nfunction importmodel!(p::MixtureModelParameters, x)::Nothing\n\nUpdate the model p using the contents of x.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"importmodelreset!","category":"page"},{"location":"api/#NMRSignalSimulator.importmodelreset!","page":"Public API","title":"NMRSignalSimulator.importmodelreset!","text":"function importmodelreset!(p::MixtureModelParameters, x, inds::SubsetVarsIndices)::Nothing\n\nThis is a version of importmodel!(p::MixtureModelParameters, x, inds::SubsetVarsIndices)::Nothing that resets the other variables not mentioned in inds.\n\nFirst, resets the model to default parameter values; 0 for shift and phase parameters, 1 for T2 multiplier parameters. Then update the model parameters as specified by inds with the entries in x.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"exportmodel!","category":"page"},{"location":"api/#NMRSignalSimulator.exportmodel!","page":"Public API","title":"NMRSignalSimulator.exportmodel!","text":"function exportmodel!(p::MixtureModelParameters)::Nothing\n\nFlattens the model parameters into p.var_flat.\n\n\n\n\n\n","category":"function"},{"location":"api/#Subset-vars","page":"Public API","title":"Subset vars","text":"","category":"section"},{"location":"api/","page":"Public API","title":"Public API","text":"VariableSetTrait","category":"page"},{"location":"api/#NMRSignalSimulator.VariableSetTrait","page":"Public API","title":"NMRSignalSimulator.VariableSetTrait","text":"abstract type VariableSetTrait end\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"Public API","title":"Public API","text":"AllVars","category":"page"},{"location":"api/#NMRSignalSimulator.AllVars","page":"Public API","title":"NMRSignalSimulator.AllVars","text":"struct AllVars <: VariableSetTrait end\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"Public API","title":"Public API","text":"SubsetVars","category":"page"},{"location":"api/#NMRSignalSimulator.SubsetVars","page":"Public API","title":"NMRSignalSimulator.SubsetVars","text":"struct SubsetVars <: VariableSetTrait\n    indices::SubsetVarsIndices\n    active_systems::Vector{Tuple{Int,Int}}\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"Public API","title":"Public API","text":"SubsetVarsIndices","category":"page"},{"location":"api/#NMRSignalSimulator.SubsetVarsIndices","page":"Public API","title":"NMRSignalSimulator.SubsetVarsIndices","text":"struct SubsetVarsIndices\n    shift::Vector{Int}\n    phase::Vector{Int}\n    T2::Vector{Int}\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"Public API","title":"Public API","text":"exportindices","category":"page"},{"location":"api/#NMRSignalSimulator.exportindices","page":"Public API","title":"NMRSignalSimulator.exportindices","text":"function exportindices(A::SubsetVarsIndices)::Vector{Int}\n\nReturns [A.shift; A.phase; A.T2]\n\n\n\n\n\nfunction exportindices(A::SubsetVars)::Vector{Int}     Returns index in A.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"setupsubsetvars","category":"page"},{"location":"api/#NMRSignalSimulator.setupsubsetvars","page":"Public API","title":"NMRSignalSimulator.setupsubsetvars","text":"function setupsubsetvars(\n    indicies_input::Vector{Int},\n    mapping::ParamsMapping;\n)::SubsetVars\n\nConvinence constructor for SubsetVars.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"getindices","category":"page"},{"location":"api/#NMRSignalSimulator.getindices","page":"Public API","title":"NMRSignalSimulator.getindices","text":"function getindices(A::AllVars)::AllVars     returns A.\n\n\n\n\n\nfunction getindices(A::SubsetVars)::SubsetVarsIndices     returns A.indices.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"getactivesystems","category":"page"},{"location":"api/#NMRSignalSimulator.getactivesystems","page":"Public API","title":"NMRSignalSimulator.getactivesystems","text":"function getactivesystems(A::AllVars)::AllVars     returns A.\n\n\n\n\n\nfunction getactivesystems(A::SubsetVars)::Vector{Tuple{Int,Int}}     returns A.active_systems.\n\n\n\n\n\n","category":"function"},{"location":"api/#Parameter-range","page":"Public API","title":"Parameter range","text":"","category":"section"},{"location":"api/","page":"Public API","title":"Public API","text":"OperationRange","category":"page"},{"location":"api/#NMRSignalSimulator.OperationRange","page":"Public API","title":"NMRSignalSimulator.OperationRange","text":"struct OperationRange{T}\n    u_min::T # in Hz.\n    u_max::T # in Hz.\n    \n    # d_max = ppm2hzfunc.(Δcs_max) .- ppm2hzfunc(zero(T))\n    d_max::Vector{T} # in Hz. ζ is in rad.\n\n    κ_λ_lb::T\n    κ_λ_ub::T\n    Δr::T\n    Δκ_λ::T\nend\n\nOperation range for the spline surrogate q(r,λ) of a resonance group. Bounds on r for groups in the i-th spin system:\n\nr_min = 2*π*(u_min - d_max[i])\nr_max = 2*π*(u_max + d_max[i])\n\nLower and upper bounds on λ are κ_λ_lb and κ_λ_ub, respectively. Able to serialize/deserialize.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"Public API","title":"Public API","text":"fetchbounds","category":"page"},{"location":"api/#NMRSignalSimulator.fetchbounds","page":"Public API","title":"NMRSignalSimulator.fetchbounds","text":"function fetchbounds(\n    p::MixtureModelParameters,\n    Bs::Vector{MoleculeType{T, SST}};\n    shift_proportion = 0.9, # between 0 to 1. Control the returned shift bounds as a proportion of the maximum allowed shift bounds used when the surrogates were created.\n    phase_lb = convert(T, -π),\n    phase_ub = convert(T, π),\n)::Tuple{Vector{T},Vector{T}} where {T,SST}\n\nReturns a Vector for lower bound and Vector for upper bound for each parameter variable. The order (first elements to last elements) of the returned vectors are: Shift (ζ), phase (κβ), then T2 (κλ) parameters.\n\nInput:\n\np – model parameters. Used to determine the size of the output.\nBs – Surrogate model. This function uses the OperationRange field from each element of Bs.\n\nOptional inputs:\n\nshift_proportion – between 0 to 1. Control the returned shift bounds as a proportion of the maximum allowed shift bounds used when the surrogates were created.\nphase_lb – the fill value for the lower bound of phase parameter (κ_β).\nphase_ub – the fill value for the lower bound of phase parameter (κ_β).\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"MoleculeParamsMapping","category":"page"},{"location":"api/#NMRSignalSimulator.MoleculeParamsMapping","page":"Public API","title":"NMRSignalSimulator.MoleculeParamsMapping","text":"struct MoleculeParamsMapping\n    st::Vector{Vector{Int}}\n    fin::Vector{Vector{Int}}\nend\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"Public API","title":"Public API","text":"ParamsMapping","category":"page"},{"location":"api/#NMRSignalSimulator.ParamsMapping","page":"Public API","title":"NMRSignalSimulator.ParamsMapping","text":"struct ParamsMapping\n    shift::MoleculeParamsMapping\n    phase::MoleculeParamsMapping\n    T2::MoleculeParamsMapping\nend\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"Public API","title":"Public API","text":"getshiftrange","category":"page"},{"location":"api/#NMRSignalSimulator.getshiftrange","page":"Public API","title":"NMRSignalSimulator.getshiftrange","text":"function getshiftrange(A::ParamsMapping)\n\nReturns the range for the shift parameters (ζ).\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"getphaserange","category":"page"},{"location":"api/#NMRSignalSimulator.getphaserange","page":"Public API","title":"NMRSignalSimulator.getphaserange","text":"function getphaserange(A::ParamsMapping)\n\nReturns the range for the phase parameters (κ_β).\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"getT2range","category":"page"},{"location":"api/#NMRSignalSimulator.getT2range","page":"Public API","title":"NMRSignalSimulator.getT2range","text":"function getT2range(A::ParamsMapping)\n\nReturns the range for the T2 parameters (κ_λ).\n\n\n\n\n\n","category":"function"},{"location":"api/#Utilies","page":"Public API","title":"Utilies","text":"","category":"section"},{"location":"api/","page":"Public API","title":"Public API","text":"getNentries","category":"page"},{"location":"api/#NMRSignalSimulator.getNentries","page":"Public API","title":"NMRSignalSimulator.getNentries","text":"function getNentries(MSS::MixtureSpinSys)::Int\n\nReturns the number of molecule entries associated with MSS\n\n\n\n\n\nfunction getNentries(model_params::MixtureModelParameters)::Int\n\nReturn the number of molecule entries associated with model_params\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"exportphysicalparams","category":"page"},{"location":"api/#NMRSignalSimulator.exportphysicalparams","page":"Public API","title":"NMRSignalSimulator.exportphysicalparams","text":"exportphysicalparams(\n    MSS::MixtureSpinSys,\n    Phys::Vector{HAM.PhysicalParamsType{T}}, # template.\n    fs::T,\n    SW::T,\n    ν_0ppm::T,\n) where T\n\nOutputs: outshift, outphase, out_T2.\n\nThe fields cs_sys and cs_singlets for out_shift, out_phase, and out_T2 does not contain chemical shift. It contains the shift difference (in ppm), phase variable κ_β (in radians), and the T2 variable ξ (a multiplier, dimensionless), respectively.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"createtablecolumns","category":"page"},{"location":"api/#NMRSignalSimulator.createtablecolumns","page":"Public API","title":"NMRSignalSimulator.createtablecolumns","text":"function createtablecolumns(Phys::Vector{NMRHamiltonian.PhysicalParamsType{T}})::Tuple{Vector{Vector{Int}},Vector{Tuple{Int,Int}},Vector{T}} where T\n\nExtract information from Phys to produce three 1-D arrays that is useful for constructing tables for visualization purposes.\n\nOutputs:\n\nnucleus_ID_set – nucleus ID numbers, as recorded from Phys[n].H_IDs, i.e., the original not the re-labelled IDs.\nlocation_set – pairs of the format (n,i), n is the compound entry index, i is the spin system index.\ncs_set – chemical shift numbers for each nucleus, in units ppm.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"ζ2Δcs","category":"page"},{"location":"api/#NMRSignalSimulator.ζ2Δcs","page":"Public API","title":"NMRSignalSimulator.ζ2Δcs","text":"function ζ2Δcs(ζ::T, ν_0ppm::T, hz2ppmfunc)::T where T\n\nConvert radial frequency to ppm. '''\n\n# test.\n\na = Δcs2ζ(0.1, ppm2hzfunc) @show ζ2Δcs(a, ν_0ppm, hz2ppmfunc) '''\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"Δcs2ζ","category":"page"},{"location":"api/#NMRSignalSimulator.Δcs2ζ","page":"Public API","title":"NMRSignalSimulator.Δcs2ζ","text":"function Δcs2ζ(Δcs::T, ppm2hzfunc)::T where T\n\nConvert ppm to radial frequency. '''\n\n# test.\n\na = Δcs2ζ(0.1, ppm2hzfunc) @show ζ2Δcs(a, ν_0ppm, hz2ppmfunc) '''\n\n\n\n\n\n","category":"function"},{"location":"api/#Serialize/Deserialize","page":"Public API","title":"Serialize/Deserialize","text":"","category":"section"},{"location":"api/","page":"Public API","title":"Public API","text":"serializclproxies","category":"page"},{"location":"api/#NMRSignalSimulator.serializclproxies","page":"Public API","title":"NMRSignalSimulator.serializclproxies","text":"function serializclproxies(\n    Bs::Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T}}},\n) where {T,ST,PT,T2T}\n\nReturns a dictionary containing the surrogate model.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"serializitpsamples","category":"page"},{"location":"api/#NMRSignalSimulator.serializitpsamples","page":"Public API","title":"NMRSignalSimulator.serializitpsamples","text":"function serializitpsamples(\n    itp_samps_set::Vector{Vector{InterpolationSamples{T}}},\n) where T\n\nReturns a dictionary of the interpolation samples used to fit the surrogate.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"deserializclproxies","category":"page"},{"location":"api/#NMRSignalSimulator.deserializclproxies","page":"Public API","title":"NMRSignalSimulator.deserializclproxies","text":"function deserializclproxies(W)\n\nW is of a data type that can be addressed via a key (of type Symbol) and returns a blue. Examples: type Dict{Symbol, Any} or JSON.Object.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"deserializitpsamples","category":"page"},{"location":"api/#NMRSignalSimulator.deserializitpsamples","page":"Public API","title":"NMRSignalSimulator.deserializitpsamples","text":"function deserializitpsamples(W)\n\nW is of a data type that can be addressed via a key (of type Symbol). Examples: type Dict{Symbol, Any} or JSON.Object.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"Public API","title":"Public API","text":"recoverclproxies","category":"page"},{"location":"api/#NMRSignalSimulator.recoverclproxies","page":"Public API","title":"NMRSignalSimulator.recoverclproxies","text":"function recoverclproxies(\n    itp_samps_set::Vector{Vector{InterpolationSamples{T}}},\n    ss_params_set::Vector{SpinSysParams{ST,PT,T2T}},\n    op_range_set::Vector{OperationRange{T}},\n    As::Vector{HAM.SHType{T}},\n    λ0::T,\n    )::Tuple{Vector{MoleculeType{T,SpinSysParams{ST,PT,T2T}}},\n    MixtureSpinSys{T,ST,PT,T2T}} where {T,ST,PT,T2T}\n\nTakes in deserialized quantities to output a surrogate model and its corresponding mixture spin system data structure.\n\n\n\n\n\n","category":"function"},{"location":"demo_code/#Demo:-code-walk-through","page":"Demo: code walk-through","title":"Demo: code walk-through","text":"","category":"section"},{"location":"demo_code/","page":"Demo: code walk-through","title":"Demo: code walk-through","text":"This section is under construction. Please see /examples/surrogate.jl in the repository folder for now.","category":"page"},{"location":"#NMRSignalSimulator.jl","page":"Overview","title":"NMRSignalSimulator.jl","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"NMRSignalSimulator is a library for simulating the 1D 1H NMR spectrum, using the output from NMRHamiltonian.","category":"page"},{"location":"#Table-of-contents","page":"Overview","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"Pages = [\"index.md\"]\nDepth = 5","category":"page"},{"location":"#Install","page":"Overview","title":"Install","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"NMRSignalSimulator is hosted on a custom Julia registry. We need to add that before installing NMRHamiltonian: ","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"using Pkg\nPkg.Registry.add(RegistrySpec(url = \"https://github.com/AI4DBiological-Systems/PublicJuliaRegistry\"))\nPkg.add(\"NMRSignalSimulator\")","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"To update this package once it is installed, do","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"using Pkg\nPkg.update(\"NMRSignalSimulator\")","category":"page"},{"location":"#Important-exported-functions","page":"Overview","title":"Important exported functions","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"The following functions are the focus of NMRSignalSimulator:","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"fitclproxies() assembles the simulation configuration container variable. This is a constructor to a composite data type (i.e., this data type is like the C programming language's struct).\nevalclmixture() evaluate the preliminary model.\nevalclproxymixture() evaluate the surrogate model.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"See the demo page for an example walk-through of how to use this library with NMRHamiltonian.","category":"page"},{"location":"#Information-to-get-started","page":"Overview","title":"Information to get started","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"See the Julia Basics, A compiled list of physical chemistry parameters, and Terminology from the NMRHamiltonian documentation.","category":"page"},{"location":"#Terminology","page":"Overview","title":"Terminology","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"The data structure for storing the surrogate model in NMRSignalSimulator is represented in index range orders by nested arrays. Let T be your choice of floating point data type, e.g., Float64. In the demo, the main data structure variables are named:","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"Bs, of type Vector{MoleculeType{T, SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}}},\nMSS  of type MixtureSpinSys{T, CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}\nmodel_param of type MixtureModelParameters{T, CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"NMRSignalSimulator provide basic evaluation routines for the model fit objective function, but leave the actual model optimization and formulation of the constants of the objective function to downstream packages.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"The approached used in NMRSignalSimulator to evaluate objective functions given a parameter is to parse a flat 1-D array into a data structure variable, update all structure fields accordingly, then evaluate the objective (also referred to as cost) function that uses that variable as input. These data structure variables together with pre-allocated buffers allow us to work with anonymous functions of the defined objective functions.","category":"page"},{"location":"#License","page":"Overview","title":"License","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"NMRSignalSimulator has the Mozilla Public License Version 2.0.","category":"page"},{"location":"#Authors","page":"Overview","title":"Authors","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"Code author:","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"Roy Chih Chung Wang","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"Supervisors:","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"Dave Campbell (Carleton University, Bank of Canada)\nMiroslava Čuperlović-Culf (National Research Council of Canada)","category":"page"},{"location":"#Funding","page":"Overview","title":"Funding","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"This projected was funded by the AI-for-Design Challenge Program from the National Research Council of Canada.","category":"page"}]
}
