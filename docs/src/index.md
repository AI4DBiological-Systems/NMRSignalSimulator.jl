# NMRSignalSimulator.jl
NMRSignalSimulator is a library for simulating the 1D 1H NMR spectrum, using the output from [NMRHamiltonian](https://github.com/AI4DBiological-Systems/NMRHamiltonian.jl).

## Table of contents
```@contents
Pages = ["index.md"]
Depth = 5
```

## Install
NMRSignalSimulator is hosted on a custom Julia registry. We need to add that before installing NMRHamiltonian: 
``` julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/AI4DBiological-Systems/PublicJuliaRegistry"))
Pkg.add("NMRSignalSimulator")
```

To update this package once it is installed, do
``` julia
using Pkg
Pkg.update("NMRSignalSimulator")
```

# Important exported functions
The following functions are the focus of NMRSignalSimulator:

* `fitclproxies()` assembles the simulation configuration container variable. This is a constructor to a composite data type (i.e., this data type is like the C programming language's `struct`).

* `evalclmixture()` evaluate the preliminary model.

* `evalclproxymixture()` evaluate the surrogate model.

See the demo page for an example walk-through of how to use this library with NMRHamiltonian.

# Information to get started

See the *Julia Basics*, *A compiled list of physical chemistry parameters*, and *Terminology* from the [NMRHamiltonian.jl](https://github.com/AI4DBiological-Systems/NMRHamiltonian.jl) documentation website.


# Terminology
The data structure for storing the surrogate model in NMRSignalSimulator is represented in index range orders by nested arrays. Let `T` be your choice of floating point data type, e.g., `Float64`. In the demo, the main data structure variables are named:

- `Bs`, of type `Vector{MoleculeType{T, SpinSysParams{CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}, OT}}`. `OT <: OperationRange` is `FIDOperationRange` for a free-induction decay time-domain model or `CLOperationRange` for a complex Lorentzian frequency-domain model.

- `MSS`  of type `CLMixtureSpinSys{T, CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}`

- `model_param` of type `MixtureModelParameters{T, CoherenceShift{T}, CoherencePhase{T}, SharedT2{T}}`.

NMRSignalSimulator provide basic evaluation routines for the model fit objective function, but leave the actual model optimization and formulation of the constants of the objective function to downstream packages.

The approached used in NMRSignalSimulator to evaluate objective functions given a parameter is to parse a *flat* 1-D array into a data structure variable, update all structure fields accordingly, then evaluate the objective (also referred to as *cost*) function that uses that variable as input. These data structure variables together with pre-allocated buffers allow us to work with anonymous functions of the defined objective functions.

# License
NMRSignalSimulator has the Mozilla Public License Version 2.0.

# Authors
Code author:
- Roy Chih Chung Wang

Supervisors:
- Dave Campbell (Carleton University)
- Miroslava Čuperlović-Culf (National Research Council of Canada)

# Funding
This projected was funded by the AI-for-Design Challenge Program from the National Research Council of Canada.