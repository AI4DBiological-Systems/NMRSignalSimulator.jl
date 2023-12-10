# NMRSignalSimulator.jl
Constructs the surrogate and preliminary models for the frequency domain 1D 1H NMR for pulse sequences that are equivalent to the 90 degrees y pulse.

# Documentation
[https://AI4DBiological-Systems.github.io/NMRSignalSimulator.jl/](https://AI4DBiological-Systems.github.io/NMRSignalSimulator.jl/)

# Install
Add the custom registries for dependencies, and then add the package.
```
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/RoyCCWang/RWPublicJuliaRegistry"))
Pkg.Registry.add(RegistrySpec(url = "https://github.com/AI4DBiological-Systems/PublicJuliaRegistry"))

pkg"add NMRSignalSimulator"
```
# Citation
Our work is undergoing peer review. Please cite our ChemRxiv version if you use this software.
[https://doi.org/10.26434/chemrxiv-2023-0s196](https://doi.org/10.26434/chemrxiv-2023-0s196)

# License
See LICENSE.md for the licence.