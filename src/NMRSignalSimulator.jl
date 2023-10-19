module NMRSignalSimulator

# Write your package code here.
using LinearAlgebra

import Interpolations, OffsetArrays

import NMRHamiltonian # consider removing this dependency, or split the data structure into a separate package.
HAM = NMRHamiltonian

# inherit dependencies.
JSON3 = NMRHamiltonian.JSON3

include("types.jl")
include("parameter_types.jl")
include("subset_types.jl")
include("serialize.jl")

include("./parameters/updates.jl")
include("./parameters/models.jl")

include("utils.jl")

include("./core/query/lorentzian.jl")
include("./core/query/proxy.jl")
include("./core/engine.jl")
include("./core/IO.jl")


include("analysis.jl") # view the parameters of the surrogate.


#include("test_helpers.jl")
export 

CLSurrogateConfig,

# surrogate-related
fitclproxies,
evalclmixture,
evalclproxymixture,

# constructors
getSpinSysParamsdatatype,
MixtureModelParameters,
ParamsMapping,

# model updates
importmodel!,
exportmodel!,
importmodelreset!,
#updatebuffer!,

# subset vars
VariableSetTrait,
AllVars,
SubsetVars,

SubsetVarsIndices,
exportindices,
setupsubsetvars,

getindices,
getactivesystems,

# parameter range
OperationRange,
fetchbounds,
MoleculeParamsMapping,
ParamsMapping,
getshiftrange,
getphaserange,
getT2range,

#utilities
getNentries,
exportphysicalparams,
createtablecolumns,
ζ2Δcs,
Δcs2ζ,

# serialize
serializclproxies,
serializitpsamples,
deserializclproxies,
deserializitpsamples,
recoverclproxies


end
