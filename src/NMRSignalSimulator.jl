module NMRSignalSimulator

# Write your package code here.
using LinearAlgebra

import JSON3, Interpolations, OffsetArrays



import NMRHamiltonian # consider removing this dependency, or split the data structure into a separate package.

include("types.jl")
include("./parameters/flatten_types.jl")
include("./lorentzian/model_types.jl")
include("./parameters/flat_model.jl") # inaccurate. paused development.

include("utils.jl")

include("./lorentzian/lorentzian.jl")
include("./lorentzian/setup_itp.jl")
include("./lorentzian/itp_evals.jl")

include("./lorentzian/derivatives.jl")

#include("./FID/FID.jl")
#include("./FID/setup_FID_itp.jl")
#include("./FID/FID_itp_evals.jl")

include("./parameters/updates.jl")

#include("test_helpers.jl")

end
