module NMRSignalSimulator

# Write your package code here.
using LinearAlgebra

import JSON3, Interpolations, OffsetArrays



import NMRHamiltonian # consider removing this dependency, or split the data structure into a separate package.


include("../src/types.jl")
include("../src/utils.jl")

include("./lorentzian/lorentzian.jl")
include("./lorentzian/setup_itp.jl")
include("./lorentzian/itp_evals.jl")

include("./FID/FID.jl")
include("./FID/setup_FID_itp.jl")
include("./FID/FID_itp_evals.jl")

end
