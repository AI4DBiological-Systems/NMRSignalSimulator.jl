module NMRSignalSimulator

# Write your package code here.
using Distributed, LinearAlgebra, FFTW, SharedArrays, OffsetArrays, Statistics

import JSON, Interpolations
#import Kronecker, NearestNeighbors, Graphs


import NMRHamiltonian # https://github.com/AI4DBiological-Systems/NMRHamiltonian.jl

include("../src/types.jl")
include("../src/utils.jl")

include("./lorentzian/lorentzian.jl")
include("./lorentzian/setup_itp.jl")
include("./lorentzian/itp_evals.jl")

include("./FID/FID.jl")
#include("./FID/setup_itp.jl")
#include("./FID/itp_evals.jl")

end
