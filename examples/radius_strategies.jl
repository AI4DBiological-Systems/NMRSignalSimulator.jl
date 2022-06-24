
import NMRHamiltonian

include("../src/NMRSignalSimulator.jl")
import .NMRSignalSimulator

using LinearAlgebra
using FFTW
import PyPlot
import JSON

#import Clustering
import Statistics

include("../examples/helpers/plot_partition.jl")
include("../examples/helpers/utils.jl")


import Random
Random.seed!(25)

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])




# try different Δc_bar grouping strategies.
println("name = ", molecule_names[1])
println("length.(A.αs)  = ", length.(A.αs))
Cs = collect( array2matrix(A.Δc_bar[i]) for i = 1:length(A.Δc_bar) )
displaymatrix(Cs[1])

dummy = 1
