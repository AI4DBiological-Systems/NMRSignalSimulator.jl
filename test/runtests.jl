using NMRSignalSimulator
using Test

#include("./helpers/updates.jl")

@testset "NMRSignalSimulator.jl" begin
    # save the MEs to JLD, and load X, Bs, As, etc.
    # in progress. taken from /examples/model.jl

    discrepancy, p = NMRSignalSimulator.testupdate!(MSS, MS)

    @show discrepancy
end
