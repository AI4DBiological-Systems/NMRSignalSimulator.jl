using NMRSignalSimulator
using Test

#include("./helpers/updates.jl")

@testset "NMRSignalSimulator.jl" begin
    # save the MEs to JLD, and load X, Bs, As, etc.
    # in progress. taken from /examples/model.jl

    p_test, maping_test = testupdateparameters1!(shifts, phases, T2s)

    discrepancy = testupdateparameters2(p_test, shifts, phases, T2s)
    @show discrepancy
end
