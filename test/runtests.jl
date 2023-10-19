
using Test

using DataDeps, Tar, CodecZlib
#import FiniteDifferences

using LinearAlgebra

import NMRSignalSimulator
SIG = NMRSignalSimulator
HAM = NMRSignalSimulator.NMRHamiltonian

import BSON

NMRHamiltonian = NMRSignalSimulator.NMRHamiltonian
JSON3 = NMRSignalSimulator.NMRHamiltonian.JSON3

serializetype = NMRSignalSimulator.serializetype
deserializetype = NMRSignalSimulator.deserializetype

serializclproxies = NMRSignalSimulator.serializclproxies
deserializclproxies = NMRSignalSimulator.deserializclproxies

serializitpsamples = NMRSignalSimulator.serializitpsamples
deserializitpsamples = NMRSignalSimulator.deserializitpsamples


include("../examples/helpers/utils.jl")
include("../examples/helpers/data.jl")

include("./helpers/setup.jl")
include("./helpers/utils.jl")

import Random
Random.seed!(25)

@testset "surrogate, serialization" begin

    Ts = [Float32; Float64]

    for T in Ts
        # # Set up.
        molecule_entries = [
            "alpha-D-Glucose";
            "Ethanol";
            "DSS";
            "Singlet - 4.7 ppm";
            "Uridine";
            "beta-D-Glucose";
        ]
        w_oracle = rand(T, length(molecule_entries))

        #type_SSParams = NMRSignalSimulator.getSpinSysParamsdatatype(NMRSignalSimulator.SharedShift{Float64})
        
        model_params, As, Bs, itp_samps, lbs, ubs, U_rad, f_U, MSS = testsetup(molecule_entries, w_oracle)
        
        # backup current parameters.
        #SIG.exportmodel!(model_params)
        #x0 = copy(model_params.var_flat[:])

        # # tests.

        # short test.
        N_serialization_tests = 2
        N_cost_tests = 2

        # ## test serialization and deserialization.
        @testset "test serialization/deserialization" begin
            # observations for the cost-related tests
            inds = 1:length(U_rad)
            U_test = U_rad[inds]
            w_test = rand(T, length(w_oracle))

            testserializationouter(
                MSS,
                w_test,
                U_test,
                As,
                Bs,
                itp_samps;
                BSON_filename = "test.bson",
                JSON_filename = "test.json",
                clean_up_files = true,
                shift_proportion = convert(T, 0.9), # between 0 and 1. proportion of interpolation range to use in the parameter generation. see fetchbounds().
                zero_tol = convert(T, 1e-12),
            )
        end
    end

end
