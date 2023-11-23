
using Test


using LinearAlgebra

import PublicationDatasets as DS

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

include("./helpers/setup.jl")
include("./helpers/utils.jl")

import Random
Random.seed!(25)

@testset "surrogate, serialization" begin

    entries = Vector{Vector{String}}(undef, 2)
    entries[1] = [
        "alpha-D-Glucose";
        "beta-D-Glucose";
        "Glycine";
        "L-Serine";
        "DSS";
        "Singlet - 4.9 ppm";
    ]
    entries[2] = [
        "alpha-D-Glucose";
        "Ethanol";
        "DSS";
        "Singlet - 4.7 ppm";
        "Uridine";
        "beta-D-Glucose";
    ]
    
    Ts = [Float32; Float64]

    for T in Ts
        for molecule_entries in entries
            # # Set up.

            w_oracle = rand(T, length(molecule_entries))

            #type_SSParams = NMRSignalSimulator.getSpinSysParamsdatatype(NMRSignalSimulator.SharedShift{Float64})
            
            model_params, As, Bs, itp_samps, lbs, ubs, U_rad, f_U, MSS = testsetup(
                molecule_entries,
                w_oracle,
            )
            
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

end
