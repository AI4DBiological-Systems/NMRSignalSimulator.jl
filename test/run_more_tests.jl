
# shared shift is actually not used in the publication.

@testset "mixture model shared shift" begin

    # # Set up.
    molecule_entries = [
        "alpha-D-Glucose";
        "Ethanol";
        "DSS";
        "Singlet - 4.7 ppm";
        "beta-D-Glucose";
    ]
    w_oracle = rand(length(molecule_entries))

    type_SSParams = NMRSignalSimulator.getSpinSysParamsdatatype(NMRSignalSimulator.SharedShift{Float64})
    #type_SSParams = NMRSignalSimulator.getSpinSysParamsdatatype(NMRSignalSimulator.CoherenceShift{Float64})

    model_params, As, Bs, itp_samps, lbs, ubs, U_rad, f_U, MSS = testsetup(molecule_entries, w_oracle, type_SSParams)
    
    # backup current parameters.
    NMRSignalSimulator.exportmodel!(model_params)
    x0 = copy(model_params.var_flat[:])

    # # tests.

    # short test.
    N_serialization_tests = 2
    N_cost_tests = 2

    # # long test.
    # N_serialization_tests = 50
    # N_cost_tests = 10

    # ## test serialization and deserialization.
    @testset "test serialization/deserialization" begin
        # observations for the cost-related tests
        inds = 1:length(U_rad)
        U_test = U_rad[inds]
        y_test = copy(f_U[inds] + randn(length(inds)))
        w_test = rand(length(w_oracle))

        testserializationouter(
            MSS,
            w_test,
            y_test,
            U_test,
            As,
            Bs,
            itp_samps;
            BSON_filename = "test.bson",
            JSON_filename = "test.json",
            clean_up_files = true,
            shift_proportion = 0.9, # between 0 and 1. proportion of interpolation range to use in the parameter generation. see fetchbounds().
            N_tests = N_serialization_tests, # number of tests per file serialization libraries. currently have two file libraries: JSON3.jl and BSON.jl.
            zero_tol = 1e-12,
        )
    end

    # ## import/export model parameters.
    @testset "test import/export" begin
        N_import_tests = 10
        testimportexport!(model_params, lbs, ubs, N_import_tests, U_rad)

        
        batch_size = 10
        N_simulator_gradient_tests = 2
        testgradient!(U_rad, model_params, lbs, ubs, batch_size, N_simulator_gradient_tests; atol = 1e-5)

        N_x_tests = 10
        testbatchgradient!(
            U_rad,
            model_params,
            lbs,
            ubs,
            N_x_tests,
        )
    end

    # ## test cost function and gradient
    # vs. numerical differentiation.
    @testset "test cost function, gradient vs. numerical differentiation" begin

        batch_size = 1000
        inds = rand(1:length(U_rad), batch_size) # for speed.
        #inds = 1:length(U_rad) # slow, but full frequency interval.
        U_test = U_rad[inds]
        #y = copy(f_U[inds])
        y = copy(f_U[inds] + randn(length(inds)))
        
        dis, ps = verifycostfunc(model_params, lbs, ubs, y, U_test, N_cost_tests)

        ND_zero_tol = 1e-5
        @show maximum(dis)
        @test maximum(dis) < ND_zero_tol
    end
end