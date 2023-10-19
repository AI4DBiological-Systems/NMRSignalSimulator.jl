

# taken from /examples/lorentzian.jl. Default test at 700 MHz spectrometer spin Hamiltonian simulation.
function testsetup(
    molecule_entries::Vector{String},
    w_oracle::Vector{T};
    rel_discrepancy_tol::T = convert(T, 1e-1),
    zero_tol::T = convert(T, 1e-12),
    fs::T = convert(T, 14005.602240896402),
    SW::T = convert(T, 20.0041938620844),
    ν_0ppm::T = convert(T, 10656.011933076665),
    shift_proportion::T = convert(T, 0.9),
    ) where T

    @assert length(w_oracle) == length(molecule_entries)

    hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
    ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

    root_data_path = getdatapath() # coupling values data repository root path

    H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.
    
    molecule_mapping_root_path = joinpath(
        root_data_path,
        "molecule_name_mapping",
    )
    molecule_mapping_file_path = joinpath(
        molecule_mapping_root_path,
        "select_molecules.json",
    )
    
    config = HAM.SHConfig{T}(
        coherence_tol = convert(T, 0.01),
        relative_α_threshold = convert(T, 0.005),
        tol_radius_1D = convert(T, 0.1), # strictly between 0 and 1. The lower, the better the approximation, but would a larger partition (i.e. more resonance groups).
        nuc_factor = convert(T, 1.5),
    )
    unique_cs_atol = convert(T, 1e-6)

    # Spin Hamiltonian simulation.
    Phys, As, MSPs = HAM.loadandsimulate(
        T,
        "700",
        molecule_entries,
        H_params_path,
        molecule_mapping_file_path;
        config = config,
        unique_cs_atol = unique_cs_atol
    )
    
    
    ###
    proxy_config = SIG.CLSurrogateConfig{T}(
        λ0 = convert(T, 3.4),
        Δr = convert(T, 1.0), # radial frequency resolution: smaller means slower to build surrogate, but more accurate.
        Δκ_λ = convert(T, 0.05), # T2 multiplier resolution. smaller means slower to build surrogate, but more accurate.
        Δcs_max_scalar = convert(T, 0.2), # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.
        κ_λ_lb = convert(T, 0.5), # lower limit for κ_λ for which the surrogate is made from.
        κ_λ_ub = convert(T, 2.5), # upper limit for κ_λ for which the surrogate is made from.
        ppm_padding = convert(T , 0.5),
    )
    
    Bs, MSS, itp_samps = SIG.fitclproxies(As, proxy_config)
    
    # test.
    w_oracle = ones(T, length(molecule_entries))

    ΩS_ppm = collect( hz2ppmfunc.( SIG.combinevectors(A.Ωs) ./ convert(T, 2*π) ) for A in As )
    ΩS_ppm_flat = SIG.combinevectors(ΩS_ppm)
    P_max = maximum(ΩS_ppm_flat) + convert(T, 0.5)
    P_min = minimum(ΩS_ppm_flat) - convert(T, 0.5)

    P = LinRange(P_min, P_max, 80000)
    U = ppm2hzfunc.(P)
    U_rad = U .* convert(T, 2*π)

    f = uu->SIG.evalclmixture(uu, As, Bs; w = w_oracle)
    q = uu->SIG.evalclproxymixture(uu, As, Bs; w = w_oracle)

    f_U = f.(U_rad)
    q_U = q.(U_rad)

    discrepancy = abs.(f_U-q_U)

    @test norm(discrepancy)/norm(f_U) < rel_discrepancy_tol


    # prepare return quantities.
    model_params = SIG.MixtureModelParameters(MSS; w = copy(w_oracle))

    lbs, ubs = SIG.fetchbounds(model_params, Bs; shift_proportion = shift_proportion)

    return model_params, As, Bs, itp_samps, lbs, ubs, U_rad, f_U, MSS
end


###################### serialization.


function testserializationouter(
    MSS,
    w::Vector{T},
    U_rad,
    As,
    Bs,
    itp_samps;
    BSON_filename = "test.bson",
    JSON_filename = "test.json",
    clean_up_files = true,
    shift_proportion::T = convert(T, 0.9), # between 0 and 1. proportion of interpolation range to use in the parameter generation. see fetchbounds().
    zero_tol::T = convert(T, 1e-12),
    ) where T

    # # Set up
    model_params = SIG.MixtureModelParameters(MSS; w = w)
    lbs, ubs = SIG.fetchbounds(model_params, Bs; shift_proportion = shift_proportion)

    # fill in non-default values.
    p_test = SIG.generateparameters(lbs, ubs)
    SIG.importmodel!(model_params, p_test)

    # # Serialize As
    S_As = HAM.serializemixture(As)
    dict_As_BSON = saveloadBSON(BSON_filename, S_As)
    #dict_As_JSON = saveloadJSON(JSON_filename, S_As)

    # # Serialize Bs.

    # ## ss_paras, op_range.
    S_Bs = serializclproxies(Bs)

    dict_Bs_BSON = saveloadBSON(BSON_filename, S_Bs)
    #dict_Bs_JSON = saveloadJSON(JSON_filename, S_Bs)

    # ## itp samps.
    S_itp = serializitpsamples(itp_samps)
    dict_itp = saveloadBSON(BSON_filename, S_itp) # only support BSON loading for this.

    ## clean up.
    cleanup(BSON_filename, JSON_filename, clean_up_files)
    
    # # test serialization using BSON.jl
    testserialization(
        U_rad,
        As,
        Bs,
        dict_As_BSON, # serialized As.
        dict_Bs_BSON,
        dict_itp,
        w;
        zero_tol = zero_tol,
    )

    # # test serialization using JSON3.jl

    # fill in non-default values.
    p_test = SIG.generateparameters(lbs, ubs)
    SIG.importmodel!(model_params, p_test)

    # # Serialize As
    S_As = HAM.serializemixture(As)
    #dict_As_BSON = saveloadBSON(BSON_filename, S_As)
    dict_As_JSON = saveloadJSON(JSON_filename, S_As)

    # # Serialize Bs.

    # ## ss_paras, op_range.
    S_Bs = serializclproxies(Bs)

    #dict_Bs_BSON = saveloadBSON(BSON_filename, S_Bs)
    dict_Bs_JSON = saveloadJSON(JSON_filename, S_Bs)

    # ## itp samps.
    S_itp = serializitpsamples(itp_samps)
    dict_itp = saveloadBSON(BSON_filename, S_itp) # only support BSON loading for this.

    ## clean up.
    cleanup(BSON_filename, JSON_filename, clean_up_files)

    testserialization(
        U_rad,
        As,
        Bs,
        dict_As_JSON, # serialized As.
        dict_Bs_JSON,
        dict_itp,
        w;
        zero_tol = zero_tol,
    )

    return nothing
end

# based on `/examples/IO.jl`.
function testserialization(
    U_rad,
    As,
    Bs,
    dict_As, # serialized As.
    dict_Bs,
    dict_itp,
    w::Vector{T};
    zero_tol::T = convert(T, 1e-12),
    ) where T

    As2 = HAM.deserializemixture(dict_As)
    
    ss_params_set2, op_range_set2, λ02 = deserializclproxies(dict_Bs)
    itp_samps2 = deserializitpsamples(dict_itp)

    Bs2, MSS2 = SIG.recoverclproxies(
        itp_samps2,
        ss_params_set2,
        op_range_set2,
        As2,
        λ02,
    )

    # verify forward eval using As, Bs.
    q = uu->SIG.evalclproxymixture(uu, As, Bs; w = w)
    q2 = uu->SIG.evalclproxymixture(uu, As2, Bs2; w = w)

    q_U = q.(U_rad)
    q2_U = q2.(U_rad)
    @test norm(q_U - q2_U) < zero_tol

    return nothing
end