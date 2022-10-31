function runSH(molecule_entries)
    root_data_path = getdatapath() # coupling values data repository root path

    H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.

    molecule_mapping_root_path = joinpath(root_data_path, "molecule_name_mapping")
    molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "select_molecules.json")
    #molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "GISSMO_names.json")


    println("Timing: getphysicalparameters")
    @time Phys, dict_molecule_to_filename = NMRHamiltonian.getphysicalparameters(molecule_entries,
        H_params_path,
        molecule_mapping_file_path;
        unique_cs_atol = 1e-6)

    # extract chemical shifts for spin systems and singlets.
    cs_sys_mixture, cs_singlets_mixture = NMRHamiltonian.extractcs(Phys)

    # if using default tolerances for coherence and intensity thresholding.
    mixture_sh_config = NMRHamiltonian.defaultmixtureshsconfig(cs_sys_mixture)
    mixture_parts_params = NMRHamiltonian.defaultmixturepartitionsparameters(cs_sys_mixture, 1.0)
    # # if loading from file.
    # file_filder = "./configs"
    # mixture_sh_config = NMRHamiltonian.loadmixtureshsconfig(
    #     cs_sys_mixture,
    #     joinpath(file_folder, "mixture_SH.json"),
    #     molecule_entries,
    #     Float64,
    # )
    # mixture_parts_params = NMRHamiltonian.loadmixturepartitionsparameters(
    #     cs_sys_mixture,
    #     joinpath(file_folder, "mixture_partition.json"),
    #     molecule_entries,
    #     2.0
    # )

    constantknnfunc, constantradiusfunc, θs,
        γs = NMRHamiltonian.setupconstantparameteroptions(molecule_entries, mixture_parts_params)

    # getgraphconfigfunc can be one of the following.:
    #   - defaultknnsearchconfig
    #   - defaultknnconfig
    #   - defaultradiusconfig
    #   - defaultradiussearchconfig
    #   - constantknnfunc
    #   - constantradiusfunc

    part_algs = NMRHamiltonian.generatemixturepartitionalgorithm(
        molecule_entries,
        θs,
        γs,
        Phys;
        #getgraphconfigfunc = constantradiusfunc,
        #getgraphconfigfunc = NMRHamiltonian.defaultknnsearchconfig,
        #getsearchθconfigfunc = NMRHamiltonian.createsearchθconfigs,
        # getsearchγconfigfunc = NMRHamiltonian.createsearchγconfigs,
        #runcgsolverlib;
        #store_trace = false,
        report_cost = true,
        verbose_kernel = true
        )

    println("Timing: setupmixtureproxies()")
    @time As, Rs = NMRHamiltonian.setupmixtureSH(part_algs,
        molecule_entries,
        fs, SW, ν_0ppm,
        Phys,
        mixture_sh_config)

    return As, Rs
end