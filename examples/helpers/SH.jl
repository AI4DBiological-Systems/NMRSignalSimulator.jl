function runSH(
    H_params_path,
    molecule_mapping_file_path,
    fs::T,
    SW::T,
    ν_0ppm::T,
    molecule_entries::Vector{String},
    max_partition_size_offset::Int;
    search_θ = true,
    θ_default = 0.0,
    γ_base = 0.1,
    γ_rate = 1.05,
    max_iters_γ = 100,
    fully_connected_convex_clustering = false,
    max_connected_components_offset = -1,
    start_knn = 60,
    ) where T <: AbstractFloat

    Phys, dict_molecule_to_filename = NMRHamiltonian.getphysicalparameters(
        molecule_entries,
        H_params_path,
        molecule_mapping_file_path;
        unique_cs_atol = 1e-6,
    )

    # extract chemical shifts for spin systems and singlets.
    cs_sys_mixture, cs_singlets_mixture = NMRHamiltonian.extractcs(Phys)

    # if using default tolerances for coherence and intensity thresholding.
    mixture_sh_config = NMRHamiltonian.defaultmixtureshsconfig(cs_sys_mixture)
    mixture_parts_params = NMRHamiltonian.defaultmixturepartitionsparameters(
        cs_sys_mixture,
        θ_default
    )
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



    searchknnconfigfunc = (nn, ii, cc, aa)->NMRHamiltonian.defaultknnsearchconfig(
        nn, ii, cc, aa;
        verbose = true,
        start_knn = start_knn,
        max_connected_components_offset  = max_connected_components_offset,
    )
    if fully_connected_convex_clustering
        searchknnconfigfunc = (nn, ii, cc, aa)->NMRHamiltonian.defaultknnsearchconfig(
            nn, ii, cc, aa;
            verbose = true,
            start_knn = length(cc),
        )
    end

    getsearchθconfigfunc = NMRHamiltonian.disablesearch
    if search_θ
        getsearchθconfigfunc = NMRHamiltonian.createsearchθconfigs
    end
    searchγconfigfunc = (nn, ii, cc, aa)->NMRHamiltonian.createsearchγconfigs(
        nn, ii, cc, aa;
        max_partition_size_offset = max_partition_size_offset,
        γ_base = γ_base,
        γ_rate = γ_rate,
        max_iters_γ = max_iters_γ,
    )
    part_algs = NMRHamiltonian.generatemixturepartitionalgorithm(
        molecule_entries,
        θs,
        γs,
        Phys;
        #getgraphconfigfunc = constantradiusfunc,
        getgraphconfigfunc = NMRHamiltonian.defaultknnsearchconfig,
        getsearchθconfigfunc = getsearchθconfigfunc,
        getsearchγconfigfunc = searchγconfigfunc,
        report_cost = true,
        verbose_kernel = true
        )

    As, Rs = NMRHamiltonian.setupmixtureSH(part_algs,
        molecule_entries,
        fs, SW, ν_0ppm,
        Phys,
        mixture_sh_config)

    return As, Rs
end