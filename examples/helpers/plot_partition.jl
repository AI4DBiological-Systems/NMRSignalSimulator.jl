
"""
save resonance groupings of a compound.
Choices for `f`` are: `real()`, `imag()`, or `abs()`. These corresponds to real part, imaginary part, and magnitude spectrum, respectively.
"""
function plotgroups(title_string::String,
    P,
    U,
    q,
    qs,
    q_singlets,
    f::Function,
    return_val_type::T;
    canvas_size::Tuple{Int,Int} = (1000, 400)) where T

    U_rad = U .* (2*π)

    q_U = q.(U_rad)
    plot_obj = Plots.plot( P,
        f.(q_U),
        title = title_string,
        label = "sum of sub-models",
        seriestype = :line,
        ticks = :native,
        xlims = (P[1],P[end]),
        hover = P,
        linewidth = 4,
        xflip = true,
        size = canvas_size)

    qs_U = Vector{Vector{Vector{Complex{T}}}}(undef, length(qs))
    for i = 1:length(qs)

        qs_U[i] = Vector{Vector{Complex{T}}}(undef, length(qs[i]))
        for k = 1:length(qs[i])

            qs_U[i][k] = qs[i][k].(U_rad)

            Plots.plot!(plot_obj, P, f.(qs_U[i][k]), label = "sub-group ($(i),$(k))",
                seriestype = :line,
                linestyle = :dot,
                xflip = true,
                linewidth = 4)
        end
    end

    q_singlets_U = q_singlets.(U_rad)
    Plots.plot!(plot_obj, P, f.(q_singlets_U), label = "group of all singlets",
        seriestype = :line,
        linestyle = :dot,
        xflip = true,
        linewidth = 4)
    #
    #Plots.plot!(plot_obj, P, f.(q_U),
    #markershape = :circle,
    #seriestype = :scatter,
    #xflip = true)

    return plot_obj, q_U, qs_U, q_singlets_U
end

"""
only plots the first compound in molecule_names.
"""
function plotgroupsfullscript(plot_title, molecule_names,
    H_params_path, project_path, tol_coherence, α_relative_lower_threshold,
    Δc_partition_radius, Δcs_max, κ_λ_lb, κ_λ_ub,
    fs, SW, ν_0ppm, λ_0ppm, save_folder;
    display_flag = false,
    prune_low_signal_for_display_flag::Bool = false,
    display_reduction_factor = 100,
    display_threshold_factor =  α_relative_lower_threshold/10,
    canvas_size = (1600, 900),
    reset_plot_dir_flag = true,
    plot_imag_and_mag_flag = false,
    save_plot_flag = true)

    load_path = joinpath(project_path, "experiment.bson")
    save_folder = joinpath(project_path, "plots")
    isdir(save_folder) || mkpath(save_folder)

    ### load block.
    dic = BSON.load(load_path)
    fs = dic[:fs]
    SW = dic[:SW]
    ν_0ppm = dic[:ν_0ppm]
    λ_0ppm = dic[:λ_0ppm]

    hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
    ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)


    Δcs_max_mixture = collect( Δcs_max for i = 1:length(molecule_names))

    # get a surrogate where K_{n,i} is encouraged to be no larger than `early_exit_part_size`.
    # println("Timing: setupmixtureproxies()")
    # @time mixture_params = NMRSignalSimulator.setupmixtureproxies(molecule_names,
    mixture_params = NMRSignalSimulator.setupmixtureproxies(molecule_names,
        H_params_path, ppm2hzfunc, fs, SW,
        λ_0ppm, ν_0ppm, dummy_SSFID;
        tol_coherence = tol_coherence,
        α_relative_lower_threshold = α_relative_lower_threshold,
        Δc_partition_radius = Δc_partition_radius)
    As = mixture_params
    A = As[1]



    ## frequency locations. For plotting.
    ΩS_ppm = NMRSignalSimulator.getPsnospininfo(mixture_params, hz2ppmfunc)
    ΩS_ppm_sorted = sort(NMRSignalSimulator.combinevectors(ΩS_ppm))

    u_offset = 0.5 # in units ppm.
    u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - u_offset)
    u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + u_offset)

    P = LinRange(hz2ppmfunc(u_min), hz2ppmfunc(u_max), 50000)
    U = ppm2hzfunc.(P)
    U_rad = U .* (2*π)

    ## end frequency locations.

    # println("Timing: fitclproxies()")
    # @time NMRSignalSimulator.fitproxies!(As;
    NMRSignalSimulator.fitproxies!(As;
        κ_λ_lb = κ_λ_lb,
        κ_λ_ub = κ_λ_ub,
        u_min = u_min,
        u_max = u_max,
        Δr = 1.0,
        Δκ_λ = 0.05)

    # store fitted surrogate for future use.
    #save_simulation_path = joinpath(project_path, "simulation.bson")

    #NMRSignalSimulator.removeauxinfo!(mixture_params)
    #BSON.bson(save_simulation_path, mixture_params = mixture_params)

    # create the functions for each resonance group.
    qs = collect( collect( ωω->A.qs[i][k](ωω-A.ss_params.d[i], A.ss_params.κs_λ[i]) for k = 1:length(A.qs[i]) ) for i = 1:length(A.qs) )
    q_singlets = ωω->NMRSignalSimulator.evalsinglets(ωω, A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets)

    # create the function for the entire compound.
    q = uu->NMRSignalSimulator.evalclproxymixture(uu, As[1:1])

    # evaluate at the plotting positions.
    q_U = q.(U_rad)

    qs_U = collect( collect( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
    q_singlets_U = q_singlets.(U_rad)

    # sanity check.
    q_check_U = q_singlets_U
    if !isempty(qs) # some compounds only have singlets.
        q_check_U += sum( sum( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
    end
    discrepancy = norm(q_check_U- q_U)

    #println("sanity check. This should be numerically zero: ", discrepancy)
    @assert discrepancy < 1e-14


    # plot.
    P_display = P
    U_display = U
    if prune_low_signal_for_display_flag

        inds, _ = NMRSignalSimulator.prunelowsignalentries(q_U, display_threshold_factor, display_reduction_factor)
        P_display = P[inds]
        U_display = U[inds]
    end

    #println("length(P) = ", length(P))
    #println("length(P_display) = ", length(P_display))

    # should be clear the plotting directory before saving plots?
    if reset_plot_dir_flag && save_plot_flag
        paths_to_delete = readdir(save_folder, join = true)
        rm.(paths_to_delete)
    end

    plots_save_path = joinpath(save_folder, "groups_real.html")
    title_string = "$(plot_title), $(molecule_names), real spectrum"
    plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P_display, U_display, q, qs, q_singlets, real, P[1]; canvas_size = canvas_size)
    if save_plot_flag
        Plots.savefig(plot_obj, plots_save_path)
    end
    if display_flag
        display(plot_obj)
    end
    #println("length(q_U) = ", length(q_U))
    #@assert 1==44

    if plot_imag_and_mag_flag

        plots_save_path = joinpath(save_folder, "groups_imag.html")
        title_string = "$(plot_title), $(molecule_names), imaginary spectrum"
        plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P_display, U_display, q, qs, q_singlets, imag, P[1]; canvas_size = canvas_size)
        if save_plot_flag
            Plots.savefig(plot_obj, plots_save_path)
        end
        if display_flag
            display(plot_obj)
        end

        plots_save_path = joinpath(save_folder, "groups_magnitude.html")
        title_string = "$(plot_title), $(molecule_names), magnitude spectrum"
        plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P_display, U_display, q, qs, q_singlets, abs, P[1]; canvas_size = canvas_size)
        if save_plot_flag
            Plots.savefig(plot_obj, plots_save_path)
        end
        if display_flag
            display(plot_obj)
        end
    end
end


"""
root_folder should contain only folders of 1D 1H NMR compound standard experiments, with the compound name being the subfolder name.
I.e., a folder containing an experiment of L-Serine and L-Histidine should have two folders in it with those names that contain the corresponding experimental data.
"""
function batchplotgroups(plot_title, root_path,
    H_params_path, tol_coherence,
    α_relative_lower_threshold, Δc_partition_radius, Δcs_max, κ_λ_lb, κ_λ_ub;
    display_flag = false,
    prune_low_signal_for_display_flag::Bool = false,
    display_reduction_factor = 100,
    display_threshold_factor =  α_relative_lower_threshold/10,
    canvas_size = (1600, 900),
    reset_plot_dir_flag = true,
    plot_imag_and_mag_flag = false,
    save_plot_flag = true)

    # get a list of experiment paths and compounds names. Assume the folder name of the experiment contains the compound.
    tmp = readdir(root_path, join = true)
    inds = findall(xx->isdir(xx), tmp)

    project_paths = tmp[inds]
    experiment_names = readdir(root_path)[inds]
    project_names = collect( "$(experiment_names[i])" for i = 1:length(experiment_names))

    for i = 1:length(project_names)

        project_path = project_paths[i]
        compound_name = experiment_names[i]

        println("Now on ", project_path)

        # make sure this compound is in our library.
        records = GISSMOReader.getGISSMOentriesall()
        record_names = collect( records[i].molecule_name for i = 1:length(records) )

        k = findfirst(xx->xx==compound_name, record_names)

        # decide if we should simulate and plot.
        if typeof(k) != Nothing
            plotgroupsfullscript(plot_title, [compound_name;],
                H_params_path, project_path, tol_coherence, α_relative_lower_threshold,
                Δc_partition_radius, Δcs_max, κ_λ_lb, κ_λ_ub;
                display_flag = display_flag,
                prune_low_signal_for_display_flag = prune_low_signal_for_display_flag,
                display_reduction_factor = display_reduction_factor,
                display_threshold_factor = display_threshold_factor,
                canvas_size = canvas_size,
                reset_plot_dir_flag = reset_plot_dir_flag,
                plot_imag_and_mag_flag = plot_imag_and_mag_flag,
                save_plot_flag = save_plot_flag)
        else
            println("Compound not in the local GISSMO library. Skip.")
        end
    end

end


"""
prunelowsignalentries(x, threshold_factor::T, reduction_factor::Int)::Vector{Int} where T

keep the high signal portion and the reduced low signal portion.
keep_threshold = maximum(abs_x) * threshold_factor.
"""
function prunelowsignalentries(x, threshold_factor::T, reduction_factor::Int) where T

    abs_x = abs.(x)

    # find indices that have a low signal.
    keep_threshold = maximum(abs_x) * threshold_factor
    inds_for_reducing = findall(xx->xx<keep_threshold, abs_x)
    reduced_inds = inds_for_reducing[1:reduction_factor:length(inds_for_reducing)]

    # find indices that have a high signal.
    keep_inds = findall(xx->xx>keep_threshold, abs_x)

    # keep the high signal portion and the reduced low signal portion.
    inds = sort([keep_inds; reduced_inds])

    return inds, keep_inds, inds_for_reducing
end
