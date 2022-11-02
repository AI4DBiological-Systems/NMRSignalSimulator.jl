
# don't use this.
function plotNMRspectrum(
    #P,
    #ys,
    qs::Vector{Function},
    label_texts::Vector{String},
    plot_lbs::Vector{T},
    plot_ubs::Vector{T};


    fig_size = fig_size,
    fig_title_string = "$(molecule_entries[1]) at $(spectrometer_freq) MHz",
    plot_title_string = "$(molecule_entries[1])\n$(spectrometer_freq) MHz",
    subplot_spacing = subplot_spacing,
    line_styles = line_styles,
    legend_space_per_row = legend_space_per_row,
    N_cols_legend = N_cols_legend,
    legend_font = legend_font,
    x_tick_font_size = x_tick_font_size,
    line_width = line_width,



    title_font_size = 30,
    grid_visible = false,
    line_width_list = ones(length(qs)) .* 3.0,
    save_folder_path = "./output",
    save_name = "figure.svg",
    pt_per_inch = 72,
    font_size = 12, # default font size for text that isn't specified.
    ,
    
    ) where T

    line_style_list = collect( :solid for _ = 1:length(qs))


    # I am here. need to add black into the sum line using the MakiePlots color func.
end
