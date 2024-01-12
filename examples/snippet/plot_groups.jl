

x_plot, y_plots, line_styles, name_strings, title_string, x_label, y_label, x_autorange = ARGS #hide

# title_string = "my plot" #hide
# y_label = "real part of spectrum" #hide
# x_label = "ppm" #hide
# x_autorange = "reversed" # or true. #hide

plot_data_set = collect(
    PLY.Config(
        x = x_plot,
        y = y_plots[m],
        name = name_strings[m],
        type = "line-chart",
        line = PLY.Config(dash = line_styles[m]),
    ) for m in eachindex(y_plots)
)
plot_data = vcat(plot_data_set...)

layout = PLY.Config(
    #type = "line-chart",
    title = title_string,
    height = 400,
    width = 800,
    margin = 10,
    font = PLY.Config(
        family = "Lato",
        size = 16,
    ),
    xaxis = PLY.Config(
        title = x_label,
        autorange = x_autorange,
        titlefont = PLY.Config(
            color = "black",
            size = 12,
        )
    ),
    yaxis = PLY.Config(
        title = y_label,        
        titlefont = PLY.Config(
            color = "black",
            size = 12,
        )
    ),
)
ph = PLY.Plot(plot_data, layout);
# display(ph) # shows up in browser. #hide
# # html_plot = PlotlyDocumenter.to_documenter(ph) #this takes a while in Literate.jl #hide

# save to HTML with Cobweb.jl #hide
plot_counter += 1 # hide
Cobweb.save(Cobweb.Page(ph), "../docs/src/generated/plot_$plot_counter.html") #hide