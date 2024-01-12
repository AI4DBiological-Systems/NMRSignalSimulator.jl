using Literate


function replace_includes(str)

    included = ["./snippet/plot_groups.jl";]
    path = "./"

    for ex in included
        content = read(path*ex, String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end
    return str
end

Literate.markdown(
    "demo.jl",
    "../docs/src/generated/";
    execute = true,
    name = "glutamine",
    preprocess = replace_includes,
)
