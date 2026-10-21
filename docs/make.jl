using Documenter
using CommToolkit

makedocs(;
    modules = [
        CommToolkit,
        CommToolkit.Modulation,
        CommToolkit.Channel,
        CommToolkit.Coding,
        CommToolkit.Equalization,
        CommToolkit.OFDM,
        CommToolkit.MIMO,
        CommToolkit.Visualization,
        CommToolkit.Metrics,
    ],
    sitename = "CommToolkit.jl",
    authors = "刘恒宇",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://longee142857.github.io/CommToolkit.jl/",
    ),
    pages = [
        "API Reference" => "index.md",
    ],
    checkdocs = :exports,
    doctest = false,
)

if get(ENV, "DOCUMENTER_KEY", "") != ""
    deploydocs(;
        repo = "github.com/longee142857/CommToolkit.jl.git",
        devbranch = "main",
        push_preview = true,
    )
end
