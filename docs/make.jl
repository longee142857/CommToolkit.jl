using Documenter
using CommToolkit

makedocs(;
    modules = [
        CommToolkit,
        CommToolkit.Modulation,
        CommToolkit.Channel,
        CommToolkit.Pulse,
        CommToolkit.Coding,
        CommToolkit.Interleaving,
        CommToolkit.Link,
        CommToolkit.Equalization,
        CommToolkit.OFDM,
        CommToolkit.MIMO,
        CommToolkit.Visualization,
        CommToolkit.Scrambling,
        CommToolkit.Sync,
        CommToolkit.Metrics,
    ],
    sitename = "CommToolkit.jl",
    authors = "刘恒宇",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://longee142857.github.io/CommToolkit.jl/",
        # One API page. 0.8.4 measured about 274 KiB; Clarke/Jakes and
        # Nakagami docstrings share that page.
        size_threshold = 450 * 1024,
        size_threshold_warn = 320 * 1024,
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
