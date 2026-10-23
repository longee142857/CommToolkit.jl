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
        # One API page. 0.8.3 (SOVA, turbo, BCH) and 0.8.4 (sync) share it.
        # The rendered HTML is about 274 KiB, above Documenter's 200 KiB default.
        size_threshold = 400 * 1024,
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
