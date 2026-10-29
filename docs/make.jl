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
        CommToolkit.BICM,
        CommToolkit.EXIT,
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
    # Sidebar selector. Inventory header uses the same string, which is
    # `Project.toml`'s version. SemVer: 1.x only adds names.
    version = string(pkgversion(CommToolkit)),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://longee142857.github.io/CommToolkit.jl/",
        inventory_version = string(pkgversion(CommToolkit)),
        # One API page. 0.8.5 was under the previous 450 KiB cap; BICM-ID
        # docstrings share that same page.
        size_threshold = 520 * 1024,
        size_threshold_warn = 400 * 1024,
    ),
    pages = [
        "API Reference" => "index.md",
        "边界" => "boundaries.md",
        "公开接口" => "api.md",
        "验证" => "validation.md",
        "已知限制" => "limitations.md",
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
