# The exported name set is the stability surface in docs/src/api.md.
# Submodule exports must match it, plus the one intentional qualified name.

@testset "Public API contract" begin
    @test pkgversion(CommToolkit) == v"1.2.0"
    root = pkgdir(CommToolkit)
    api = read(joinpath(root, "docs", "src", "api.md"), String)
    @test occursin("## v1.0 承诺", api)
    section = match(r"## 稳定导出\n(.*?)(?=\n## |\z)"s, api)
    @test section !== nothing
    listed = String[]
    for line in split(section.captures[1], '\n')
        item = match(r"^- `([^`]+)`\s*$", line)
        item === nothing && continue
        push!(listed, item.captures[1])
    end
    # `names` includes the module's own name. That binding is not an algorithm export.
    surface = setdiff(names(CommToolkit), (nameof(CommToolkit),))
    exported = sort!(string.(surface))
    @test sort(listed) == exported
    @test length(listed) == length(unique(listed))

    @test isdefined(CommToolkit, :modulation)
    @test !(:modulation in surface)
    @test occursin("CommToolkit.modulation", api)

    submodules = (
        :Metrics, :Scrambling, :Modulation, :Channel, :Pulse, :Coding,
        :Interleaving, :BICM, :EXIT, :Link, :XApp, :Equalization, :OFDM, :MIMO,
        :Visualization, :Sync,
    )
    union_exports = Set{Symbol}()
    for name in submodules
        @test isdefined(CommToolkit, name)
        for exported_name in names(getproperty(CommToolkit, name))
            push!(union_exports, exported_name)
        end
    end
    # Each submodule's `names` also contains that submodule. `modulation` stays qualified.
    extra = setdiff(union_exports, surface)
    @test extra == Set([submodules..., :modulation])
    @test Set(surface) == setdiff(union_exports, extra)

    reference = read(joinpath(root, "docs", "src", "index.md"), String)
    documented = Set{String}()
    for block in eachmatch(r"```@docs(.*?)```"s, reference)
        for line in split(block.captures[1], '\n')
            line = strip(line)
            isempty(line) && continue
            push!(documented, line)
        end
    end
    for name in exported
        @test name in documented
    end
    @test "CommToolkit" in documented
    @test "CommToolkit.modulation" in documented

    validation = read(joinpath(root, "docs", "src", "validation.md"), String)
    for file in readdir(joinpath(root, "test"))
        endswith(file, ".jl") || continue
        @test occursin(file, validation)
    end

    bounds = read(joinpath(root, "docs", "src", "boundaries.md"), String)
    @test occursin("3GPP", bounds)
    @test occursin("O-RAN", bounds)
    @test occursin("38.901", bounds)
    @test occursin("CDL", bounds)
    @test occursin("调度", bounds)
end
