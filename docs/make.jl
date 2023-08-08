using CARS
using Documenter

DocMeta.setdocmeta!(CARS, :DocTestSetup, :(using CARS); recursive=true)

makedocs(;
    modules=[CARS],
    authors="Max Greifenstein, TU Darmstadt RSM",
    repo="https://github.com/tudarsm/CARS.jl/blob/{commit}{path}#{line}",
    sitename="CARS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tudarsm.github.io/CARS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Example Usage" => "examples.md",
        "How to..." => "howto.md",
        "Function Reference" => "allfunctions.md",
    ],
)

deploydocs(; repo = "github.com/tudarsm/CARS.jl")
