using Ellasy
using Documenter

DocMeta.setdocmeta!(Ellasy, :DocTestSetup, :(using Ellasy); recursive=true)

makedocs(;
    modules=[Ellasy],
    authors="Stephan Scholz",
    sitename="Ellasy.jl",
    format=Documenter.HTML(;
        canonical="https://stephans3.github.io/Ellasy.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/stephans3/Ellasy.jl",
    devbranch="main",
)
