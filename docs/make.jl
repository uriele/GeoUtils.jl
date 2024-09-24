using GeoUtils
using Documenter

DocMeta.setdocmeta!(GeoUtils, :DocTestSetup, :(using GeoUtils); recursive=true)

makedocs(;
    modules=[GeoUtils],
    authors="uriele <menarini.marco@gmail.com> and contributors",
    sitename="GeoUtils.jl",
    format=Documenter.HTML(;
        canonical="https://isaccnrbo.github.io/GeoUtils.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/isaccnrbo/GeoUtils.jl",
    devbranch="master",
)
