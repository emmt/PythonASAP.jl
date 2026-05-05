using PythonASAP
using Documenter

DocMeta.setdocmeta!(PythonASAP, :DocTestSetup, :(using PythonASAP); recursive=true)

makedocs(;
    modules=[PythonASAP],
    authors="Éric Thiébaut <eric.thiebaut@univ-lyon1.fr> and contributors",
    sitename="PythonASAP.jl",
    format=Documenter.HTML(;
        canonical="https://emmt.github.io/PythonASAP.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/emmt/PythonASAP.jl",
    devbranch="main",
)
