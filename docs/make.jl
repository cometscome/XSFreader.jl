using XSFreader
using Documenter

DocMeta.setdocmeta!(XSFreader, :DocTestSetup, :(using XSFreader); recursive=true)

makedocs(;
    modules=[XSFreader],
    authors="Yuki Nagai <cometscome@gmail.com> and contributors",
    repo="https://github.com/cometscome/XSFreader.jl/blob/{commit}{path}#{line}",
    sitename="XSFreader.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cometscome.github.io/XSFreader.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cometscome/XSFreader.jl",
    devbranch="main",
)
