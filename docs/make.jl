push!(LOAD_PATH,"../src/")

using Documenter, Oxygen

makedocs(
    sitename="Oxygen.jl",
    authors="Rushil Mallarapu",
    modules = [Oxygen],
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Home" => "index.md",
        "API" => [
            "api/atoms.md",
            "api/mols.md",
            "api/smiles.md",
            "api/molops.md",
        ],
        "TODO" => "todo.md"
    ]
)

deploydocs(
    repo = "github.com/sudo-rushil/Oxygen.jl.git",
)
