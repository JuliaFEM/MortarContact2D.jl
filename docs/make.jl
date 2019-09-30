# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2D.jl/blob/master/LICENSE

using Documenter, MortarContact2D

makedocs(modules=[MortarContact2D],
         format = Documenter.HTML(),
         checkdocs = :all,
         sitename = "MortarContact2D.jl",
         pages = ["index.md"]
        )
