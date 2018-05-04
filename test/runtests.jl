# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2D.jl/blob/master/LICENSE

using MortarContact2D
using Base.Test

@testset "Mortar coupling, example 1" begin include("test_mortar2d_ex1.jl") end
