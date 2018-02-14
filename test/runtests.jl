# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2D.jl/blob/master/LICENSE

using FEMBase
using FEMBase.Test
using MortarContact2D

@testset "MortarContact2D.jl" begin
    include("test_meshtie2d.jl")
end
