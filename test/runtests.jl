# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2D.jl/blob/master/LICENSE

using MortarContact2D
using Base.Test

@testset "test MortarContact2D.jl" begin

    @testset "Projecting vertices between surfaces" begin
        include("test_mortar2d_calculate_projection.jl")
    end

    @testset "Mortar coupling, example 1" begin
        include("test_mortar2d_ex1.jl")
    end

    @testset "Mortar coupling, example 2" begin
        include("test_mortar2d_ex2.jl")
    end

    @testset "Contact basic test" begin
        include("test_contact2d.jl")
    end

end
