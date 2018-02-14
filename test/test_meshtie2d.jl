# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2D.jl/blob/master/LICENSE

using FEMBase
using FEMBase.Test
using MortarContact2D

@testset "MeshTie2D" begin
    Xs = Dict(1 => [0.0, 1.0], 2 => [5/4, 1.0], 3 => [2.0, 1.0])
    Xm = Dict(4 => [0.0, 1.0], 5 => [1.0, 1.0], 6 => [2.0, 1.0])
    X = merge(Xm , Xs)
    es1 = Element(Seg2, [1, 2])
    es2 = Element(Seg2, [2, 3])
    em1 = Element(Seg2, [4, 5])
    em2 = Element(Seg2, [5, 6])
    elements = [es1, es2, em1, em2]
    update!(elements, "geometry", X)
    problem = Problem(MeshTie2D, "test interface", 1, "u")
    @test_throws ErrorException add_elements!(problem, [es1, es2])
    add_slave_elements!(problem, [es1, es2])
    add_master_elements!(problem, [em1, em2])
    assemble!(problem, 0.0)
    s = get_slave_dofs(problem)
    m = get_master_dofs(problem)
    @test s == [1, 2, 3]
    @test m == [4, 5, 6]
    P = get_projection_matrix(problem)
    u = zeros(6)
    u[m] = ones(3)
    u[s] = P*u[m]
    @test isapprox(u[s], u[m])
end
