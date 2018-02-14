# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2D.jl/blob/master/LICENSE

"""
    2D mortar contact mechanics for JuliaFEM.

# Problem types

- MeshTie2D
- Contact2DSS
- Contact2DFSAD

"""
module MortarContact2D

using FEMBase
using Mortar2D: calculate_mortar_assembly

import FEMBase: get_unknown_field_name, assemble_elements!, add_elements!

type MeshTie2D <: BoundaryProblem
    master_elements :: Vector{Element}
end

function MeshTie2D()
    return MeshTie2D([])
end

function add_elements!(::Problem{MeshTie2D}, ::Any)
    error("use `add_slave_elements!` and `add_master_elements!` to add ",
          "elements to the problem.")
end

function add_slave_elements!(problem::Problem{MeshTie2D}, elements)
    for element in elements
        push!(problem.elements, element)
    end
end

function add_master_elements!(problem::Problem{MeshTie2D}, elements)
    for element in elements
        push!(problem.properties.master_elements, element)
    end
end

function get_slave_elements(problem::Problem{MeshTie2D})
    return problem.elements
end

function get_master_elements(problem::Problem{MeshTie2D})
    return problem.properties.master_elements
end

function get_slave_dofs(problem::Problem{MeshTie2D})
    dofs = Int64[]
    for element in get_slave_elements(problem)
        append!(dofs, get_gdofs(problem, element))
    end
    return sort(unique(dofs))
end

function get_master_dofs(problem::Problem{MeshTie2D})
    dofs = Int64[]
    for element in get_master_elements(problem)
        append!(dofs, get_gdofs(problem, element))
    end
    return sort(unique(dofs))
end

# This is already implemented in Mortar2D.jl
function assemble_elements!(problem::Problem{MeshTie2D}, assembly::Assembly,
                            slave_elements::Vector{Element{Seg2}}, time::Float64)

    coords = Dict()
    element_conn = Dict()
    element_types = Dict()
    slave_element_ids = []
    master_element_ids = []
    master_elements = get_master_elements(problem)

    function collect_info!(elements, element_id_list)
        for element in elements
            X_el = element("geometry", time)
            elcon = get_connectivity(element)
            for (ci, Xi) in zip(elcon, X_el) coords[ci] = Xi end
            while element.id == -1 element.id = rand(Int) end
            element_types[element.id] = :Seg2
            element_conn[element.id] = elcon
            push!(element_id_list, element.id)
        end
    end
    collect_info!(slave_elements, slave_element_ids)
    collect_info!(master_elements, master_element_ids)

    s, m, D, M = calculate_mortar_assembly(element_conn, element_types, coords,
                                           slave_element_ids, master_element_ids)

    append!(assembly.C1, convert(SparseMatrixCOO, D'))
    append!(assembly.C1, convert(SparseMatrixCOO, M'))
    append!(assembly.C2, convert(SparseMatrixCOO, D))
    append!(assembly.C2, convert(SparseMatrixCOO, M))

    return nothing
end

function get_projection_matrix(problem::Problem{MeshTie2D})
    m = get_master_dofs(problem)
    s = get_slave_dofs(problem)
    C2 = sparse(problem.assembly.C2)
    D_ = C2[s,s]
    M_ = C2[s,m]
    P = ldltfact(1/2*(D_ + D_')) \ M_
    return P
end

export MeshTie2D
export add_slave_elements!, add_master_elements!
export get_slave_elements, get_master_elements
export get_slave_dofs, get_master_dofs
export get_projection_matrix

end
