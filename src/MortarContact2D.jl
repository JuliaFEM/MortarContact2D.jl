# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2D.jl/blob/master/LICENSE

"""
    2D mortar contact mechanics for JuliaFEM.

# Problem types

- `Mortar2D` to couple non-conforming meshes (mesh tying).
- `MortarContact2D` to add contact constraints.

"""
module MortarContact2D

using Reexport
@reexport using FEMBase

include("mortar2d.jl")

export Mortar2D
export add_slave_elements!, add_master_elements!
export get_slave_elements, get_master_elements
export get_slave_dofs, get_master_dofs
export get_mortar_matrix_D, get_mortar_matrix_M, get_mortar_matrix_P

include("contact2d.jl")

export Contact2D

end
