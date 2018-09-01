# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2D.jl/blob/master/LICENSE

using MortarContact2D
using FEMBase.Test

D = [1.0e-5, 0.0]

X = Dict(1 => [-2.0, 0.0], 2 => [ 0.0, 0.0],
         3 => [ 2.0, 0.0], 4 => [-2.0, 0.0],
         5 => [ 0.0, 0.0])

for j in [4, 5]
    X[j] += D
end

contact = Problem(Contact2D, "contact", 2, "displacement")
contact.properties.contact_state_in_first_iteration = :ACTIVE
element1 = Element(Seg2, [1, 2])
element2 = Element(Seg2, [2, 3])
element3 = Element(Seg2, [5, 4])
add_slave_elements!(contact, [element1, element2])
add_master_elements!(contact, [element3])
all_elements = [element1, element2, element3]
update!(all_elements, "geometry", X)

assemble!(contact, 0.0)

C1 = full(contact.assembly.C1, 6, 10)
C2 = full(contact.assembly.C2, 6, 10)
D = full(contact.assembly.D, 6, 10)
g = full(contact.assembly.g, 6, 1)
@test isapprox(C1[1:6,1:6], diagm([1.0, 1.0, 2.0, 2.0, 1.0, 1.0]))
