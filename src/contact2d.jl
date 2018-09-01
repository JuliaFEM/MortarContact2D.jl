# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/MortarContact2D.jl/blob/master/LICENSE

mutable struct Contact2D <: BoundaryProblem
    master_elements :: Vector{Element}
    rotate_normals :: Bool
    dual_basis :: Bool
    max_distance :: Float64
    iteration :: Int
    contact_state_in_first_iteration :: Symbol
    always_active :: Set{Int64}
    always_inactive :: Set{Int64}
    always_stick :: Set{Int64}
    always_slip :: Set{Int64}
    active_nodes :: Set{Int64}
    inactive_nodes :: Set{Int64}
    stick_nodes :: Set{Int64}
    slip_nodes :: Set{Int64}
    store_fields :: Vector{String}
    use_scaling :: Bool
    alpha :: Dict{Int64, Float64}
    beta :: Dict{Int64, Float64}
    nadj_nodes :: Dict{Int64, Float64}
    scaling_factors :: Dict{Int64, Float64}
end

function Contact2D()
    return Contact2D([], false, true, Inf, 0, :AUTO,
                    Set(), Set(), Set(), Set(),
                    Set(), Set(), Set(), Set(),
                    [], true, Dict(), Dict(), Dict(), Dict())
end

function FEMBase.add_elements!(::Problem{Contact2D}, ::Any)
    error("use `add_slave_elements!` and `add_master_elements!` to add ",
          "elements to the Contact2D problem.")
end

function FEMBase.get_elements(problem::Problem{Contact2D})
    return [problem.elements; problem.properties.master_elements]
end

function FEMBase.add_slave_elements!(problem::Problem{Contact2D}, elements)
    for element in elements
        push!(problem.elements, element)
    end
end

function FEMBase.add_master_elements!(problem::Problem{Contact2D}, elements)
    for element in elements
        push!(problem.properties.master_elements, element)
    end
end

function FEMBase.get_slave_elements(problem::Problem{Contact2D})
    return problem.elements
end

function FEMBase.get_master_elements(problem::Problem{Contact2D})
    return problem.properties.master_elements
end

function get_slave_dofs(problem::Problem{Contact2D})
    dofs = Int64[]
    for element in get_slave_elements(problem)
        append!(dofs, get_gdofs(problem, element))
    end
    return sort(unique(dofs))
end

function get_slave_nodes(problem::Problem{Contact2D})
    slave_nodes = Set{Int64}()
    for element in get_slave_elements(problem)
        for j in get_connectivity(element)
            push!(slave_nodes, j)
        end
    end
    return slave_nodes
end

function get_master_dofs(problem::Problem{Contact2D})
    dofs = Int64[]
    for element in get_master_elements(problem)
        append!(dofs, get_gdofs(problem, element))
    end
    return sort(unique(dofs))
end

function create_rotation_matrix(element::Element{Seg2}, time::Float64)
    n = element("normal", time)
    R = [0.0 -1.0; 1.0 0.0]
    t1 = R'*n[1]
    t2 = R'*n[2]
    Q1 = [n[1] t1]
    Q2 = [n[2] t2]
    Z = zeros(2, 2)
    Q = [Q1 Z; Z Q2]
    return Q
end

function create_contact_segmentation(problem::Problem{Contact2D}, slave_element::Element{Seg2}, master_elements::Vector, time::Float64; deformed=false)

    result = []

    max_distance = problem.properties.max_distance

    x1 = slave_element("geometry", time)

    if deformed && haskey(slave_element, "displacement")
        x1 = map(+, x1, slave_element("displacement", time))
    end

    slave_midpoint = 1/2*(x1[1] + x1[2])
    slave_length = norm(x1[2]-x1[1])

    for master_element in master_elements

        x2 = master_element("geometry", time)

        if deformed && haskey(master_element, "displacement")
            x2 = map(+, x2, master_element("displacement", time))
        end

        master_midpoint = 1/2*(x2[1] + x2[2])
        master_length = norm(x2[2]-x2[1])
        # charasteristic length
        dist = norm(slave_midpoint - master_midpoint)
        cl = dist/max(slave_length, master_length)
        if cl > max_distance
            continue
        end

        # 3.1 calculate segmentation
        x21, x22 = x2
        xi1a = project_from_master_to_slave(slave_element, x21, time)
        xi1b = project_from_master_to_slave(slave_element, x22, time)
        xi1 = clamp.([xi1a; xi1b], -1.0, 1.0)
        l = 1/2*abs(xi1[2]-xi1[1])
        if isapprox(l, 0.0)
            continue # no contribution in this master element
        end
        push!(result, (master_element, xi1, l))
    end
    return result
end

function FEMBase.assemble_elements!(problem::Problem{Contact2D}, assembly::Assembly,
                                    elements::Vector{Element{Seg2}}, time::Float64)

    props = problem.properties
    props.iteration += 1

    slave_elements = get_slave_elements(problem)

    field_dim = get_unknown_field_dimension(problem)
    field_name = get_parent_field_name(problem)

    # 1. calculate nodal normals for slave element nodes j ∈ S
    normals = calculate_normals([slave_elements...], time; rotate_normals=props.rotate_normals)
    tangents = Dict(j => [t2, -t1] for (j, (t1,t2)) in normals)
    update!(slave_elements, "normal", time => normals)
    update!(slave_elements, "tangent", time => tangents)

    Rn = 0.0

    alpha = empty!(problem.properties.alpha)
    beta = empty!(problem.properties.beta)
    nadj_nodes = empty!(problem.properties.nadj_nodes)

    # starting situation is that all slave nodes in contact interface are active
    # except nodes which are marked to be "always inactive"
    S = get_slave_nodes(problem) # slave element nodes
    props.active_nodes = setdiff(S, props.always_inactive)

    # 2. loop all slave elements
    for slave_element in slave_elements
        slave_connectivity = get_connectivity(slave_element)
        if all(j in props.always_inactive for j in slave_connectivity)
            continue
        end
        sdofs = get_gdofs(problem, slave_element)
        nsl = length(slave_element)
        X1 = slave_element("geometry", time)
        x1 = slave_element("geometry", time)
        if haskey(slave_element, "displacement")
            u1 = slave_element("displacement", time)
            x1 = map(+, x1, u1)
        end
        la1 = slave_element("lambda", time)
        n1 = slave_element("normal", time)
        t1 = slave_element("tangent", time)

        contact_area = 0.0
        contact_error = 0.0
        Q2 = create_rotation_matrix(slave_element, time)

        master_elements = get_master_elements(problem)
        segmentation = create_contact_segmentation(problem, slave_element, master_elements, time)
        if length(segmentation) == 0 # no overlapping in master and slave surfaces with this slave element
            setdiff!(props.active_nodes, get_connectivity(slave_element))
            continue
        end

        ae = zeros(nsl)
        be = zeros(nsl)

        for ip in get_integration_points(slave_element, 3)
            detJ = slave_element(ip, time, Val{:detJ})
            w = ip.weight * detJ
            N1 = slave_element(ip, time)
            ae += w*vec(N1)/detJ
        end

        Ae = eye(nsl)
        if props.dual_basis
            De = zeros(nsl, nsl)
            Me = zeros(nsl, nsl)
            for (master_element, xi1, l) in segmentation
                for ip in get_integration_points(slave_element, 3)
                    detJ = slave_element(ip, time, Val{:detJ})
                    w = ip.weight*detJ*l
                    xi = ip.coords[1]
                    xi_s = dot([1/2*(1-xi); 1/2*(1+xi)], xi1)
                    N1 = vec(get_basis(slave_element, xi_s, time))
                    De += w*diagm(N1)
                    Me += w*N1*N1'
                    be += w*vec(N1)/detJ
                end
            end
            Ae = De*inv(Me)
            update!(slave_element, "dual basis coefficients", time => Ae)
        end

        for (i, j) in enumerate(get_connectivity(slave_element))
            alpha[j] = get(alpha, j, 0.0) + ae[i]
            beta[j] = get(beta, j, 0.0) + be[i]
            nadj_nodes[j] = get(nadj_nodes, j, 0) + 1
        end

        # loop all segments
        for (master_element, xi1, l) in segmentation

            nm = length(master_element)
            X2 = master_element("geometry", time)
            x2 = master_element("geometry", time)
            if haskey(master_element, "displacement")
                u2 = master_element("displacement", time)
                x2 = map(+, x2, u2)
            end

            # 3.3. loop integration points of one integration segment and calculate
            # local mortar matrices
            De = zeros(nsl, nsl)
            Me = zeros(nsl, nsl)
            Ne = zeros(nsl, 2*nsl)
            Te = zeros(nsl, 2*nsl)
            He = zeros(nsl, 2*nsl)
            ce = zeros(nsl)
            ge = zeros(nsl)
            for ip in get_integration_points(slave_element, 3)
                detJ = slave_element(ip, time, Val{:detJ})
                w = ip.weight*detJ*l
                xi = ip.coords[1]
                xi_s = dot([1/2*(1-xi); 1/2*(1+xi)], xi1)
                N1 = vec(get_basis(slave_element, xi_s, time))
                Phi = Ae*N1

                # project gauss point from slave element to master element in direction n_s
                X_s = interpolate(N1, X1) # coordinate in gauss point
                n_s = interpolate(N1, n1) # normal direction in gauss point
                t_s = interpolate(N1, t1) # tangent condition in gauss point
                n_s /= norm(n_s)
                t_s /= norm(t_s)
                xi_m = project_from_slave_to_master(master_element, X_s, n_s, time)
                N2 = vec(get_basis(master_element, xi_m, time))
                X_m = interpolate(N2, X2)

                u_s = interpolate(N1, u1)
                u_m = interpolate(N2, u2)
                x_s = map(+, X_s, u_s)
                x_m = map(+, X_m, u_m)
                la_s = interpolate(Phi, la1)

                # contact virtual work
                De += w*Phi*N1'
                Me += w*Phi*N2'

                # contact constraints
                Ne += w*reshape(kron(N1, n_s, Phi), 2, 4)
                Te += w*reshape(kron(N2, n_s, Phi), 2, 4)
                He += w*reshape(kron(N1, t_s, Phi), 2, 4)
                ge += w*Phi*dot(n_s, x_m-x_s)
                Rn += w*dot(n_s, la_s)

                contact_area += w
                contact_error += 1/2*w*dot(n_s, x_s-x_m)^2
            end

            sdofs = get_gdofs(problem, slave_element)
            mdofs = get_gdofs(problem, master_element)

            # add contribution to contact virtual work
            for i=1:field_dim
                lsdofs = sdofs[i:field_dim:end]
                lmdofs = mdofs[i:field_dim:end]
                add!(problem.assembly.C1, lsdofs, lsdofs, De)
                add!(problem.assembly.C1, lsdofs, lmdofs, -Me)
            end

            # add contribution to contact constraints
            add!(assembly.C2, sdofs[1:field_dim:end], sdofs, Ne)
            add!(assembly.C2, sdofs[1:field_dim:end], mdofs, -Te)
            add!(assembly.D, sdofs[2:field_dim:end], sdofs, He)
            add!(assembly.g, sdofs[1:field_dim:end], ge)

        end # master elements done

        if "contact area" in props.store_fields
            update!(slave_element, "contact area", time => contact_area)
        end

        if "contact error" in props.store_fields
            update!(slave_element, "contact error", time => contact_error)
        end

    end # slave elements done, contact virtual work ready

    # Next, determine which nodes are active, which are inactive, which are
    # sticking and which are slipping.

    union!(props.active_nodes, props.always_active)
    if isempty(props.active_nodes) # No active nodes at all, we're done!
        empty!(problem.assembly)
        return nothing
    end

    C1 = sparse(assembly.C1)
    ndofs = maximum(size(C1))
    C2 = sparse(assembly.C2, ndofs, ndofs)
    D = sparse(assembly.D, ndofs, ndofs)
    g = full(assembly.g, ndofs, 1)

    # 1. if we are in first iteration, check that do we have planar surface.
    # If have, calculate average and mean gap and set contact active if
    # surfaces are "close enough". This helps with convergence issues if one
    # of the bodies is missing Dirichlet boundary condition.
    gaps = [g[2*(j-1)+1] for j in props.active_nodes]
    avg_gap = mean(gaps)
    std_gap = std(gaps)
    if props.iteration == 1 &&
            props.contact_state_in_first_iteration == :AUTO &&
            !isempty(props.active_nodes) &&
            avg_gap < 1.0e-2 &&
            std_gap < 1.0e-3
        info("Average weighted gap = $avg_gap, std gap = $std_gap.")
        info("Looks that contact surfaces are planar and close each other.")
        info("$(length(props.active_nodes)) nodes active.")
    else # 2. Otherwise, check contact condition based on PDASS
        la = problem("lambda", time)
        for j in props.active_nodes
            lan = dot(normals[j], get(la, j, zeros(2)))
            if lan-g[2*(j-1)+1] < 0
                setdiff!(props.active_nodes, j)
            end
        end
    end

    # 3. Determine tangential direction for the active nodes
    props.slip_nodes = setdiff(props.active_nodes, props.always_stick)

    info("active nodes: $(props.active_nodes)")

    # solve variational inequality

    props.inactive_nodes = setdiff(S, props.active_nodes)
    props.stick_nodes = setdiff(props.active_nodes, props.slip_nodes)

    # remove all inactive nodes from assembly
    for j in setdiff(S, props.active_nodes)
        dofs = [2*(j-1)+1, 2*(j-1)+2]
        info("Node $j is inactive, removing dofs $dofs")
        C1[dofs,:] = 0.0
        C2[dofs,:] = 0.0
        D[dofs,:] = 0.0
        g[dofs,:] = 0.0
    end

    for j in props.slip_nodes
        dofs = [2*(j-1)+1, 2*(j-1)+2]
        info("Node $j is slipping, removing tangential constraint $(dofs[2])")
        C2[dofs[2],:] = 0.0
        g[dofs[2]] = 0.0
        D[dofs[2], dofs] = tangents[j]
    end

    if problem.properties.use_scaling
        scaling = empty!(problem.properties.scaling_factors)
        for j in props.active_nodes
            dofs = [2*(j-1)+1, 2*(j-1)+2]
            isapprox(beta[j], 0.0) && continue
            scaling[j] = alpha[j] / beta[j]
            C1[dofs,:] *= scaling[j]
            C2[dofs[1],:] *= scaling[j]
            g[dofs[1]] *= scaling[j]
        end
    end

    problem.assembly.C1 = C1
    problem.assembly.C2 = C2
    problem.assembly.D = D
    problem.assembly.g = g

end
