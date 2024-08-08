# Working:
function statics_assembly(mesh::Mesh, materials::Vector{Material})
    
    assembly = Assembly()

    # DEFINING PARAMETERS
    # mesh/material/other constants definition
    nelem = mesh.nelem
    nnel = size(mesh.IEN,1)
    
    # Define integration rules for non-sing
    rules = define_rules(mesh)

    # Matrices initialization
    max_DOF = maximum(mesh.LM)
    assembly.H = zeros(max_DOF,max_DOF)
    assembly.G = zeros(max_DOF,max_DOF)

    # variable to hold information about integration and indices

    # Loop over elements

    p = Progress(nelem,1, "Computing static G and H...", 50)

    for body in mesh.bodies
        # Threads.@threads for e in body
        for e in body

            nodesidx = mesh.IEN[:,e]
            nodes = mesh.nodes[nodesidx,2:end]

            pointsidx = mesh.IEN_geo[:,abs(mesh.EEN[e])]

            if sign(mesh.EEN[e]) == -1
                points = pointsidx[end:-1:1]
            end
            
            points = mesh.points[pointsidx,2:end]

            # Loop over source node
            for se in body
                for n in 1:nnel
                    s = mesh.IEN[n,se]

                    source = mesh.nodes[s,2:end]

                    # Calculating matrix entries
                    # Non-singular integration
                    if s ∉ nodesidx
                        # Integrate
                        HELEM, GELEM = integrate_nonsing(source, points, rules, materials[1])
                        # @infiltrate
                    else
                    # Singular integration
                        # Local index of source node on the element
                        idx = findfirst(nodesidx.==s)
                        # Integrate
                        HELEM, GELEM = integrate_sing(source, points, rules, materials[1], idx)
                    end

                    #Assembly on matrix
                    assembly.H[mesh.LM[3*(n-1)+1:3*n,e], mesh.LM[:,e]] = HELEM
                    assembly.G[mesh.LM[3*(n-1)+1:3*n,e], mesh.LM[:,e]] = GELEM

                end
                # End loop over elements
            end
            next!(p)

        # End loop over source
        end
    end

    # RIGID BODY MOTION STRATEGY
    integrate_rigid_body!(assembly.H,mesh)

    return assembly

end


function dynamics_assembly!(mesh,problem,materials,assembly,freq)

    # DEFINING PARAMETERS
    # mesh/material/other constants definition
    nelem = mesh.nelem
    
    # Define integration constants
    C_dyn = cmplx_consts(materials[1],freq)

    # Define integration rules for non-sing
    rules = define_rules(mesh)

    # Matrices initialization
    max_DOF = maximum(mesh.ID)
    assembly.zH = zeros(ComplexF64,max_DOF,max_DOF)
    assembly.zG = zeros(ComplexF64,max_DOF,max_DOF)

    # Loop over elements

    p = Progress(nelem,1, "Computing static G and H...", 50)

    Threads.@threads for e in 1:nelem
    # for e in 1:mesh.nelem

        nodesidx = mesh.IEN[:,e]
        points = mesh.points[mesh.IEN_geo[:,e],2:end]
        nodes = mesh.nodes[nodesidx,2:end]

        # Loop over source node
        for s in 1:mesh.nnodes

            source = mesh.nodes[s,2:end]

            # Calculating matrix entries
            # Non-singular integration
            if s ∉ nodesidx
                # Integrate
                zHELEM, zGELEM = integrate_nonsing_dyn(source, points, rules, C_dyn)
            else
            # Singular integration
                # Local index of source node on the element
                idx = findfirst(nodesidx.==s)
                # Integrate
                zHELEM, zGELEM = integrate_sing_dyn(source, points, rules, materials[1],C_dyn, idx)
            end

            #Assembly on matrix
            assembly.zH[mesh.ID[:,s], mesh.LM[:,e]] = zHELEM
            assembly.zG[mesh.ID[:,s], mesh.LM[:,e]] = zGELEM

        # End loop over elements
        end
        next!(p)

    # End loop over source
    end

    # Sum static part again
    for n in 1:mesh.nnodes
        assembly.zH[mesh.ID[:,n],mesh.ID[:,n]] += complex.(assembly.H[mesh.ID[:,n],mesh.ID[:,n]])
    end


end