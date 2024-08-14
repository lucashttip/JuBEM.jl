# Working:
function statics_assembly(mesh::Mesh, materials::Vector{Material},problem::Problem)
    
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

    # for body in mesh.bodies
    for body in mesh.bodies
        Threads.@threads for e in body
            tag = mesh.tag[e]
            matidx = problem.taginfo[tag,3]
            nodesidx = mesh.IEN[:,e]
            nodes = mesh.nodes[nodesidx,2:end]

            pointsidx = mesh.IEN_geo[:,abs(mesh.EEN[e])]
            points = mesh.points[pointsidx,2:end]
            normalsign = sign(mesh.EEN[e])            

            # Loop over source node
            for se in body
                for n in 1:nnel
                    s = mesh.IEN[n,se]

                    source = mesh.nodes[s,2:end]

                    # Calculating matrix entries
                    # Non-singular integration
                    if s ∉ nodesidx
                        # Integrate
                        HELEM, GELEM = integrate_nonsing(source, points, rules, materials[matidx],normalsign)
                        # @infiltrate
                    else
                    # Singular integration
                        # Local index of source node on the element
                        # Integrate
                        HELEM, GELEM = integrate_sing(source, points, rules, materials[matidx],normalsign, n)
                    end

                    #Assembly on matrix
                    assembly.H[mesh.LM[3*(n-1)+1:3*n,se], mesh.LM[:,e]] = HELEM
                    assembly.G[mesh.LM[3*(n-1)+1:3*n,se], mesh.LM[:,e]] = GELEM

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


function dynamics_assembly!(mesh::Mesh,problem::Problem,materials::Vector{Material},assembly::Assembly,freq)

    # DEFINING PARAMETERS
    # mesh/material/other constants definition
    nelem = mesh.nelem
    nnel = size(mesh.IEN,1)
    
    # Define integration constants
    C_dyn = cmplx_consts[]
    for m in eachindex(materials)
        push!(C_dyn,cmplx_consts(materials[m],freq))
    end

    # Define integration rules for non-sing
    rules = define_rules(mesh)

    # Matrices initialization
    max_DOF = maximum(mesh.LM)
    assembly.zH = zeros(ComplexF64,max_DOF,max_DOF)
    assembly.zG = zeros(ComplexF64,max_DOF,max_DOF)

    # Loop over elements

    p = Progress(nelem,1, "Computing static G and H...", 50)

    for body in mesh.bodies
        Threads.@threads for e in body
        # for e in body
            tag = mesh.tag[e]
            matidx = problem.taginfo[tag,3]
            nodesidx = mesh.IEN[:,e]
            nodes = mesh.nodes[nodesidx,2:end]

            pointsidx = mesh.IEN_geo[:,abs(mesh.EEN[e])]
            points = mesh.points[pointsidx,2:end]
            normalsign = sign(mesh.EEN[e])            

            # Loop over source node
            for se in body
                for n in 1:nnel
                    s = mesh.IEN[n,se]

                    source = mesh.nodes[s,2:end]

                    # Calculating matrix entries
                    # Non-singular integration
                    if s ∉ nodesidx
                        # Integrate
                        zHELEM, zGELEM = integrate_nonsing_dyn(source, points, rules, C_dyn[matidx],normalsign)
                    else
                    # Singular integration
                        # Local index of source node on the element
                        idx = findfirst(nodesidx.==s)
                        # Integrate
                        zHELEM, zGELEM = integrate_sing_dyn(source, points, rules, materials[matidx],C_dyn[matidx],normalsign, idx)
                    end

                    #Assembly on matrix
                    assembly.zH[mesh.LM[3*(n-1)+1:3*n,se], mesh.LM[:,e]] = zHELEM
                    assembly.zG[mesh.LM[3*(n-1)+1:3*n,se], mesh.LM[:,e]] = zGELEM

                end
            end
            next!(p)

        # End loop over source
        end
    end

    # Sum static part again
    for e in 1:mesh.nelem
        for n in 1:nnel
            idxs = mesh.LM[3*(n-1)+1:3*n,e]
            assembly.zH[idxs,idxs] += complex.(assembly.H[idxs,idxs])
        end
    end


end