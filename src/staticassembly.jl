# Working:
function statics_assembly!(mesh::Mesh, material::Vector{Material}, assembly::Assembly)
    
    # DEFINING PARAMETERS
    # mesh/material/other constants definition
    nelem = mesh.nelem
    delta = I(3)
    nnel = (mesh.eltype+1)^2

    C_stat = calc_static_constants(material[1])

    # normal integration constants
    csis_cont, csis_descont, rules  = calc_nonsing_consts(mesh)
    N, dNc, dNe, Nd = calc_N_nonsing(csis_cont, csis_descont, rules)

    # singular integration definitions
    N_sing, dNc_sing, dNe_sing, Nd_sing, Jb_sing, omega_sing = calc_N_sing(mesh, assembly, csis_cont, csis_descont, rules.npg_sing)

    #matrices initialization
    max_GL = maximum(mesh.ID)
    assembly.H = zeros(max_GL,max_GL)
    assembly.G = zeros(max_GL,max_GL)

    p = Progress(nelem,1, "Computing static G and H...", 50)
    
    # FIELD ELEMENT LOOP
    Threads.@threads for fe in 1:nelem
    # for fe in 1:nelem

        # FIELD ELEMENT PARAMETER FOR INTEGRATION
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]
        normal, J, gauss_points = calc_nJgp(N, dNc, dNe, rules.gp, field_points)
        
        # SOURCE ELEMENT LOOP
        for se in 1:nelem
            # SOURCE NODE LOOP
            for n = 1:nnel
                # SOURCE NODE DEFINITION
                sn = mesh.IEN[n,se]
                source_node = mesh.nodes[sn,2:end]

                # NON-SINGULAR INTEGRATION
                if fe != se
                    r, c, d = calc_dist(source_node, field_points, rules.dists,csis_cont)
                    if r > 0
                        # NORMAL INTEGRATION
                        HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points[r],Nd[r],normal[r],J[r], rules.gp[r].omega, delta, C_stat)
                    else
                        # NEAR INTEGRATION
                        gauss_points_near, Nd_near, J_near, normal_near, weights_near = calc_nearpoints(rules.gp_near.csi, rules.gp_near.omega,c, d, csis_cont, csis_descont, field_points)
                        HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points_near,Nd_near,normal_near,J_near, weights_near, delta, C_stat)
                    end

                # SINGULAR INTEGRATION
                else
                    normal_sing, gauss_points_sing, weights_sing = calc_points_sing(N_sing,dNc_sing, dNe_sing,Nd_sing, field_points, Jb_sing, nnel, n, omega_sing)

                    HELEM, GELEM = integrate_sing_static(source_node, gauss_points_sing, normal_sing, delta, C_stat, weights_sing, n)
                end

                # INSERT INTO GLOBAL MATRICES
                assembly.H[mesh.ID[:,sn], mesh.LM[:,fe]] = HELEM
                assembly.G[mesh.ID[:,sn], mesh.LM[:,fe]] = GELEM
            end
            
        end
        next!(p)
    end

    # RIGID BODY MOTION STRATEGY
    for n in 1:mesh.nnodes
        assembly.H[mesh.ID[:,n],mesh.ID[:,n]] .= 0
    end
    integrate_rigid_body!(assembly.H,mesh)

    return assembly

end

function statics_assembly2!(mesh::Mesh, material::Vector{Material}, assembly::Assembly)
    
    # DEFINING PARAMETERS
    # mesh/material/other constants definition
    nelem = mesh.nelem
    delta = I(3)
    nnel = (mesh.eltype+1)^2
    
    # Define integration constants


    # Loop over elements

    for e in 1:mesh.nelem

        nodesidx = mesh.IEN[:,e]

        # Loop over source node
        for s in 1:mesh.nnodes

            # Calculating matrix entries
            # Non-singular integration
            if s ∉ nodesidx

            end
                
            # Singular integration


            #Assembly on matrix

        # End loop over elements

    # End loop over source







    C_stat = calc_static_constants(material[1])

    # normal integration constants
    csis_cont, csis_descont, rules  = calc_nonsing_consts(mesh)
    N, dNc, dNe, Nd = calc_N_nonsing(csis_cont, csis_descont, rules)

    # singular integration definitions
    N_sing, dNc_sing, dNe_sing, Nd_sing, Jb_sing, omega_sing = calc_N_sing(mesh, assembly, csis_cont, csis_descont, rules.npg_sing)

    #matrices initialization
    max_GL = maximum(mesh.ID)
    assembly.H = zeros(max_GL,max_GL)
    assembly.G = zeros(max_GL,max_GL)

    p = Progress(nelem,1, "Computing static G and H...", 50)
    
    # FIELD ELEMENT LOOP
    Threads.@threads for fe in 1:nelem
    # for fe in 1:nelem

        # FIELD ELEMENT PARAMETER FOR INTEGRATION
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]
        normal, J, gauss_points = calc_nJgp(N, dNc, dNe, rules.gp, field_points)
        
        # SOURCE ELEMENT LOOP
        for se in 1:nelem
            # SOURCE NODE LOOP
            for n = 1:nnel
                # SOURCE NODE DEFINITION
                sn = mesh.IEN[n,se]
                source_node = mesh.nodes[sn,2:end]

                # NON-SINGULAR INTEGRATION
                if fe != se
                    r, c, d = calc_dist(source_node, field_points, rules.dists,csis_cont)
                    if r > 0
                        # NORMAL INTEGRATION
                        HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points[r],Nd[r],normal[r],J[r], rules.gp[r].omega, delta, C_stat)
                    else
                        # NEAR INTEGRATION
                        gauss_points_near, Nd_near, J_near, normal_near, weights_near = calc_nearpoints(rules.gp_near.csi, rules.gp_near.omega,c, d, csis_cont, csis_descont, field_points)
                        HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points_near,Nd_near,normal_near,J_near, weights_near, delta, C_stat)
                    end

                # SINGULAR INTEGRATION
                else
                    normal_sing, gauss_points_sing, weights_sing = calc_points_sing(N_sing,dNc_sing, dNe_sing,Nd_sing, field_points, Jb_sing, nnel, n, omega_sing)

                    HELEM, GELEM = integrate_sing_static(source_node, gauss_points_sing, normal_sing, delta, C_stat, weights_sing, n)
                end

                # INSERT INTO GLOBAL MATRICES
                assembly.H[mesh.ID[:,sn], mesh.LM[:,fe]] = HELEM
                assembly.G[mesh.ID[:,sn], mesh.LM[:,fe]] = GELEM
            end
            
        end
        next!(p)
    end

    # RIGID BODY MOTION STRATEGY
    for n in 1:mesh.nnodes
        assembly.H[mesh.ID[:,n],mesh.ID[:,n]] .= 0
    end
    integrate_rigid_body!(assembly.H,mesh)

    return assembly

end