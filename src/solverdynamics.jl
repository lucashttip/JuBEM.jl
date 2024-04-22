function calc_GH_dynamic!(mesh::Mesh, material::Vector{Material}, solver_var::Svar, frequency::Float64)


    # DEFINING PARAMETERS
    # mesh/material/other constants definition
    nelem = mesh.nelem
    delta = I(3)
    nnel = (mesh.eltype+1)^2
    C_stat = calc_static_constants(material[1])
    zconsts = cmplx_consts(material[1],frequency)


    # normal integration constants
    csis_cont, csis_descont, rules  = calc_nonsing_consts(mesh)
    N, dNc, dNe, Nd = calc_N_nonsing(csis_cont, csis_descont, rules)

    # singular integration definitions
    N_sing, dNc_sing, dNe_sing, Nd_sing, Jb_sing, omega_sing = calc_N_sing(mesh, solver_var, csis_cont, csis_descont, rules.npg_sing)

    #matrices initialization
    max_GL = maximum(mesh.ID)
    solver_var.zH = zeros(ComplexF64,max_GL,max_GL)
    solver_var.zG = zeros(ComplexF64,max_GL,max_GL)

    p = Progress(nelem,1, "Computing dynamic zG and zH...", 50)


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
                    r, c, d = calc_dist2(source_node, field_points, rules.dists,csis_cont)
                    if r > 0
                        # NORMAL INTEGRATION
                        zHELEM, zGELEM = integrate_nonsing_dynamic(source_node,gauss_points[r],Nd[r],normal[r],J[r], rules.gp[r].omega, delta, zconsts)
                    else
                        # NEAR INTEGRATION
                        @infiltrate
                        gauss_points_near, Nd_near, J_near, normal_near, weights_near = calc_nearpoints(rules.gp_near.csi, rules.gp_near.omega,c, d, csis_cont, csis_descont, field_points)
                        zHELEM, zGELEM = integrate_nonsing_dynamic(source_node,gauss_points_near,Nd_near,normal_near,J_near, weights_near, delta, zconsts)
                    end

                # SINGULAR INTEGRATION
                else
                    normal_sing, gauss_points_sing, weights_sing = calc_points_sing(N_sing,dNc_sing, dNe_sing,Nd_sing, field_points, Jb_sing, nnel, n, omega_sing)

                    zHELEM1, zGELEM = integrate_sing_dynamic(source_node, gauss_points_sing, normal_sing, delta, zconsts, weights_sing, n)
                    zHELEM2 = integrate_sing_dynamic2(source_node,gauss_points_sing,normal_sing, delta, zconsts,C_stat,weights_sing,n)

                    zHELEM = zHELEM1 + zHELEM2                
                end

                # INSERT INTO GLOBAL MATRICES
                solver_var.zH[mesh.ID[:,sn], mesh.LM[:,fe]] = zHELEM
                solver_var.zG[mesh.ID[:,sn], mesh.LM[:,fe]] = zGELEM
            end

        end
        next!(p)
    end

    for n in 1:mesh.nnodes
        solver_var.zH[mesh.ID[:,n],mesh.ID[:,n]] += solver_var.H[mesh.ID[:,n],mesh.ID[:,n]]
    end

    # @infiltrate

    return solver_var

end



function calc_GH_dynamic_gpu!(mesh::Mesh, material::Vector{Material}, solver_var::Svar, frequency::Float64)


    # DEFINING PARAMETERS
    # mesh/material/other constants definition
    nelem = mesh.nelem
    delta = I(3)
    nnel = (mesh.eltype+1)^2
    C_stat = calc_static_constants(material[1])
    zconsts = cmplx_consts(material[1],frequency)


    # normal integration constants
    csis_cont, csis_descont, rules  = calc_nonsing_consts(mesh)
    N, dNc, dNe, Nd = calc_N_nonsing(csis_cont, csis_descont, rules)

    # singular integration definitions
    N_sing, dNc_sing, dNe_sing, Nd_sing, Jb_sing, omega_sing = calc_N_sing(mesh, solver_var, csis_cont, csis_descont, rules.npg_sing)

    #matrices initialization
    max_GL = maximum(mesh.ID)
    solver_var.zH = zeros(ComplexF64,max_GL,max_GL)
    solver_var.zG = zeros(ComplexF64,max_GL,max_GL)

    p = Progress(nelem,1, "Computing dynamic zG and zH...", 50)


    # FIELD ELEMENT LOOP
    Threads.@threads for fe in 1:nelem
    # for fe in 1:nelem

        # FIELD ELEMENT PARAMETER FOR INTEGRATION
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]
        normal, J, gauss_points = calc_nJgp(N, dNc, dNe, rules.gp, field_points)
        
        fns = mesh.nodes[mesh.IEN[:,fe],2:end]
        # DEFINITION OF NON-SINGULAR NODES
        nonsing_nodes = setdiff(1:mesh.nnodes,fns)

        for sn in nonsing_nodes
            source_node = mesh.nodes[sn,2:end]

            # NON-SINGULAR INTEGRATION
            r, c, d = calc_dist2(source_node, field_points, rules.dists,csis_cont)
            if r > 0
                # NORMAL INTEGRATION
                zHELEM, zGELEM = integrate_nonsing_dynamic(source_node,gauss_points[r],Nd[r],normal[r],J[r], rules.gp[r].omega, delta, zconsts)
            else
                # NEAR INTEGRATION
                gauss_points_near, Nd_near, J_near, normal_near, weights_near = calc_nearpoints(rules.gp_near.csi, rules.gp_near.omega,c, d, csis_cont, csis_descont, field_points)
                zHELEM, zGELEM = integrate_nonsing_dynamic(source_node,gauss_points_near,Nd_near,normal_near,J_near, weights_near, delta, zconsts)
            end

            # INSERT INTO GLOBAL MATRICES
            solver_var.zH[mesh.ID[:,sn], mesh.LM[:,fe]] = zHELEM
            solver_var.zG[mesh.ID[:,sn], mesh.LM[:,fe]] = zGELEM
        end

        # SINGULAR INTEGRATION
        for n in 1:4
            sn = fns[n]
            source_node = mesh.nodes[sn,2:end]

            normal_sing, gauss_points_sing, weights_sing = calc_points_sing(N_sing,dNc_sing, dNe_sing,Nd_sing, field_points, Jb_sing, nnel, n, omega_sing)

            zHELEM1, zGELEM = integrate_sing_dynamic(source_node, gauss_points_sing, normal_sing, delta, zconsts, weights_sing, n)
            zHELEM2 = integrate_sing_dynamic2(source_node,gauss_points_sing,normal_sing, delta, zconsts,C_stat,weights_sing,n)

            zHELEM = zHELEM1 + zHELEM2       
            
            # INSERT INTO GLOBAL MATRICES
            solver_var.zH[mesh.ID[:,sn], mesh.LM[:,fe]] = zHELEM
            solver_var.zG[mesh.ID[:,sn], mesh.LM[:,fe]] = zGELEM
        end

        next!(p)
    end

    for n in 1:mesh.nnodes
        solver_var.zH[mesh.ID[:,n],mesh.ID[:,n]] += solver_var.H[mesh.ID[:,n],mesh.ID[:,n]]
    end

    return solver_var

end