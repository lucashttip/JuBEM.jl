function calc_GH_static_non_const!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type)
    
    # gp = rules.gp
    nelem = mesh.nelem
    m = 1
    C_stat = calc_static_constants(material[m])
    delta = I(3) 
    nnel = (mesh.eltype+1)^2

    max_GL = maximum(mesh.ID)
    solver_var.H = zeros(max_GL,max_GL)
    solver_var.G = zeros(max_GL,max_GL)

    csis_cont = range(-1,1,length = mesh.eltype+1)
    csis_descont = range(-1+mesh.offset,1 - mesh.offset,length=mesh.eltype+1)
    omegas = calc_omegas(solver_var.omega)

    # Cálculos para elementos não singulares
    gp, dists = calc_points_weights()

    Nc_nonsing = []
    dNcdcsi_nonsing = []
    dNcdeta_nonsing = []
    Nd_nonsing = []

    for i in eachindex(gp)
        N = calc_N_matrix(csis_cont,gp[i].csi)
        Nd = calc_N_matrix(csis_descont,gp[i].csi)
        Ncsi = calc_dNdcsi_matrix(csis_cont,gp[i].csi)
        Neta = calc_dNdeta_matrix(csis_cont,gp[i].csi)
        push!(Nc_nonsing,N)
        push!(dNcdcsi_nonsing,Ncsi)
        push!(dNcdeta_nonsing,Neta)
        push!(Nd_nonsing,Nd)
    end


    # Cálculos para elementos singulares
    csis = calc_csis_grid(solver_var.csi)
    csis_cont_lin = range(-1,1,length = 2)
    Nc_lin = calc_N_matrix(csis_cont_lin,csis)
    dNcdcsi_lin = calc_dNdcsi_matrix(csis_cont_lin,csis)
    dNcdeta_lin = calc_dNdeta_matrix(csis_cont_lin,csis)

    if mesh.eltype == 2
        nperm = 3
    elseif mesh.eltype < 2
        nperm = 1
    end

    csi_sing, Jb_sing = csis_sing(mesh.offset, Nc_lin, dNcdcsi_lin,dNcdeta_lin,mesh.eltype)
    npg_sing = size(csi_sing,1)
    Nc_sing = zeros(npg_sing,nnel,nperm)
    dNcdcsi_sing = zeros(npg_sing,nnel,nperm)
    dNcdeta_sing = zeros(npg_sing,nnel,nperm)
    Nd_sing = zeros(npg_sing,nnel,nperm)

    for i in 1:nperm
        Nc_sing[:,:,i] = calc_N_matrix(csis_cont,csi_sing[:,:,i])
        dNcdcsi_sing[:,:,i] = calc_dNdcsi_matrix(csis_cont,csi_sing[:,:,i])
        dNcdeta_sing[:,:,i] = calc_dNdeta_matrix(csis_cont,csi_sing[:,:,i])
        Nd_sing[:,:,i] = calc_N_matrix(csis_descont,csi_sing[:,:,i])
    end
    omega_sing = repeat(omegas,4)

    csis_telles, omegas_telles =  gausslegendre(12)


    # @infiltrate
    # Field loop:
    # Threads.@threads for fe in 1:nelem
    for fe in 1:nelem
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]
        normal = []
        J = []
        gauss_points = []
        for i in eachindex(gp)
            normal2, J2 = calc_n_J_matrix(dNcdcsi_nonsing[i], dNcdeta_nonsing[i], field_points)
            gauss_points2 = Nc_nonsing[i]*field_points
            push!(normal,normal2)
            push!(J,J2)
            push!(gauss_points,gauss_points2)
        end

        # source loop:
        for se in 1:nelem
            for n = 1:nnel
                sn = mesh.IEN[n,se]
                source_node = mesh.nodes[sn,2:end]

                if fe != se

                    r, c, d = calc_dist(source_node, field_points, dists,csis_cont)
                    if r > 0
                        # normal integration
                        HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points[r],Nd_nonsing[r],normal[r],J[r], gp[r].omega, delta, C_stat)
                    else
                        # near integration
                        gp_local_near, weights_near = pontos_pesos_local_telles(csis_telles, omegas_telles, c[1:2],d)
                        # gp_local_near, weights_near = pontos_pesos_local_subelem(c[1], c[2], Nc_lin, dNcdcsi_lin, dNcdeta_lin, omegas)
                        N_near = calc_N_matrix(csis_cont,gp_local_near)
                        dNc_near = calc_dNdcsi_matrix(csis_cont,gp_local_near)
                        dNe_near = calc_dNdeta_matrix(csis_cont,gp_local_near)
                        Nd_near = calc_N_matrix(csis_descont,gp_local_near)
                        gp_near = N_near*field_points
                        normal_near,J_near = calc_n_J_matrix(dNc_near, dNe_near, field_points)
                        HELEM, GELEM = integrate_nonsing_static(source_node,gp_near,Nd_near,normal_near,J_near, weights_near, delta, C_stat)
                    end
                else

                    idx_ff, np = calc_idx_permutation(nnel,n)

                    normal_sing, J_sing = calc_n_J_matrix(dNcdcsi_sing[:,idx_ff,np], dNcdeta_sing[:,idx_ff,np], field_points)

                    # @infiltrate
                    J_sing = J_sing.*Jb_sing[:,np]
                    pesos = zeros(length(omega_sing),nnel)

                    for kn in 1:nnel
                        pesos[:,kn] = J_sing.*omega_sing.*Nd_sing[:,idx_ff[kn],np]
                    end
                    # @infiltrate
                    gauss_points_sing = Nc_sing[:,idx_ff,np]*field_points
                    # @infiltrate
                    # GELEM1 = integrate_nonsing_static2(source_node,gauss_points[end],Nd_nonsing[end],normal[end],J[end], gp[end].omega, delta, C_stat,n)

                    HELEM, GELEM = integrate_sing_static(source_node, gauss_points_sing, normal_sing, delta, C_stat, pesos, n)
                    # GELEM = GELEM + GELEM1
                    # @infiltrate

                    # HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points,Nd,normal,J, omegas, delta, C_stat)
                end

                solver_var.H[mesh.ID[:,sn], mesh.LM[:,fe]] = HELEM
                solver_var.G[mesh.ID[:,sn], mesh.LM[:,fe]] = GELEM
            end

        end
    
    end

    for n in 1:mesh.nnodes
        solver_var.H[mesh.ID[:,n],mesh.ID[:,n]] .= 0
    end

    integrate_rigid_body!(solver_var.H,mesh)

    
    return solver_var

end

function calc_GH_static_const!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type)

    nelem = mesh.nelem
    m = 1
    C_stat = calc_static_constants(material[m])
    delta = I(3) 
    nnel = (mesh.eltype+1)^2

    max_GL = maximum(mesh.ID)
    solver_var.H = zeros(max_GL,max_GL)
    solver_var.G = zeros(max_GL,max_GL)

    csis_cont = range(-1,1,length = 2)
    omegas = calc_omegas(solver_var.omega)

    csis = calc_csis_grid(solver_var.csi)
    # Nc = calc_N_matrix(csis_cont,csis)
    # dNcdcsi = calc_dNdcsi_matrix(csis_cont,csis)
    # dNcdeta = calc_dNdeta_matrix(csis_cont,csis)


     # Cálculos para elementos não singulares
     gp, dists = calc_points_weights()
 
    Nc = []
    dNcdcsi = []
    dNcdeta = []
 
    for i in eachindex(gp)
        N = calc_N_matrix(csis_cont,gp[i].csi)
        Ncsi = calc_dNdcsi_matrix(csis_cont,gp[i].csi)
        Neta = calc_dNdeta_matrix(csis_cont,gp[i].csi)
        push!(Nc,N)
        push!(dNcdcsi,Ncsi)
        push!(dNcdeta,Neta)
    end
    
    csis_telles, omegas_telles =  gausslegendre(12)

    # Field loop:
    # Threads.@threads for fe in 1:nelem
    for fe in 1:nelem
    
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]

        normal = []
        J = []
        gauss_points = []
        for i in eachindex(gp)
            normal2, J2 = calc_n_J_matrix(dNcdcsi[i], dNcdeta[i], field_points)
            gauss_points2 = Nc[i]*field_points
            push!(normal,normal2)
            push!(J,J2)
            push!(gauss_points,gauss_points2)
        end

        # source loop:
        for se in 1:nelem
            for n = 1:nnel
                sn = mesh.IEN[n,se]
                source_node = mesh.nodes[sn,2:end]

                if fe != se

                    r, c, d = calc_dist(source_node, field_points, dists,csis_cont)

                    if r > 0
                        # HELEM, GELEM = integration_const_raw(source_node, field_points, normal[1,:], solver_var.csi, solver_var.omega, delta, C_stat)
                        HELEM, GELEM = integrate_const_static(source_node,gauss_points[r],normal[r],J[r], gp[r].omega, delta, C_stat)
                    else
                        # near integration
                        gp_local_near, weights_near = pontos_pesos_local_telles(csis_telles, omegas_telles, c[1:2],d)
                        N_near = calc_N_matrix(csis_cont,gp_local_near)
                        dNc_near = calc_dNdcsi_matrix(csis_cont,gp_local_near)
                        dNe_near = calc_dNdeta_matrix(csis_cont,gp_local_near)
                        gp_near = N_near*field_points
                        normal_near,J_near = calc_n_J_matrix(dNc_near, dNe_near, field_points)
                        HELEM, GELEM = integrate_const_static(source_node,gp_near,normal_near,J_near,weights_near, delta, C_stat)
                        # @infiltrate
                    end
                else
                    GELEM = zeros(3,3)
                    HELEM = zeros(3,3)
                    points2 = zeros(4,3)
                    points2[3,:] = source_node
                    points2[4,:] = source_node

                    for k in 1:4
                        if k ==1
                            points2[1,:] = field_points[1,:]
                            points2[2,:] = field_points[2,:]
                        else 
                            if k == 2
                                points2[1,:] = field_points[2,:]
                                points2[2,:] = field_points[3,:]
                            else
                                if k ==3
                                    points2[1,:] = field_points[3,:]
                                    points2[2,:] = field_points[4,:]
                                else
                                    if k ==4
                                        points2[1,:] = field_points[4,:]
                                        points2[2,:] = field_points[1,:]
                                    end
                                end
                            end
                        end
                    

                        normal2, J2= calc_n_J_matrix(dNcdcsi[end], dNcdeta[end], points2)
                        gauss_points2 = Nc[end]*points2

                        _, GELEM2 = integrate_const_static(source_node, gauss_points2, normal2, J2, gp[end].omega, delta, C_stat)
                        # _, GELEM2 = integration_const_raw(source_node, points2, normal2[1,:], solver_var.csi, solver_var.omega, delta, C_stat)
                        # @infiltrate
                        GELEM = GELEM + GELEM2
                    end
                end

                solver_var.H[mesh.ID[:,sn], mesh.LM[:,fe]] = HELEM
                solver_var.G[mesh.ID[:,sn], mesh.LM[:,fe]] = GELEM
            end

        end
    
    end


    
    for n in 1:mesh.nnodes
        solver_var.H[mesh.ID[:,n],mesh.ID[:,n]] .= 0
    end

    integrate_rigid_body!(solver_var.H,mesh)
    
    return solver_var  

end

# Working:
function calc_GH_static!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type)
    
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
    N_sing, dNc_sing, dNe_sing, Nd_sing, Jb_sing, omega_sing = calc_N_sing(mesh, solver_var, csis_cont, csis_descont, rules.gp[end].omega)

    #matrices initialization
    max_GL = maximum(mesh.ID)
    solver_var.H = zeros(max_GL,max_GL)
    solver_var.G = zeros(max_GL,max_GL)

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
                solver_var.H[mesh.ID[:,sn], mesh.LM[:,fe]] = HELEM
                solver_var.G[mesh.ID[:,sn], mesh.LM[:,fe]] = GELEM
            end
            
        end
        next!(p)
    end

    # RIGID BODY MOTION STRATEGY
    for n in 1:mesh.nnodes
        solver_var.H[mesh.ID[:,n],mesh.ID[:,n]] .= 0
    end
    integrate_rigid_body!(solver_var.H,mesh)

    return solver_var

end

function calc_GH_static_non_const_azim!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type)

    nelem = mesh.nelem
    m = 1
    C_stat = calc_static_constants(material[m])
    delta = I(3) 
    nnel = (mesh.eltype+1)^2

    max_GL = maximum(mesh.ID)
    solver_var.H = zeros(max_GL,max_GL)
    solver_var.G = zeros(max_GL,max_GL)

    csis_cont = range(-1,1,length = mesh.eltype+1)
    csis_descont = range(-1+mesh.offset,1 - mesh.offset,length=mesh.eltype+1)
    omegas = calc_omegas(solver_var.omega)

    # Cálculos para elementos não singulares
    gp, dists = calc_points_weights()

    Nc_nonsing = []
    dNcdcsi_nonsing = []
    dNcdeta_nonsing = []
    Nd_nonsing = []

    for i in eachindex(gp)
        N = calc_N_matrix(csis_cont,gp[i].csi)
        Nd = calc_N_matrix(csis_descont,gp[i].csi)
        Ncsi = calc_dNdcsi_matrix(csis_cont,gp[i].csi)
        Neta = calc_dNdeta_matrix(csis_cont,gp[i].csi)
        push!(Nc_nonsing,N)
        push!(dNcdcsi_nonsing,Ncsi)
        push!(dNcdeta_nonsing,Neta)
        push!(Nd_nonsing,Nd)
    end


    # Cálculos para elementos singulares
    csis = calc_csis_grid(solver_var.csi)
    csis_cont_lin = range(-1,1,length = 2)
    Nc_lin = calc_N_matrix(csis_cont_lin,csis)
    dNcdcsi_lin = calc_dNdcsi_matrix(csis_cont_lin,csis)
    dNcdeta_lin = calc_dNdeta_matrix(csis_cont_lin,csis)

    if mesh.eltype == 2
        nperm = 3
    elseif mesh.eltype < 2
        nperm = 1
    end

    csi_sing, Jb_sing = csis_sing(mesh.offset, Nc_lin, dNcdcsi_lin,dNcdeta_lin,mesh.eltype)
    npg_sing = size(csi_sing,1)
    Nc_sing = zeros(npg_sing,nnel,nperm)
    dNcdcsi_sing = zeros(npg_sing,nnel,nperm)
    dNcdeta_sing = zeros(npg_sing,nnel,nperm)
    Nd_sing = zeros(npg_sing,nnel,nperm)

    for i in 1:nperm
        Nc_sing[:,:,i] = calc_N_matrix(csis_cont,csi_sing[:,:,i])
        dNcdcsi_sing[:,:,i] = calc_dNdcsi_matrix(csis_cont,csi_sing[:,:,i])
        dNcdeta_sing[:,:,i] = calc_dNdeta_matrix(csis_cont,csi_sing[:,:,i])
        Nd_sing[:,:,i] = calc_N_matrix(csis_descont,csi_sing[:,:,i])
    end
    omega_sing = repeat(omegas,4)

    csis_telles, omegas_telles =  gausslegendre(12)


    # @infiltrate
    # Field loop:
    Threads.@threads for fe in 1:nelem
    # for fe in 1:nelem
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]
        normal = []
        J = []
        gauss_points = []
        for i in eachindex(gp)
            normal2, J2 = calc_n_J_matrix(dNcdcsi_nonsing[i], dNcdeta_nonsing[i], field_points)
            gauss_points2 = Nc_nonsing[i]*field_points
            push!(normal,normal2)
            push!(J,J2)
            push!(gauss_points,gauss_points2)
        end

        # source loop:
        for se in 1:nelem
            for n = 1:nnel
                sn = mesh.IEN[n,se]
                source_node = mesh.nodes[sn,2:end]

                if fe != se

                    r, c, d = calc_dist(source_node, field_points, dists,csis_cont)
                    if r > 0
                        # normal integration
                        HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points[r],Nd_nonsing[r],normal[r],J[r], gp[r].omega, delta, C_stat)
                    else
                        # near integration
                        gp_local_near, weights_near = pontos_pesos_local_telles(csis_telles, omegas_telles, c[1:2],d)
                        # gp_local_near, weights_near = pontos_pesos_local_subelem(c[1], c[2], Nc_lin, dNcdcsi_lin, dNcdeta_lin, omegas)
                        N_near = calc_N_matrix(csis_cont,gp_local_near)
                        dNc_near = calc_dNdcsi_matrix(csis_cont,gp_local_near)
                        dNe_near = calc_dNdeta_matrix(csis_cont,gp_local_near)
                        Nd_near = calc_N_matrix(csis_descont,gp_local_near)
                        gp_near = N_near*field_points
                        normal_near,J_near = calc_n_J_matrix(dNc_near, dNe_near, field_points)
                        HELEM, GELEM = integrate_nonsing_static(source_node,gp_near,Nd_near,normal_near,J_near, weights_near, delta, C_stat)
                    end
                else

                    idx_ff, np = calc_idx_permutation(nnel,n)

                    normal_sing, J_sing = calc_n_J_matrix(dNcdcsi_sing[:,idx_ff,np], dNcdeta_sing[:,idx_ff,np], field_points)

                    # @infiltrate
                    J_sing = J_sing.*Jb_sing[:,np]
                    pesos = zeros(length(omega_sing),nnel)

                    for kn in 1:nnel
                        pesos[:,kn] = J_sing.*omega_sing.*Nd_sing[:,idx_ff[kn],np]
                    end
                    # @infiltrate
                    gauss_points_sing = Nc_sing[:,idx_ff,np]*field_points
                    # @infiltrate
                    # GELEM1 = integrate_nonsing_static2(source_node,gauss_points[end],Nd_nonsing[end],normal[end],J[end], gp[end].omega, delta, C_stat,n)

                    HELEM, GELEM = integrate_sing_static(source_node, gauss_points_sing, normal_sing, delta, C_stat, pesos, n)
                    # GELEM = GELEM + GELEM1
                    # @infiltrate

                    # HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points,Nd,normal,J, omegas, delta, C_stat)
                end

                solver_var.H[mesh.ID[:,sn], mesh.LM[:,fe]] = HELEM
                solver_var.G[mesh.ID[:,sn], mesh.LM[:,fe]] = GELEM
            end

        end
    
    end

    for n in 1:mesh.nnodes
        solver_var.H[mesh.ID[:,n],mesh.ID[:,n]] .= 0
    end

    integrate_rigid_body!(solver_var.H,mesh)

    for n in 1:mesh.nnodes
        solver_var.H[mesh.ID[:,n],mesh.ID[:,n]] += delta./2
    end
    
    return solver_var

end

function calc_GH_static_const_azim!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type)

    nelem = mesh.nelem
    m = 1
    C_stat = calc_static_constants(material[m])
    delta = I(3) 
    nnel = (mesh.eltype+1)^2

    max_GL = maximum(mesh.ID)
    solver_var.H = zeros(max_GL,max_GL)
    solver_var.G = zeros(max_GL,max_GL)

    csis_cont = range(-1,1,length = 2)
    omegas = calc_omegas(solver_var.omega)

    csis = calc_csis_grid(solver_var.csi)
    # Nc = calc_N_matrix(csis_cont,csis)
    # dNcdcsi = calc_dNdcsi_matrix(csis_cont,csis)
    # dNcdeta = calc_dNdeta_matrix(csis_cont,csis)


     # Cálculos para elementos não singulares
     gp, dists = calc_points_weights()
 
    Nc = []
    dNcdcsi = []
    dNcdeta = []
 
    for i in eachindex(gp)
        N = calc_N_matrix(csis_cont,gp[i].csi)
        Ncsi = calc_dNdcsi_matrix(csis_cont,gp[i].csi)
        Neta = calc_dNdeta_matrix(csis_cont,gp[i].csi)
        push!(Nc,N)
        push!(dNcdcsi,Ncsi)
        push!(dNcdeta,Neta)
    end
    
    csis_telles, omegas_telles =  gausslegendre(12)

    # Field loop:
    # Threads.@threads for fe in 1:nelem
    for fe in 1:nelem
    
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]

        normal = []
        J = []
        gauss_points = []
        for i in eachindex(gp)
            normal2, J2 = calc_n_J_matrix(dNcdcsi[i], dNcdeta[i], field_points)
            gauss_points2 = Nc[i]*field_points
            push!(normal,normal2)
            push!(J,J2)
            push!(gauss_points,gauss_points2)
        end

        # source loop:
        for se in 1:nelem
            for n = 1:nnel
                sn = mesh.IEN[n,se]
                source_node = mesh.nodes[sn,2:end]

                if fe != se

                    r, c, d = calc_dist(source_node, field_points, dists,csis_cont)

                    if r > 0
                        # HELEM, GELEM = integration_const_raw(source_node, field_points, normal[1,:], solver_var.csi, solver_var.omega, delta, C_stat)
                        HELEM, GELEM = integrate_const_static(source_node,gauss_points[r],normal[r],J[r], gp[r].omega, delta, C_stat)
                    else
                        # near integration
                        gp_local_near, weights_near = pontos_pesos_local_telles(csis_telles, omegas_telles, c[1:2],d)
                        N_near = calc_N_matrix(csis_cont,gp_local_near)
                        dNc_near = calc_dNdcsi_matrix(csis_cont,gp_local_near)
                        dNe_near = calc_dNdeta_matrix(csis_cont,gp_local_near)
                        gp_near = N_near*field_points
                        normal_near,J_near = calc_n_J_matrix(dNc_near, dNe_near, field_points)
                        HELEM, GELEM = integrate_const_static(source_node,gp_near,normal_near,J_near,weights_near, delta, C_stat)
                        # @infiltrate
                    end
                else
                    GELEM = zeros(3,3)
                    HELEM = zeros(3,3)
                    points2 = zeros(4,3)
                    points2[3,:] = source_node
                    points2[4,:] = source_node

                    for k in 1:4
                        if k ==1
                            points2[1,:] = field_points[1,:]
                            points2[2,:] = field_points[2,:]
                        else 
                            if k == 2
                                points2[1,:] = field_points[2,:]
                                points2[2,:] = field_points[3,:]
                            else
                                if k ==3
                                    points2[1,:] = field_points[3,:]
                                    points2[2,:] = field_points[4,:]
                                else
                                    if k ==4
                                        points2[1,:] = field_points[4,:]
                                        points2[2,:] = field_points[1,:]
                                    end
                                end
                            end
                        end
                    

                        normal2, J2= calc_n_J_matrix(dNcdcsi[end], dNcdeta[end], points2)
                        gauss_points2 = Nc[end]*points2

                        _, GELEM2 = integrate_const_static(source_node, gauss_points2, normal2, J2, gp[end].omega, delta, C_stat)
                        # _, GELEM2 = integration_const_raw(source_node, points2, normal2[1,:], solver_var.csi, solver_var.omega, delta, C_stat)
                        # @infiltrate
                        GELEM = GELEM + GELEM2
                    end
                end

                solver_var.H[mesh.ID[:,sn], mesh.LM[:,fe]] = HELEM
                solver_var.G[mesh.ID[:,sn], mesh.LM[:,fe]] = GELEM
            end

        end
    
    end


    
    for n in 1:mesh.nnodes
        solver_var.H[mesh.ID[:,n],mesh.ID[:,n]] .= 0
    end

    integrate_rigid_body!(solver_var.H,mesh)

    for n in 1:mesh.nnodes
        solver_var.H[mesh.ID[:,n],mesh.ID[:,n]] += delta./2
    end
    
    return solver_var  

end
