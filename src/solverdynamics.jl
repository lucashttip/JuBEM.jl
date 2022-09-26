function calc_GH_dynamic_non_const!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type, frequency::Float64)

    nelem = mesh.nelem
    m = 1
    C_stat = calc_static_constants(material[m])
    zconsts = cmplx_consts(material[m],frequency)
    delta = I(3) 
    nnel = (mesh.eltype+1)^2

    max_GL = maximum(mesh.ID)
    solver_var.zH = zeros(ComplexF64,max_GL,max_GL)
    solver_var.zG = zeros(ComplexF64,max_GL,max_GL)

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
                        zHELEM, zGELEM = integrate_nonsing_dynamic(source_node,gauss_points[r],Nd_nonsing[r],normal[r],J[r], gp[r].omega, delta, zconsts)
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
                        zHELEM, zGELEM = integrate_nonsing_dynamic(source_node,gp_near,Nd_near,normal_near,J_near, weights_near, delta, zconsts)
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
                    nk = 5
                    zGELEM1 = integrate_nonsing_dynamic2(source_node,gauss_points[nk],Nd_nonsing[nk],normal[nk],J[nk], gp[nk].omega, delta, zconsts,n)
                    zHELEM1, zGELEM2 = integrate_sing_dynamic(source_node, gauss_points_sing, normal_sing, delta, zconsts, pesos, n)
                    zHELEM2 = integrate_sing_dynamic2(source_node,gauss_points_sing,normal_sing, delta, zconsts,C_stat,pesos,n)

                    zGELEM = zGELEM1 + zGELEM2
                    zHELEM = zHELEM1 + zHELEM2
                    # @infiltrate

                    # HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points,Nd,normal,J, omegas, delta, C_stat)
                end

                solver_var.zH[mesh.ID[:,sn], mesh.LM[:,fe]] = zHELEM
                solver_var.zG[mesh.ID[:,sn], mesh.LM[:,fe]] = zGELEM
            end

        end
    
    end

    for n in 1:mesh.nnodes
        solver_var.zH[mesh.ID[:,n],mesh.ID[:,n]] += solver_var.H[mesh.ID[:,n],mesh.ID[:,n]]
    end
    
    return solver_var

end


function calc_GH_dynamic_const!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type, frequency::Float64)

    nelem = mesh.nelem
    m = 1
    C_stat = calc_static_constants(material[m])
    zconsts = cmplx_consts(material[m],frequency)
    delta = I(3) 
    nnel = (mesh.eltype+1)^2

    max_GL = maximum(mesh.ID)
    solver_var.zH = zeros(ComplexF64,max_GL,max_GL)
    solver_var.zG = zeros(ComplexF64,max_GL,max_GL)

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
    Threads.@threads for fe in 1:nelem
    # for fe in 1:nelem
    
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
                    # @infiltrate
                    r, c, d = calc_dist(source_node, field_points, dists,csis_cont)

                    if r > 0
                        zHELEM, zGELEM = integrate_const_dynamic(source_node,gauss_points[r],normal[r],J[r], gp[r].omega, delta, zconsts)
                    else
                        # near integration
                        gp_local_near, weights_near = pontos_pesos_local_telles(csis_telles, omegas_telles, c[1:2],d)
                        N_near = calc_N_matrix(csis_cont,gp_local_near)
                        dNc_near = calc_dNdcsi_matrix(csis_cont,gp_local_near)
                        dNe_near = calc_dNdeta_matrix(csis_cont,gp_local_near)
                        gp_near = N_near*field_points
                        normal_near,J_near = calc_n_J_matrix(dNc_near, dNe_near, field_points)
                        zHELEM, zGELEM = integrate_const_dynamic(source_node,gp_near,normal_near,J_near,weights_near, delta, zconsts)
                    end

                    
                else
                    zGELEM = zeros(ComplexF64,3,3)
                    zHELEM = zeros(ComplexF64,3,3)
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


                        _, zGELEM2 = integrate_const_dynamic(source_node, gauss_points2, normal2, J2, gp[end].omega, delta, zconsts)
                        # _, GELEM2 = integration_const_raw(source_node, points2, normal2[1,:], solver_var.csi, solver_var.omega, delta, C_stat)
                        # @infiltrate
                        zGELEM = zGELEM + zGELEM2
                        zHELEM = zHELEM + integrate_const_sing_dynamic(source_node,gauss_points2,normal2,J2,gp[end].omega, delta, zconsts,C_stat)
                    end
                    # zHELEM = integrate_const_sing_dynamic(source_node,gauss_points,normal,J,omegas, delta, zconsts,C_stat)
                end

                solver_var.zH[mesh.ID[:,sn], mesh.LM[:,fe]] = zHELEM
                solver_var.zG[mesh.ID[:,sn], mesh.LM[:,fe]] = zGELEM
            end

        end
    
    end


    
    for n in 1:mesh.nnodes
        solver_var.zH[mesh.ID[:,n],mesh.ID[:,n]] += solver_var.H[mesh.ID[:,n],mesh.ID[:,n]]
    end
    
    return solver_var


end