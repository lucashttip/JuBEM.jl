function calc_GH_static_non_const!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type)

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

    k = calc_k(nnel)
    csis = calc_csis_grid(solver_var.csi)
    Nc = calc_N_matrix(csis_cont,csis)
    dNcdcsi = calc_dNdcsi_matrix(csis_cont,csis)
    dNcdeta = calc_dNdeta_matrix(csis_cont,csis)

    G = calc_G(csis_cont,csis_descont,k)

    Nd = Nc*G
    # @infiltrate
    
    csi_sing, Jb_sing = csis_sing2(mesh.offset, Nc, dNcdcsi,dNcdeta)
    Nc_sing2 = calc_N_matrix(csis_cont,csi_sing)
    dNcdcsi_sing2 = calc_dNdcsi_matrix(csis_cont,csi_sing)
    dNcdeta_sing2 = calc_dNdeta_matrix(csis_cont,csi_sing)
    Nd_sing2 = calc_N_matrix(csis_descont,csi_sing)
    omega_sing = repeat(omegas,4)

    # c_sing, IEN_sing = divide_elem(mesh.offset)
    # npg2 = size(Nc,1)
    # nsubelem = size(IEN_sing,2)
    # Nc_sing = zeros(nsubelem*npg2,nnel)
    # dNcdcsi_sing = zeros(nsubelem*npg2,nnel)
    # dNcdeta_sing = zeros(nsubelem*npg2,nnel)
    # omega_sing = repeat(omegas,nsubelem)

    # for e in 1:nsubelem
    #     tmp_N = calc_N_matrix(csis_cont,c_sing[IEN_sing[:,e],:])
    #     Nc_sing[(e-1)*npg2+1:e*npg2,:] = Nc*tmp_N
    #     dNcdcsi_sing[(e-1)*npg2+1:e*npg2,:] = dNcdcsi*tmp_N
    #     dNcdeta_sing[(e-1)*npg2+1:e*npg2,:] = dNcdeta*tmp_N
    # end

    # Field loop:
    # Threads.@threads for fe in 1:nelem
    for fe in 1:nelem
    
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]

        normal, J = calc_n_J_matrix(dNcdcsi, dNcdeta, field_points)
        gauss_points = Nc*field_points
        # @infiltrate
        # source loop:
        for se in 1:nelem
            for n = 1:nnel
                sn = mesh.IEN[n,se]
                source_node = mesh.nodes[sn,2:end]

                if fe != se
                    HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points,Nd,normal,J, omegas, delta, C_stat)
                else

                    idx_ff, idx_points = calc_idx_permutation(nnel,n)

                    # normal_sing, J_sing, omega_sing, gauss_points_sing = divide_elem(source_node,field_points[idx_points,:],Nc,dNcdcsi,dNcdeta,omegas)

                    normal_sing, J_sing = calc_n_J_matrix(dNcdcsi_sing2, dNcdeta_sing2, field_points[idx_points,:])
                    # @infiltrate
                    J_sing = J_sing.*Jb_sing
                    # @infiltrate
                    gauss_points_sing = Nc_sing2[:,idx_ff]*field_points
                    # @infiltrate
                    # GELEM1 = integrate_nonsing_static2(source_node,gauss_points,Nd,normal,J, omegas, delta, C_stat,n)
                    HELEM, GELEM = integrate_sing_static(source_node, gauss_points_sing, Nd_sing2[:,idx_ff], normal_sing, J_sing, omega_sing, delta, C_stat, n)
                    # GELEM = GELEM + GELEM1

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
    Nc = calc_N_matrix(csis_cont,csis)
    dNcdcsi = calc_dNdcsi_matrix(csis_cont,csis)
    dNcdeta = calc_dNdeta_matrix(csis_cont,csis)

    points2 = zeros(4,3)
    
    # Field loop:
    Threads.@threads for fe in 1:nelem
    # for fe in 1:nelem
    
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]

        normal, J = calc_n_J_matrix(dNcdcsi, dNcdeta, field_points)
        gauss_points = Nc*field_points

        # source loop:
        for se in 1:nelem
            for n = 1:nnel
                sn = mesh.IEN[n,se]
                source_node = mesh.nodes[sn,2:end]

                if fe != se
                    # HELEM, GELEM = integration_const_raw(source_node, field_points, normal[1,:], solver_var.csi, solver_var.omega, delta, C_stat)
                    HELEM, GELEM = integrate_const_static(source_node,gauss_points,normal,J, omegas, delta, C_stat)
                else
                    GELEM = zeros(3,3)
                    HELEM = zeros(3,3)
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
                    end

                    normal2, J2= calc_n_J_matrix(dNcdcsi, dNcdeta, field_points)
                    gauss_points2 = Nc*field_points

                    _, GELEM2 = integrate_const_static(source_node, gauss_points2, normal2, J2, omegas, delta, C_stat)
                    # _, GELEM2 = integration_const_raw(source_node, points2, normal[1,:], solver_var.csi, solver_var.omega, delta, C_stat)

                    GELEM = GELEM + GELEM2
                end

                solver_var.H[mesh.ID[:,sn], mesh.LM[:,fe]] = HELEM
                solver_var.G[mesh.ID[:,sn], mesh.LM[:,fe]] = GELEM
            end

        end
    
    end

    integrate_rigid_body!(solver_var.H,mesh)
    
    return solver_var

end