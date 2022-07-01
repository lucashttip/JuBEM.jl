function calc_GH_static_non_const2!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type)

    nelem = mesh.nelem
    m = 1
    C_stat = calc_static_constants(material[m])
    delta = I(3)
    delta2 = delta./2
    z = zeros(3,3)
    zHselem = complex(delta)./2.0
    

    max_GL = maximum(mesh.ID)
    solver_var.H = zeros(max_GL,max_GL)
    solver_var.G = zeros(max_GL,max_GL)

    nnel = (mesh.eltype+1)^2

    csis_cont = range(-1,1,length = mesh.eltype+1)
    csis_descont = range(-1+mesh.offset,1 - mesh.offset,length=mesh.eltype+1)
    omegas = calc_omegas(solver_var.omega)

    Nc, Nd, dNdcsi, dNdeta, dNcdcsi, dNcdeta = calc_Ns(csis_cont, csis_descont, solver_var.csi,solver_var.csi)

    csi_sing, omega_sing = csis_sing_3(mesh.offset, Nc, omegas)
    Nc_sing, Nd_sing, dNdcsi_sing, dNdeta_sing = calc_Ns_sing(csis_cont, csis_descont, csi_sing)
    
    # Field loop:
    for fe in 1:nelem
    
        field_nodes = mesh.nodes[mesh.IEN[:,fe],2:end]
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]

        normal, J = calc_n_J_matrix(dNdcsi, dNdeta, field_nodes)
        gauss_points = generate_points_in_elem(Nd,field_nodes)
        # source loop:
        for se in 1:nelem
            for n = 1:nnel
                sn = mesh.IEN[n,se]
                source_node = mesh.nodes[sn,2:end]

                if fe != se
                    HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points,Nd,normal,J, omegas, delta, C_stat)
                else

                    idx = collect((1:nnel) .+ (n-1))
                    idx[idx.>nnel] = idx[idx.>nnel] .-nnel
                    idx = [1,2,3,4]
                    if n == 2
                        idx = [4,1,2,3]
                    elseif n == 3
                        idx = [3,4,1,2]
                    elseif n==4
                        idx = [2,3,4,1]
                    end
                    
                    Nd_sing2 = Nd_sing[:,idx]
                    dNdcsi_sing2 = dNdcsi_sing[:,idx]
                    dNdeta_sing2 = dNdeta_sing[:,idx]
                    dNcdcsi2 = dNcdcsi[:,idx]
                    dNcdeta2 = dNcdeta[:,idx]

                    # normal_sing, J_sing = calc_n_J_matrix(dNdcsi_sing2, dNdeta_sing2, field_nodes)
                    normal_sing, J_sing = calc_n_J_matrix_sing(dNcdcsi2, dNcdeta2, field_points,field_nodes[n,:])
                    gauss_points_sing = generate_points_in_elem(Nd_sing2,field_nodes)
                    # Plots.scatter(gauss_points_sing[:,2],gauss_points_sing[:,3])
                    # @infiltrate
                    HELEM, GELEM = integrate_sing_static_1(source_node, gauss_points_sing, Nd_sing2, normal_sing, J_sing, omega_sing, delta, C_stat, n)
                    # HELEM = zeros(3,3*nnel)
                    # GELEM = zeros(3,3*nnel)

                    # HELEM = HELEM + [repeat(z,n-1);delta2;repeat(z,nnel-n)]'
                end

                solver_var.H[mesh.ID[:,sn], mesh.LM[:,fe]] = HELEM
                solver_var.G[mesh.ID[:,sn], mesh.LM[:,fe]] = GELEM
            end

        end
    
        integrate_rigid_body2!(solver_var.H,nnel)

    end
    return solver_var

end

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
    # dNddcsi = dNcdcsi*G
    # dNddeta = dNcdeta*G
    
    csi_sing, omega_sing = csis_sing_3(mesh.offset, Nc, omegas)
    Nc_sing = calc_N_matrix(csis_cont,csi_sing)
    Nd_sing = Nc_sing*G
    
    # Field loop:
    for fe in 1:nelem
    
        field_nodes = mesh.nodes[mesh.IEN[:,fe],2:end]
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]

        normal, J = calc_n_J_matrix(dNcdcsi, dNcdeta, field_points)
        gauss_points = Nc*field_points
        
        # source loop:
        for se in 1:nelem
            for n = 1:nnel
                sn = mesh.IEN[n,se]
                source_node = mesh.nodes[sn,2:end]

                if fe != se
                    HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points,Nd,normal,J, omegas, delta, C_stat)
                else

                    idx = collect((1:nnel) .- (n-1))
                    idx[idx.<1] = idx[idx.<1] .+nnel

                    gauss_points_sing = Nc_sing[:,idx]*field_points
                    normal_sing,J_sing = divide_elem(source_node,field_points,dNcdcsi,dNcdeta)

                    HELEM, GELEM = integrate_sing_static_1(source_node, gauss_points_sing, Nd_sing[:,idx], normal_sing, J_sing, omega_sing, delta, C_stat, n)
                    # HELEM = zeros(3,3*nnel)
                    # GELEM = zeros(3,3*nnel)
                end

                solver_var.H[mesh.ID[:,sn], mesh.LM[:,fe]] = HELEM
                solver_var.G[mesh.ID[:,sn], mesh.LM[:,fe]] = GELEM
            end

        end
    
        # integrate_rigid_body2!(solver_var.H,nnel)
        
        # for i in 1:size(solver_var.H)

    end
    integrate_rigid_body2!(solver_var.H,nnel)
    return solver_var

end