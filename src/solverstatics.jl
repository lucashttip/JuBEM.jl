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
    
    csi_sing, _ = csis_sing(mesh.offset, Nc, omegas)
    Nc_sing = calc_N_matrix(csis_cont,csi_sing)
    Nd_sing = Nc_sing*G
    

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
                    HELEM, GELEM = integrate_nonsing_static(source_node,gauss_points,Nd,normal,J, omegas, delta, C_stat)
                else

                    idx_ff = collect((1:nnel) .- (n-1))
                    idx_ff[idx_ff.<1] = idx_ff[idx_ff.<1] .+nnel

                    idx_points = collect((1:nnel) .+ (n-1))
                    idx_points[idx_points.>nnel] = idx_points[idx_points.>nnel] .-nnel

                    normal_sing, J_sing, omega_sing, gauss_points_sing = divide_elem(source_node,field_points[idx_points,:],Nc,dNcdcsi,dNcdeta,omegas)
                    GELEM1 = integrate_nonsing_static2(source_node,gauss_points,Nd,normal,J, omegas, delta, C_stat,n)
                    HELEM, GELEM2 = integrate_sing_static(source_node, gauss_points_sing, Nd_sing[:,idx_ff], normal_sing, J_sing, omega_sing, delta, C_stat, n)

                    GELEM = GELEM1 + GELEM2
                end

                solver_var.H[mesh.ID[:,sn], mesh.LM[:,fe]] = HELEM
                solver_var.G[mesh.ID[:,sn], mesh.LM[:,fe]] = GELEM
            end

        end
    
    end

    integrate_rigid_body!(solver_var.H,mesh)
    
    return solver_var

end