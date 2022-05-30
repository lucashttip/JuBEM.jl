function calc_GH_dynamic_non_const!(mesh::mesh_type, material::Vector{material_table_type}, solver_var::solver_var_type, frequency::Float64)

    nelem = mesh.nelem
    m = 1
    C_stat = calc_static_constants(material[m])
    zconsts = cmplx_consts(material[m],frequency)
    delta = I(3)
    zHselem = complex(delta)./2.0

    max_GL = maximum(mesh.ID)
    solver_var.H = complex(zeros(max_GL,max_GL))
    solver_var.G = complex(zeros(max_GL,max_GL))

    nnel = (mesh.eltype+1)^2

    csis_cont = range(-1,1,length = mesh.eltype+1)
    csis_descont = range(-1+mesh.offset,1 - mesh.offset,length=mesh.eltype+1)

    Nc = calc_N_matrix(csis_cont, solver_var.csi)
    dNcdcsi = calc_dNdcsi_matrix(csis_cont, solver_var.csi)
    dNcdeta = calc_dNdeta_matrix(csis_cont,solver_var.csi)
    G = calc_G(csis_cont, csis_descont)
    Nd = calc_N_matrix_descont(Nc,G)
    dNdcsi = calc_N_matrix_descont(dNcdcsi,G)
    dNdeta = calc_N_matrix_descont(dNcdeta,G)

    # Field loop:
    for fe in 1:nelem
    
        field_nodes = mesh.nodes[mesh.IEN[:,fe],2:end]

        normal, J = calc_n_J_matrix(dNdcsi, dNdeta, field_nodes)
        gauss_points = generate_points_in_elem(Nd,field_nodes)
        # source loop:
        for se in 1:nelem
            for n = 1:nnel
                sn = mesh.IEN[n,se]
                source_node = mesh.nodes[sn,2:end]

                if fe != se
                    HELEM, GELEM = calc_nonsing(source_node,gauss_points,Nd,normal,J, solver_var.omega, delta, zconsts)
                else
                    # HELEM, GELEM = calc_sing()
                    HELEM = zeros(ComplexF64,3,3*nnel)
                    GELEM = zeros(ComplexF64,3,3*nnel)
                end

                solver_var.H[mesh.ID[:,sn], mesh.LM[:,fe]] = HELEM
                solver_var.G[mesh.ID[:,sn], mesh.LM[:,fe]] = GELEM
            end

        end
    
    end
    return solver_var
end