function calc_N_nonsing(csis_cont, csis_descont, rules)

    N = []
    dNc = []
    dNe = []
    Nd = []

    for i in eachindex(rules.gp)
        N_t = calc_N_matrix(csis_cont,rules.gp[i].csi)
        dNc_t = calc_dNdcsi_matrix(csis_cont,rules.gp[i].csi)
        dNe_t = calc_dNdeta_matrix(csis_cont,rules.gp[i].csi)
        Nd_t = calc_N_matrix(csis_descont,rules.gp[i].csi)
        push!(N,N_t)
        push!(dNc,dNc_t)
        push!(dNe,dNe_t)
        push!(Nd,Nd_t)
    end

    return N, dNc, dNe, Nd
end

function calc_N_sing(mesh, solver_var, csis_cont, csis_descont,npg_sing)

    # Cálculos para elementos singulares
    csi,omega = gausslegendre(npg_sing)
    csis = calc_csis_grid(csi)
    omegas = calc_omegas(omega)
    csis_cont_lin = range(-1,1,length = 2)
    N_lin = calc_N_matrix(csis_cont_lin,csis)
    dNc_lin = calc_dNdcsi_matrix(csis_cont_lin,csis)
    dNe_lin = calc_dNdeta_matrix(csis_cont_lin,csis)

    eltype_cont = mesh.eltype
    if mesh.eltype == 0
        eltype_cont = 1
    end

    if mesh.eltype == 2
        nperm = 3
    elseif mesh.eltype < 2
        nperm = 1
    end

    npel = (eltype_cont+1)^2
    nnel = (mesh.eltype+1)^2

    csi_sing, Jb_sing = csis_sing(mesh.offset, N_lin, dNc_lin,dNe_lin,mesh.eltype)
    npg_sing = size(csi_sing,1)
    
    N_sing = zeros(npg_sing,npel,nperm)
    dNc_sing = zeros(npg_sing,npel,nperm)
    dNe_sing = zeros(npg_sing,npel,nperm)
    Nd_sing = zeros(npg_sing,nnel,nperm)

    for i in 1:nperm
        N_sing[:,:,i] = calc_N_matrix(csis_cont,csi_sing[:,:,i])
        dNc_sing[:,:,i] = calc_dNdcsi_matrix(csis_cont,csi_sing[:,:,i])
        dNe_sing[:,:,i] = calc_dNdeta_matrix(csis_cont,csi_sing[:,:,i])
        Nd_sing[:,:,i] = calc_N_matrix(csis_descont,csi_sing[:,:,i])
    end
    nregs = Int(length(Jb_sing)/length(omegas))
    omega_sing = repeat(omegas,nregs)

    return N_sing, dNc_sing, dNe_sing, Nd_sing, Jb_sing, omega_sing

end

function calc_integration_rules!(rules)
    
    pw = gausslegendre.(rules.npgs)

    for i in eachindex(rules.npgs)
        csis = calc_csis_grid(pw[i][1])
        omegas = calc_omegas(pw[i][2])
        push!(rules.gp,gauss_points_type(csis,omegas))
    end

    pnear = gausslegendre(rules.npg_near)
    csis = calc_csis_grid(pnear[1])
    omegas = calc_omegas(pnear[2])
    rules.gp_near = gauss_points_type(csis,omegas)

    return rules
end

function calc_nonsing_consts(mesh)

    npgs = Int[4,5,6,8]
    dists = [4,2,0.5,0.3]
    npg_near = 12
    npg_sing = 5
    # npgs = Int[4,4,4,8]
    # dists = [4,2,0.5,0.2]
    # npg_near = 12

    rules = integration_rules_type(npgs,dists,npg_near,npg_sing)
    calc_integration_rules!(rules)



    eltype_cont = mesh.eltype
    offset = mesh.offset
    if mesh.eltype == 0
        eltype_cont = 1
        offset = 1.0
    end

    csis_cont = range(-1,1,length = eltype_cont+1)
    csis_descont = range(-1+offset,1 - offset,length=mesh.eltype+1)

    return csis_cont, csis_descont, rules

end

function calc_nJgp(N, dNc, dNe, gp, field_points)
    normal = []
    J = []
    gauss_points = []
    for i in eachindex(gp)
        normal2, J2 = calc_n_J_matrix(dNc[i], dNe[i], field_points)
        gauss_points2 = N[i]*field_points
        push!(normal,normal2)
        push!(J,J2)
        push!(gauss_points,gauss_points2)
    end
    return normal, J, gauss_points
end

function calc_nearpoints(csis, omegas,c, d, csis_cont, csis_descont,field_points)

    gp_local_near, weights_near = pontos_pesos_local_telles(csis, omegas, c[1:2],d)
    # gp_local_near, weights_near = points_weights_local_near_combined(csis, omegas, c,minimum(d))

    # gp_local_near, weights_near = pontos_pesos_local_subelem(c[1], c[2], Nc_lin, dNcdcsi_lin, dNcdeta_lin, omegas)
    N_near = calc_N_matrix(csis_cont,gp_local_near)
    dNc_near = calc_dNdcsi_matrix(csis_cont,gp_local_near)
    dNe_near = calc_dNdeta_matrix(csis_cont,gp_local_near)
    Nd_near = calc_N_matrix(csis_descont,gp_local_near)
    gauss_points_near = N_near*field_points
    normal_near,J_near = calc_n_J_matrix(dNc_near, dNe_near, field_points)


    return gauss_points_near, Nd_near, J_near, normal_near, weights_near
end

function calc_points_sing(Nc_sing,dNc_sing, dNe_sing,Nd_sing, field_points, Jb_sing, nnel, n, omega_sing)
    
    idx_cont, idx_descont, np = calc_idx_permutation2(nnel,n)
                    
    normal_sing, J_sing = calc_n_J_matrix(dNc_sing[:,idx_cont,np], dNe_sing[:,idx_cont,np], field_points)
    J_sing = J_sing.*Jb_sing[:,np]
    pesos = zeros(length(omega_sing),nnel)
    # @infiltrate
    for kn in 1:nnel
        pesos[:,kn] = J_sing.*omega_sing.*Nd_sing[:,idx_descont[kn],np]
    end

    gauss_points_sing = Nc_sing[:,idx_cont,np]*field_points

    return normal_sing, gauss_points_sing, pesos

end