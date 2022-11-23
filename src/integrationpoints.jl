function calc_N_nonsing(mesh)

    eltype_fis = mesh.eltype
    eltype_geom = mesh.eltype
    offset = mesh.offset

    if mesh.eltype == 0
        offset = 1.0
        eltype_geom = 1
    end

    csis_cont = range(-1,1,length = eltype_geom+1)
    csis_descont = range(-1+offset,1 - offset,length=eltype_fis+1)

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

    return Nc_nonsing, dNcdcsi_nonsing, dNcdeta_nonsing, Nd_nonsing, gp, dists
end


function calc_N_sing(mesh, solver_var)

    eltype_fis = mesh.eltype
    eltype_geom = mesh.eltype
    offset = mesh.offset

    if mesh.eltype == 0
        offset = 1.0
        eltype_geom = 1
    end

    csis_cont = range(-1,1,length = eltype_geom+1)
    csis_descont = range(-1+offset,1 - offset,length=eltype_fis+1)
    omegas = calc_omegas(solver_var.omega)


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

end