function calc_CD_discont(mesh, solver_var, rbcenter, rbidx)

    rbelem = findall(mesh.bc[:,1].==rbidx)

    nnel = (mesh.eltype+1)^2
    omegas = calc_omegas(solver_var.omega)
    csi_grid = calc_csis_grid(solver_var.csi)
    if mesh.eltype == 0
        csis_cont = range(-1,1,length = 2)
        Nd = ones(size(csi_grid,1))
    else
        csis_cont = range(-1,1,length = mesh.eltype+1)
        csis_descont = range(-1+mesh.offset,1 - mesh.offset,length=mesh.eltype+1)
        Nd = calc_N_matrix(csis_descont,csi_grid)
    end
    Nc = calc_N_matrix(csis_cont,csi_grid)
    dNcdcsi = calc_dNdcsi_matrix(csis_cont, csi_grid)
    dNcdeta = calc_dNdeta_matrix(csis_cont, csi_grid)

    nrbelem = length(rbelem)

    nrbnodes = nrbelem*nnel
    C = zeros(3*nrbnodes,6)
    D = zeros(6,3*nrbnodes)

    for e in 1:nrbelem
        # Calculating C
        rbe = rbelem[e]
        for n in 1:nnel
            node = mesh.nodes[mesh.IEN[n,rbe],2:end]
            r = node - rbcenter
            idx1 = 3*(nnel*(e-1)+n-1)+1
            idx2 = 3*(nnel*(e-1)+n)
            C[idx1:idx2,:] = [
                1 0 0 0 r[3] -r[2]
                0 1 0 -r[3] 0 r[1]
                0 0 1 r[2] -r[1] 0
            ]
        end

        # Calculating D
        points = mesh.points[mesh.IEN_geo[:,rbe],2:end]
        _,J = calc_n_J_matrix(dNcdcsi, dNcdeta, points)
        gauss_points = Nc*points
        De = view(D,:,3*nnel*(e-1)+1:3*nnel*(e))
        for g in eachindex(J)
            r = gauss_points[g,:] - rbcenter

            Dg = [
                1 0 0
                0 1 0
                0 0 1
                0 -r[3] r[2]
                r[3] 0 -r[1]
                -r[2] r[1] 0
            ]

            for k in 1:nnel
                De[:,3*(k-1)+1:3*k] = De[:,3*(k-1)+1:3*k] + Dg.*J[g].*Nd[g,k]*omegas[g]
            end
        end
    end
    return C, D
end