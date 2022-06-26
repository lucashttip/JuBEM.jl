function calc_N_matrix(csisij,csis)
    n = size(csis,1)
    ordem = length(csisij)
    nnel = ordem^2

    N = zeros(n,nnel)
    kij = calc_k(nnel)

    for i in 1:n
        for k in 1:nnel
            csi = csis[i,1]
            eta = csis[i,2]
            N[i,k] = calc_N(csisij,ordem,kij[k][1], kij[k][2], csi, eta)
        end
    end

    return N
end

function calc_N_matrix(csisij,csis,etas)

    nl = length(csis)
    
    csis_grid = [repeat(csis',nl)[:] repeat(etas,nl)]

    N = calc_N_matrix(csisij,csis_grid)

    return N
end

function calc_csis_grid(csis)
    nl = length(csis)
    
    csis_grid = [repeat(csis',nl)[:] repeat(csis,nl)]
    return csis_grid
end

function calc_dNdcsi_matrix(csisij,csis)

    n = size(csis,1)
    ordem = length(csisij)
    nnel = ordem^2

    dNdcsi = zeros(n,nnel)
    kij = calc_k(nnel)

    for i in 1:n
            csi = csis[i,1]
            eta = csis[i,2]
            dNdcsi[i,:] = calc_dNdcsi(csisij,ordem,kij, csi, eta)
            # dNdcsi[i,:] = calc_dNdcsi2(eta)
    end

    return dNdcsi
end

function calc_dNdcsi_matrix(csisij,csis, etas)

    nl = length(csis)
    
    csis_grid = [repeat(csis',nl)[:] repeat(etas,nl)]

    dNdcsi = calc_dNdcsi_matrix(csisij,csis_grid)

    return dNdcsi

end

function calc_dNdeta_matrix(csisij,csis)
    n = size(csis,1)
    ordem = length(csisij)
    nnel = ordem^2

    dNdeta = zeros(n,nnel)
    kij = calc_k(nnel)

    for i in 1:n
            csi = csis[i,1]
            eta = csis[i,2]
            # dNdeta[i,:] = calc_dNdeta(csisij,ordem,kij, csi, eta)
            dNdeta[i,:] = calc_dNdeta2(csi)
    end

    return dNdeta
end

function calc_dNdeta_matrix(csisij,csis,etas)

    nl = length(csis)
    
    csis_grid = [repeat(csis',nl)[:] repeat(etas,nl)]

    dNdeta = calc_dNdeta_matrix(csisij,csis_grid)

    return dNdeta
end

function calc_N(csisij,nnel,i, j, csi, eta)
    N = 1

    for l in 1:nnel
        if l != i
            N = N*(csi - csisij[l])/(csisij[i] - csisij[l])
        end
    end
    
    for l in 1:nnel
        if l != j
            N = N*(eta - csisij[l])/(csisij[j] - csisij[l])
        end
    end
    return N
end

function calc_dNdcsi(csisij,ordem,kij,csi,eta)

    nnel = ordem^2

    dNdcsi = zeros(nnel)

    for k in 1:nnel
        dNdcsi[k] = 0
        i = kij[k][1]
        j = kij[k][2]
        
        for m in 1:ordem
            dN = 1
            for l in 1:ordem
                if l != i && l != m
                    dN = dN*(csi - csisij[l])/(csisij[i] - csisij[l])
                end
            end
            if m != i
                dNdcsi[k] = dNdcsi[k] + (dN/(csisij[i] - csisij[m]))
            end
        end

        
        for l in 1:ordem
            if l != j
                dNdcsi[k] = dNdcsi[k]*(eta - csisij[l])/(csisij[j] - csisij[l])
            end
        end
    end
    return dNdcsi
end

function calc_dNdcsi2(eta)
    dNdcsi = (1/4)*[-(1-eta), (1-eta), (1+eta), -(1+eta)]
    return dNdcsi
end

function calc_dNdeta2(csi)
    dNdeta = (1/4)*[-(1-csi), -(1+csi), (1+csi), (1-csi)]
    return dNdeta
end

function calc_dNdeta(csisij,ordem,kij,csi,eta)

    nnel = ordem^2

    dNdeta = zeros(nnel)

    for k in 1:nnel
        i = kij[k][1]
        j = kij[k][2]
        dNdeta[k] = 0

        for m in 1:ordem
            dN = 1
            
            for l in 1:ordem
                if l != j && l != m
                    dN = dN*(eta - csisij[l])/(csisij[j] - csisij[l])
                end
            end
            if m != j
                dNdeta[k] = dNdeta[k] + (dN/(csisij[j] - csisij[m]))
            end
        end

        for l in 1:ordem
            if l != i
                dNdeta[k] = dNdeta[k]*(csi - csisij[l])/(csisij[i] - csisij[l])
            end
        end
        
    end
    
    return dNdeta
end

function calc_k(nnel)

    if nnel == 4
        k = (
            [1,1],
            [2,1],
            [2,2],
            [1,2]
        )
    elseif nnel == 9
        k = (
            [1,1],
            [1,3],
            [3,3],
            [3,1],
            [1,2],
            [2,3],
            [3,2],
            [2,1],
            [2,2]
        )
    else
        error("Element type not supported")
    end
end

function calc_G(csis_cont, csis_descont, k)

    nl = sqrt(length(k))

    idx = [Int(nl*(k[i][1]-1)+k[i][2]) for i in 1:length(k)]

    L = calc_N_matrix(csis_cont, csis_descont, csis_descont)

    L = L[idx,:]

    G = inv(L)

    return G

end

function calc_N_matrix_descont(N,G)
    
    Nd = N*G

    return Nd

end

function calc_Ns(csis_cont, csis_descont, csis, etas)

    nnel = length(csis_cont)^2

    k = calc_k(nnel)

    Nc = calc_N_matrix(csis_cont, csis, etas)
    dNcdcsi = calc_dNdcsi_matrix(csis_cont, csis, etas)
    dNcdeta = calc_dNdeta_matrix(csis_cont, csis, etas)
    G = calc_G(csis_cont, csis_descont,k)
    Nd = Nc*G
    dNdcsi = dNcdcsi*G
    dNdeta = dNcdeta*G

    return Nc, Nd, dNdcsi, dNdeta, dNcdcsi, dNcdeta
end

function calc_Ns_sing(csis_cont, csis_descont, csis_sing)

    nnel = length(csis_cont)^2

    k = calc_k(nnel)

    Nc = calc_N_matrix2(csis_cont, csis_sing)
    dNcdcsi = calc_dNdcsi_matrix2(csis_cont, csis_sing)
    dNcdeta = calc_dNdeta_matrix2(csis_cont, csis_sing)
    G = calc_G(csis_cont, csis_descont,k)
    Nd = calc_N_matrix_descont(Nc,G)
    dNdcsi = calc_N_matrix_descont(dNcdcsi,G)
    dNdeta = calc_N_matrix_descont(dNcdeta,G)

    return Nc, Nd, dNdcsi, dNdeta
end

function calc_omegas(omega)
    nGP = length(omega)
    omegas = zeros(nGP*nGP)
    for i in 1:nGP
        for j in 1:nGP
            omegas[nGP*(i-1)+j] = omega[i]*omega[j]
        end
    end
    return omegas
end