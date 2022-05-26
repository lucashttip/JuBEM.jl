function calc_N_matrix(csisij,csis)
    nGP = length(csis)
    ordem = length(csisij)
    nnel = ordem^2

    N = zeros(nGP,nGP,nnel)
    kij = calc_k(nnel)
    
    for i in 1:nGP
        for j in 1:nGP
            for k in 1:nnel
                N[i,j,k] = calc_N(csisij,ordem,kij[k][1], kij[k][2], csis[i], csis[j])
            end
        end
    end

    return N
end

function calc_dNdcsi_matrix(csisij,csis)
    nGP = length(csis)
    ordem = length(csisij)
    nnel = ordem^2

    dNdcsi = zeros(nGP,nGP,nnel)
    kij = calc_k(nnel)

    for i in 1:nGP
        for j in 1:nGP
            for k in 1:nnel
                dNdcsi[i,j,k] = calc_dNdcsi(csisij,ordem,kij[k][1], kij[k][2], csis[i], csis[j])
            end
        end
    end

    return dNdcsi
end

function calc_dNdeta_matrix(csisij,csis)
    nGP = length(csis)
    ordem = length(csisij)
    nnel = ordem^2

    dNdcsi = zeros(nGP,nGP,nnel)
    kij = calc_k(nnel)

    for i in 1:nGP
        for j in 1:nGP
            for k in 1:nnel
                dNdcsi[i,j,k] = calc_dNdeta(csisij,ordem,kij[k][1], kij[k][2], csis[i], csis[j])
            end
        end
    end

    return dNdcsi
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

function calc_dNdcsi(csisij,ordem,i,j,csi,eta)

    dNdcsi = 0
    
    if ordem > 2
        for m in 1:ordem
            dN = 1
            for l in 1:ordem
                if l != i && l != m
                    dN = dN*(csi - csisij[l])/(csisij[i] - csisij[l])
                end
            end
            if m != i
                dNdcsi = dNdcsi + (dN/(csisij[i] - csisij[m]))
            end
        end
    else
        dNdcsi =1
    end
    
    for l in 1:ordem
        if l != j
            dNdcsi = dNdcsi*(eta - csisij[l])/(csisij[j] - csisij[l])
        end
    end
    return dNdcsi
end

function calc_dNdeta(csisij,ordem,i,j,csi,eta)

    dNdeta = 1
    

    for l in 1:ordem
        if l != i
            dNdeta = dNdeta*(csi - csisij[l])/(csisij[i] - csisij[l])
        end
    end


    if ordem >2
        for m in 1:ordem
            dN = 1
            
            for l in 1:ordem
                if l != j && l != m
                    dN = dN*(eta - csisij[l])/(csisij[j] - csisij[l])
                end
            end
            if m != j
                dNdeta = dNdeta + (dN/(csisij[j] - csisij[m]))
            end
        end
    end
    
    return dNdeta
end

function calc_u_in_N(N, u)
    lx = size(N,1)
    ly = size(N,2)
    u2 = zeros(lx,ly)

    for i in 1:lx
        for j in 1:ly
            u2[i,j] = N[i,j,:]'*u
        end
    end

    return u2
end

function calc_k(nnel)

    if nnel == 4
        k = (
            [1,1],
            [1,2],
            [2,2],
            [2,1]
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

function calc_G(csis_cont, csis_descont)
    N = calc_N_matrix(csis_cont, csis_descont)

    nnel = length(csis_descont)^2

    L = remap_N(N,nnel)

    G = inv(L)

    return G

end

function calc_N_matrix_descont(N,G)
    Nd = zeros(size(N))

    for i in 1:size(N,1)
        for j in 1:size(N,2)
            Nd[i,j,:] = N[i,j,:]'*G
        end
    end

    return Nd

end

function remap_N(N,nnel)
    Nn = zeros(nnel,nnel)
    
    a,b,c = size(N)
    kij = calc_k(nnel)

    for k in 1:nnel
        for j in 1:nnel
            Nn[j,k] = N[kij[k][1], kij[k][2],j]  
        end
    end
    return Nn
end

