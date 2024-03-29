function calc_csis_grid(csis)
    nl = length(csis)
    
    csis_grid = [repeat(csis',nl)[:] repeat(csis,nl)]
    return csis_grid
end

function calc_k(nnel)
    if nnel == 1
        k = ([1,1],)
    elseif nnel == 4
        k = (
            [1,1],
            [2,1],
            [2,2],
            [1,2]
        )
    elseif nnel == 9
        k = (
            [1,1],
            [3,1],
            [3,3],
            [1,3],
            [2,1],
            [3,2],
            [2,3],
            [1,2],
            [2,2]
        )
    else
        error("Element type not supported")
    end
end

function calc_N_matrix(csisij,csis)
    n = size(csis,1)
    ordem = length(csisij)
    nnel = ordem^2

    N = zeros(n,nnel)
    kij = calc_k(nnel)

    for i in 1:n
        csi = csis[i,1]
        eta = csis[i,2]
        N[i,:] = calc_N(csisij,ordem,kij, csi, eta)
    end

    return N
end

function calc_N_gen(csisij,csis;dg = :N)

    if dg == :N
        f = calc_N
    elseif dg == :dNdc
        f = calc_dNdcsi
    elseif dg == :dNde
        f = calc_dNdeta
    elseif dg == :d2Ndc2
        f = calc_dNdcsicsi
    elseif dg == :d2Ndce
        f = calc_dNdcsieta
    else
        error("dg not defined. Defined symbols are :N, :dNdc, :dNde, :d2Ndc2 and :d2Ndce. Symbol passed: " ,dg)
    end

    n = size(csis,1)
    ordem = length(csisij)
    nnel = ordem^2

    N = zeros(n,nnel)
    kij = calc_k(nnel)

    for i in 1:n
        csi = csis[i,1]
        eta = csis[i,2]
        N[i,:] = f(csisij,ordem,kij, csi, eta)
    end

    return N


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
    end

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
            dNdeta[i,:] = calc_dNdeta(csisij,ordem,kij, csi, eta)
    end

    return dNdeta
end

function calc_N(csisij,ordem,kij, csi, eta)
    nnel = ordem^2

    N = ones(nnel)

    for k in 1:nnel
        i = kij[k][1]
        j = kij[k][2]

        for l in 1:ordem
            if l != i
                N[k] = N[k]*(csi - csisij[l])/(csisij[i] - csisij[l])
            end
        end
        
        for l in 1:ordem
            if l != j
                N[k] = N[k]*(eta - csisij[l])/(csisij[j] - csisij[l])
            end
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

function calc_dNdcsicsi(csisij,ordem,kij,csi,eta)

    nnel = ordem^2

    dNdcsi = zeros(nnel)

    for k in 1:nnel
        dNdcsi[k] = 0
        i = kij[k][1]
        j = kij[k][2]
        


        for n in 1:ordem
            dN1 = 0
            for m in 1:ordem
                dN2 = 1
                for l in 1:ordem
                    if l != i && l != m && l != n
                        dN2 = dN2*(csi - csisij[l])/(csisij[i] - csisij[l])
                    end
                end
                if m != i && m != n
                    dN1 = dN1 + (dN2/(csisij[i] - csisij[m]))
                end
            end

            if n != i
                dNdcsi[k] = dNdcsi[k] + (dN1/(csisij[i] - csisij[n]))
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

function calc_dNdcsieta(csisij,ordem,kij,csi,eta)

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


        dN1 = 0
        for m in 1:ordem
            dN = 1
            for l in 1:ordem
                if l != j && l != m
                    dN = dN*(csi - csisij[l])/(csisij[j] - csisij[l])
                end
            end
            if m != j
                dN1 = dN1 + (dN/(csisij[j] - csisij[m]))
            end
        end
        dNdcsi[k] = dNdcsi[k]*dN1

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

function calc_G(csis_cont, csis_descont, k)

    nl = sqrt(length(k))

    idx = [Int(nl*(k[i][1]-1)+k[i][2]) for i in 1:length(k)]

    csis_vec_decont = calc_csis_grid(csis_descont)

    L = calc_N_matrix(csis_cont, csis_vec_decont)

    L = L[idx,:]

    G = inv(L)

    return G

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