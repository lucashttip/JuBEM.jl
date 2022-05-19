function calc_N_matriz(csisij,csis)
    nGP = length(csis)
    ordem = length(csisij)
    nnel = ordem^2

    N = zeros(nGP,nGP,nnel)
    kij = calc_k(nnel)

    println(nnel)
    println(kij)
    println(N)
    for i in 1:nGP
        for j in 1:nGP
            for k in 1:nnel
                N[i,j,k] = calc_N(csisij,ordem,kij[k][1], kij[k][2], csis[i], csis[j])
            end
        end
    end

    return N
end

function calc_dNdcsi_matriz(csisij,csis)
    nGP = length(csis)
    ordem = length(csisij)
    nnel = ordem^2

    dNdcsi = zeros(nGP,nGP,nnel)
    kij = calc_k(nnel)

    println(nnel)
    println(kij)
    println(N)
    for i in 1:nGP
        for j in 1:nGP
            for k in 1:nnel
                dNdcsi[i,j,k] = calc_dNdcsi(csisij,ordem,kij[k][1], kij[k][2], csis[i], csis[j])
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
    
    for l in 1:ordem
        if l != j
            dNdcsi = dNdcsi*(eta - csisij[l])/(csisij[j] - csisij[l])
        end
    end
    return dNdcsi
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

using FastGaussQuadrature

csis, omega = gausslegendre(8)

csisij = [-1.0,0.0, 1.0]

N = calc_N_matriz(csisij, csis)
dNdcsi = calc_dNdcsi_matriz(csisij,csis)
