"""
    csis2d, (omegas2d) = calc_csis2d(csis1d,(omegas1d))

This function takes a set of gaussian quadrature points from -1 to 1 (csis1d) and their weights(omegas1d) and maps them to two dimensions returning a matrix os the positions and a vector of the new weights.
the returned points are csis2d = [ξi, ηj], and the returned weights are omegas2d = [omegai*omegaj]

it goes 
csis_grid = [ξ1, η1; ξ1, η2; ...; ξ1, ηn; ξ2, η1; ξ2, η2; ...; ξn, ηn]

"""
function calc_csis2d(csis1d,omegas1d)
    # nl = length(csis1d)
    
    # csis_grid = [repeat(csis1d',nl)[:] repeat(csis1d,nl)]

    nGP = length(omegas1d)
    omegas2d = zeros(nGP*nGP)
    csis2d = zeros(nGP*nGP,2)
    for i in 1:nGP
        for j in 1:nGP
            csis2d[nGP*(i-1)+j,:] = [csis1d[i],csis1d[j]]
            omegas2d[nGP*(i-1)+j] = omegas1d[i]*omegas1d[j]
        end
    end

    return csis2d, omegas2d
end

function calc_csis2d(csis1d)
    nGP = length(csis1d)
    csis2d = zeros(nGP*nGP,2)
    for i in 1:nGP
        for j in 1:nGP
            csis2d[nGP*(i-1)+j,:] = [csis1d[i],csis1d[j]]
        end
    end

    return csis2d
end


function map_k(nnel)
    if nnel == 1
        k = ([1,1],)
    elseif nnel == 4
        # k = [(i,j)] mapping
        #        i1  i2
        # j = 2: 4 - 3
        #        |   |
        # j = 1: 1 - 2 
        k = (
            [1,1],
            [2,1],
            [2,2],
            [1,2]
        )
    elseif nnel == 9
        # k = [(i,j)] mapping
        #        i1  i2  i3
        # j = 3: 4 - 7 - 3
        #        |       |
        # j = 2: 8   9   6
        #        |       |
        # j = 1: 1 - 5 - 2 
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

function calc_N_gen(nodal_csis,csis_grid;dg = :N)

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

    n = size(csis_grid,1)
    ordem = length(nodal_csis)
    nnel = ordem^2

    N = zeros(n,nnel)
    kij = map_k(nnel)

    for i in 1:n
        csi = csis_grid[i,1]
        eta = csis_grid[i,2]
        N[i,:] = f(nodal_csis,ordem,kij, csi, eta)
    end

    return N


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