function csis_sing(offset,Nc,dNcdcsi,dNcdeta,eltype)

    if eltype <= 1
        c, IEN = divide_elem_lin(offset)
        n = 1
    elseif eltype == 2
        c, IEN = divide_elem_quad(offset)
        n = 3
    end

    nsubelem = size(IEN,2)
    npg = size(Nc,1)
    nnel = (eltype+1)^2
    csi_sing = zeros(nsubelem*npg,2,n)
    Jb_sing = zeros(nsubelem*npg,n)

    for j in 1:n       
        for i in 1:nsubelem
            csi_sing[npg*(i-1)+1:npg*i,:,j] = generate_points_in_elem(Nc,c[j][IEN[:,i],:])
            _, Jb_sing[npg*(i-1)+1:npg*i,j] = calc_n_J_matrix(dNcdcsi, dNcdeta, [c[j][IEN[:,i],:] zeros(4)])
        end
    end

    return csi_sing, Jb_sing

end

function divide_elem_lin(e)
    c = [[
        -1 -1
        1 -1
        1 1
        -1 1
        -1+e -1+e
    ]]

    IEN = [
        1 2 3 4
        2 3 4 1
        5 5 5 5
        5 5 5 5
    ]

    return c, IEN
end

function divide_elem_quad(e)
    c = [
        [-1 -1
        1 -1
        1 1
        -1 1
        -1+e -1+e],
        [-1 -1
        1 -1
        1 1
        -1 1 
        0 -1+e],
        [-1 -1
        1 -1
        1 1
        -1 1 
        0 0]
    ]

    IEN = [
        1 2 3 4
        2 3 4 1
        5 5 5 5
        5 5 5 5
    ]
    return c, IEN
end

function divide_elem_quasising()

    c = [[
        -1 -1
        1 -1
        1 1
        -1 1
        0 0
        0 -1
        1 0
        0 1
        -1 0
    ]]

    IEN = [
        1 6 5 9
        6 2 7 5
        5 7 3 8
        9 5 8 4
    ]

    return c, IEN
end

function divide_elem(s)

    c = [
        -1 -1
        1 -1
        1 1
        -1 1
        s[1] s[2]
    ]        

    IEN = [
        1 2 3 4
        2 3 4 1
        5 5 5 5
        5 5 5 5
    ]


    if all(abs.(s) .==  1.0)
        if s == c[1,:]
            IEN = IEN[:,[2,3]]
        elseif s == c[2,:]
            IEN = IEN[:,[3,4]]
        elseif s == c[3,:]
            IEN = IEN[:,[4,1]]
        elseif s == c[4,:]
            IEN = IEN[:,[1,2]]
        end

    elseif any(abs.(s) .==  1.0)

        if s[1] == -1
            IEN = IEN[:,[1,2,3]]
        elseif s[1] == 1
            IEN = IEN[:,[1,3,4]]
        elseif s[2] == -1
            IEN = IEN[:,[2,3,4]]
        elseif s[2] == 1
            IEN = IEN[:,[1,2,4]]
        end
    end
    return c, IEN
end

function csis_quasi_sing(Nc)
    c, IEN = divide_elem_quasising()

    nsubelem = size(IEN,2)
    npg = size(Nc,1)
    csi = zeros(nsubelem*npg,2)

    for i in 1:nsubelem
        csi[npg*(i-1)+1:npg*i,:] = generate_points_in_elem(Nc,c[1][IEN[:,i],:])
    end


    return csi
end

function csis_quasi_sing(Nc, nsubx, nsuby)
    p0 = [-1 -1]
    d = 2
    dx = d/nsubx
    dy = d/nsuby

    npg = size(Nc,1)
    csi = zeros(nsubx*nsuby*npg,2)
    
    k = 1
    for i in 1:nsubx
        for j in 1:nsuby
            p1 = p0+[(i-1)*dx (j-1)*dy]
            p2 = p0+[(i)*dx (j-1)*dy]            
            p3 = p0+[(i)*dx (j)*dy]            
            p4 = p0+[(i-1)*dx (j)*dy]            

            p = [p1;p2;p3;p4]
            csi[npg*(k-1)+1:npg*k,:] = generate_points_in_elem(Nc,p)
            k +=1
        end
    end


    return csi
end

function divide_elem2(e)
    c = [
        -1 -1
        1 -1
        1 1
        -1 1
        -1+e -1+e
        -1+2e -1
        1 -1+2e
        -1+2e 1
        -1 -1+2e
        -1+2e -1+2e
    ]

    IEN = [
        1 6 10 9 6 10 9
        6 10 9 1 2 7 10
        5 5 5 5 7 3 8
        5 5 5 5 10 8 4
    ]

    return c, IEN
end

function calc_subelems(source_node,points,Nc,dNcdcsi,dNcdeta,omegas)

    npg2 = size(dNcdcsi,1)
    At = calc_area(points)

    normal_sing = zeros(4*npg2,3)
    J_sing = zeros(4*npg2)
    omega_sing = zeros(4*npg2)
    gauss_points_sing = zeros(4*npg2,3)


    for i in 1:4

        idx = i != 4 ? [i,i+1] : [i,1]

        ps = [points[idx,:]; repeat(source_node',2)]
        a = calc_area(ps)
        n1, J1 = calc_n_J_matrix(dNcdcsi, dNcdeta, ps)

        gauss_points_sing[(i-1)*npg2+1:i*npg2,:] = Nc*ps
        # omega_sing[(i-1)*npg2+1:i*npg2] = omegas.*a./At
        omega_sing[(i-1)*npg2+1:i*npg2] = omegas

        normal_sing[(i-1)*npg2+1:i*npg2,:] = n1
        J_sing[(i-1)*npg2+1:i*npg2] = J1

    end
    return normal_sing, J_sing, omega_sing, gauss_points_sing
end

function calc_idx_permutation2(nnel,n)

    if nnel == 1
        idx_cont = [1,2,3,4]
        idx_descont = [1]
        nperm = 1

    elseif nnel == 4
        idx_cont = collect((1:nnel) .- (n-1))
        idx_cont[idx_cont.<1] = idx_cont[idx_cont.<1] .+nnel
        idx_descont = idx_cont

        nperm = 1

    elseif nnel == 9
        idx = [0,1,2,3,0,1,2,3,0]
        idx_ff1 = collect((1:4) .- idx[n])
        idx_ff1[idx_ff1.<1] = idx_ff1[idx_ff1.<1] .+4
        idx_ff2 = collect((5:8) .- idx[n])
        idx_ff2[idx_ff2.<5] = idx_ff2[idx_ff2.<5] .+4
        idx_cont = [idx_ff1; idx_ff2; 9]
        idx_descont = idx_cont
        
        if n < 5
            nperm = 1
        elseif n>=5 && n<=8
            nperm = 2
        elseif n == 9
            nperm = 3
        end
    end
    return idx_cont,idx_descont, nperm
    
end

function calc_idx_permutation(nnel,n)

    if nnel == 1
        idx_ff = [1]
        nperm = 1

    elseif nnel == 4
        idx_ff = collect((1:nnel) .- (n-1))
        idx_ff[idx_ff.<1] = idx_ff[idx_ff.<1] .+nnel

        nperm = 1

    elseif nnel == 9
        idx = [0,1,2,3,0,1,2,3,0]
        idx_ff1 = collect((1:4) .- idx[n])
        idx_ff1[idx_ff1.<1] = idx_ff1[idx_ff1.<1] .+4
        idx_ff2 = collect((5:8) .- idx[n])
        idx_ff2[idx_ff2.<5] = idx_ff2[idx_ff2.<5] .+4
        idx_ff = [idx_ff1; idx_ff2; 9]
        
        if n < 5
            nperm = 1
        elseif n>=5 && n<=8
            nperm = 2
        elseif n == 9
            nperm = 3
        end
    end
    return idx_ff, nperm
    
end

function calc_points_weights()

    subelem = [(1,1),(1,1),(1,1),(1,1)]
    npg = [4,5,6,8]
    dists = [4,2,0.5,0.2]
    pw = gausslegendre.(npg)

    gp = gauss_points_type[]

    for i in eachindex(npg)
        csi = calc_csis_grid(pw[i][1])
        Nc = calc_N_matrix([-1,1],csi)
        csis = csis_quasi_sing(Nc,subelem[i][1],subelem[i][2])
        fact = subelem[i][1]*subelem[i][2]
        omegas = repeat(calc_omegas(pw[i][2])./fact,fact)
        push!(gp,gauss_points_type(csis,omegas))
    end
    return gp, dists
end

function pontos_pesos_local_subelem(csi_source, eta_source,Nc_lin,dNcdcsi_lin,dNcdeta_lin,omegas)
    points = [-1 -1
        1 -1
        1 1
        -1 1]
    source = [csi_source eta_source]

    npg = length(omegas)
    pontos_gauss = zeros(4*npg,2)
    pesos = zeros(4*npg)

    p1p2 = [
        1 2
        2 3
        3 4
        4 1
    ]

    for t in 1:4
        points2 = [points[p1p2[t,:],:];source; source]
        i1 = npg*(t-1) +1
        i2 = npg*t
        pontos_gauss[i1:i2,:] = Nc_lin*points2
        _,J = calc_n_J_matrix(dNcdcsi_lin, dNcdeta_lin, [points2 zeros(4)])
        pesos[i1:i2] = J.*omegas

    end
    return pontos_gauss, pesos
end