function csis_sing(offset,Nc,dNcdcsi,dNcdeta,eltype)

    if eltype == 1
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
            csi_sing[npg*(i-1)+1:npg*i,:,n] = generate_points_in_elem(Nc,c[j][IEN[:,i],:])
            _, Jb_sing[npg*(i-1)+1:npg*i,n] = calc_n_J_matrix(dNcdcsi, dNcdeta, [c[j][IEN[:,i],:] zeros(4)])
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
        -1+e 0],
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

function calc_idx_permutation(nnel,n)
    if nnel == 4
        idx_ff = collect((1:nnel) .- (n-1))
        idx_ff[idx_ff.<1] = idx_ff[idx_ff.<1] .+nnel

        idx_points = collect((1:nnel) .+ (n-1))
        idx_points[idx_points.>nnel] = idx_points[idx_points.>nnel] .-nnel
    end

    return idx_ff, idx_points
end