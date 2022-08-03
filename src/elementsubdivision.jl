function csis_sing(offset,Nc,omegas)


    p1 = [-1 -1]
    p2 = [1 -1]
    p3 = [1 1]
    p4 = [-1 1]
    p5 = [-1+offset -1+offset]

    q1 = [p1;p2;p5;p5]
    q2 = [p2;p3;p5;p5]
    q3 = [p3;p4;p5;p5]
    q4 = [p4;p1;p5;p5]

    csi_sing = generate_points_in_elem(Nc,q1)
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q2)]
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q3)]
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q4)]

    return csi_sing

end

function csis_sing2(offset,Nc,dNcdcsi,dNcdeta)

    p1 = [-1 -1]
    p2 = [1 -1]
    p3 = [1 1]
    p4 = [-1 1]
    p5 = [-1+offset -1+offset]

    q1 = [p1;p2;p5;p5]
    q2 = [p2;p3;p5;p5]
    q3 = [p3;p4;p5;p5]
    q4 = [p4;p1;p5;p5]

    csi_sing = generate_points_in_elem(Nc,q1)
    # @infiltrate
    _, J = calc_n_J_matrix(dNcdcsi, dNcdeta, [q1 zeros(4)])
    Jb_sing = J
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q2)]
    _, J = calc_n_J_matrix(dNcdcsi, dNcdeta, [q2 zeros(4)])
    Jb_sing = [Jb_sing; J]
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q3)]
    _, J = calc_n_J_matrix(dNcdcsi, dNcdeta, [q3 zeros(4)])
    Jb_sing = [Jb_sing; J]
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q4)]
    _, J = calc_n_J_matrix(dNcdcsi, dNcdeta, [q4 zeros(4)])
    Jb_sing = [Jb_sing; J]

    return csi_sing, Jb_sing

end

function csis_sing3(offset,Nc,dNcdcsi,dNcdeta)

    c, IEN = divide_elem2(offset)

    nsubelem = size(IEN,2)

    csi_sing = generate_points_in_elem(Nc,c[IEN[:,1],:])
    _, J = calc_n_J_matrix(dNcdcsi, dNcdeta, [c[IEN[:,1],:] zeros(4)])
    Jb_sing = J
    
    for i in 2:nsubelem
        csi_sing = [csi_sing;generate_points_in_elem(Nc,c[IEN[:,i],:])]
        _, J = calc_n_J_matrix(dNcdcsi, dNcdeta, [c[IEN[:,i],:] zeros(4)])
        Jb_sing = [Jb_sing; J]
    end
    return csi_sing, Jb_sing

end

function divide_elem(source_node,points,Nc,dNcdcsi,dNcdeta,omegas)

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


function divide_elem(e)
    c = [
        -1 -1
        1 -1
        1 1
        -1 1
        -1+e -1+e
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

function calc_idx_permutation(nnel,n)
    if nnel == 4
        idx_ff = collect((1:nnel) .- (n-1))
        idx_ff[idx_ff.<1] = idx_ff[idx_ff.<1] .+nnel

        idx_points = collect((1:nnel) .+ (n-1))
        idx_points[idx_points.>nnel] = idx_points[idx_points.>nnel] .-nnel
    end

    return idx_ff, idx_points
end