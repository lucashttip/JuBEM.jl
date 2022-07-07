function csis_sing_1(offset , Nc, omegas)
    l1 = 2*offset

    p1 = [-1 -1]
    p2 = [1 -1]
    p3 = [1 1]
    p4 = [-1 1]
    p5 = [-1+l1 -1]
    p6 = [1 -1+l1]
    p7 = [-1+l1 1]
    p8 = [-1 -1+l1]
    p9 = [-1+l1 -1+l1]

    q1 = [p1;p5;p9;p8]
    q2 = [p5;p2;p6;p9]
    q3 = [p9;p6;p3;p7]
    q4 = [p8;p9;p7;p4]

    a2 = calc_area(q2)
    a3 = calc_area(q3)
    a4 = calc_area(q4)
    
    csi_sing = generate_points_in_elem(Nc,q2)
    omegas_sing = omegas.*a2./4
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q3)]
    omegas_sing = [omegas_sing; omegas.*a3./4]
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q4)]
    omegas_sing = [omegas_sing; omegas.*a4./4]


    return csi_sing, omegas_sing
end

function csis_sing_2(offset , Nc, omegas)
    l1 = 2*offset

    p1 = [-1 -1]
    p2 = [1 -1]
    p3 = [1 1]
    p4 = [-1 1]
    p5 = [-1+l1 -1]
    p6 = [1 -1+l1]
    p7 = [-1+l1 1]
    p8 = [-1 -1+l1]
    p9 = [-1+l1 -1+l1]
    p10 = [-1+offset -1+offset]

    q2 = [p5;p2;p6;p9]
    q3 = [p9;p6;p3;p7]
    q4 = [p8;p9;p7;p4]

    q11 = [p1; p5;p10;p10]
    q12 = [p5;p9;p10;p10]
    q13 = [p9;p8;p10;p10]
    q14 = [p8;p1;p10;p10]


    a2 = calc_area(q2)
    a3 = calc_area(q3)
    a4 = calc_area(q4)

    a11 = calc_area(q11)
    a12 = calc_area(q12)
    a13 = calc_area(q13)
    a14 = calc_area(q14)
    
    csi_sing = generate_points_in_elem(Nc,q2)
    omegas_sing = omegas.*a2./4
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q3)]
    omegas_sing = [omegas_sing; omegas.*a3./4]
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q4)]
    omegas_sing = [omegas_sing; omegas.*a4./4]

    csi_sing = [csi_sing;generate_points_in_elem(Nc,q11)]
    omegas_sing = [omegas_sing; omegas.*a11./4]
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q12)]
    omegas_sing = [omegas_sing; omegas.*a12./4]
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q13)]
    omegas_sing = [omegas_sing; omegas.*a13./4]
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q14)]
    omegas_sing = [omegas_sing; omegas.*a14./4]


    return csi_sing, omegas_sing
end

function csis_sing_3(offset,Nc,omegas)

    l1 = 2*offset

    p1 = [-1 -1]
    p2 = [1 -1]
    p3 = [1 1]
    p4 = [-1 1]
    p5 = [-1+offset -1+offset]

    q1 = [p1;p2;p5;p5]
    q2 = [p2;p3;p5;p5]
    q3 = [p3;p4;p5;p5]
    q4 = [p4;p1;p5;p5]

    a1 = calc_area(q1)
    a2 = calc_area(q2)
    a3 = calc_area(q3)
    a4 = calc_area(q4)

    csi_sing = generate_points_in_elem(Nc,q1)
    omegas_sing = omegas.*a1./4
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q2)]
    omegas_sing = [omegas_sing; omegas.*a2./4]
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q3)]
    omegas_sing = [omegas_sing; omegas.*a3./4]
    csi_sing = [csi_sing;generate_points_in_elem(Nc,q4)]
    omegas_sing = [omegas_sing; omegas.*a4./4]

    return csi_sing, omegas_sing


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


function divide_elem2(e)
    c = [
        -1      -1
        1       -1
        1       1
        -1      1
        -1+2e   -1
        1       -1+2e
        -1+2e   1
        -1      -1+2e
        -1+2e   -1+2e
        -1+e    -1+e
    ]    

    IEN = [
        1 5 9 8 5 9 8
        5 9 8 1 2 6 9
        10 10 10 10 6 3 7
        10 10 10 10 9 7 4
    ]

    return c, IEN
end