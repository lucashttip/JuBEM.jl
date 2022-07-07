function csis_sing(offset,Nc,omegas)

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