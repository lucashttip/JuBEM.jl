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

