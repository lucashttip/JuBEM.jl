function bilinear_strat1(nodes, offset,N, csi, omega)
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

    csi_sing = []

    push!(csi_sing,generate_points_in_elem(N,q1))
    push!(csi_sing,generate_points_in_elem(N,q2))
    push!(csi_sing,generate_points_in_elem(N,q3))
    push!(csi_sing,generate_points_in_elem(N,q4))

    gauss_points = generate_points_in_elem(Nd,field_nodes)


    dNdcsi_sing = 
    dNdeta_sing = 

    return csi_sing
end