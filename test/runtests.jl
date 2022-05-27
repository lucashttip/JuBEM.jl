using JuBEM
# using Test

# @testset "JuBEM.jl" begin
#     # Write your tests here.
# end

begin
    import CairoMakie
    using Meshes
    using MeshViz
    using JuBEM

    function visualize_mesh(points, elem)
        conn = connect.([Tuple(elem[i,:]) for i in 1:size(elem,1)])
        ps = [Meshes.Point3(points[i,:]) for i in 1:size(points,1)]
    
        m = SimpleMesh(ps, conn)
    
        return viz(m,showfacets=true)
    end 

    nodes = [
    -1.0 -1.0 0.0
    0.0 -1.0 0.0
    1.0 -1.0 0.0
    1.0 1.0 0.0
    0.0 1.0 0.0
    -1.0 1.0 0.0
    ]

    elem = [
        1 2 5 6
        2 3 4 5
    ]

    csis_cont = [-1.0, 1.0]
    csis_descont = [-0.5, 0.5]

    N = calc_N_matrix(csis_cont,csis_descont)
    k = calc_k(4)

    pd = generate_nodes_in_elem(N,nodes[elem[1,:],:],k)

    # scene = visualize_mesh(nodes, elem)
end

begin
    using JuBEM
    import Plots
    nodes = [
    -1.0 -1.0 0.0
    0.0 -1.0 0.0
    1.0 -1.0 0.0
    1.0 1.0 0.0
    0.0 1.0 0.0
    -1.0 1.0 0.0
    ]

    elem = [
        1 2 5 6
        2 3 4 5
    ]

    csis_cont = [-1.0, 1.0]
    csis_descont = [-0.5, 0.5]

    N = calc_N_matrix(csis_cont,csis_descont)
    k = calc_k(4)

    pd = generate_nodes_in_elem(N,nodes[elem[1,:],:],k)
    Plots.scatter(nodes[elem[1,:],1], nodes[elem[1,:],2], label="elem continuo")
    Plots.scatter!(pd[:,1], pd[:,2], label="elem descontinuo")

end

begin   # Para verificar localização dos pontos de gauss e normais
    fe = 2

    csis_cont = range(-1,1,length = mesh.eltype+1)
    csis_descont = range(-1+mesh.offset,1 - mesh.offset,length=mesh.eltype+1)

    Nc = calc_N_matrix(csis_cont, solver_var.csi)
    dNcdcsi = calc_dNdcsi_matrix(csis_cont, solver_var.csi)
    dNcdeta = calc_dNdeta_matrix(csis_cont,solver_var.csi)
    G = calc_G(csis_cont, csis_descont)
    Nd = calc_N_matrix_descont(Nc,G)
    dNdcsi = calc_N_matrix_descont(dNcdcsi,G)
    dNdeta = calc_N_matrix_descont(dNcdeta,G)

    nnel = (mesh.eltype+1)^2
    k = calc_k(nnel)
    
    nodes = mesh.nodes[mesh.IEN[:,fe],2:end]
    points = mesh.points[mesh.IEN_geo[:,fe],2:end]
    gauss_points = generate_points_in_elem(Nd,nodes)
    normal, J = calc_n_J_matrix(dNdcsi, dNdeta, nodes)
    # normal = reshape(normal,Int64(solver_var.nGP^2),3)


    x = gauss_points[:,1]
    y = gauss_points[:,2]
    z = gauss_points[:,3]


    Plots.scatter(nodes[:,1], nodes[:,2], label="elem descont")
    Plots.scatter!(gauss_points[:,1], gauss_points[:,2], label="gauss points")
    plt1 = Plots.scatter!(points[:,1], points[:,2], label="elem continuo")

    plt2 = Plots.quiver(x,y,z,quiver=(normal[:,1], normal[:,2], normal[:,3]))


end

begin
    using JuBEM, Plots

    mesh,material,problem,solver_var = read_msh("mesh.msh")

    generate_mesh!(mesh)

    plt1 = scatter(mesh.points[mesh.IEN_geo[:,1],2],mesh.points[mesh.IEN_geo[:,1],3])
    plt1 = scatter!(plt1,mesh.nodes[mesh.IEN[:,1],2],mesh.nodes[mesh.IEN[:,1],3])

    plt2 = scatter(mesh.points[:,2],mesh.points[:,3])
    plt2 = scatter!(plt2,mesh.nodes[:,2],mesh.nodes[:,3])
end

begin
    using JuBEM
    using FastGaussQuadrature
    import Plots
    csis, omega = gausslegendre(8)

    csisij = [-1.0, 1.0]

    N = calc_N_matrix(csisij, csis)
    dN = calc_dNdcsi_matrix(csisij,csis)
    dNdeta = calc_dNdeta_matrix(csisij,csis)

    
end

begin
    i = 2
    plt1 = Plots.surface(csis,csis,N[:,:,i],xlabel="csi", ylabel = "eta")
    plt2 = Plots.surface(csis,csis,dN[:,:,i],xlabel="csi", ylabel = "eta")
    plt3 = Plots.surface(csis,csis,dNdeta[:,:,i],xlabel="csi", ylabel = "eta")

    display(plt1)
    display(plt2)
    display(plt3)

end

