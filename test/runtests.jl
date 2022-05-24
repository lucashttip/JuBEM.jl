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

    N = calc_N_matriz(csis_cont,csis_descont)
    k = calc_k(9)

    pd = generate_nodes_in_elem(N,nodes[elem[1,:],:],k)

    # scene = visualize_mesh(nodes, elem)
end

begin
    using JuBEM

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

    csisij = [-1.0,0.0, 1.0]

    N = calc_N_matriz(csisij, csis)
    dN = calc_dNdcsi_matriz(csisij,csis)
end

begin
    plt1 = Plots.surface(csis,csis,N[:,:,9])
    plt2 = Plots.surface(csis,csis,dN[:,:,9])
    display(plt1)
    display(plt2)
end

