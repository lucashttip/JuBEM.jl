import CairoMakie
# import GLMakie
using Meshes
using MeshViz

function visualize_mesh(mesh)
    points = mesh.points[:,2:end]
    elem = mesh.IEN_geo'
    conn = connect.([Tuple(elem[i,:]) for i in 1:size(elem,1)])
    ps = [Meshes.Point3(points[i,:]) for i in 1:size(points,1)]

    m = SimpleMesh(ps, conn)

    viz(m,showfacets=true)

    # return m
end 