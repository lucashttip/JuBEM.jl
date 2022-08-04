# using CairoMakie
import GLMakie
import GeometryBasics
const mk = GLMakie
const gb = GeometryBasics
# const mk = CairoMakie
# using Meshes
# using MeshViz


# function visualize_mesh(mesh)
#     points = mesh.points[:,2:end]
#     elem = mesh.IEN_geo'

#     visualize_mesh_raw(points, elem)
#     # return m
# end 

# function visualize_mesh_raw(points, elem)
#     conn = connect.([Tuple(elem[i,:]) for i in 1:size(elem,1)])
#     ps = [Meshes.Point3(points[i,:]) for i in 1:size(points,1)]

#     m = SimpleMesh(ps, conn)

#     viz(m,showfacets=true)
# end


# function visualize_result(mesh, u, idx)
    


#     points = mesh.points[:,2:end]
#     elem = mesh.IEN_geo'
#     conn = connect.([Tuple(elem[i,:]) for i in 1:size(elem,1)])
#     ps = [Meshes.Point3(points[i,:]) for i in 1:size(points,1)]

#     m = SimpleMesh(ps, conn)

#     # return CairoMakie.mesh(m,color=u_points[:,idx],showfacets=true)
#     return CairoMakie.mesh(m)
    

# end

function view_mesh(mesh)
    # fig = mk.Figure()
    # ax = mk.Axis(fig[1, 1])
    fig,ax,plt = mk.scatter(mesh.points[:,2],mesh.points[:,3],mesh.points[:,4], markersize=80, color = :blue)
    # mk.scatter!(ax,mesh.points[:,2],mesh.points[:,3],mesh.points[:,4], markersize=70, color = :blue)

    mk.scatter!(ax,mesh.nodes[:,2],mesh.nodes[:,3],mesh.nodes[:,4], markersize=80, color = :red,marker=:cross)

    for e in 1:mesh.nelem
        p1, p2, p3, p4 = mesh.IEN_geo[1:4,e]
        px = mesh.points[[p1,p2,p3,p4,p1],2]
        py = mesh.points[[p1,p2,p3,p4,p1],3]
        pz = mesh.points[[p1,p2,p3,p4,p1],4]

        pxt = [mesh.points[[p1,p2,p3],2];mesh.points[[p3,p4,p1],2]]
        pyt = [mesh.points[[p1,p2,p3],3];mesh.points[[p3,p4,p1],3]]
        pzt = [mesh.points[[p1,p2,p3],4];mesh.points[[p3,p4,p1],4]]

        xyz = reshape([pxt[:] pyt[:] pzt[:]]', :)
        # mk.poly!(ax,gb.connect(xyz, gb.Point{3}), gb.connect(1:length(pxt), gb.TriangleFace), color = :white, strokecolor = :black, strokewidth = 0)
        mk.lines!(ax,px,py,pz,color=:black)
    end
    return fig
end

