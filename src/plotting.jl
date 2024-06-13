# using CairoMakie
import GLMakie
import GeometryBasics
const mk = GLMakie
const gb = GeometryBasics
# const mk = CairoMakie
# using Meshes
# using MeshViz

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

function view_res(mesh,u)
    fig,ax,plt = mk.scatter(mesh.points[:,2],mesh.points[:,3],mesh.points[:,4], markersize=80, color = :blue)
    # mk.scatter!(ax,mesh.points[:,2],mesh.points[:,3],mesh.points[:,4], markersize=70, color = :blue)

    mk.scatter!(ax,mesh.nodes[:,2],mesh.nodes[:,3],mesh.nodes[:,4], markersize=80, color = :red,marker=:cross)

    csi_cont = range(-1,1,length = mesh.eltype+1)
    csi_descont = range(-1+mesh.offset,1-mesh.offset,length = mesh.eltype+1)

    csis = calc_csis_grid(csi_cont)
    N = calc_N_gen(csi_descont,csis)
    

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

function animate_res_freq(mesh,u,freq;dir=0,frac = 0.5, filename = "anim.mp4",res = (800, 600),t = 0)

    u1,_ = calc_utpoints(mesh,u,u).*frac

    points = mesh.points[:,2:end] .+ real.(u1)
    elem = mesh.IEN_geo'
    conn = connect.([Tuple(elem[i,:]) for i in axes(elem,1)])
    faces = [gb.QuadFace(elem[n,:]) for n in axes(elem,1)]
    ps = [gb.Point3f0(points[n,:]) for n in axes(points,1)]

    m = gb.normal_mesh(ps, faces)

    if dir == 0
        c = [norm(real.(u1[n,:])) for n in axes(u1,1)]
    else
        c = real.(u1[:,dir])
    end

    fig = mk.Figure(resolution = res)
    ax = mk.Axis3(fig[1, 1],aspect=:data)
    # fig,ax,plt = mk.mesh(m,color = c)
    plt = mk.mesh!(ax,m,color = c)
    plt2 = mk.wireframe!(ax, m, color=(:black, 0.5), linewidth=2, transparency=true)

    # animation settings
    ti = 0
    if t <= 0
        tf = 2*pi/freq
    else
        tf = t
    end
    framerate = 30
    nframes = Int(round(30*(tf-ti)))
    time_iterator = range(ti, tf, length=nframes)

    mk.record(fig, filename, time_iterator; framerate = framerate) do time
        change = exp(freq*time*im)
        u2 = u1.*change
        points = mesh.points[:,2:end] .+ real.(u2)
        ps = [gb.Point3f0(points[n,:]) for n in axes(points,1)]
        m = gb.normal_mesh(ps, faces)

        if dir == 0
            c = [norm(real(u2[n,:])) for n in axes(u2,1)]
        else
            c = real.(u2[:,dir])
        end

        plt.input_args[1][] = m
        plt.attributes.color[] = c
        plt2.input_args[1][] = m
    end


end