import PlotlyJS

function triangulate(mesh::Mesh)
    triangles = zeros(Int64,2*mesh.nelem,3)

    for e in axes(mesh.IEN,2)
        triidx = [2*(e-1)+1,2*e]

        ijk = [
            1 2 4
            2 3 4 
            ]
        triangles[triidx,:] = mesh.IEN[ijk,e].-1
    end
    return triangles
end

function map_plocs(csis)
    nnel = length(csis)^2
    k = JuBEM.map_k(nnel)
    plocs = zeros(nnel,2)
    for i in eachindex(k)
        plocs[i,:] = [csis[k[i][1]], csis[k[i][2]]]
    end
    return plocs
end

function plot_disp(mesh::Mesh,sol::Solution,dim)
    
    u = sol.u

    order = mesh.eltype
    offset = mesh.offset
    csis_points = range(-1,1,length=order+1) 
    csis_nodes = range(-1+offset,1-offset,length=order+1) 
    plocs = map_plocs(csis_points)
    N = calc_N_gen(csis_nodes,plocs)
    upoints = zeros(size(u))
    xpoints = zeros(size(u))

    tri = triangulate(mesh)

    for e in axes(mesh.IEN,2)
        pidx = mesh.IEN_geo[:,e]
        nidx = mesh.IEN[:,e]
        nodes = mesh.nodes[nidx,2:end]




        unodes = sol.u[nidx,:]
        upoints[nidx,:] = N*unodes
        xpoints[nidx,:] = N*nodes

    end

    t = PlotlyJS.mesh3d(x=xpoints[:,1],y=xpoints[:,2],z=xpoints[:,3],i=tri[:,1],j=tri[:,2],k=tri[:,3], intensity = upoints[:,dim],showlegend=true)

    max_x = maximum(abs.(xpoints[:,1]))
    max_y = maximum(abs.(xpoints[:,2]))
    max_z = maximum(abs.(xpoints[:,3]))

    max_dim = maximum([max_x,max_y,max_z])

    asp_x = 1.5*max_x/max_dim
    asp_y = 1.5*max_y/max_dim
    asp_z = 1.5*max_z/max_dim

    layout = PlotlyJS.Layout( 
        coloraxis=PlotlyJS.attr(autocolorscale=true)    ,
        scene_aspectratio=PlotlyJS.attr(x=asp_x, y=asp_y, z=asp_z),
        margin=PlotlyJS.attr(autoexpand=true))

    return PlotlyJS.plot(t, layout)
end


function triangulate_geo(mesh::Mesh)
    triangles = zeros(Int64,2*mesh.nelem,3)

    for e in axes(mesh.IEN_geo,2)
        triidx = [2*(e-1)+1,2*e]

        ijk = [
            1 2 4
            2 3 4 
            ]
        triangles[triidx,:] = mesh.IEN_geo[ijk,e].-1
    end
    return triangles
end

function plot_meshtags(mesh::Mesh)

    tri = triangulate_geo(mesh)

    xpoints = mesh.points[:,2:end]
    cpoints = zeros(size(mesh.points,1))

    for e in axes(mesh.IEN_geo,2)

        pidx = mesh.IEN_geo[:,e]

        cpoints[pidx] .= mesh.tag[e]

    end

    colors = [
    	"rgb(255, 0, 0)",
    	"rgb(0, 255, 0)",
    	"rgb(0, 0, 255)",
    	"rgb(200, 200, 50)",
    	"rgb(230, 200, 10)",
    	"rgb(255, 140, 0)"
    ]
    facecolor = repeat(colors[mesh.tag], inner=[2])

    t = PlotlyJS.mesh3d(x=xpoints[:,1],y=xpoints[:,2],z=xpoints[:,3],i=tri[:,1],j=tri[:,2],k=tri[:,3], facecolor = facecolor,showlegend=true,colorscale=[
        [0, "rgb(255, 0, 255)"],
        [0.5, "rgb(0, 255, 0)"],
        [1, "rgb(0, 0, 255)"]
    ])

    # t = PlotlyJS.mesh3d(x=xpoints[:,1],y=xpoints[:,2],z=xpoints[:,3],i=tri[:,1],j=tri[:,2],k=tri[:,3], intensity = cpoints,showlegend=true,colorscale=[
    #     [0, "rgb(255, 0, 255)"],
    #     [0.5, "rgb(0, 255, 0)"],
    #     [1, "rgb(0, 0, 255)"]
    # ])

    # @infiltrate
    
    # layout = PlotlyJS.Layout( 
    #     coloraxis=PlotlyJS.attr(autocolorscale=true)    ,
    #     scene_aspectratio=PlotlyJS.attr(x=asp_x, y=asp_y, z=asp_z),
    #     margin=PlotlyJS.attr(autoexpand=true))

    # return PlotlyJS.plot(t, layout)
    return PlotlyJS.plot(t)

end