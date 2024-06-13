function writevtk(mesh,u,t,filename)

    if mesh.eltype > 1
        error("Element type not yet supported")
    end
    if mesh.offset > 0
        u,t = calc_utpoints(mesh,u,t)
    end

    if typeof(u[1]) == ComplexF64
        u = real.(u)
        t = real.(t)
    end

    points = mesh.points[:,2:end]'

    cells = [MeshCell(VTKCellTypes.VTK_QUAD,mesh.IEN_geo[:,e]) for e in 1:mesh.nelem]


    vtk_grid(filename, points, cells) do vtk
        vtk["disp_x"] = u[:,1]
        vtk["disp_y"] = u[:,2]
        vtk["disp_z"] = u[:,3]
        vtk["disp"] = u'
        vtk["tractions_x"] = t[:,1]
        vtk["tractions_y"] = t[:,2]
        vtk["tractions_z"] = t[:,3]
    end

end

function writevtk2(mesh, u,t, filename)

    u2,t2,points,elem = extrapolateres(mesh, u, t)


    if typeof(u2[1]) == ComplexF64
        u2 = real.(u2)
        t2 = real.(t2)
    end

    points = points[:,2:end]'

    cells = [MeshCell(VTKCellTypes.VTK_QUAD,elem[:,e]) for e in axes(elem,2)]


    vtk_grid(filename, points, cells) do vtk
        vtk["disp_x"] = u2[:,1]
        vtk["disp_y"] = u2[:,2]
        vtk["disp_z"] = u2[:,3]
        vtk["disp"] = u2'
        vtk["tractions_x"] = t2[:,1]
        vtk["tractions_y"] = t2[:,2]
        vtk["tractions_z"] = t2[:,3]
    end


end

function extrapolateres(mesh, u, t)

    u2 = zeros(typeof(u[1]),size(u))
    t2 = zeros(typeof(t[1]),size(t))
    points = zeros(typeof(mesh.nodes[1]),size(mesh.nodes))
    elem = mesh.IEN

    eltype = mesh.eltype
    offset = mesh.offset
    eltype_cont = eltype

    if eltype == 0
        offset = 1.0
        eltype_cont = 1
        mesh2 = copy(mesh)
        mesh2.eltype = 1
        mesh2.offset = 0.5
        generate_mesh!(mesh2)
        u2 = zeros(typeof(u[1]),mesh2.nnodes,3)
        t2 = zeros(typeof(t[1]),mesh2.nnodes,3)
        points = zeros(typeof(mesh2.nodes[1]),size(mesh2.nodes))
        remove_EE_mesh!(mesh2)
        elem = mesh2.IEN
    end

    csis_cont = range(-1,1,length=eltype_cont+1) 
    csis_descont = range(-1+offset,1-offset,length = eltype+1)

    grid = calc_csis_grid(csis_cont)

    N = calc_N_gen(csis_descont,grid)

    idx = [1,4,2,3]


    for e in axes(elem,2)
        nodes_in = mesh.IEN[:,e]
        nodes_out = elem[:,e]
        idxp = mesh.IEN_geo[:,e]
        u2[nodes_out[idx],:] = N*u[nodes_in,:]
        t2[nodes_out[idx],:] = N*t[nodes_in,:]

        points[nodes_out,:] = mesh.points[idxp,:]

    end

    return u2, t2, points, elem

end