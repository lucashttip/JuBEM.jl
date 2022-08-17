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