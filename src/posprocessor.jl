function calc_interior_static(mesh,material,u,t,points_int)

    nelem = mesh.nelem
    m = 1
    delta = I(3) 
    nnel = (mesh.eltype+1)^2
    npoints = size(points_int,1)
    u_int = zeros(typeof(u[1]),npoints,3)
    s_int = zeros(typeof(t[1]),npoints,3)

    eltype_cont = mesh.eltype
    offset = mesh.offset

    if eltype_cont == 0
        eltype_cont = 1
        offset = 1.0
    end

    csis_cont = range(-1,1,length = eltype_cont+1)
    csis_descont = range(-1+offset,1 - offset,length=mesh.eltype+1)
    
    # Cálculos para elementos não singulares
    gp, dists = calc_points_weights()
    
    Nc_nonsing = []
    dNcdcsi_nonsing = []
    dNcdeta_nonsing = []
    Nd_nonsing = []
    
    for i in eachindex(gp)
        N = calc_N_matrix(csis_cont,gp[i].csi)
        Nd = calc_N_matrix(csis_descont,gp[i].csi)
        Ncsi = calc_dNdcsi_matrix(csis_cont,gp[i].csi)
        Neta = calc_dNdeta_matrix(csis_cont,gp[i].csi)
        push!(Nc_nonsing,N)
        push!(dNcdcsi_nonsing,Ncsi)
        push!(dNcdeta_nonsing,Neta)
        push!(Nd_nonsing,Nd)
    end

    
    C_stat = calc_static_constants(material[m])

    # Field loop:
    # Threads.@threads for fe in 1:nelem
    for fe in 1:nelem
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]
        nodes_idx = mesh.IEN[:,fe]
        normal = []
        J = []
        gauss_points = []
        for i in eachindex(gp)
            normal2, J2 = calc_n_J_matrix(dNcdcsi_nonsing[i], dNcdeta_nonsing[i], field_points)
            gauss_points2 = Nc_nonsing[i]*field_points
            push!(normal,normal2)
            push!(J,J2)
            push!(gauss_points,gauss_points2)
        end

        # source loop:
        for sp in 1:npoints
                source_point = points_int[sp,:]

                r,c,d = calc_dist(source_point, field_points, dists,csis_cont)
                HELEM, GELEM = integrate_nonsing_static(source_point,gauss_points[r],Nd_nonsing[r],normal[r],J[r], gp[r].omega, delta, C_stat)
                u_int[sp,:] += GELEM*vec(t[nodes_idx,:]') - HELEM*vec(u[nodes_idx,:]')
        end
    
    end


    return u_int, s_int
    
end


function calc_interior_static_const(mesh,material,u,t,points_int)

    nelem = mesh.nelem
    m = 1
    delta = I(3) 
    nnel = (mesh.eltype+1)^2
    npoints = size(points_int,1)
    u_int = zeros(typeof(u[1]),npoints,3)
    s_int = zeros(typeof(t[1]),npoints,3)

    csis_cont = range(-1,1,length = 2)
    
    # Cálculos para elementos não singulares
    gp, dists = calc_points_weights()
    
    Nc_nonsing = []
    dNcdcsi_nonsing = []
    dNcdeta_nonsing = []
    
    for i in eachindex(gp)
        N = calc_N_matrix(csis_cont,gp[i].csi)
        Ncsi = calc_dNdcsi_matrix(csis_cont,gp[i].csi)
        Neta = calc_dNdeta_matrix(csis_cont,gp[i].csi)
        push!(Nc_nonsing,N)
        push!(dNcdcsi_nonsing,Ncsi)
        push!(dNcdeta_nonsing,Neta)
    end

    
    C_stat = calc_static_constants(material[m])

    # Field loop:
    # Threads.@threads for fe in 1:nelem
    for fe in 1:nelem
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]
        nodes_idx = mesh.IEN[:,fe]
        normal = []
        J = []
        gauss_points = []
        for i in eachindex(gp)
            normal2, J2 = calc_n_J_matrix(dNcdcsi_nonsing[i], dNcdeta_nonsing[i], field_points)
            gauss_points2 = Nc_nonsing[i]*field_points
            push!(normal,normal2)
            push!(J,J2)
            push!(gauss_points,gauss_points2)
        end

        # source loop:
        for sp in 1:npoints
                source_point = points_int[sp,:]

                r,c,d = calc_dist(source_point, field_points, dists,csis_cont)
                HELEM, GELEM = integrate_const_static(source_point, gauss_points[r],normal[r],J[r], gp[r].omega, delta, C_stat)
                u_int[sp,:] += GELEM*vec(t[nodes_idx,:]') - HELEM*vec(u[nodes_idx,:]')
        end
    
    end


    return u_int, s_int
    
end