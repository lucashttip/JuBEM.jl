# AS FUNÇÕES NESTE ARQUIVO ESTÃO ANTIGAS E NÃO FUNCIONAM MAIS NA API DO PROGRAMA.
# DEVEM SER MODIFICADAS EM ATUALIZAÇÕES FUTURAS


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
        N = calc_N_gen(csis_cont,gp[i].csi)
        Nd = calc_N_gen(csis_descont,gp[i].csi)
        Ncsi = calc_N_gen(csis_cont,gp[i].csi;dg=:dNdc)
        Neta = calc_N_gen(csis_cont,gp[i].csi;dg=:dNde)
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
        N = calc_N_gen(csis_cont,gp[i].csi)
        Ncsi = calc_N_gen(csis_cont,gp[i].csi;dg=:dNdc)
        Neta = calc_N_gen(csis_cont,gp[i].csi;dg=:dNde)
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


function calc_interior_dynamic(mesh,material,u,t,points_int,frequency)

    npoints = size(points_int,1)
    u_int = zeros(typeof(u[1]),npoints,3)
    s_int = zeros(typeof(t[1]),npoints,3)


    # DEFINING PARAMETERS
    # mesh/material/other constants definition
    nelem = mesh.nelem
    delta = I(3)
    nnel = (mesh.eltype+1)^2
    C_stat = calc_static_constants(material[1])
    zconsts = cmplx_consts(material[1],frequency)


    # normal integration constants
    csis_cont, csis_descont, rules  = calc_nonsing_consts(mesh)
    N, dNc, dNe, Nd = calc_N_nonsing(csis_cont, csis_descont, rules)
    
    # Field loop:
    # Threads.@threads for fe in 1:nelem
    for fe in 1:nelem

        nodes_idx = mesh.IEN[:,fe]
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]
        normal, J, gauss_points = calc_nJgp(N, dNc, dNe, rules.gp, field_points)


        # source loop:
        for sp in 1:npoints
                source_point = points_int[sp,:]

                r, c, d = calc_dist(source_point, field_points, rules.dists,csis_cont)
                if r == 0
                    r = 4
                end
                zHELEM, zGELEM = integrate_nonsing_dynamic(source_point,gauss_points[r],Nd[r],normal[r],J[r], rules.gp[r].omega, delta, zconsts)
                u_int[sp,:] += zGELEM*vec(t[nodes_idx,:]') - zHELEM*vec(u[nodes_idx,:]')
        end
    
    end


    return u_int, s_int
    
end


function calc_interior_static2(mesh,material,u,t,points_int)

    npoints = size(points_int,1)
    u_int = zeros(typeof(u[1]),npoints,3)
    s_int = zeros(typeof(t[1]),npoints,3)


    # DEFINING PARAMETERS
    # mesh/material/other constants definition
    nelem = mesh.nelem
    delta = I(3)
    nnel = (mesh.eltype+1)^2
    C_stat = calc_static_constants(material[1])

    # normal integration constants
    csis_cont, csis_descont, rules  = calc_nonsing_consts(mesh)
    N, dNc, dNe, Nd = calc_N_nonsing(csis_cont, csis_descont, rules)
    
    # Field loop:
    # Threads.@threads for fe in 1:nelem
    for fe in 1:nelem

        nodes_idx = mesh.IEN[:,fe]
        field_points = mesh.points[mesh.IEN_geo[:,fe],2:end]
        normal, J, gauss_points = calc_nJgp(N, dNc, dNe, rules.gp, field_points)


        # source loop:
        for sp in 1:npoints
                source_point = points_int[sp,:]

                r, c, d = calc_dist(source_point, field_points, rules.dists,csis_cont)
                if r == 0
                    r = 4
                end
                HELEM, GELEM = integrate_nonsing_static(source_point,gauss_points[r],Nd[r],normal[r],J[r], rules.gp[r].omega, delta, C_stat)

                u_int[sp,:] += GELEM*vec(t[nodes_idx,:]') - HELEM*vec(u[nodes_idx,:]')
        end
    
    end


    return u_int, s_int
    
end