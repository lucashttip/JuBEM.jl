function applyBC!(mesh,solver_var)

    elements = [i for i in 1:mesh.nelem]
    nnel = (mesh.eltype+1)^2

    neu = count(==(1),mesh.bc)
    net = count(==(2),mesh.bc)
    nerb = count(>=(3), mesh.bc)

    c = zeros(3*nnel*nerb,6)
    d = zeros(6, 3*nnel*nerb)
    mb = zeros(typeof(solver_var.H[1]), 3*nnel*mesh.nelem+6, 3*nnel*(neu+net)+6)
    solver_var.ma = zeros(typeof(solver_var.H[1]), 3*nnel*mesh.nelem+6, 3*nnel*mesh.nelem+6)

    elem_u = elements[mesh.bc.==1]
    elem_t = elements[mesh.bc.==2]
    elem_rb = elements[mesh.bc.==3]
    elem_nrb = [elem_u; elem_t]

    y = zeros(typeof(solver_var.H[1]), 3*nnel*(neu+net)+6)

    for i in 1:neu
        i1 = 3*nnel*(i-1)+1
        i2 = 3*nnel*i
        y[i1:i2] = repeat(mesh.bcvalue[3*(elem_u[i]-1)+1:3*(elem_u[i])],nnel)
    end
    j = 3*nnel*neu
    for i in 1:net
        i1 = j + 3*nnel*(i-1)+1
        i2 = j + 3*nnel*i
        y[i1:i2] = repeat(mesh.bcvalue[3*(elem_t[i]-1)+1:3*(elem_t[i])],nnel)
    end

    y[end-5:end] = mesh.bcvalue[end-5:end]


    y_elem = mesh.bcvalue[[elem_u;elem_t;3*mesh.nelem+1:3*mesh.nelem+6]]

    y = convert.(typeof(solver_var.H[1]),y)

    origin = [0.0,0.0,0.0]

    if nerb > 0
        c,d = calc_CeD(mesh, nerb, elem_rb, origin)
    end

    mb.=0.0
    solver_var.ma.=0.0

    i2 = 3*nerb
    i3 = 3*mesh.nelem
    if nerb > 0
    j1 = 1
    j2 = 6
    # ! Montando Hff*C em ma
        ma[1:i2, j1:j2] = solver_var.H[[mesh.LM[:,elem_rb]],[mesh.LM[:,elem_rb]]]*c

        # ! Montando Hsf*C em ma
        ma[i2+1:i3, j1:j2] = solver_var.H[[mesh.LM[:,elem_nrb]],[mesh.LM[:,elem_rb]]]*c
    
    j1 = j2+1
    j2 = j2+3*nerb
    # ! Montando -Gff em ma
    ma[1:i2, j1:j2] = -solver_var.G[[mesh.LM[:,elem_rb]],[mesh.LM[:,elem_rb]]]

    # ! Montando -Gsf em ma
    ma[i2+1:i3,j1:j2] = -solver_var.G[[mesh.LM[:,elem_nrb]],[mesh.LM[:,elem_rb]]]

    # ! Montando D em ma
    ma[i3+1:end,j1:j2] = d
    end

    j1 = j2+1
    j2 = j2+3*net
    # ! Montando Hfs2 em ma
    ma[1:i2,j1:j2] = solver_var.H[[mesh.LM[:,elem_rb]],[mesh.LM[:,elem_t]]]

    # ! Montando Hss2 em ma
    ma[i2+1:i3,j1:j2] = solver_var.H[[mesh.LM[:,elem_nrb]],[mesh.LM[:,elem_t]]]

    j1 = j2+1
    j2 = j2+3*neu
    # ! Montando -Gfs1 em ma
    ma[1:i2,j1:j2] = - solver_var.G[[mesh.LM[:,elem_rb]],[mesh.LM[:,elem_u]]]
    # ! print *, mesh.LM[:,elem_u)

    # ! Montando -Gss1 em ma
    ma[i2+1:i3,j1:j2] = - solver_var.G[[mesh.LM[:,elem_nrb]],[mesh.LM[:,elem_u]]]

    i2 = 3*nerb
    i3 = 3*mesh.nelem
    j1 = 1
    j2 = 3*neu
    # ! Montando -Hfs1 em mb
    mb[1:i2, j1:j2] = -solver_var.H[[mesh.LM[:,elem_rb]],[mesh.LM[:,elem_u]]]

    # ! Montando -Hss1 em mb
    mb[i2+1:i3, j1:j2] = -solver_var.H[[mesh.LM[:,elem_nrb]],[mesh.LM[:,elem_u]]]

    j1 = j2+1
    j2 = j2+3*net
    # ! Montando -Hfs1 em mb
    mb[1:i2, j1:j2] = solver_var.G[[mesh.LM[:,elem_rb]],[mesh.LM[:,elem_t]]]

    # ! Montando -Hss1 em mb
    mb[i2+1:i3, j1:j2] = solver_var.G[[mesh.LM[:,elem_nrb]],[mesh.LM[:,elem_t]]]

    j1 = j2+1
    j2 = j2+6
    eye = convert.(typeof(mb[1]),I(3))

    # ! Montando identidade em mb
    mb[i3+1:end,j1:j2] = eye

    mesh.zbcvalue = mb*y

    return mesh,solver_var

end

function calc_CeD(mesh,nerb,elem_rb,origin)
    nnel = (mesh.eltype+1)^2

    c = zeros(3*nnel*nerb,6)
    d = zeros(6, 3*nnel*nerb)

   
    # ! montagem da matriz rb_nodes com as coordenadas nodais do corpo rigido
    # !
    rb_nodes = mesh.nodes[elem_rb,2:end]

    # ! zerar a matriz c
    # !
    c .= 0.0

    # ! montagem da matriz de compatibilidade cinematica "C"
    # !
    for i in 1:nerb
        for n in 1:nnel
            ii=3*nnel*(i-1) + 3*(n-1)

            dx = (rb_nodes[i,1]-origin[1])
            dy = (rb_nodes[i,2]-origin[2])
            dz = (rb_nodes[i,3]-origin[3])

            c[ii+1:ii+3, 1:6] = [
                1.0 0.0 0.0 0.0 +dz -dy
                0.0 1.0 0.0 -dz 0.0 +dx
                0.0 0.0 1.0 +dy -dx 0.0
            ]
        end
    end
    # !
    # ! zerar a matriz d
    # !
    d .= 0.0
    
    # ! montagem da matriz de equilibrio "D"
    # !
    # ! a area do elemento e calculada pelo produto vetorial de dois vetoras
    # ! que formam dois lados do elemento
    # !
    for i in 1:nerb
        j=elem_rb[i]
        n1=mesh.IEN_geo[1, j]
        n2=mesh.IEN_geo[2, j]
        n3=mesh.IEN_geo[3, j]
        n4=mesh.IEN_geo[4, j]

        l1 = mesh.points[n2,2:end] - mesh.points[n1,2:end]
        l2 = mesh.points[n3,2:end] - mesh.points[n2,2:end]
        l3 = mesh.points[n4,2:end] - mesh.points[n3,2:end]
        l4 = mesh.points[n1,2:end] - mesh.points[n4,2:end]

        cross12 = cross(l1,l2)
        cross34 = cross(l3,l4)

        # cross12(1) = l1(2)*l2(3) - l1(3)*l2(2)
        # cross12(2) = l1(3)*l2(1) - l1(1)*l2(3)
        # cross12(3) = l1(1)*l2(2) - l1(2)*l1(1)

        # cross34(1) = l3(2)*l4(3) - l3(3)*l4(2)
        # cross34(2) = l3(3)*l4(1) - l3(1)*l4(3)
        # cross34(3) = l3(1)*l4(2) - l3(2)*l4(1)

        area1 = norm(cross12)/2.0

        area2 = norm(cross34)/2.0

        area= (area1 + area2)/nnel
        for n in 1:nnel
            ii=3*nnel*(i-1) + 3*(n-1)

            dx = (rb_nodes[i,1]-origin[1])
            dy = (rb_nodes[i,2]-origin[2])
            dz = (rb_nodes[i,3]-origin[3])

            d[1:6,ii+1:ii+3] = [
                1.0 0.0 0.0
                0.0 1.0 0.0
                0.0 0.0 1.0
                0.0 -dz +dy
                +dz 0.0 -dx
                -dy +dx 0.0
            ].*area
        end
    end

    return c, d
end


function applyBC_nonrb!(mesh, solver_var)
    # elements = [i for i in 1:mesh.nelem]
    elements = 1:mesh.nelem
    nnel = (mesh.eltype+1)^2

    neu = count(==(1),mesh.bc)
    net = count(==(2),mesh.bc)

    mb = zeros(typeof(solver_var.H[1]), 3*mesh.nnodes, 3*mesh.nnodes)
    solver_var.ma = zeros(typeof(solver_var.H[1]), 3*mesh.nnodes, 3*mesh.nnodes)

    elem_u = elements[mesh.bc.==1]
    elem_t = elements[mesh.bc.==2]
    
    y = zeros(typeof(solver_var.H[1]), 3*mesh.nnodes)

    for i in 1:neu
        i1 = 3*nnel*(i-1)+1
        i2 = 3*nnel*i
        y[i1:i2] = repeat(mesh.bcvalue[3*(elem_u[i]-1)+1:3*(elem_u[i])],nnel)
    end

    j = 3*nnel*neu

    for i in 1:net
        i1 = j + 3*nnel*(i-1)+1
        i2 = j + 3*nnel*i
        y[i1:i2] = repeat(mesh.bcvalue[3*(elem_t[i]-1)+1:3*(elem_t[i])],nnel)
    end

    solver_var.ma[:,1:3*nnel*neu] = - solver_var.G[:,mesh.LM[:,elem_u][:]]
    solver_var.ma[:,3*nnel*neu+1:end] = solver_var.H[:,mesh.LM[:,elem_t][:]]
    
    mb[:,1:3*nnel*neu] = - solver_var.H[:,mesh.LM[:,elem_u][:]]
    mb[:,3*nnel*neu+1:end] = solver_var.G[:,mesh.LM[:,elem_t][:]]
    # @infiltrate
    mesh.zbcvalue = mb*y

    return mesh,solver_var
end


function applyBC_nonrb2!(mesh, solver_var)

    nnel = size(mesh.IEN,1)
   
    solver_var.ma = zeros(size(solver_var.H))
    mb = zeros(size(solver_var.H))

    y = zeros(3*nnel*length(mesh.bc))

    for e in 1:mesh.nelem
        if mesh.bc[e] == 1
            solver_var.ma[:,mesh.LM[:,e]] = -solver_var.G[:,mesh.LM[:,e]]
            mb[:,mesh.LM[:,e]] = - solver_var.H[:,mesh.LM[:,e]]
        elseif mesh.bc[e] == 2
            solver_var.ma[:,mesh.LM[:,e]] = solver_var.H[:,mesh.LM[:,e]]
            mb[:,mesh.LM[:,e]] = solver_var.G[:,mesh.LM[:,e]]
        end
    end

    for e in 1:mesh.nelem
        i1 = 3*nnel*(e-1)+1
        i2 = 3*nnel*e
        y[i1:i2] = repeat(mesh.bcvalue[3*(e-1)+1:3*(e)],nnel)
    end
    
    mesh.zbcvalue = mb*y

    return mesh,solver_var
end


function returnut(mesh,x)
    nnel = (mesh.eltype+1)^2
    u = zeros(3*mesh.nnodes)
    t = zeros(3*mesh.nnodes)

    elements = [i for i in 1:mesh.nelem]
    neu = count(==(1),mesh.bc)
    net = count(==(2),mesh.bc)
    elem_u = elements[mesh.bc.==1]
    elem_t = elements[mesh.bc.==2]

    t[mesh.LM[:,elem_u]] = x[1:3*nnel*neu]
    u[mesh.LM[:,elem_t]] = x[3*nnel*neu+1:end]

    u = [u[1:3:end] u[2:3:end] u[3:3:end]]
    t = [t[1:3:end] t[2:3:end] t[3:3:end]]

    return u, t
end

function returnut2(mesh,x)
    nnel = (mesh.eltype+1)^2
    u = zeros(3*mesh.nnodes)
    t = zeros(3*mesh.nnodes)


    for e in 1:mesh.nelem
        if mesh.bc[e] == 1
            u[mesh.LM[:,e]] = repeat(mesh.bcvalue[3*(e-1)+1:3*(e)],nnel)
            t[mesh.LM[:,e]] = x[mesh.LM[:,e]]
        elseif mesh.bc[e] == 2
            t[mesh.LM[:,e]] = repeat(mesh.bcvalue[3*(e-1)+1:3*(e)],nnel)
            u[mesh.LM[:,e]] = x[mesh.LM[:,e]]
        end
    end

    u = [u[1:3:end] u[2:3:end] u[3:3:end]]
    t = [t[1:3:end] t[2:3:end] t[3:3:end]]


    return u, t
end

function calc_utpoints(mesh,u,t)
    csi_points = range(-1,1,length = mesh.eltype+1)
    csi_nodes = range(-1+mesh.offset,1-mesh.offset,length = mesh.eltype+1)
    csi_vec = calc_csis_grid(csi_points)
    N = calc_N_matrix(csi_nodes,csi_vec)

    up = zeros(mesh.npoints,3)
    tp = zeros(mesh.npoints,3)

    for e in 1:mesh.nelem
        tmpu = N*u[mesh.IEN[:,e],:]
        tmpt = N*t[mesh.IEN[:,e],:]

        points = mesh.IEN_geo[:,e]

        for p in 1:length(points)
            np = count(==(points[p]),mesh.IEN_geo)
            up[points[p],:] .= up[points[p],:] .+ tmpu[p,:]./np
            tp[points[p],:] .= tp[points[p],:] .+ tmpt[p,:]./np

        end
    end
    return up, tp
end