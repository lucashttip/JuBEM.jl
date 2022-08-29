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

function applyBC_rb(mesh, solver_var,H,G)

    nnel = (mesh.eltype+1)^2

    bcvalue = reshape(mesh.bcvalue',length(mesh.bcvalue))
    LM = zeros(Int64,nnel,size(mesh.LM,2),3)
    for i in 1:3
        LM[:,:,i] = mesh.LM[i:3:end,:]
    end

    nu = count(==(1),mesh.bc)
    nt = count(==(2),mesh.bc)
    ifu = findall(x->x==1,mesh.bc)
    ift = findall(x->x==2,mesh.bc)
    elem_F = findall(mesh.bc[:,1].<3)


    nrrb = count(>=(3),unique(mesh.bc))
    nrb = zeros(nrrb)
    elem_r = findall(mesh.bc[:,1].>=3)
    elem_ri = Array{Any,1}(undef,nrrb)
    C = Array{Any,1}(undef,nrrb)
    D = Array{Any,1}(undef,nrrb)
    for i in 1:nrrb
        nrb[i] = count(==(2+i), mesh.bc[:,1])
        rbidx = 2+i
        elem_ri[i] = sort(findall(mesh.bc[:,1].==rbidx))
        iaux = findfirst(mesh.bc[:,1].==rbidx)
        rbcenter = mesh.bcvalue[iaux,:]
        c,d = calc_CD_discont(mesh,solver_var,rbcenter,rbidx)
        C[i] = c
        D[i] = d
    end
    nerb = Int(sum(nrb))
    ngF = Int((nu+nt))


    Hrirj = Array{Any,2}(undef,nrrb,nrrb)
    HFri = Array{Any,1}(undef,nrrb)
    Grirj = Array{Any,2}(undef,nrrb,nrrb)
    GFri= Array{Any,1}(undef,nrrb)

    for i in 1:nrrb
        for j in 1:nrrb
            Hrirj[i,j] = H[mesh.LM[:,elem_ri[i]][:], mesh.LM[:,elem_ri[j]][:]]*C[j]
            Grirj[i,j] = G[mesh.LM[:,elem_ri[i]][:], mesh.LM[:,elem_ri[j]][:]]
        end
        HFri[i] = H[mesh.LM[:,elem_F][:], mesh.LM[:,elem_ri[i]][:]]*C[i]
        GFri[i] = G[mesh.LM[:,elem_F][:], mesh.LM[:,elem_ri[i]][:]]
    end

    HrFt = H[mesh.LM[:,elem_r][:], LM[:,ift][:]]
    HFFt = H[mesh.LM[:,elem_F][:], LM[:,ift][:]]
    HrFu = H[mesh.LM[:,elem_r][:], LM[:,ifu][:]]
    HFFu = H[mesh.LM[:,elem_F][:], LM[:,ifu][:]]
    GrFt = G[mesh.LM[:,elem_r][:], LM[:,ift][:]]
    GFFt = G[mesh.LM[:,elem_F][:], LM[:,ift][:]]
    GrFu = G[mesh.LM[:,elem_r][:], LM[:,ifu][:]]
    GFFu = G[mesh.LM[:,elem_F][:], LM[:,ifu][:]]

    Hrr = zeros(typeof(H[1]), 3*nnel*nerb, 6*nrrb)
    HFr = zeros(typeof(H[1]), nnel*ngF, 6*nrrb)
    Grr = zeros(typeof(H[1]), 3*nnel*nerb, 3*nnel*nerb)
    GFr = zeros(typeof(H[1]), nnel*ngF, 3*nnel*nerb)

    i1 = 1
    i2 = 0
    for i in 1:nrrb
        neri = Int(nrb[i])
        i2 = i2+3*nnel*neri
        j1 = 1
        j2 = 0
        for j in 1:nrrb
            nerj = Int(nrb[j])
            j2 = j2+3*nnel*nerj
            Hrr[i1:i2,6*(j-1)+1:6*j] = Hrirj[i,j]
            Grr[i1:i2,j1:j2] = Grirj[i,j]
            j1 = j2+1
        end
        HFr[:,6*(i-1)+1:6*i] = HFri[i]
        GFr[:,i1:i2] = GFri[i]
        i1 = i2+1
    end

    Dg = zeros(6*nrrb, 3*nnel*nerb)
    If = I(6*nrrb)

    j1 = 1
    j2 = 0
    for i in 1:nrrb
        nerj = Int(nrb[i])
        j2 = j2+3*nnel*nerj
        Dg[6*(i-1)+1:6*i,j1:j2] = D[i]
        j1 = j2+1
    end

    # mb = zeros(typeof(H[1]), 3*mesh.nnodes+6*nrrb, 3*mesh.nnodes+6*nrrb)
    # solver_var.ma = zeros(typeof(H[1]), 3*mesh.nnodes+6*nrrb, 3*mesh.nnodes+6*nrrb)

    y = zeros(typeof(H[1]), nnel*(nu+nt)+6*nrrb)

    for i in 1:nu
        i1 = nnel*(i-1)+1
        i2 = nnel*i
        y[i1:i2] .= mesh.bcvalue[ifu[i]]
    end

    j = nnel*nu

    for i in 1:nt
        i1 = j + nnel*(i-1)+1
        i2 = j + nnel*i
        y[i1:i2] .= mesh.bcvalue[ift[i]]
    end

    j = nnel*(nu+nt)

    for i in 1:nrrb
        i1 = j+6*(i-1)+1
        i2 = j+6*i
        y[i1:i2] = mesh.forces[:,i]
    end
    
    l1 = 6*nrrb
    O1 = zeros(l1,size(Hrr,2))
    O2 = zeros(l1,size(HFFt,2))
    O3 = zeros(l1,size(GFFu,2))
    O4 = zeros(size(HrFu,1),l1)
    O5 = zeros(size(HFFu,1),l1)
    O6 = zeros(l1,size(HFFu,2))
    O7 = zeros(l1,size(GFFt,2))

    solver_var.ma = [
        Hrr HrFt -Grr -GrFu
        HFr HFFt -GFr -GFFu
        O1 O2 Dg O3
    ]
    mb = [
        -HrFu GrFt O4
        -HFFu GFFt O5
        O6 O7 collect(If)
    ]

    # solver_var.ma[:,1:nnel*nu] = - G[:,LM[:,iu][:]]
    # solver_var.ma[:,nnel*nu+1:end] = H[:,LM[:,it][:]]
    
    # mb[:,1:nnel*nu] = - H[:,LM[:,iu][:]]
    # mb[:,nnel*nu+1:end] = G[:,LM[:,it][:]]
    mesh.zbcvalue = mb*y

    return mesh,solver_var, C
end

function returnut_rb(mesh,x, C=[])
    nnel = Int((mesh.eltype+1)^2.0)
    u = zeros(typeof(x[1]),3*mesh.nnodes)
    t = zeros(typeof(x[1]),3*mesh.nnodes)

    LM = zeros(Int64,nnel,size(mesh.LM,2),3)
    for i in 1:3
        LM[:,:,i] = mesh.LM[i:3:end,:]
    end

    nu = count(==(1),mesh.bc)
    nt = count(==(2),mesh.bc)
    nrbg = count(>=(3),mesh.bc)

    iu = findall(x->x==1,mesh.bc)
    it = findall(x->x==2,mesh.bc)

    nrrb = count(>=(3),unique(mesh.bc))
    nrb = zeros(Int64,nrrb)
    elem_ri = Array{Any,1}(undef,nrrb)
    urb = zeros(typeof(x[1]),6,nrrb)
    for i in 1:nrrb
        rbidx = 2+i
        nrb[i] = count(==(rbidx), mesh.bc[:,1])
        elem_ri[i] = sort(findall(mesh.bc[:,1].==rbidx))
    end

    for i in 1:nrrb
        urb[:,i] = x[6*(i-1)+i:6*i]
        u[mesh.LM[:,elem_ri[i]][:],:] = C[i]*urb[:,i]
    end

    j = 6*nrrb

    u[LM[:,it][:]] = x[j+1:j+nnel*nt]
    j = j+nnel*nt
    @infiltrate
    for i in 1:nrrb
        t[mesh.LM[:,elem_ri[i]]] = x[j+1:j+3*nnel*nrb[i]]
        j = j + 3*nnel*nrb[i]
    end

    t[LM[:,iu][:]] = x[j+1+1:end]

    for i in 1:nu
        u[LM[:,iu[i]][:]] .= mesh.bcvalue[iu[i]]
    end

    for i in 1:nt
        t[LM[:,it[i]][:]] .= mesh.bcvalue[it[i]]
    end


    u = [u[1:3:end] u[2:3:end] u[3:3:end]]
    t = [t[1:3:end] t[2:3:end] t[3:3:end]]

    return u, t, urb
end

function applyBC_nonrb!(mesh, solver_var,H,G)

    nnel = (mesh.eltype+1)^2

    bcvalue = reshape(mesh.bcvalue',length(mesh.bcvalue))
    LM = zeros(Int64,nnel,size(mesh.LM,2),3)
    for i in 1:3
        LM[:,:,i] = mesh.LM[i:3:end,:]
    end

    nu = count(==(1),mesh.bc)
    nt = count(==(2),mesh.bc)

    mb = zeros(typeof(H[1]), 3*mesh.nnodes, 3*mesh.nnodes)
    solver_var.ma = zeros(typeof(H[1]), 3*mesh.nnodes, 3*mesh.nnodes)

    iu = findall(x->x==1,mesh.bc)
    it = findall(x->x==2,mesh.bc)
    
    y = zeros(typeof(H[1]), 3*mesh.nnodes)

    for i in 1:nu
        i1 = nnel*(i-1)+1
        i2 = nnel*i
        y[i1:i2] .= mesh.bcvalue[iu[i]]
    end

    j = nnel*nu

    for i in 1:nt
        i1 = j + nnel*(i-1)+1
        i2 = j + nnel*i
        y[i1:i2] .= mesh.bcvalue[it[i]]
    end

    solver_var.ma[:,1:nnel*nu] = - G[:,LM[:,iu][:]]
    solver_var.ma[:,nnel*nu+1:end] = H[:,LM[:,it][:]]
    
    mb[:,1:nnel*nu] = - H[:,LM[:,iu][:]]
    mb[:,nnel*nu+1:end] = G[:,LM[:,it][:]]
    mesh.zbcvalue = mb*y

    return mesh,solver_var
end

function returnut(mesh,x)
    nnel = (mesh.eltype+1)^2
    u = zeros(typeof(x[1]),3*mesh.nnodes)
    t = zeros(typeof(x[1]),3*mesh.nnodes)

    bcvalue = reshape(mesh.bcvalue',length(mesh.bcvalue))
    LM = zeros(Int64,nnel,size(mesh.LM,2),3)
    for i in 1:3
        LM[:,:,i] = mesh.LM[i:3:end,:]
    end

    nu = count(==(1),mesh.bc)
    nt = count(==(2),mesh.bc)

    iu = findall(x->x==1,mesh.bc)
    it = findall(x->x==2,mesh.bc)


    u[LM[:,it][:]] = x[1:nnel*nt]
    t[LM[:,iu][:]] = x[nnel*nt+1:end]

    for i in 1:nu
        u[LM[:,iu[i]][:]] .= mesh.bcvalue[iu[i]]
    end

    for i in 1:nt
        t[LM[:,it[i]][:]] .= mesh.bcvalue[it[i]]
    end
    
    # t[mesh.LM[:,elem_u]] = x[1:3*nnel*neu]
    # u[mesh.LM[:,elem_t]] = x[3*nnel*neu+1:end]

    # for e in elem_t
    #     t[mesh.LM[:,e]] = repeat(mesh.bcvalue[3*(e-1)+1:3*(e)],nnel)
    # end

    # for e in elem_u
    #     u[mesh.LM[:,e]] = repeat(mesh.bcvalue[3*(e-1)+1:3*(e)],nnel)
    # end

    u = [u[1:3:end] u[2:3:end] u[3:3:end]]
    t = [t[1:3:end] t[2:3:end] t[3:3:end]]

    return u, t
end

function applyBC_nonrb3!(mesh, solver_var,H,G)

    nnel = size(mesh.IEN,1)
   
    solver_var.ma = zeros(typeof(H[1]),size(H))
    mb = zeros(typeof(H[1]),size(H))

    y = zeros(typeof(H[1]),3*mesh.nnodes)

    for e in 1:mesh.nelem
        for n in 1:nnel
            for i in 1:3
                if mesh.bc[e,i] == 1
                    solver_var.ma[:,mesh.ID[i,mesh.IEN[n,e]]] = - G[:,mesh.ID[i,mesh.IEN[n,e]]]
                    mb[:,mesh.ID[i,mesh.IEN[n,e]]] = - H[:,mesh.ID[i,mesh.IEN[n,e]]]
                elseif mesh.bc[e,i] == 2
                    solver_var.ma[:,mesh.ID[i,mesh.IEN[n,e]]] = H[:,mesh.ID[i,mesh.IEN[n,e]]]
                    mb[:,mesh.ID[i,mesh.IEN[n,e]]] = G[:,mesh.ID[i,mesh.IEN[n,e]]]
                end
            end
        end
    end

    for e in 1:mesh.nelem
        for n in 1:nnel
            for i in 1:3
            y[mesh.ID[i,mesh.IEN[n,e]]] = mesh.bcvalue[e,i]
            end
        end
    end
    
    mesh.zbcvalue = mb*y

    return mesh,solver_var
end

function returnut3(mesh,x)
    nnel = (mesh.eltype+1)^2
    u = zeros(typeof(x[1]),mesh.nnodes,3)
    t = zeros(typeof(x[1]),mesh.nnodes,3)

        

    for e in 1:mesh.nelem
        for i in 1:3
            for n in 1:nnel
                if mesh.bc[e,i] == 1
                    u[mesh.IEN[n,e],i] = mesh.bcvalue[e,i]
                    t[mesh.IEN[n,e],i] = x[mesh.ID[i,mesh.IEN[n,e]]]
                elseif mesh.bc[e,i] == 2
                    t[mesh.IEN[n,e],i] = mesh.bcvalue[e,i]
                    u[mesh.IEN[n,e],i] = x[mesh.ID[i,mesh.IEN[n,e]]]
                end
            end
        end
    end

    return u, t
end

function calc_utpoints(mesh,u,t)
    if mesh.eltype>0
        csi_points = range(-1,1,length = mesh.eltype+1)
        csi_nodes = range(-1+mesh.offset,1-mesh.offset,length = mesh.eltype+1)
        csi_vec = calc_csis_grid(csi_points)
        N = calc_N_matrix(csi_nodes,csi_vec)
    else
        csi_points = [-1.0,1.0]
        csi_nodes = 0.0
        csi_vec = calc_csis_grid(csi_points)
        N = calc_N_matrix(csi_nodes,csi_vec)
    end

    up = zeros(typeof(u[1]),mesh.npoints,3)
    tp = zeros(typeof(t[1]),mesh.npoints,3)

    for e in 1:mesh.nelem
        tmpu = N*u[mesh.IEN[:,e],:]
        tmpt = N*t[mesh.IEN[:,e],:]

        points = mesh.IEN_geo[:,e]

        for p in eachindex(points)
            np = count(==(points[p]),mesh.IEN_geo)
            up[points[p],:] .= up[points[p],:] .+ tmpu[p,:]./np
            tp[points[p],:] .= tp[points[p],:] .+ tmpt[p,:]./np

        end
    end
    return up, tp
end