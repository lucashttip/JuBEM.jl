function applyBC(mesh, solver_var,H,G)

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

function applyBC_simple(mesh::Mesh,problem::Problem,assembly::Assembly,H,G)

    nnel = size(mesh.IEN,1)

    nDofs = size(assembly.H,1)

    LHS = zeros(eltype(assembly.H),nDofs,nDofs)
    RHS = zeros(eltype(assembly.H),nDofs)

    RHS_mat = zeros(eltype(assembly.H),nDofs,nDofs)
    RHS_vec = zeros(eltype(assembly.H),nDofs)
    

    for e in 1:mesh.nelem

        Dofs = mesh.LM[:,e]

        tag = mesh.tag[e]
        bctag = problem.taginfo[tag,2]
        bctypes = problem.bctype[bctag,2:end]
        bcvalues = problem.bcvalue[bctag,:]

        LHS_aux = zeros(nDofs,3*nnel)
        RHS_aux = zeros(nDofs,3*nnel)

        for i in 1:3
            if bctypes[i] == 1
                LHS_aux[:,i:3:3*nnel] = -assembly.G[:,Dofs[i:3:3*nnel]]
                RHS_aux[:,i:3:3*nnel] = -assembly.H[:,Dofs[i:3:3*nnel]]
                RHS_vec[Dofs[i:3:3*nnel]] .= bcvalues[i]

            elseif bctypes[i] == 2
                LHS_aux[:,i:3:3*nnel] = assembly.H[:,Dofs[i:3:3*nnel]]
                RHS_aux[:,i:3:3*nnel] = assembly.G[:,Dofs[i:3:3*nnel]]
                RHS_vec[Dofs[i:3:3*nnel]] .= bcvalues[i]
            else
                error("not supported by this func")
            end
        end

        LHS[:,Dofs] = LHS_aux
        RHS_mat[:,Dofs] = RHS_aux
    end

    RHS = RHS_mat*RHS_vec

    return LHS, RHS
end

function returnut_simple(x,mesh::Mesh,problem::Problem)
    
    nnel = size(mesh.IEN,1)
    
    u = zeros(mesh.nnodes,3)
    t = zeros(mesh.nnodes,3)

    
    for e in 1:mesh.nelem
        Dofs = mesh.LM[:,e]

        tag = mesh.tag[e]
        bctag = problem.taginfo[tag,2]
        bctypes = problem.bctype[bctag,2:end]
        bcvalues = problem.bcvalue[bctag,:]

        u_aux = zeros(nnel,3)
        t_aux = zeros(nnel,3)

        
        for i in 1:3
            if bctypes[i] == 1
                u_aux[:,i] .= bcvalues[i]
                t_aux[:,i] = x[Dofs[i:3:3*nnel]]

            elseif bctypes[i] == 2
                u_aux[:,i] = x[Dofs[i:3:3*nnel]]
                t_aux[:,i] .= bcvalues[i]
            else
                error("not supported by this func")
            end
        end

        u[mesh.IEN[:,e],:] = u_aux
        t[mesh.IEN[:,e],:] = t_aux
    end


    return u, t
end

function returnut(mesh,x, C=[])
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
    for i in 1:nrrb
        t[mesh.LM[:,elem_ri[i]]] = x[j+1:j+3*nnel*nrb[i]]
        j = j + 3*nnel*nrb[i]
    end

    t[LM[:,iu][:]] = x[j+1:end]

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

function calc_utpoints(mesh,u,t)
    if mesh.eltype>0
        csi_points = range(-1,1,length = mesh.eltype+1)
        csi_nodes = range(-1+mesh.offset,1-mesh.offset,length = mesh.eltype+1)
        csi_vec = calc_csis_grid(csi_points)
        N = calc_N_gen(csi_nodes,csi_vec)
    else
        csi_points = [-1.0,1.0]
        csi_nodes = 0.0
        csi_vec = calc_csis_grid(csi_points)
        N = calc_N_gen(csi_nodes,csi_vec)
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