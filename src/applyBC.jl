"""
    applyBC(mesh::Mesh, assembly::Assembly,problem::Problem,H,G)

  # Description
  
  Applies boundary conditions defined in problem to the matrices H and G, returning the matrix LHS and the vector RHS. The problem is then described as
  
  ```math
  [LHS]{x} = [RHS]
  ```
  
  # Arguments
  - `n::Integer`: .
  - `mesh::Mesh`: .
  - `assembly::Assembly`: .
  - `problem::Problem`: .
  - `H::Any`: Must be a matrix.
  - `G::Any`: Must be a matrix.
"""
function applyBC_rb(mesh::Mesh, problem::Problem,H,G)

    nnel = size(mesh.IEN,1)
        
    # Definir as tuplas GDLr[i]
    nrb, elemr, GDLr, rbcenters = JuBEM.get_GDLr(mesh, problem)
    
    # Montar matrizes Hrr e Grr
    nGDLflex = size(H,1)
    nGDL =  nGDLflex + 6*nrb
    nGDLr = isempty(GDLr) ? 0 : sum(length,GDLr)
    Hr = zeros(eltype(H),nGDLflex,6*nrb)
    Gr = zeros(eltype(H),nGDLflex,nGDLr)
    D = zeros(eltype(H),6*nrb,nGDLr)

    jg1 = 1
    for r in 1:nrb

        cols = 6*(r-1)+1:6*r

        # Calcula as matrizes C e D
        jg2 = length(GDLr[r])

        Cr, Dr = calc_CD(mesh,rbcenters[r,:],elemr[r])
        
        # Pega as entradas da matrizes H com os graus de liberdade do corpo rígido r
        Hr[:,cols] = H[:,GDLr[r]]*Cr
        Gr[:,jg1:jg2] = G[:,GDLr[r]]
        D[:,jg1:jg2] = Dr

        jg1 = jg2

    end


    GDLu, GDLt,GDLi, y = find_bcidxs(mesh, problem)
    nGDLu = length(GDLu)
    nGDLt = length(GDLt)
    Hu = H[:,GDLu]
    Gu = G[:,GDLu]
    Ht = H[:,GDLt]
    Gt = G[:,GDLt]


    LHS = zeros(eltype(H), nGDL, nGDL)
    RHS = zeros(eltype(H), nGDL)

    LHS[1:nGDLflex,:] = [Hr Ht -Gr -Gu]
    j1 = 6*nrb+nGDLt+1
    j2 = 6*nrb+nGDLt+nGDLr
    LHS[nGDLflex+1:end,j1:j2] = D
   
    RHS = [
        [-Hu Gt]*y
        reshape(problem.forces,:)
    ]

    return LHS, RHS
end

"""
    GDLu, GDLt, y = find_bcidxs(mesh::Mesh, problem::Problem)

  # Description

  Neste função assume-se que as regiões que tem condição de contorno u ou t tem as 3 condições de contorno u ou t

  # Arguments
  - `n::Integer`: desc.
"""
function find_bcidxs(mesh::Mesh, problem::Problem)

    nnel = size(mesh.IEN,1)
    # bcsut = findall(problem.bctype[:,2] .== 1 .|| problem.bctype[:,2] .== 2)
    bcsut = findall([x ∈ [-1,1,2] for x in problem.bctype[:,2]])

    ngdlu = 0
    ngdlt = 0
    ngdli = 0

    for b in bcsut
        bc = problem.bctype[b,2:end]
        bcvalue = problem.bcvalue[b,:]

        ni = count(bc.==-1 .&& bcvalue[1] == 1)
        nu = count(bc.==1)
        nt = count(bc.==2)

        tags = findall(problem.taginfo[:,2].==b)
        nelem = count([x ∈ tags for x in mesh.tag])
        ngdli = ngdli + ni*nnel*nelem
        ngdlu = ngdlu + nu*nnel*nelem
        ngdlt = ngdlt + nt*nnel*nelem
    end

    GDLi = zeros(Int64,ngdli)
    GDLu = zeros(Int64,ngdlu)
    GDLt = zeros(Int64,ngdlt)
    yu = zeros(ngdlu)
    yt = zeros(ngdlt)
    y = zeros(ngdlu+ngdlt)

    ii = 1
    iu = 1
    it = 1
    for b in axes(problem.bctype,1)
        bctype = problem.bctype[b,2:end]
        bcvalue = problem.bcvalue[b,:]
        if any(bctype.==1 .|| bctype.==2 .|| bctype.==-1)
            for t in axes(problem.taginfo,1)
                if problem.taginfo[t,2]==b
                    elems = findall(mesh.tag.==t)
                    for e in elems
                        for n in axes(mesh.IEN,1)
                            nidx = mesh.IEN[n,e]
                            for i in 1:3
                                gdl = mesh.ID[i,nidx]
                                if bctype[i] == 1
                                    GDLu[iu] = gdl
                                    yu[iu] = bcvalue[i]
                                    iu = iu+1
                                elseif bctype[i] == 2
                                    GDLt[it] = gdl
                                    yt[it] = bcvalue[i]
                                    it = it+1
                                elseif bctype[i] == -1 && bcvalue[1] == 1
                                    GDLi[ii] = gdl
                                    ii = ii + 1
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    y = [yu;yt]

    return GDLu, GDLt, GDLi, y

end

"""
    nrb, nelemr, elemr, GDLr = get_GDLr(mesh::Mesh, problem::Problem)

  # Description
  
  Retorna os outputs descritos abaixo 
  
  # Output
  - `nrb::Integer`: número de corpos rígidos
  - `nelemr::Array{Integer,1}`: Vetor que contém a quantidade de elementos em cada corpo rígido
  - `elemr`: Array de vetores que contém os índices dos elementos de cada corpo rígido
  - `GDLr`: Array de vetores que contém os graus de liberdade de cada corpo rígido
"""
function get_GDLr(mesh::Mesh, problem::Problem)

    bctags = @view problem.taginfo[:,2]
    # Indices dos bctypes que são rb
    ibcrb = findall(problem.bctype[:,2].>2)
    nrb = length(ibcrb)
    elemr = Array{Array{Int64,1},1}(undef,nrb)
    GDLr = Array{Array{Int64,1},1}(undef,nrb)
    nelemr = zeros(nrb)
    rbcenters = zeros(nrb,3)

    for i in 1:nrb
        irbtag = findfirst(bctags.==ibcrb[i])
        elemr[i] = findall(mesh.tag.==irbtag)
        nodesr = vec(mesh.IEN[:,elemr[i]])
        GDLr[i] = vec(mesh.ID[:,nodesr])
        nelemr[i] = length(elemr[i])
        rbcenters[i,:] = problem.bcvalue[ibcrb[i],:]
    end 


    return nrb, elemr, GDLr, rbcenters

end


function returnut_rb(x,mesh::Mesh,problem::Problem)


    nnel = Int((mesh.eltype+1)^2.0)
    u = zeros(typeof(x[1]),mesh.nnodes,3)
    t = zeros(typeof(x[1]),mesh.nnodes,3)


    nrb, elemr, GDLr, rbcenters = JuBEM.get_GDLr(mesh, problem)
    
    ## Montar matrizes Hrr e Grr
    nGDL0 = 6*nrb
    nGDLr = isempty(GDLr) ? 0 : sum(length,GDLr)
    urb = zeros(typeof(x[1]),6,nrb)

    for r in 1:nrb

        idxr = 6*(r-1)+1:6*r
        u0 = x[idxr]
        urb[:,r] = u0

        Cr,_ = calc_CD(mesh,rbcenters[r,:],elemr[r])
        
        ur = Cr*u0 
    
        nodesr = vec(mesh.IEN[:,elemr[r]])
        u[nodesr,:] = reshape(ur,3,:)'

    end



    GDLu, GDLt, _ = find_bcidxs(mesh, problem)
    nGDLt = length(GDLt)


    ut = x[nGDL0+1:nGDL0+nGDLt]
    tr = x[nGDL0+nGDLt+1:nGDL0+nGDLt+nGDLr]
    tu = x[nGDL0+nGDLt+nGDLr+1:end]

    j = 1
    for r in 1:nrb
        nodesr = vec(mesh.IEN[:,elemr[r]])
        nnodesr = length(nodesr)
        t[nodesr,:] = reshape(tr[j:3*nnodesr],3,:)'
        j = j+nnodesr+1
    end

    iu = 1
    it = 1
    for b in axes(problem.bctype,1)
        bctype = problem.bctype[b,2:end]
        bcvalue = problem.bcvalue[b,:]
        if any(bctype.==1 .||bctype.==2)
            for j in axes(problem.taginfo,1)
                if problem.taginfo[j,2]==b
                    elems = findall(mesh.tag.==j)
                    for e in elems
                        for n in axes(mesh.IEN,1)
                            nidx = mesh.IEN[n,e]
                            for i in 1:3
                                if bctype[i] == 1
                                    t[nidx,i] = tu[iu]
                                    iu = iu+1
                                elseif bctype[i] == 2
                                    u[nidx,i] = ut[it]
                                    it = it+1
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return u, t, urb
end

function calc_utpoints(mesh::Mesh,u,t)
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

# TODO: Verificar essa formulação com mais calma e estruturar melhor os graus de liberdade de entrada
function calc_CD(mesh::Mesh, rbcenter, rbelem)

    nnel = (mesh.eltype+1)^2

    csis1d, omegas1d = gausslegendre(4)

    csis2d, omegas2d = calc_csis2d(csis1d,omegas1d)

    if mesh.eltype == 0
        csis_cont = range(-1,1,length = 2)
        Nd = ones(size(csis2d,1))
    else
        csis_cont = range(-1,1,length = mesh.eltype+1)
        csis_descont = range(-1+mesh.offset,1 - mesh.offset,length=mesh.eltype+1)
        Nd = calc_N_gen(csis_descont,csis2d)
    end
    Nc = calc_N_gen(csis_cont,csis2d)
    dNcdcsi = calc_N_gen(csis_cont, csis2d;dg=:dNdc)
    dNcdeta = calc_N_gen(csis_cont, csis2d;dg=:dNde)

    nrbelem = length(rbelem)

    nrbnodes = nrbelem*nnel
    C = zeros(3*nrbnodes,6)
    D = zeros(6,3*nrbnodes)

    for e in 1:nrbelem

        # Calculating C
        rbe = rbelem[e]

        for n in 1:nnel
            idxnode = mesh.IEN[n,rbe]
            node_pos = mesh.nodes[idxnode,2:end]
            r = node_pos - rbcenter
            idxno1rb = nnel*(e-1)+n
            idx1 = 3*(idxno1rb-1)+1
            idx2 = 3*idxno1rb
            C[idx1:idx2,:] = [
                1 0 0   0       r[3]    -r[2]
                0 1 0   -r[3]   0       r[1]
                0 0 1   r[2]    -r[1]   0
            ]
        end

        # Calculating D
        rbelem_geo = abs(mesh.EEN[rbe])
        points = mesh.points[mesh.IEN_geo[:,rbelem_geo],2:end]
        _,J = calc_n_J_matrix(dNcdcsi, dNcdeta, points)
        gauss_points = Nc*points
        De = view(D,:,3*nnel*(e-1)+1:3*nnel*(e))
        for g in eachindex(J)
            r = gauss_points[g,:] - rbcenter

            Dg = [
                1 0 0
                0 1 0
                0 0 1
                0 -r[3] r[2]
                r[3] 0 -r[1]
                -r[2] r[1] 0
            ]

            for k in 1:nnel
                De[:,3*(k-1)+1:3*k] = De[:,3*(k-1)+1:3*k] + Dg.*J[g].*Nd[g,k]*omegas2d[g]
            end
        end
    end
    return C, D
end

function applyBC_simple(mesh::Mesh,problem::Problem,H,G)

    nnel = size(mesh.IEN,1)

    nDofs = size(H,1)

    LHS = zeros(eltype(H),nDofs,nDofs)
    RHS = zeros(eltype(H),nDofs)

    RHS_mat = zeros(eltype(H),nDofs,nDofs)
    RHS_vec = zeros(eltype(H),nDofs)
    

    for e in 1:mesh.nelem

        Dofs = mesh.LM[:,e]

        tag = mesh.tag[e]
        bctag = problem.taginfo[tag,2]
        bctypes = problem.bctype[bctag,2:end]
        bcvalues = problem.bcvalue[bctag,:]

        LHS_aux = zeros(eltype(H),nDofs,3*nnel)
        RHS_aux = zeros(eltype(H),nDofs,3*nnel)

        for i in 1:3
            if bctypes[i] == 1
                LHS_aux[:,i:3:3*nnel] = -G[:,Dofs[i:3:3*nnel]]
                RHS_aux[:,i:3:3*nnel] = -H[:,Dofs[i:3:3*nnel]]
                RHS_vec[Dofs[i:3:3*nnel]] .= bcvalues[i]

            elseif bctypes[i] == 2
                LHS_aux[:,i:3:3*nnel] = H[:,Dofs[i:3:3*nnel]]
                RHS_aux[:,i:3:3*nnel] = G[:,Dofs[i:3:3*nnel]]
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
    
    u = zeros(eltype(x),mesh.nnodes,3)
    t = zeros(eltype(x),mesh.nnodes,3)

    
    for e in 1:mesh.nelem
        Dofs = mesh.LM[:,e]

        tag = mesh.tag[e]
        bctag = problem.taginfo[tag,2]
        bctypes = problem.bctype[bctag,2:end]
        bcvalues = problem.bcvalue[bctag,:]

        u_aux = zeros(eltype(x),nnel,3)
        t_aux = zeros(eltype(x),nnel,3)

        
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


function applyBC_multi_simple(mesh::Mesh,problem::Problem,H,G)

    nnel = size(mesh.IEN,1)

    nDofs = size(H,1)

    LHS = zeros(eltype(H),nDofs,nDofs)
    RHS = zeros(eltype(H),nDofs)
    
    # Primeiro colapsar as colunas de H e de G

    H2, G2 = reorganizeHG_multi(mesh,problem,H,G)

    # Reorganiza colunas de u e t

    for e in 1:mesh.nelem

        tagidx = mesh.tag[e]
        bcidx = problem.taginfo[tagidx,2]
        bctypes = problem.bctype[bcidx,2:end]
        bcvalues = problem.bcvalue[bcidx,:]

        if bctypes[1] != -1

            for n in 1:nnel
                nodeidx = mesh.IEN[n,e]
                Dofs = mesh.ID[:,nodeidx]

                for i in 1:3

                    if bctypes[i] == 1
                        RHS = RHS -H2[:,Dofs[i]].*bcvalues[i]
                        H2[:,Dofs[i]] = -G2[:,Dofs[i]]
        
                    elseif bctypes[i] == 2
                        RHS = RHS + G2[:,Dofs[i]].*bcvalues[i]
                    else                    
                        error("not supported by this func")
                    end
                end
            end

        end

    end

    # Find all columns of intersection of G2, put them on the end of H2 and put zeros on them back
    idxbcinterfaces = findall(problem.bctype[:,2].==-1)
    idxtaginterfaces = findall([x ∈ idxbcinterfaces for x in problem.taginfo[:,2]])
    idxelemsinterfaces = findall([x ∈ idxtaginterfaces for x in mesh.tag])
    nodesidxinterfaces = unique(vec(mesh.IEN[:,idxelemsinterfaces]))
    gdlinterfaces = vec(mesh.ID[:,nodesidxinterfaces])

    LHS = [H2 -G2[:,gdlinterfaces]]
    # G2[:,gdlinterfaces] .= 0

    # RHS = G2*RHS_vec

    return LHS, RHS

end

function reorganizeHG_multi(mesh::Mesh,problem::Problem,H,G)

    nnel = size(mesh.IEN,1)
    nDofsi = maximum(mesh.LM)
    nDofsj = maximum(mesh.ID)

    Haux = zeros(eltype(H),nDofsi,nDofsj)
    Gaux = zeros(eltype(G),nDofsi,nDofsj)

    for e in 1:mesh.nelem

        for n in 1:nnel
            nodeidx = mesh.IEN[n,e]
            sinal = 1

            tagidx = mesh.tag[e]
            bcidx = problem.taginfo[tagidx,2]
            bctype = problem.bctype[bcidx,2]

            if bctype == -1 
                sinal = problem.bcvalue[bcidx,1]
            end
            Haux[:,mesh.ID[:,nodeidx]] = Haux[:,mesh.ID[:,nodeidx]] + H[:,mesh.LM[3*(n-1)+1:3*n,e]]
            Gaux[:,mesh.ID[:,nodeidx]] = Gaux[:,mesh.ID[:,nodeidx]] + sinal.*G[:,mesh.LM[3*(n-1)+1:3*n,e]]
            
        end

    end

    return Haux, Gaux
end

function returnut_multi_simple(x,mesh::Mesh,problem::Problem)

    nnel = size(mesh.IEN,1)
    
    nnodeDofs = maximum(mesh.ID)

    nnodes = maximum(mesh.IEN)

    u = zeros(eltype(x),nnodes,3)
    t = zeros(eltype(x),nnodes,3)

    
    for e in 1:mesh.nelem
        for n in 1:nnel
            nodeidx = mesh.IEN[n,e]
            Dofs = mesh.ID[:,nodeidx]

            tag = mesh.tag[e]
            bctag = problem.taginfo[tag,2]
            bctypes = problem.bctype[bctag,2:end]
            bcvalues = problem.bcvalue[bctag,:]
            
            for i in 1:3
                if bctypes[i] == 1
                    u[nodeidx,i] = bcvalues[i]
                    t[nodeidx,i] = x[Dofs[i]]

                elseif bctypes[i] == 2
                    u[nodeidx,i] = x[Dofs[i]]
                    t[nodeidx,i] = bcvalues[i]
                elseif bctypes[i] == -1
                    u[nodeidx,1] = x[Dofs[i]]
                end
            end
        end
    end

    idxbcinterfaces = findall(problem.bctype[:,2].==-1)
    idxtaginterfaces = findall([x ∈ idxbcinterfaces for x in problem.taginfo[:,2]])
    idxelemsinterfaces = findall([x ∈ idxtaginterfaces for x in mesh.tag])
    nodesidxinterfaces = unique(vec(mesh.IEN[:,idxelemsinterfaces]))
    gdlinterfaces = vec(mesh.ID[:,nodesidxinterfaces])

    j = 1
    for n in nodesidxinterfaces
        for i in 1:3
            t[n,i] = x[nnodeDofs+j]    
            j = j+1
        end
    end

    return u, t
end

function applyBC_rb_multi(mesh::Mesh,problem::Problem,H,G)

    nnel = size(mesh.IEN,1)
    

    nnel = size(mesh.IEN,1)

    nDofs = size(H,1)

    # LHS = zeros(eltype(H),nDofs,nDofs)
    # RHS = zeros(eltype(H),nDofs)
    
    # Primeiro colapsar as colunas de H e de G

    HID, GID = reorganizeHG_multi(mesh,problem,H,G)


    # Definir as tuplas GDLr[i]
    nrb, elemr, GDLr, rbcenters = JuBEM.get_GDLr(mesh, problem)
    
    # Montar matrizes Hrr e Grr
    nGDLflex = size(H,1)
    nGDL =  nGDLflex + 6*nrb
    nGDLr = isempty(GDLr) ? 0 : sum(length,GDLr)
    Hr = zeros(eltype(H),nGDLflex,6*nrb)
    Gr = zeros(eltype(H),nGDLflex,nGDLr)
    D = zeros(eltype(H),6*nrb,nGDLr)

    jg1 = 1
    for r in 1:nrb

        cols = 6*(r-1)+1:6*r

        # Calcula as matrizes C e D
        jg2 = length(GDLr[r])

        Cr, Dr = calc_CD(mesh,rbcenters[r,:],elemr[r])
        
        # Pega as entradas da matrizes H com os graus de liberdade do corpo rígido r
        Hr[:,cols] = HID[:,GDLr[r]]*Cr
        Gr[:,jg1:jg1+jg2-1] = GID[:,GDLr[r]]
        D[:,jg1:jg1+jg2-1] = Dr

        jg1 = jg2

    end


    GDLu, GDLt,GDLi, y = find_bcidxs(mesh, problem)
    nGDLu = length(GDLu)
    nGDLt = length(GDLt)
    nGDLi = length(GDLi)
    Hu = @view HID[:,GDLu]
    Gu = @view GID[:,GDLu]
    Ht = @view HID[:,GDLt]
    Gt = @view GID[:,GDLt]
    Hi = @view HID[:,GDLi]
    Gi = @view GID[:,GDLi]


    LHS = zeros(eltype(H), nGDL, nGDL)
    RHS = zeros(eltype(H), nGDL)
    # @infiltrate
    LHS[1:nGDLflex,:] = [Hr -Gr Ht -Gu Hi -Gi]
    j1 = 6*nrb+1
    j2 = 6*nrb+nGDLr
    LHS[nGDLflex+1:end,j1:j2] = D
   
    RHS = [
        [-Hu Gt]*y
        reshape(problem.forces,:)
    ]

    return LHS, RHS
end

function returnut_rb_multi(x,mesh::Mesh,problem::Problem)
    
    nnel = Int((mesh.eltype+1)^2.0)
    nnodes = maximum(mesh.IEN)
    u = zeros(typeof(x[1]),nnodes,3)
    t = zeros(typeof(x[1]),nnodes,3)


    nrb, elemr, GDLr, rbcenters = JuBEM.get_GDLr(mesh, problem)
    
    # Montar matrizes Hrr e Grr
    nGDL0 = 6*nrb
    nGDLr = isempty(GDLr) ? 0 : sum(length,GDLr)
    urb = zeros(typeof(x[1]),6,nrb)

    for r in 1:nrb

        idxr = 6*(r-1)+1:6*r
        u0 = x[idxr]
        urb[:,r] = u0

        Cr,_ = calc_CD(mesh,rbcenters[r,:],elemr[r])
        
        ur = Cr*u0 
    
        nodesr = vec(mesh.IEN[:,elemr[r]])
        u[nodesr,:] = reshape(ur,3,:)'

    end



    GDLu, GDLt,GDLi, _ = find_bcidxs(mesh, problem)
    nGDLt = length(GDLt)
    nGDLi = length(GDLi)
    nGDLu = length(GDLu)

    tr = x[nGDL0+1:nGDL0+nGDLr]
    ut = x[nGDL0+nGDLr+1:nGDL0+nGDLr+nGDLt]
    tu = x[nGDL0+nGDLt+nGDLr+1:nGDL0+nGDLt+nGDLr+nGDLu]
    ui = x[nGDL0+nGDLt+nGDLr+nGDLu+1:nGDL0+nGDLt+nGDLr+nGDLu+nGDLi]
    ti = x[nGDL0+nGDLt+nGDLr+nGDLu+nGDLi+1:end]

    j = 1
    for r in 1:nrb
        nodesr = vec(mesh.IEN[:,elemr[r]])
        nnodesr = length(nodesr)
        t[nodesr,:] = reshape(tr[j:3*nnodesr],3,:)'
        j = j+nnodesr+1
    end

    iu = 1
    it = 1
    ii = 1
    for b in axes(problem.bctype,1)
        bctype = problem.bctype[b,2:end]
        bcvalue = problem.bcvalue[b,:]
        if any(bctype.==1 .||bctype.==2)
            for j in axes(problem.taginfo,1)
                if problem.taginfo[j,2]==b
                    elems = findall(mesh.tag.==j)
                    for e in elems
                        for n in axes(mesh.IEN,1)
                            nidx = mesh.IEN[n,e]
                            for i in 1:3
                                if bctype[i] == 1
                                    t[nidx,i] = tu[iu]
                                    iu = iu+1
                                elseif bctype[i] == 2
                                    u[nidx,i] = ut[it]
                                    it = it+1
                                elseif bctype[i] == -1 && bcvalue[i] == 1
                                    u[nidx,i] = ui[ii]
                                    t[nidx,i] = ti[ii]
                                    ii = ii+1
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return u, t, urb
end