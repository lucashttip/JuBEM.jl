function remove_EE!(mesh::Mesh,assembly::Assembly,problem::Problem)

    bctypeidxee = findall(x->x == 0 , problem.bctype[:,2])
    tagidxee = findall(x->x ∈ bctypeidxee , problem.taginfo[:,2])
    elemidxee = findall(x->x ∈ tagidxee, mesh.tag)
    elemidxnonee = setdiff(1:mesh.nelem, elemidxee)

    nnel = size(mesh.IEN,1)
    nelem_new = length(elemidxnonee)
    nnodes_new = nnel*length(elemidxnonee)

    ID_new = zeros(Int32,size(mesh.ID))
    LM_new = zeros(Int32,3*nnel,nelem_new)
    IEN_new = zeros(Int32,nnel,nelem_new)
    IEN_geo_new = zeros(Int32,nnel,nelem_new)
    tag_new = zeros(Int32,nelem_new)

    
    # Remover as colunas de LM, IEN, IEN_geo
    mesh.nelem = nelem_new

    pos = 1

    for e in 1:nelem_new
        IEN_new[:,e] = mesh.IEN[:,elemidxnonee[e]]
        IEN_geo_new[:,e] = mesh.IEN_geo[:,elemidxnonee[e]]
        for i in 1:nnel
            ID_new[:,IEN_new[i,e]] = [pos, pos+1, pos+2]
            pos = pos+3
        end
        LM_new[:,e] = ID_new[:,IEN_new[:,e]][:]
    end

    tag_new = mesh.tag[elemidxnonee]


    G = zeros(3*nnodes_new, 3*nnodes_new)
    H = zeros(3*nnodes_new, 3*nnodes_new)
    # Remover as colunas de G e de H


    G[LM_new[:,elemidxnonee][:],LM_new[:,elemidxnonee][:]] = assembly.G[mesh.LM[:,elemidxnonee][:],mesh.LM[:,elemidxnonee][:]]
    H[LM_new[:,elemidxnonee][:],LM_new[:,elemidxnonee][:]] = assembly.H[mesh.LM[:,elemidxnonee][:],mesh.LM[:,elemidxnonee][:]]

    println(size(assembly.G[mesh.LM[:,elemidxnonee][:],mesh.LM[:,elemidxnonee][:]]))

    assembly.G = G
    assembly.H = H

    # Update mesh
    mesh.nelem = nelem_new
    mesh.nnodes = nnodes_new
    mesh.ID = ID_new
    mesh.LM = LM_new
    mesh.IEN = IEN_new
    mesh.IEN_geo = IEN_geo_new
    mesh.tag = tag_new

    
    return mesh, assembly
end



function remove_EE_mesh!(mesh)

    nonee = findall(x->x != 0 , mesh.bc[:,1])
    
    # Remover as colunas de LM, IEN, IEN_geo
    mesh.nelem = length(nonee)
    mesh.LM = mesh.LM[:,nonee]
    mesh.IEN = mesh.IEN[:,nonee]
    mesh.IEN_geo = mesh.IEN_geo[:,nonee]
    mesh.bc = mesh.bc[nonee,:]
    mesh.bcvalue = mesh.bcvalue[nonee,:]

    # Remover as colunas de ID, nodes, bc e bcvalue
    mesh.nnodes = length(unique(mesh.IEN[:]))  #Apenas para descontínuo
    mesh.ID = mesh.ID[:,unique(sort(mesh.IEN[:]))]
    mesh.nodes = mesh.nodes[unique(sort(mesh.IEN[:])),:]

    return mesh
end