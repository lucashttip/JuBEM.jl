function remove_EE!(mesh::Mesh,assembly::Assembly)

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

    # Remover as colunas de G e de H
    assembly::Assembly.G = assembly::Assembly.G[mesh.LM[:],mesh.LM[:]]
    assembly::Assembly.H = assembly::Assembly.H[mesh.LM[:],mesh.LM[:]]


    return mesh, solver_var
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