# incomplete
# Inspirado no pacote MeshIO.jl: https://github.com/JuliaIO/MeshIO.jl/blob/master/src/io/msh.jl
@enum MSHBlockType MSHFormatBlock MSHPhysicalNamesBlock MSHNodesBlock MSHElementsBlock MSHUnknownBlock MSHMaterialBlock MSHFrequenciesBlock MSHMeshTypeBlock MSHForcesBlock MSHEntitiesBlock

function read_msh(inp_file)
    
    io = open(inp_file,"r")

    material = material_table_type[]
    mesh = mesh_type()
    problem = problem_type()
    solver_var = solver_var_type()

    F = zeros(6)

    vars=[]

    while !eof(io)
        BlockType = parse_blocktype!(io)
        if BlockType == MSHMaterialBlock
            parse_materials!(io, material)
        elseif BlockType == MSHFrequenciesBlock
            parse_frequencies!(io, problem)
        elseif BlockType == MSHMeshTypeBlock
            parse_meshtype!(io, mesh, solver_var)
        elseif BlockType == MSHForcesBlock
            parse_forces!(io,F)
        elseif BlockType == MSHPhysicalNamesBlock
            bc, bcvalue, mat = parse_physicalnames(io)
            push!(vars,bc)
            push!(vars,bcvalue)
            push!(vars,mat)
        elseif BlockType == MSHEntitiesBlock
            s_entities = parse_entities(io)
            push!(vars,s_entities)
            # println(s_entities)
        elseif BlockType == MSHNodesBlock
            parse_nodes!(io, mesh)
        elseif BlockType == MSHElementsBlock
            parse_elements!(io, mesh, vars[4], vars[1], vars[2], vars[3])
        else
            skip_block!(io)
        end
    end

    return mesh,material,problem,solver_var
end

function parse_blocktype!(io)
    header = readline(io)
    if header == "\$Material"
        return MSHMaterialBlock
    elseif header == "\$Frequencies"
        return MSHFrequenciesBlock
    elseif header == "\$MeshType"
        return MSHMeshTypeBlock
    elseif header == "\$Forces"
        return MSHForcesBlock
    elseif header == "\$MeshFormat"
        return MSHFormatBlock
    elseif header == "\$PhysicalNames"
        return MSHPhysicalNamesBlock
    elseif header == "\$Entities"
        return MSHEntitiesBlock
    elseif header == "\$Nodes"
        return MSHNodesBlock
    elseif header == "\$Elements"
        return MSHElementsBlock
    else
        return MSHUnknownBlock
    end
end

function parse_materials!(io, material)
    nmat = parse(Int,readline(io))
    for i in 1:nmat
        Ge, Nu, Dam, Rho = parse.(Float64, split(readline(io)))
        push!(material,material_table_type(Ge,Nu,Dam,Rho))
    end
    endblock = readline(io)
    if endblock != "\$EndMaterial"
        error("expected end block tag, got $endblock")
    end
    return material
end

function parse_frequencies!(io, problem)
    nfr_range = parse(Int,readline(io))
    for i in 1:nfr_range
        fi, ff, fr_range = parse.(Float64, split(readline(io)))
        if i == 1
            problem.fr_range = [fi; ff]
        else
            problem.fr_range = [problem.fr_range; ff]
        end
        problem.nFr = [problem.nFr; fr_range]
    end
    endblock = readline(io)
    if endblock != "\$EndFrequencies"
        error("expected end block tag, got $endblock")
    end
    return problem
end

function parse_meshtype!(io, mesh, solver_var)
    solver_var.nGP = parse(Int16, readline(io))
    mesh.offset = parse(Float64, readline(io))
    mesh.eltype = parse(Int8, readline(io))
    endblock = readline(io)
    if endblock != "\$EndMeshType"
        error("expected end block tag, got $endblock")
    end
    return mesh, solver_var
end

function parse_forces!(io, F)
    for i in 1:6
        F[i] = parse(Float64, readline(io))
    end
    endblock = readline(io)
    if endblock != "\$EndForces"
        error("expected end block tag, got $endblock")
    end
    return F
end

function parse_physicalnames(io)
    
    nphys = parse(Int,readline(io))

    bc = zeros(nphys,4)
    bcvalue = zeros(nphys,3)
    mat = zeros(nphys,2)

    for i in 1:nphys
        data = split(replace(readline(io), "\"" =>""))
        bc[i,1] = parse(Float64, data[2])
        mat[i,:] = [parse(Float64, data[2]); parse(Float64, data[end])]
        bcvalue[i,:] = parse.(Float64,data[[4,6,8]])
        for j in 1:3
            if data[2*j+1] == "u"
                bc[i,j+1] = 1
            elseif data[2*j+1] == "t"
                bc[i,j+1] = 2
            elseif data[2*j+1] == "rb"
                bc[i,j+1] = 3
            end
        end
    end
    endblock = readline(io)

    if endblock != "\$EndPhysicalNames"
        error("expected end block tag, got $endblock")
    end

    return bc, bcvalue, mat
end

function parse_entities(io)
    npoints, ncurves, nsurfaces = parse.(Int,split(readline(io))[1:end-1])

    s_entities = zeros(nsurfaces,2)

    for i in 1:(npoints+ncurves)
        readline(io)
    end

    for i in 1:nsurfaces
        s_entities[i,:] = parse.(Float64,split(readline(io))[[1,9]])
    end
    endblock = readline(io)
    if endblock != "\$EndEntities"
        error("expected end block tag, got $endblock")
    end
    return s_entities
end

function parse_nodes!(io,mesh)
    entity_blocks, num_nodes, min_node_tag, max_node_tag = parse.(Int, split(readline(io)))
    mesh.npoints = num_nodes
    mesh.points = zeros(num_nodes,4)

    for index_entity in 1:entity_blocks
        dim, tag, parametric, nodes_in_block = parse.(Int, split(readline(io)))
        node_tags = zeros(Int,nodes_in_block)
        for i in 1:nodes_in_block
            node_tags[i] = parse(Int, readline(io))
        end
        mesh.points[node_tags,1] = node_tags
        for i in 1:nodes_in_block
            xyz = parse.(Float64, split(readline(io)))
            mesh.points[node_tags[i],2:end] = xyz
        end
    end
    endblock = readline(io)
    if endblock != "\$EndNodes"
        error("expected end block tag, got $endblock")
    end
    return mesh
end

function parse_elements!(io, mesh, s_entities, bc, bcvalue, mat)
    num_elements = parse.(Int, split(readline(io)))
    
    num_entity_blocks, num_elements, min_element_tag, max_element_tag = num_elements
    mesh.nelem = num_elements
    
    mesh.IEN_geo = zeros(Int32, 4, num_elements)
    mesh.bc = zeros(Int16,num_elements,3)
    mesh.bcvalue = zeros(Float64,num_elements,3)
    mesh.material = zeros(Int16,num_elements)

    for index_entity in 1:num_entity_blocks

        dim, tag, element_type, elements_in_block = parse.(Int, split(readline(io)))
        s = findfirst(s_entities[:,1].==tag)
        ptag = Int(s_entities[s,2])

        if element_type == 3 # Quadrangles
            for i in 1:elements_in_block
                e, n1, n2, n3, n4 = parse.(Int, split(readline(io)))
                mesh.IEN_geo[:,e] = [n1, n2, n3, n4]
                mesh.bc[e,:] = bc[ptag,2:end]
                mesh.bcvalue[e,:] = bcvalue[ptag,:]
                mesh.material[e] = mat[ptag,2]
            end
        else
            # for now we ignore all other elements (points, lines, hedrons, etc)
            for i in 1:elements_in_block
                readline(io)
            end
        end
    end
    endblock = readline(io)
    if endblock != "\$EndElements"
        error("expected end block tag, got $endblock")
    end
    return mesh
end

function skip_block!(io)
    while true
        line = readline(io)
        if length(line) < 4
            continue
        end
        if line[1:4] == "\$End"
            break
        end
    end
    return nothing
end