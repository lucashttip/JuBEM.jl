# incomplete
# Inspirado no pacote MeshIO.jl: https://github.com/JuliaIO/MeshIO.jl/blob/master/src/io/msh.jl
@enum MSHBlockType MSHFormatBlock MSHPhysicalNamesBlock MSHNodesBlock MSHElementsBlock MSHUnknownBlock MSHMaterialBlock MSHFrequenciesBlock MSHMeshTypeBlock MSHForcesBlock MSHEntitiesBlock MSHBoundaryConditionsBlock MSHTagInformationBlock MSHProblemNameBlock


"""
mesh = read_msh(msh_file)
"""
function read_msh(msh_file)

    io = open(msh_file,"r")

    mesh = Mesh()

    s_entities=[]

    while !eof(io)
        BlockType = parse_blocktype!(io)
        if BlockType == MSHEntitiesBlock
            s_entities = parse_entities(io)
        elseif BlockType == MSHPhysicalNamesBlock
            parse_physicalnames!(io, mesh)
        elseif BlockType == MSHNodesBlock
            parse_nodes!(io, mesh)
        elseif BlockType == MSHElementsBlock
            parse_elements!(io, mesh, s_entities)
        else
            skip_block!(io)
        end
    end
    
    close(io)

    return mesh


end

"""
problem,materials = read_problem(prob_file,mesh)
"""
function read_problem(prob_file,mesh)
    
    io = open(prob_file,"r")

    problem = Problem()
    materials = Material[]

    while !eof(io)
        BlockType = parse_blocktype!(io)
        if BlockType == MSHMaterialBlock
            parse_materials!(io,materials)
        elseif BlockType == MSHFrequenciesBlock
            parse_frequencies!(io, problem)
        elseif BlockType == MSHMeshTypeBlock
            parse_meshtype!(io, mesh)
        elseif BlockType == MSHForcesBlock
            parse_forces!(io, problem)
        elseif BlockType == MSHBoundaryConditionsBlock
            parse_BC!(io, problem)
        elseif BlockType == MSHTagInformationBlock
            parse_taginfo!(io, problem)
        elseif BlockType == MSHProblemNameBlock
            parse_problemname!(io, problem)
        else
            skip_block!(io)
        end
    end
    
    close(io)

    return problem,materials
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
    elseif header == "\$BoundaryConditions"
        return MSHBoundaryConditionsBlock
    elseif header == "\$TagInformation"
        return MSHTagInformationBlock
    elseif header == "\$ProblemName"
        return MSHProblemNameBlock
    else
        return MSHUnknownBlock
    end
end

function parse_problemname!(io,problem)

    probname = readline(io)
    problem.name = probname
    endblock = readline(io)
    if endblock != "\$EndProblemName"
        error("expected end block tag, got $endblock")
    end
    return problem
end

function parse_materials!(io, materials)
    nmat = parse(Int,readline(io))
    for i in 1:nmat
        Ge, Nu, Dam, Rho = parse.(Float64, split(readline(io)))
        push!(materials,Material(Ge,Nu,Dam,Rho))
    end
    endblock = readline(io)
    if endblock != "\$EndMaterial"
        error("expected end block tag, got $endblock")
    end
    return materials
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

function parse_meshtype!(io, mesh)
    mesh.eltype = parse(Int8, readline(io))
    mesh.offset = parse(Float64, readline(io))
    endblock = readline(io)
    if endblock != "\$EndMeshType"
        error("expected end block tag, got $endblock")
    end
    return mesh
end

function parse_forces!(io, problem)
    nrb = parse(Int64, readline(io))
    problem.forces = zeros(6,nrb)
    for i in 1:nrb
        data = split(readline(io))
        problem.forces[:,i] = parse.(Float64, data)
    end
    endblock = readline(io)
    if endblock != "\$EndForces"
        error("expected end block tag, got $endblock")
    end
    return problem
end

function parse_physicalnames!(io,mesh)
    
    nphys = parse(Int,readline(io))

    for i in 1:nphys
        data = split(readline(io),"\"")
        push!(mesh.tagnames,data[2])
    end
    endblock = readline(io)

    if endblock != "\$EndPhysicalNames"
        error("expected end block tag, got $endblock")
    end

    return mesh
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

function parse_elements!(io, mesh, s_entities)
        
    num_entity_blocks, num_elements, min_element_tag, max_element_tag = parse.(Int, split(readline(io)))

    mesh.nelem = num_elements
    
    pos = position(io)

    dim, tag, element_type, elements_in_block = parse.(Int, split(readline(io)))

    seek(io,pos)

    if element_type == 3
        npel = 4
    elseif element_type == 10
        npel = 9
    else
        error("Element type not supported by JuBEM")
    end

    mesh.IEN_geo = zeros(Int32, npel, num_elements)
    mesh.tag = zeros(Int16, num_elements)

    for index_entity in 1:num_entity_blocks

        dim, tag, element_type, elements_in_block = parse.(Int, split(readline(io)))
        s = findfirst(s_entities[:,1].==tag)
        # @infiltrate
        ptag = Int(s_entities[s,2])

        if element_type == 3 || element_type == 10
            for i in 1:elements_in_block
                e, n... = parse.(Int, split(readline(io)))
                mesh.IEN_geo[:,e] = n
                mesh.tag[e] = ptag
            end
        else
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

function parse_BC!(io,problem)
    nbc = parse(Int,readline(io))

    bc = zeros(nbc,4)
    bcvalue = zeros(nbc,3)
    mat = zeros(nbc,2)

    rbidx = 3
    rbidx2 = 3

    for i in 1:nbc
        data = split(readline(io))
        bc[i,1] = parse(Float64, data[1])
        bcvalue[i,:] = parse.(Float64,data[[3,5,7]])
        for j in 1:3
            if data[2*j] == "u"
                bc[i,j+1] = 1
            elseif data[2*j] == "t"
                bc[i,j+1] = 2
            elseif data[2*j] == "rb"
                bc[i,j+1] = rbidx
                rbidx2 = rbidx+1
            elseif data[2*j] == "ee"
                bc[i,j+1] = 0
            elseif data[2*j] == "i"
                bc[i,j+1] = -1
            end
        end

        if rbidx != rbidx2
            rbidx = rbidx2
        end
    end

    problem.bctype = bc
    problem.bcvalue = bcvalue

    endblock = readline(io)

    if endblock != "\$EndBoundaryConditions"
        error("expected end block tag, got $endblock")
    end

    return problem
end

function parse_taginfo!(io,problem)
    ntags = parse(Int64, readline(io))

    problem.taginfo = zeros(Int16,ntags,4)

    for i in 1:ntags
        data = split(readline(io))
        problem.taginfo[i,:] = parse.(Float64, data)
    end
    endblock = readline(io)
    if endblock != "\$EndTagInformation"
        error("expected end block tag, got $endblock")
    end
    return problem
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