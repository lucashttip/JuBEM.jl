function derive_data!(mesh::Mesh,material::Vector{Material}, problem::Problem)

      
        
        # ! Escreve propriedades complexas dos materiais
        for i in 1:length(material)
            material[i].zGe = complex(material[i].Ge, 2.0*material[i].Ge*material[i].Dam)
            material[i].zSwv = sqrt(material[i].zGe/material[i].Rho)
            material[i].zPwv = material[i].zSwv*sqrt((2.0-2.0*material[i].Nu)/(1.0-2.0*material[i].Nu))
            
            C_stat = zeros(4)
            C_stat[1]=1.0/(16.0*pi*material[i].Ge*(1.0-material[i].Nu))
            C_stat[2]=3.0-(4.0*material[i].Nu)
            C_stat[3]=-1.0/(8.0*pi*(1.0-material[i].Nu))
            C_stat[4]=1.0-(2.0*material[i].Nu)
            material[i].C_stat = C_stat
        end


        if !isempty(problem.nFr)
            calc_frequencies!(problem)
        end
        
    # Definicao de bodies
    nbodies = maximum(problem.taginfo[:,4])

    for b in 1:nbodies
        tagidxs = findall(problem.taginfo[:,4].==b)
        elems = findall([x ∈ tagidxs for x in mesh.tag])
        push!(mesh.bodies,elems)
    end


end

function calc_frequencies!(problem)
    # ! Calcula as frequencias
    k = 1
        
    problem.frequencies = zeros(sum(problem.nFr) - length(problem.nFr)+1)

    problem.frequencies[1] = problem.fr_range[1]
    for i in 1:length(problem.nFr)
        fr_step = (problem.fr_range[i+1] - problem.fr_range[i])/float(problem.nFr[i]-1)
        for j in 1:problem.nFr[i]-1
            problem.frequencies[k+j] = problem.fr_range[i] + j*fr_step
        end
        k = k + problem.nFr[i] - 1 
    end

    return problem
end

function calc_static_constants(material)

    C_stat = zeros(4)
    
    C_stat[1]=1.0/(16.0*pi*material.Ge*(1.0-material.Nu))
    C_stat[2]=3.0-(4.0*material.Nu)
    C_stat[3]=-1.0/(8.0*pi*(1.0-material.Nu))
    C_stat[4]=1.0-(2.0*material.Nu)

    return C_stat

end