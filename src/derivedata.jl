function derive_data!(material::Vector{Material}, problem::Problem)

      
        
        # ! Escreve propriedades complexas dos materiais
        for i in 1:length(material)
            material[i].zGe = complex(material[i].Ge, 2.0*material[i].Ge*material[i].Dam)
            material[i].zSwv = sqrt(material[i].zGe/material[i].Rho)
            material[i].zPwv = material[i].zSwv*sqrt((2.0-2.0*material[i].Nu)/(1.0-2.0*material[i].Nu))
        end

        # ! Calcula os pontos de gauss
        assembly = Assembly()
        assembly.nGP = 8
        assembly.csi, assembly.omega = gausslegendre(assembly.nGP)
 

        if !isempty(problem.nFr)
            calc_frequencies!(problem)
        end
        
        return material, problem, assembly

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