function derive_data(mesh::mesh_type, material::material_table_type, problem::problem_type, solver_var::solver_var_type)

      
        
        # ! Escreve propriedades complexas dos materiais
        for i in 1:length(material)
            material[i].zGe = complex(material[i].Ge, 2.0*material[i].Ge*material[i].Dam)
            material[i].zSwv = sqrt(material[i].zGe/material[i].Rho)
            material[i].zPwv = material[i].zSwv*sqrt((2.0-2.0*material[i].Nu)/(1.0-2.0*material[i].Nu))
        end

        # ! Calcula os pontos de gauss
        solver_var.csi, solver_var.omega = gausslegendre(solver_var.nGP)
 
        # ! Calcula as frequencias
        # k = 1
        # associate(frequencies => problem%frequencies, nFr => problem%nFr, fr_range=> problem%fr_range)
        #     frequencies(1) = fr_range(1)
        #     do i = 1,size(nFr)
        #         fr_step = (fr_range(i+1) - fr_range(i))/real(nFr(i)-1,kind=dp)
        #         do j = 1, nFr(i)-1
        #             frequencies(k+j) = fr_range(i) + j*fr_step
        #         end do
        #         k = k + nFr(i) - 1 
        #     end do
        # end associate


end