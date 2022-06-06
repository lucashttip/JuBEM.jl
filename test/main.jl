# !> @brief     Main program to calculate the Dynamic Soil response with BEM
# !> @author    Lucas Agatti Pacheco
# !> @author    Ronaldo Carrion
# !> @author    Euclides Mesquita Neto
# !> @author    Amauri Ferraz
# !> @version   j_0.1
# !> @date      2022
#     ! =======================================================
#     ! =======================================================
#     ! ================ Modules and variables ================
#     ! =======================================================
#     ! =======================================================
using Revise
using JuBEM

    inp_file = "bar.msh"
#     ! =================================================
#     ! =================================================
#     ! ===================== Input =====================
#     ! =================================================
#     ! =================================================

    # ! ## Read from input
    mesh, material, problem, solver_var = read_msh(inp_file)

    # Calculate material, problem and solver variables:
    derive_data!(material, problem, solver_var)

    # ! Generate physical mesh:
    generate_mesh!(mesh)

    # ! Create outputs
    # call create_output (out_file, mesh, material(1), solver_var%nGP, problem%nFr,problem%fr_range, out_fid,&
    #     ux_fid, uy_fid, uz_fid, phix_fid, phiy_fid, phiz_fid)
    
    # ! ================================================
    # ! ================================================
    # ! ==================== Solver ====================
    # ! ================================================
    # ! ================================================


    # do i = 1, size(problem%frequencies)
    # problem.frequency = problem.frequencies[1]
    # println("Rodando o para a frequencia 1: ", string(problem.frequency))
        # problem%frequency = problem%frequencies(i)
        # write(*,"(A, X, I3, 2X, A, F10.4)") 'rodando para frequencia ',i , 'Frequencia: ', problem%frequency
        
        # ! ## Calculate G and H
        calc_GH!(mesh, material, problem, solver_var)

        # ! ## Apply BC, arranging Ax = b
        applyBC_nonrb!(mesh, solver_var)

        # ! ## Solve system Ax = b
        # call zcgesv (3*mesh%nelem+6, 1, solver_var%zma, 3*mesh%nelem+6, IPIV, mesh%zbcvalue, 3*mesh%nelem+6, &
        # solver_var%zvetsol, 3*mesh%nelem+6, WORK, SWORK, RWORK, itersolve, infosolve)
        # print *, itersolve
        # print *, infosolve

        # ! output
        # call output(mesh,mesh%zbcvalue,solver_var%zvetsol,problem%frequency,out_fid,ux_fid, uy_fid, uz_fid, &
        # phix_fid, phiy_fid, phiz_fid)

        # ! ## Retrieve u and t
        

    # end do
    
    # ! ===================================================
    # ! ===================================================
    # ! ================= Post processing =================
    # ! ===================================================
    # ! ===================================================

    ## calculate surface tensions


    ## Calculate responses on internal points


    ## Save output


    ## Free spaces and close files