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

#     ! # Variable Declaration
#=
    #     ! ## IO Variables
    #     character(len=1024), parameter :: inp_file = "../data/mesh.msh" 
    #     character(len=1024), parameter :: out_file = "../out/output.dat" 
    #     integer :: out_fid,ux_fid, uy_fid, uz_fid, phix_fid, phiy_fid, phiz_fid

    #     ! ## Mesh Variables
    #     type(mesh_type) :: mesh

    #     ! ## Material Variables
    #     type(material_table_type), dimension(:), allocatable :: material

    #     ! ## Problem variables
    #     type(problem_type) :: problem

    #     ! ## Solver variables
    #     type(solver_var_type) :: solver_var

    #     ! ## Auxiliary variables
    #     integer :: i, itersolve, infosolve
    #     integer, dimension(:), allocatable :: IPIV
    #     complex(kind=dp), dimension(:), allocatable :: WORK          
    #     complex(kind=dp), dimension(:), allocatable :: SWORK
    #     real(kind=dp), dimension(:), allocatable :: RWORK
=#
    inp_file = "mesh.msh"
#     ! =================================================
#     ! =================================================
#     ! ===================== Input =====================
#     ! =================================================
#     ! =================================================

    # ! ## Read from input
    mesh, material, problem, solver_var = readmsh(inp_file)

    # ! ## Allocate data as needed 
    # !    and calculate material, problem and solver variables:
    # call derive_data (mesh, material, problem, solver_var)

    # ! Generate physical mesh:
    # call generate_mesh(mesh)

    # ! Create outputs
    # call create_output (out_file, mesh, material(1), solver_var%nGP, problem%nFr,problem%fr_range, out_fid,&
    #     ux_fid, uy_fid, uz_fid, phix_fid, phiy_fid, phiz_fid)
    
    # ! ================================================
    # ! ================================================
    # ! ==================== Solver ====================
    # ! ================================================
    # ! ================================================


    # do i = 1, size(problem%frequencies)
        # problem%frequency = problem%frequencies(i)
        # write(*,"(A, X, I3, 2X, A, F10.4)") 'rodando para frequencia ',i , 'Frequencia: ', problem%frequency
        
        # ! ## Calculate G and H
        # call calc_GH(mesh, material, solver_var, problem)

        # ! ## Apply BC, arranging Ax = b
        # call apply_bc (mesh, solver_var)

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

