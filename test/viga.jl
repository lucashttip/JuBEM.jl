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

    inp_file = "meshes/vigas/viga_2_20.msh"
    # ref = -3.2
    G = 4e3
    v = 0.25
    E = 2*G*(1+v)
    L = 10
    q = 1
    l = 1
    I = l^4/12
    A = l^2
    ref = -(q*l^4)/(24*E*I) # Carga distribuída ao longo de x
    # ref = -(q*l^3)/(3*E*I) # Carga concentrada em x = L
    # ref2 = (q*A)*L/(E*A)

    mesh, material, problem, solver_var, u, t = solvestatic(inp_file)

    # applyBC_nonrb!(mesh, solver_var)

    # Post-processing
    ut = u[mesh.IEN[:,mesh.bc[:,2].==2][:],:]

    neu = sum(mesh.bc[:,2].==1)

    ud = sum(ut[1:4*neu,2])./(4*neu)
    erro = ((ud-ref)/ref)*100

    up,tp = calc_utpoints(mesh,u,t)

    writevtk(mesh,up,tp,"vis")
    
    erro, ud

#=
    # ref = 1
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
        # applyBC_nonrb2!(mesh, solver_var)


        x = solver_var.ma \ mesh.zbcvalue

        # u,t = returnut2(mesh,x)
        u,t = returnut(mesh,x)
       
    
        

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
 =#
    