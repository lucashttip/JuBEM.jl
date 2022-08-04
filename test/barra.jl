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

    inp_file = "meshes/bars/bar_32.msh"
     # ref = -3.2
     E = 10
     L = 10
     q = 1
     l = 1
     A = l^2
     I = l^4/12
     ref = (q*A)*L/(E*A)

    mesh, material, problem, solver_var, u, t = solvestatic(inp_file)

    # Post-processing
    ut = u[mesh.IEN[:,mesh.bc[:,1].==2][:],:]

    neu = sum(mesh.bc[:,1].==1)

    nnel = (mesh.eltype+1)^2

    # ud = sum(ut[1:4*neu,1])/(4*neu)
    ud = sum(ut[1:nnel*neu,1])/(nnel*neu)

    erro = ((ud-ref)/ref)*100

    up,tp = calc_utpoints(mesh,u,t)

    writevtk(mesh,up,tp,"vis")
    erro, ud

    # maximum(ut, dims=1)
