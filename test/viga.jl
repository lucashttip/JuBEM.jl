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

    inp_file = "meshes/static/vigas_reg/viga_4_40.msh"
    # ref = -3.2
    G = 5e3
    v = 0.0
    E = 2*G*(1+v)
    L = 10
    q = 1
    l = 1
    I = l^4/12
    A = l^2
    # ref = -(3*q*L^4)/(24*E*I) # Carga distribuída ao longo de x
    ref = -(q*L^3)/(3*E*I) # Carga concentrada em x = L
    # ref2 = (q*A)*L/(E*A)

    solve(inp_file)
    mesh,material,problem,solver_var = readvars_out("output")
    u,t = getfreqres_out("output",0)

    # applyBC_nonrb!(mesh, solver_var)

    # Post-processing
    ut = u[mesh.IEN[:,mesh.bc[:,2].==2][:],:]

    neu = sum(mesh.bc[:,2].==1)
    nnel = (mesh.eltype+1)^2
    ud = sum(ut[1:nnel*neu,2])./(nnel*neu)
    erro = ((ud-ref)/ref)*100

    up,tp = calc_utpoints(mesh,u,t)

    # writevtk(mesh,up,tp,"vis")
    
    erro, ud, ref
