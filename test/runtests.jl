using Revise
using JuBEM
using Test
using LinearAlgebra

function teste_bar_io(inp_file)

    file_out = "teste"
    solve(inp_file;file_out=file_out)
    mesh,material,problem,solver_var = readvars_out(file_out)    
    u,t = getfreqres_out(file_out,0)

    rm(string(file_out,".h5"))

    ut = u[mesh.IEN[:,mesh.bc[:,1].==2][:],:]
    neu = sum(mesh.bc[:,1].==1)
    nnel = (mesh.eltype+1)^2
    ud = sum(ut[1:nnel*neu,1])/(nnel*neu)

    return  ud

end

function teste_bar_static(inp_file, t=0)


    mesh, material, problem, solver_var = read_msh(inp_file)
    mesh.eltype = t
    derive_data!(material, problem, solver_var)
    generate_mesh!(mesh)
    calc_GH!(mesh, material, solver_var,-1.0)
    mesh, solver_var, C = JuBEM.applyBC(mesh, solver_var,solver_var.H,solver_var.G)
    solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
    u,t,urb = JuBEM.returnut(mesh,solver_var.zvetsol, C)

    ut = u[mesh.IEN[:,mesh.bc[:,1].==2][:],:]
    neu = sum(mesh.bc[:,1].==1)
    nnel = (mesh.eltype+1)^2
    ud = sum(ut[1:nnel*neu,1])/(nnel*neu)

    return  ud

end

function teste_bar_static2(inp_file, t=0)

    mesh, material, problem, solver_var = read_msh(inp_file)
    mesh.eltype = t
    derive_data!(material, problem, solver_var)
    generate_mesh!(mesh)
    JuBEM.calc_GH_static!(mesh, material, solver_var)
    mesh, solver_var, C = JuBEM.applyBC(mesh, solver_var,solver_var.H,solver_var.G)
    solver_var.zvetsol = solver_var.ma \ mesh.zbcvalue
    u,t,urb = JuBEM.returnut(mesh,solver_var.zvetsol, C)

    ut = u[mesh.IEN[:,mesh.bc[:,1].==2][:],:]
    neu = sum(mesh.bc[:,1].==1)
    nnel = (mesh.eltype+1)^2
    ud = sum(ut[1:nnel*neu,1])/(nnel*neu)

    return  ud

end

function teste_soil_dyn(inp_file)


    solve_flex_dyn(inp_file;file_out = "test_output")
    N, freqs = getflex_out("test_output")
    rm("test_output.h5", force=true)


    return  N, freqs

end

function teste_soil_dyn2(inp_file)


    JuBEM.solve_flex_dyn2(inp_file;file_out = "test_output")
    N, freqs = getflex_out("test_output")
    rm("test_output.h5", force=true)


    return  N, freqs

end


@testset "IO" begin
    inp_file = "../meshes/static/bars/bar_2_3.msh"
    ud = teste_bar_io(inp_file)

    @test ud > 0.5

end

@testset "bar_static" begin

    inp_file = "../meshes/static/bars/bar_2_3.msh"
    ud_const = teste_bar_static(inp_file, 0)
    ud_lin = teste_bar_static(inp_file, 1)

    @test abs(ud_const-1) < 0.5
    @test abs(ud_lin-1) < 0.02

end

# @testset "soil_dyn" begin
begin

    inp_file = "./meshes/dynamic/soils/soilEE_109_rb.msh"
    ref_file = "./refres/test_output"

    N2,freqs = teste_soil_dyn2(inp_file)
    N,freqs = getflex_out(ref_file)

    @test norm(N - N2,1) < 1e-7
end