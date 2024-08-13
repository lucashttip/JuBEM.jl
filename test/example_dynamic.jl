using JuBEM

# mesh_file = "./input/meshes/bar_1_10.msh"
# mesh_file = "./input/meshes/soilEE_109.msh"
mesh_file = "./input/meshes/multibar_ne=1x(2+2)_l=(5+5)x1x1_eo=1.msh"
# problem_file = "./input/problems/bar_dynamic.prob"
# problem_file = "./input/problems/soilrb_EE_dynamic.prob"
problem_file = "./input/problems/multi_bar_dynamic.prob"

output_file = "test_out.h5"

mesh = read_msh(mesh_file)
problem, materials = read_problem(problem_file,mesh)


generate_mesh!(mesh)
derive_data!(mesh,materials,problem)

output_vars(output_file, mesh, problem, materials)

assembly = statics_assembly(mesh,materials,problem)
JuBEM.remove_EE!(mesh,assembly,problem)


for freq in problem.frequencies
    println("Running for freq $freq")
    dynamics_assembly!(mesh,problem,materials,assembly,freq)
    # LHS, RHS = JuBEM.applyBC_rb(mesh::Mesh,problem::Problem,assembly.zH,assembly.zG)
    LHS, RHS = JuBEM.applyBC_multi_simple(mesh::Mesh,problem::Problem,assembly.zH,assembly.zG)
    x = LHS\RHS
    # u,t,urb = JuBEM.returnut_rb(x,mesh,problem)
    u,t = JuBEM.returnut_multi_simple(x,mesh,problem)
    sol = Solution(u,t,freq)
    output_solution(output_file,sol)
end

## Pos-processing

e = findfirst(mesh.tag.==4)
n = mesh.IEN[1,e]

u,t,freq = getnoderes_out(output_file,n)

## Pos-processing flex

urb,freqs = geturb_out(output_file,3)


##
L = 10
G = 10
nu = 0.0
rho = 1.0
E = 2*G*(1+nu)
c = sqrt(E/rho)
nwn = 3
wn = [(2n-1)*pi*c/(2*L) for n in 1:nwn]

