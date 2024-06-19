using JuBEM

mesh_file = "./input/meshes/bar_1_10.msh"
problem_file = "./input/problems/bar_dynamic.prob"
output_file = "test_out.h5"


# function solve(mesh,problem,materials)

#     output_vars_h5(file_out, mesh, problem, material)

#     static_ass = statics_assembly(mesh,problem,materials)

#     remove_ee!(mesh,problem,static_ass)

#     for freq in problem.frequencies

#         assembly = dynamic_assembly(mesh,problem,materials,static_ass,freq)
#         solution = solve(assembly)


#         output_freq_h5(file_out,freq,solution)

#     end

# end

##

mesh = read_msh(mesh_file)
problem, materials = read_problem(problem_file,mesh)

generate_mesh!(mesh)
assembly = derive_data!(materials, problem)

output_vars(output_file, mesh, problem, materials)

assembly = statics_assembly(mesh,materials)

for freq in problem.frequencies
    println("Running for freq $freq")
    dynamics_assembly!(mesh,problem,materials,assembly,freq)
    LHS, RHS = JuBEM.applyBC_simple(mesh::Mesh,problem::Problem,assembly::Assembly,assembly.zH,assembly.zG)
    x = LHS\RHS
    u,t = JuBEM.returnut_simple(x,mesh,problem)
    sol = Solution(u,t,ComplexF64[],0.0,freq)
    output_solution(output_file,sol)
end

## Pos-processing

e = findfirst(mesh.tag.==2)
n = mesh.IEN[1,e]

u,t,freq = getnoderes_out(output_file,n)

##
L = 10
G = 10
nu = 0.0
rho = 1.0
E = 2*G*(1+nu)
c = sqrt(E/rho)
nwn = 3
wn = [(2n-1)*pi*c/(2*L) for n in 1:nwn]

