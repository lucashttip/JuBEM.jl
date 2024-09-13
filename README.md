# JuBEM.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://lucashttip.github.io/JuBEM.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lucashttip.github.io/JuBEM.jl/dev/)


JuBEM.jl is a Boundary Element Method package written in Julia. It currently supports 3D elastostatics and elastodynamics problems with quadrilateral linear discontinuous elements.
Infinite domains are supported.

Rigid-body boundary condition is also supported.

JuBEM has been implemented with the purpose of solving Soil-Structure Interaction Problems.

Supported meshes are generated with [Gmsh](https://gmsh.info/) and created using the [JuBEMeshes.jl](https://github.com/lucashttip/JuBEMeshes.jl) package.

To-dos:

- [x] Update/Improve Type definitions
- [x] Update input to be separate files for mesh and problem
- [x] Update/Clean functions
    - [x] Update shape functions
    - [x] Update generate_mesh! and derive_data!
    - [x] Create statics_assembly from solver_statics
    - [x] Verify statics_assembly
    - [x] Create new apply BC_simple
    - [x] Solver problems and verify statics
    - [x] Update output functions
    - [x] Create new plotting function (lighter)
    - [x] Create dynamics_assembly
    - [x] Update remove_ee
    - [x] Update apply BC rb
    - [x] Verify apply BC rb
    - [x] Juntar integratiorules e integrationconstants
    - [x] Remover integrationstatics e integrationdynamics?  
    - [x] Rethink/Remove solver, solverstatics, solverdynamics
- [x] Add multiregions
    - [x] Update generatemesh with needed 
        - [x] Include EEN to link geometric and physical elements
    - [x] Update assemblies to be over bodies
    - [x] Make applybc_multi_simple and returnut_multi_simple
    - [x] Update remove_EE!
    - [x] Make applybc_rb_multi and returnut_rb_multi
- [ ] Improve multiregion to generate banded matrices
- [ ] Improve API for running code/make it more concise
- [ ] Add documentations
- [ ] Eventually improve package to contain CI, Docs and Codecov
- [ ] Improve distance calculation
- [ ] Improve integrations
- [ ] Revisit and report element subdivision