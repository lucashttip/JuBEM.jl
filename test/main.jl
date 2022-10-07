using Revise
using JuBEM
using Plots

inp_file = "meshes/dynamic/soils/SoBEE.msh"

file_out = "SoBEE_lin_Fz"

solve(inp_file;file_out=file_out)