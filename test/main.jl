using Revise
using JuBEM
using Plots

inp_file = "./meshes/dynamic/soils/soil_rb_6_4_fs_10_FxFz.msh"

# inp_file = "./meshes/static/bars/bar_2_3.msh"

file_out = "teste"

solve(inp_file;file_out=file_out)