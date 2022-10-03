using HDF5
function read_mat_h5s(filename)
    fid = h5open(filename,"r")
    
    # Real part
    Ar = read(fid,"real")

    # Imag part
    Ai = read(fid,"imag")

    return Ar .+ Ai.*im
    
end


H = read_mat_h5s("./test/H.h5")
G = read_mat_h5s("./test/G.h5")
zma = read_mat_h5s("./test/zma.h5")
zbcvalue = read_mat_h5s("./test/zbcvalue.h5")