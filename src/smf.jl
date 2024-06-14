

function load_smf_data(dir)
    files = glob("*.posbin", dir)
    @show files
    PM = load_pos_matrix.(files)
end


function get_mod_data(readindex, modification, PM)
    mi = PM.modifications[modification]
    off = (mi - 1)*2
    rind = off .+ 3 .+ (1:2)
    pind = PM.readindex[rind[1], readindex]:PM.readindex[rind[2], readindex]
    d = decompress(Int32, PM.PM[pind])
    reshape(d, 2, div(length(d), 2))
end