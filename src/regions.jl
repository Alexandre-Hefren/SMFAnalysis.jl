intervals(table; start=:start, stop=:stop, off = 0) = IntervalCollection(Interval.(table.chrom, table[!, start] .+ off, table[!, stop], '.', 1:size(table, 1)))


function markintersection!(tableA, tableB, label=:Inter)
    ivA = intervals(tableA)::IntervalCollection{Int64} 
    ivB = intervals(tableB)::IntervalCollection{Int64} 
    tableA[!, label] = falses(size(tableA, 1))
    for (ia, ib) in eachoverlap(ivA, ivB)
        tableA[GenomicFeatures.metadata(ia), label] = true
    end
    tableA
end

function loadblacklist(file="hg19-blacklist.v2.bed.gz", projdir=getprojectdir())
    filepath = joinpath(projdir, "data", file)

    blacklist = CSV.read(filepath, DataFrame, header=[:chrom, :start, :stop, :name]);
    sort!(blacklist, [:chrom, :start, :stop]);
    blacklist
end

function load_tf_data(dir, blacklist=loadblacklist())
    files = glob("*.tsv", dir)
    tfs = getindex.(split.(basename.(files), r"[_.]"), 2)
    peakmotif = Dict(tf => @subset(CSV.read(f, DataFrame), :totalmotifs .> 0) for (tf, f) in zip(tfs, files))

    for tf in tfs
        markintersection!(peakmotif[tf], blacklist, :Blacklist);
        # markintersection!(peakmotif[tf], targets, :Target);
        peakmotif[tf] = @subset(peakmotif[tf], .!:Blacklist)
    end
    
    peakmotif
end