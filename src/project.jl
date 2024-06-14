using DataFrames, DataFramesMeta, CSV
using Glob
using SMFTools, Blosc

using Statistics
using Plots, StatsPlots, Measures
theme(:wong2)

using GenomicFeatures, GenomeFragments

include("smf.jl")
include("regions.jl")

function getprojectdir()
    d = pwd()
    if basename(d) == "notebooks"
        return dirname(d)
    else
        return d
    end
end