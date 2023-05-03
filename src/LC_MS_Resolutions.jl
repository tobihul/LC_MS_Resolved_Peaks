__precompile__()
module LC_MS_Resolutions

#Packages
using Statistics 
using SAFD
using CSV 
using DataFrames 
using LoopVectorization
using StatsPlots
using Distributions
using MS_Import
using LinearAlgebra 
using CompCreate
using JLD2
using VoronoiCells
using GeometryBasics
using Random

try 
    using MS_Import
catch
    using Pkg
    Pkg.add(PackageSpec(url="https://bitbucket.org/SSamanipour/ms_import.jl/src/master/"))
    using MS_Import
end
try 
        using CompCreate
    catch
        using Pkg
        Pkg.add(PackageSpec(url="https://bitbucket.org/SSamanipour/compcreate.jl/src/master/"))
        using CompCreate
end
    try 
        using SAFD
    catch
        using Pkg
        Pkg.add(PackageSpec(url="https://bitbucket.org/SSamanipour/safd.jl/src/master/"))
        using SAFD
end

include("Functionalized and sped up code.jl")

export mass_align, resolutions, unresolved_per_window_Rt_and_MS, window_split_Rt, Peaks_p_window,
        plot_heatmap, gradient_curve, normalize_lc_ms, surface_voronoi, load_and_prep_data
end




