
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

#python
#setup conda
if Sys.iswindows()
        #move miniconda installer in case it is in the Conda/3/ folder
        pathConda = pathof(Conda)
        if contains(pathConda, "\\")
            pathConda = split(pathConda, "\\packages")[1]
        end
        pathConda = joinpath(pathConda, "conda", "3","installer.exe")
    
    
    
        if isfile(pathConda)
            mv(pathConda, pathConda[1:end-15]*pathConda[end-12:end], force = true)
        else
            Downloads.download("https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe",pathConda)
        end
    end
    
    #
    ENV["PYTHON"]=""
    Conda.pip_interop(true, Conda.ROOTENV)
    using PyCall
    println("Installing scipy")
    if Sys.isapple()
        Conda.add("scipy",Conda.ROOTENV)
    else
        Conda.pip("install","scipy",Conda.ROOTENV)
    end
    #rebuild PyCall
    # Pkg.build("PyCall")
    const PyStat = PyNULL()
    function __init__()
        copy!(PyStat,pyimport("scipy.stats"))
    end
include("Functionalized and sped up code.jl")

export mass_align, resolutions, unresolved_per_window_Rt_and_MS, window_split_Rt, Peaks_p_window,
        plot_heatmap, gradient_curve, normalize_lc_ms, surface_voronoi, load_and_prep_data, Resolved_peaks_algorithm
end




