__precompile__()
module LC_MS_Resolutions

using Statistics, SAFD, CSV, DataFrames, LoopVectorization, StatsPlots, Distributions
using MS_Import, LinearAlgebra, CompCreate, JLD2, VoronoiCells, GeometryBasics, Random
include("Functionalized and sped up code.jl")

    export mass_align, resolutions, unresolved_per_window_Rt_and_MS, window_split_Rt, Peaks_p_window,
        plot_heatmap, gradient_curve, normalize_lc_ms, surface_voronoi, load_and_prep_data
end




