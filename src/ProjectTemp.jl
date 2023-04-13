module ProjectTemp2

include("/Users/tobias/.julia/dev/ProjectTemp2/src/Functionalized and sped up code.jl")

#For loading in data
#For PC
    path="C:\\Users\\tehul\\Downloads"
    SAFD_output = CSV.read("C:\\Users\\tehul\\OneDrive - UvA\\Julia Research project\\PestMix1-8_1000ug-L_Tea_1-10dil_1ul_AllIon_pos_18_report.csv", DataFrame)
    SAFD_output_100ppb = CSV.read("C:\\Users\\tehul\\OneDrive - UvA\\Julia Research project\\210705_ToF_JO_004-100ppb Drugs_report.csv", DataFrame)
    SAFD_output_Stef = CSV.read("C:\\Users\\tehul\\OneDrive - UvA\\Julia Research project\\Stef paper data_report.csv", DataFrame)
#For MAC
    path = "/Users/tobias/Downloads"    
    SAFD_output = CSV.read("/Users/tobias/Downloads/PestMix1-8_1000ug-L_Tea_1-10dil_1ul_AllIon_pos_18_report.csv", DataFrame)
    SAFD_output_100ppb = CSV.read("/Users/tobias/Downloads/210705_ToF_JO_004-100ppb Drugs_report.csv", DataFrame)
    SAFD_output_Stef = CSV.read("/Users/tobias/Downloads/Stef paper data_report.csv", DataFrame)
    filenames = ["pes_mix9_m1.mzXML"]
#Set filename
    filenames = ["Stef paper data.mzXML"]
    filenames = ["PestMix1-8_1000ug-L_Tea_1-10dil_1ul_AllIon_pos_18.mzXML"]
    filenames = ["210705_ToF_JO_004-100ppb Drugs.mzXML"]
    mz_thresh = [0, 0]

#Import MS data
    GC.gc()

    @time mz_val, mz_int, t0, t_end, m, path, msModel, msIonisation, msManufacturer,
    polarity, Rt = import_files_MS1(path, filenames, mz_thresh)

FileName = m[1]
#Run SAFD 3D
max_numb_iter = 50
max_t_peak_w=300 # or 20
res=20000
min_ms_w=0.02
r_thresh=0.9
min_int=2000
sig_inc_thresh=5
S2N=3

min_peak_w_s=3
    GC.gc()
    @time rep_table, SAFD_output = safd_s3D(mz_val, mz_int, Rt, FileName, path, max_numb_iter,
    max_t_peak_w, res, min_ms_w, r_thresh, min_int, sig_inc_thresh, S2N, min_peak_w_s)

#Componentization of features
    # Parameters for CompCreate
    pathin = "/Users/tobias/Downloads" 
    path2features= "/Users/tobias/Downloads/PestMix1-8_1000ug-L_Tea_1-10dil_1ul_AllIon_pos_18_report.csv"
    path2features= "/Users/tobias/Downloads/pes_mix9_m1_report.csv"
    mass_win_per=0.8
    ret_win_per=0.5
    r_thresh=0.9
    delta_mass=0.004
    min_int=300 # Not needed for comp_ms1()
    
    
    chrom=import_files(pathin,filenames,mz_thresh)
    
    
    ## For only MS1
    SAFD_output = comp_ms1(chrom,path2features,mass_win_per,ret_win_per,r_thresh,delta_mass, min_int)
    
#Align masses and run resolutions algorithm
    unique_mz_values, plot_matrix = mass_align(Rt, mz_val, mz_int)
    gradient_pest = CSV.read("/Users/tobias/Library/CloudStorage/OneDrive-UvA/Gradient pesticides.csv", DataFrame)
    results, colors, df_1, gradient = @time unresolved_per_window_Rt_and_MS(Rt, SAFD_output, 12, 1.5, gradient_pest)
    results
    colors[37] = 4
    df_1
    gradient
#Plot the heatmap with features and windows
    plot_heatmap(SAFD_output, Rt, unique_mz_values, plot_matrix, 12, gradient)

#Figuring out the optimal number of windows to use
    unresolved_peaks = Vector{Int32}(zeros(150))
    for i = 1:150
        result, colors, df = unresolved_per_window_Rt_and_MS(Rt, SAFD_output, i, 1.5)
        unresolved_peaks[i] = sum(result[:,3])

    end
    unresolved_peaks
    windows = collect(0:5:150)
    windows_all = collect(1:150)
    ratio = unresolved_peaks./windows_all
    scatter(ratio, size = (1280,720), xticks = windows, legend = false)
    plot!(ratio, size = (1280,720), xticks = windows, legend = false)

    scatter(unresolved_peaks, size = (1920,1080), title = "Optimal nr of windows-Drug data", xlabel = "Nr of windows", ylabel = "Nr of unresolved peaks",left_margin=15Plots.mm, top_margin=7.5Plots.mm,
    bottom_margin=8.5Plots.mm, grid = false, legend = false, xticks = windows)
    plot!(unresolved_peaks)


    savefig("/Users/tobias/Library/CloudStorage/OneDrive-UvA/Julia Research project/Figures/Gradient/Heatmap with gradient overlay.png")
SAFD_output
    export mass_align, resolutions, unresolved_per_window_Rt_and_MS, window_split_Rt, Peaks_p_window,
        plot_heatmap
end
