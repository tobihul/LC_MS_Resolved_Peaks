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
    GC.gc()
    @time rep_table, SAFD_output_Stef = safd_s3D(mz_val, mz_int, Rt, FileName, path, max_numb_iter,
    max_t_peak_w, res, min_ms_w, r_thresh, min_int, sig_inc_thresh, S2N, min_peak_w_s)

#Componentization of features
    # Parameters for CompCreate
    pathin = "/Users/tobias/Downloads" 
    path2features= "/Users/tobias/Downloads/PestMix1-8_1000ug-L_Tea_1-10dil_1ul_AllIon_pos_18_report.csv"
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
    results, colors, df_1 = @time unresolved_per_window_Rt_and_MS(Rt, SAFD_output, 12, 1.5)
    results
    colors
    df_1
#Plot the heatmap with features and windows
    plot_heatmap(SAFD_output, Rt, unique_mz_values, plot_matrix, 12)


#Calculate the coverage of the measurement
    coverage, p = calculate_percentage_coverage(plot_matrix,10000,23.5,900)
    coverage
    p
#Split the data into a 9 part grid
    A1, A2, A3, B1, B2, B3, C1, C2, C3, h_line_1, h_line_2, v_line_1, v_line_2, p  = mat_split(plot_matrix, 900, 23.6)
    p

#Calculate average coverage and standard deviation
    mean_c, std_c = @time calc_coverage_grid(A1, A2, A3, B1, B2, B3, C1, C2, C3, 10)
    mean_c
    std_c



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

#Now with percentage B mobile phase 

    Rt
    L_step_1 = Vector{Float32}(ones(length(Rt[1:findfirst(x->x >=1,Rt)])).*5)
    
    G_step_1 = collect(5:(90/(length(Rt[(findfirst(x->x >=1,Rt))+1:findfirst(x->x >=21,Rt)-1]))):95)

    L_step_2 = Vector{Float32}(ones(length(Rt[findfirst(x->x >=21,Rt)+1:findfirst(x->x >=23,Rt)])).*95)

    G_step_2 = collect(95:-(90/(length(Rt[(findfirst(x->x >=23,Rt))+1:findfirst(x->x >=23.1,Rt)-1]))):5)

    L_step_3 = Vector{Float32}(ones(length(Rt[findfirst(x->x >=23.1,Rt)+1:end])).*5)

    B_gradient_final = [L_step_1; G_step_1; L_step_2; G_step_2; L_step_3]
    
    plot(Rt, B_gradient_final)

#Also delta % B Now per window
    Rt
    B_gradient_final
    split
    delta = Vector{Float64}(zeros((length(split))-1))
    # first define pos at the first split for the gradient 
    pos = findfirst(x->x>=split[1], Rt)
    @time for i = 1:length(split)-1
        # save the index of the i + 1 location in temp to not overwrite the pos variable
        tmp = findfirst(x->x>=split[i+1], Rt)
        delta[i] = (B_gradient_final[tmp] - B_gradient_final[pos])/(split[i+1]-split[i])
        # overwrite pos since the index in tmp is equal to the next position that we want to use
        pos = tmp
    end
    delta
    
    
    


filtered_features = comp_list

split = window_split_Rt(Rt, 12)
    heatmap(Rt, unique_mz_values, plot_matrix',
        #c = cgrad([:white,:navy,:indigo,:teal,:green,:yellow,:red],[0,0.04,1]),
        color=:plasma,
        clims=(25000, 80000),
        #ylims = (50,600),
        size=(1280, 720),
        xlabel="Rt",
        ylabel="m/z",
        title="Heat map of pest mix",
        left_margin=5Plots.mm, right_margin=7.5Plots.mm,
        bottom_margin=8.5Plots.mm,
        colorbar = false,
        xticks = (round.(split; digits = 1)),
        yticks = (0:100:1200)

    )
    # Create a scatter plot using the x and y coordinates and the colors and symbols vectors
    paint = paint = [:Red, :Orange, :Green, :Yellow]
    colors_final = paint[colors]
    mapping = Dict(1 => "Unresolved in RT and MS", 2 => "Resolved in MS only", 3 => "Resolved in RT only", 4 => "Fully resolved")
    labels_final = map(x -> mapping[x], colors)
    p2 = scatter!(SAFD_output[:, 4], SAFD_output[:, 8],
        #series_annotations = text.(1:length(SAFD_output[:,4]),size = 1),
        markershape=:xcross,
        color=colors_final,
        group=labels_final,
        legend=true,
        markersize=2.5,
        title="Pest mix, 2000 iterations, S/N = 3, r = 0.9, accepted_res = 1.5, With componetization (866 features)",
        left_margin=5Plots.mm, right_margin=7.5Plots.mm,
        bottom_margin=8.5Plots.mm,
    )

    for i = 1:length(split)
        @show i
        p2 = plot!(ones(Int32(ceil(maximum(unique_mz_values)))) .* (split[i]), collect(1:1:Int32(ceil(maximum(unique_mz_values)))), color=:red, legend=true, label=false)
        display(p2)
    end
    plot!(twinx(), Rt, B_gradient_final, yticks = (0:5:100), label = ("gradient"), ylabel = ("%B"), linewidth = 5, linestyle = :dot)

    savefig("/Users/tobias/Library/CloudStorage/OneDrive-UvA/Julia Research project/Figures/Gradient/Heatmap with gradient overlay.png")

    


    

    export mass_align, resolutions, unresolved_per_window_Rt_and_MS, window_split_Rt, Peaks_p_window,
        plot_heatmap, calculate_percentage_coverage
end
