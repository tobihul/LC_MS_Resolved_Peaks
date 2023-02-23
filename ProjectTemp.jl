module ProjectTemp

include("functionalized and sped up code.jl")

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
    
filenames = ["Stef paper data.mzXML"]
filenames = ["PestMix1-8_1000ug-L_Tea_1-10dil_1ul_AllIon_pos_18.mzXML"]
filenames = ["210705_ToF_JO_004-100ppb Drugs.mzXML"]
mz_thresh = [0, 1200]

#Import MS data
GC.gc()

@time mz_val, mz_int, t0, t_end, m, path, msModel, msIonisation, msManufacturer,
polarity, Rt = import_files_MS1(path, filenames, mz_thresh)

FileName = m[1]
#Run SAFD 3D
    GC.gc()
    @time rep_table, SAFD_output_Stef = safd_s3D(mz_val, mz_int, Rt, FileName, path, max_numb_iter,
    max_t_peak_w, res, min_ms_w, r_thresh, min_int, sig_inc_thresh, S2N, min_peak_w_s)


#Align masses and run resolutions algorithm
    unique_mz_values, plot_matrix = mass_align(Rt, mz_val, mz_int)
    results, colors, df_1 = @time unresolved_per_window_Rt_and_MS(Rt, SAFD_output, 12, 1.5)
    results
#Plot the heatmap with features and windows
    plot_heatmap(SAFD_output, 1200, Rt, unique_mz_values, plot_matrix, 12)


#Calculate the coverage of the measurement
    coverage, p = calculate_percentage_coverage(plot_matrix,10000,23.5,900)
    coverage
    p



unresolved_peaks = Vector{Int32}(zeros(150))
for i = 1:150
    result, colors, df = unresolved_per_window_Rt_and_MS_new(Rt, SAFD_output, i, 1.5)
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



M = rand(Int(trunc(rand() * 99 + 1)), Int(trunc(rand() * 99 + 1)))







function mat_split(M::Matrix{Float64})

    M = plot_matrix'
    m,n = size(M)

    sub_matrix_row = Int(trunc(m/3))    
    sub_matrix_col = Int(trunc(n/3))

    A1 = M[1:Int((sub_matrix_row)*1),1:Int((sub_matrix_col)*1)]
    A2 = M[Int((sub_matrix_row)*1):Int((sub_matrix_row)*2), 1:Int((sub_matrix_col)*1)]
    A3 = M[Int((sub_matrix_row)*2):Int((sub_matrix_row)*3), 1:Int((sub_matrix_col)*1)]

    B1 = M[1:Int((sub_matrix_row)*1), Int((sub_matrix_col)*1):Int((sub_matrix_col)*2)]
    B2 = M[Int((sub_matrix_row)*1):Int((sub_matrix_row)*2), Int((sub_matrix_col)*1):Int((sub_matrix_col)*2)]
    B3 = M[Int((sub_matrix_row)*2):Int((sub_matrix_row)*3), Int((sub_matrix_col)*1):Int((sub_matrix_col)*2)]

    C1 = M[1:Int((sub_matrix_row)*1), Int((sub_matrix_col)*2):Int((sub_matrix_col)*3)]
    C2 = M[Int((sub_matrix_row)*1):Int((sub_matrix_row)*2), Int((sub_matrix_col)*2):Int((sub_matrix_col)*3)]
    C3 = M[Int((sub_matrix_row)*2):Int((sub_matrix_row)*3), Int((sub_matrix_col)*2):Int((sub_matrix_col)*3)]

    #For plotting

    h_line_1 = ones(Int(ceil(Rt[end]))).* unique_mz_values[sub_matrix_row*1]
    h_line_2 = ones(Int(ceil(Rt[end]))).* unique_mz_values[sub_matrix_row*2]
    v_line_1 = ones(Int(ceil(unique_mz_values[end]))).* Rt[sub_matrix_col*1]
    v_line_2 = ones(Int(ceil(unique_mz_values[end]))).* Rt[sub_matrix_col*2]
    return A1, A2, A3, B1, B2, B3, C1, C2, C3, h_line_1, h_line_2, v_line_1, v_line_2
end
Rt
A1, A2, A3, B1, B2, B3, C1, C2, C3, h_line_1, h_line_2, v_line_1, v_line_2  = mat_split(plot_matrix)

heatmap(Rt, unique_mz_values, plot_matrix',
color=:plasma,
        clims=(25000, 80000),
        #ylims = (50,600),
        size=(1280, 720),
        xlabel="Rt",
        ylabel="m/z",
        title="Heat map of pest mix",
        left_margin=5Plots.mm, right_margin=7.5Plots.mm,
        bottom_margin=8.5Plots.mm)

plot!(h_line_1, c =:red)
plot!(h_line_2, c =:red)
plot!(v_line_1, c =:red)
plot!(v_line_2, c =:red)



round(50.5)

plot!(h_line_1, c =:red)





A1
A2
A3
B1
B2
B3
C1
C2
C3

size(plot_matrix)





#Better implementation??
function mat_split_new(M::Matrix{Float64})
    m,n = size(M)

    sub_matrix_row = Int(ceil(m/3))
    sub_matrix_col = Int(ceil(n/3))

    A1 = M[1:min(sub_matrix_row, m), 1:min(sub_matrix_col, n)]
    A2 = M[(sub_matrix_row+1):min(2*sub_matrix_row, m), 1:min(sub_matrix_col, n)]
    A3 = M[(2*sub_matrix_row+1):min(3*sub_matrix_row, m), 1:min(sub_matrix_col, n)]

    B1 = M[1:min(sub_matrix_row, m), (sub_matrix_col+1):min(2*sub_matrix_col, n)]
    B2 = M[(sub_matrix_row+1):min(2*sub_matrix_row, m), (sub_matrix_col+1):min(2*sub_matrix_col, n)]
    B3 = M[(2*sub_matrix_row+1):min(3*sub_matrix_row, m), (sub_matrix_col+1):min(2*sub_matrix_col, n)]

    C1 = M[1:min(sub_matrix_row, m), (2*sub_matrix_col+1):min(3*sub_matrix_col, n)]
    C2 = M[(sub_matrix_row+1):min(2*sub_matrix_row, m), (2*sub_matrix_col+1):min(3*sub_matrix_col, n)]
    C3 = M[(2*sub_matrix_row+1):min(3*sub_matrix_row, m), (2*sub_matrix_col+1):min(3*sub_matrix_col, n)]

    return A1, A2, A3, B1, B2, B3, C1, C2, C3
end

A1, A2, A3, B1, B2, B3, C1, C2, C3 = mat_split_new(plot_matrix)

    export mass_align, resolutions_new, unresolved_per_window_Rt_and_MS_new, window_split_Rt, Peaks_p_window_new,
        plot_heatmap, calculate_percentage_coverage
end

