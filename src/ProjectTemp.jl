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


function mat_split(M::Matrix{Float32}, max_mz::Int64, Gradient_end::Float64)
    #Trimming for max_mz and end_gradient
    matrix_grad = M[1:(findfirst(x->Gradient_end<x, Rt)),:]
    matrix_mz = matrix_grad[:,1:(findfirst(x->x>max_mz, unique_mz_values))]

    M_final = matrix_mz'
    m,n = size(M_final)

    sub_matrix_row = Int(trunc(m/3))    
    sub_matrix_col = Int(trunc(n/3))

    A1 = M_final[1:min(sub_matrix_row, m), 1:min(sub_matrix_col, n)]
    A2 = M_final[(sub_matrix_row+1):min(2*sub_matrix_row, m), 1:min(sub_matrix_col, n)]
    A3 = M_final[(2*sub_matrix_row+1):min(3*sub_matrix_row, m), 1:min(sub_matrix_col, n)]

    B1 = M_final[1:min(sub_matrix_row, m), (sub_matrix_col+1):min(2*sub_matrix_col, n)]
    B2 = M_final[(sub_matrix_row+1):min(2*sub_matrix_row, m), (sub_matrix_col+1):min(2*sub_matrix_col, n)]
    B3 = M_final[(2*sub_matrix_row+1):min(3*sub_matrix_row, m), (sub_matrix_col+1):min(2*sub_matrix_col, n)]

    C1 = M_final[1:min(sub_matrix_row, m), (2*sub_matrix_col+1):min(3*sub_matrix_col, n)]
    C2 = M_final[(sub_matrix_row+1):min(2*sub_matrix_row, m), (2*sub_matrix_col+1):min(3*sub_matrix_col, n)]
    C3 = M_final[(2*sub_matrix_row+1):min(3*sub_matrix_row, m), (2*sub_matrix_col+1):min(3*sub_matrix_col, n)]


    #For plotting

    h_line_1 = ones(n).* unique_mz_values[sub_matrix_row*1]
    h_line_2 = ones(n).* unique_mz_values[sub_matrix_row*2]
    v_line_1 = ones(m).* Rt[sub_matrix_col*1]
    v_line_2 = ones(m).* Rt[sub_matrix_col*2]

    heatmap(Rt[1:length(matrix_grad[:,1])], unique_mz_values[1:length(matrix_mz[1,:])], matrix_mz',
        color=:plasma,
        clims=(25000, 80000),
        legend = false,
        size=(1280, 720),
        xlabel="Rt",
        ylabel="m/z",
        title="Heat map of pest mix",
        left_margin=5Plots.mm, right_margin=7.5Plots.mm,
        bottom_margin=8.5Plots.mm)


    plot!(Rt[1:length(matrix_grad[:,1])], h_line_1, c =:red, linestyle = :dash)
    plot!(Rt[1:length(matrix_grad[:,1])], h_line_2, c =:red, linestyle = :dash)
    plot!(v_line_1, unique_mz_values[1:length(matrix_mz[1,:])], c =:red, linestyle = :dash)
    p_f = plot!(v_line_2, unique_mz_values[1:length(matrix_mz[1,:])], c =:red, linestyle = :dash)


   
    return A1, A2, A3, B1, B2, B3, C1, C2, C3, h_line_1, h_line_2, v_line_1, v_line_2, p_f
end

A1, A2, A3, B1, B2, B3, C1, C2, C3, h_line_1, h_line_2, v_line_1, v_line_2, p  = mat_split(plot_matrix, 900, 23.6)
p

A1
A2
A3
B1
B2
B3
C1
C2
C3

sub = sum(A1) + sum(A2) + sum(A3) + sum(B1) + sum(B2) + sum(B3) + sum(C1) + sum(C2) + sum(C3)
all = sum(matrix_mz)

function calc_coverage_grid(A1, A2, A3, B1, B2, B3, C1, C2, C3, threshold)

    A1_cov = calculate_percentage_coverage(A1, threshold)
    A2_cov = calculate_percentage_coverage(A2, threshold)
    A3_cov = calculate_percentage_coverage(A3, threshold)
    B1_cov = calculate_percentage_coverage(B1, threshold)
    B2_cov = calculate_percentage_coverage(B2, threshold)
    B3_cov = calculate_percentage_coverage(B3, threshold)
    C1_cov = calculate_percentage_coverage(C1, threshold)
    C2_cov = calculate_percentage_coverage(C2, threshold)
    C3_cov = calculate_percentage_coverage(C3, threshold)

    mean_cov = mean([A1_cov, A2_cov, A3_cov, B1_cov, B2_cov, B3_cov,
                     C1_cov, C2_cov, C3_cov])
    std_cov = std([A1_cov, A2_cov, A3_cov, B1_cov, B2_cov, B3_cov,
                   C1_cov, C2_cov, C3_cov])

    return mean_cov, std_cov
end

mean_c, std_c = @time calc_coverage_grid(A1, A2, A3, B1, B2, B3, C1, C2, C3, 100000)
mean_c
std_c





#Old Implementation
function mat_split_old(M::Matrix{Float64})
    m,n = size(M)

    sub_matrix_row = Int(trunc(m/3))
    sub_matrix_col = Int(trunc(n/3))

    A1 = M_final[1:Int((sub_matrix_row)*1),1:Int((sub_matrix_col)*1)]
    A2 = M_final[Int((sub_matrix_row)*1):Int((sub_matrix_row)*2), 1:Int((sub_matrix_col)*1)]
    A3 = M_final[Int((sub_matrix_row)*2):Int((sub_matrix_row)*3), 1:Int((sub_matrix_col)*1)]

    B1 = M_final[1:Int((sub_matrix_row)*1), Int((sub_matrix_col)*1):Int((sub_matrix_col)*2)]
    B2 = M_final[Int((sub_matrix_row)*1):Int((sub_matrix_row)*2), Int((sub_matrix_col)*1):Int((sub_matrix_col)*2)]
    B3 = M_final[Int((sub_matrix_row)*2):Int((sub_matrix_row)*3), Int((sub_matrix_col)*1):Int((sub_matrix_col)*2)]

    C1 = M_final[1:Int((sub_matrix_row)*1), Int((sub_matrix_col)*2):Int((sub_matrix_col)*3)]
    C2 = M_final[Int((sub_matrix_row)*1):Int((sub_matrix_row)*2), Int((sub_matrix_col)*2):Int((sub_matrix_col)*3)]
    C3 = M_final[Int((sub_matrix_row)*2):Int((sub_matrix_row)*3), Int((sub_matrix_col)*2):Int((sub_matrix_col)*3)]



    return A1, A2, A3, B1, B2, B3, C1, C2, C3
end

export mass_align, resolutions, unresolved_per_window_Rt_and_MS, window_split_Rt, Peaks_p_window,
        plot_heatmap, calculate_percentage_coverage
end
