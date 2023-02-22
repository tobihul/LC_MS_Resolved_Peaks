module ProjectTemp

include("functionalized and sped up code.jl")
# Write your package code here.
 
export mass_align, resolutions_new, unresolved_per_window_Rt_and_MS_new, window_split_Rt, Peaks_p_window_new,
       plot_heatmap, calculate_percentage_coverage
end

j
max_numb_iter = 500
max_t_peak_w = 300 # or 20
res = 20000
min_ms_w = 0.02
r_thresh = 0.9
min_int = 2000
sig_inc_thresh = 5
S2N = 3
min_peak_w_s = 3


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

GC.gc()

@time mz_val, mz_int, t0, t_end, m, path, msModel, msIonisation, msManufacturer,
polarity, Rt = import_files_MS1(path, filenames, mz_thresh)

FileName = m[1]

GC.gc()

#@time rep_table,final_table=safd(mz_vals,mz_int,t0,t_end,FileName,path,max_numb_iter,
#    max_t_peak_w,res,min_ms_w,r_tresh,min_int,sig_inc_thresh,S2N,min_peak_w_s)

@time rep_table, SAFD_output_Stef = safd_s3D(mz_val, mz_int, Rt, FileName, path, max_numb_iter,
    max_t_peak_w, res, min_ms_w, r_thresh, min_int, sig_inc_thresh, S2N, min_peak_w_s)

@time unique_mz_values, plot_matrix = mass_align(Rt, mz_val, mz_int)

unique_mz_values, plot_matrix = mass_align(Rt, mz_val, mz_int)
results, colors, df_1 = @time unresolved_per_window_Rt_and_MS_new(Rt, SAFD_output, 10, 1.5)
results

plot_heatmap(SAFD_output_Stef,Rt,unique_mz_values,plot_matrix,10)



results, colors, df_1 = @time unresolved_per_window_Rt_and_MS_new(Rt, SAFD_output, 10, 1.5)
results
colors

results_old, colors_old, df = @time unresolved_per_window_Rt_and_MS(Rt, SAFD_output, 10, 1.5)
results_old
colors_old
colors == colors_old
results == results_old

coverage, p = calculate_percentage_coverage(plot_matrix,10000,23.5,900)
coverage
p
@code_warntype mass_align(Rt, mz_val, mz_int_1)
@code_warntype unresolved_per_window_Rt_and_MS_new(Rt, SAFD_output, 10, 1.5)
@code_warntype Peaks_p_window(10, Rt, SAFD_output)

@time resolutions_new(SAFD_output)
@time resolutions(SAFD_output)

unresolved_peaks = Vector{Int32}(zeros(150))
for i = 1:150
    result, colors, df = unresolved_per_window_Rt_and_MS_new(Rt, SAFD_output, i, 1.5)
    unresolved_peaks[i] = sum(result[:,3])
    @show i
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


function simple_moving_average(vec::Array{Float64,1}, window::Int)
    n = length(vec)
    avg = similar(vec)
    for i in window:n
        avg[i] = sum(vec[i-window+1:i]) / window
    end
    return avg
end

tic = vec(sum(mz_int, dims=2))
sma = simple_moving_average(tic,3)
p_tic = plot(Rt, sma, size=(1920, 1080), xlabel="Rt",
        ylabel="Intensity", left_margin=15Plots.mm, right_margin=7.5Plots.mm,
        bottom_margin=7.5Plots.mm, title="TIC of pest mix data", legend=false, c=:black)