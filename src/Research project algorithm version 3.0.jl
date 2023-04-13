using Statistics, SAFD, CSV, DataFrames, LoopVectorization, StatsPlots, Distributions
using MS_Import, Colors, LinearAlgebra

#Import data and run SAFD-3D
path = "/Users/tobias/Downloads"
filenames = ["Stef paper data.mzXML"]
filenames = ["PestMix1-8_1000ug-L_Tea_1-10dil_1ul_AllIon_pos_18.mzXML"]
filenames = ["210705_ToF_JO_004-100ppb Drugs.mzXML"]
mz_thresh = [0, 0]
mz_val
#Or just load some old data

SAFD_output = CSV.read("/Users/tobias/Downloads/PestMix1-8_1000ug-L_Tea_1-10dil_1ul_AllIon_pos_18_report.csv", DataFrame)
SAFD_output_100ppb = CSV.read("/Users/tobias/Downloads/210705_ToF_JO_004-100ppb Drugs_report.csv", DataFrame)
SAFD_output_Stef = CSV.read("/Users/tobias/Downloads/Stef paper data_report.csv", DataFrame)
# Feature detection parameters
max_numb_iter = 2000
max_t_peak_w = 300 # or 20
res = 20000
min_ms_w = 0.02
r_thresh = 0.9
min_int = 2000
sig_inc_thresh = 5
S2N = 3

min_peak_w_s = 3

GC.gc()

@time mz_vals, mz_int, t0, t_end, m, pathin, msModel, msIonisation, msManufacturer,
polarity, Rt = import_files_MS1(path, filenames, mz_thresh)

FileName = m[1]

GC.gc()

#@time rep_table,final_table=safd(mz_vals,mz_int,t0,t_end,FileName,path,max_numb_iter,
#    max_t_peak_w,res,min_ms_w,r_tresh,min_int,sig_inc_thresh,S2N,min_peak_w_s)

@time rep_table, SAFD_output = safd_s3D(mz_vals, mz_int, Rt, FileName, path, max_numb_iter,
    max_t_peak_w, res, min_ms_w, r_thresh, min_int, sig_inc_thresh, S2N, min_peak_w_s)


#Align the masses 
function mass_align(Rt, Mz_values, Mz_intensity)
    # Round the MS-values to the nearest integer
    Mz_values = round.(Mz_values .* 1) ./ 1
    # Get the unique masses (that were just rounded to the nearest integer)
    unique_mz_values = sort(unique(Mz_values))
    # Define the intensity matrix
    plot_matrix = zeros(length(Rt), length(unique_mz_values))
    # Make variables to save memory allocations
    inner_mz_values = Vector{Float64}(undef, size(Mz_values, 2))
    inner_mz_intensity = Vector{Float64}(undef, size(Mz_values, 2))
    inner_plot_matrix = Vector{Float64}(undef, size(plot_matrix, 2))
    i = 1
    for i in 1:length(Rt)
        # Saving a view of the matrix (less memory allocation)
        inner_mz_values .= view(Mz_values, i, :)
        inner_mz_intensity .= view(Mz_intensity, i, :)
        inner_plot_matrix .= view(plot_matrix, i, :)
        # Set pos to 1, start of Julia counting in arrays
        pos = 1
        #last_pos = 1 # removed code
        # Make a turbo loop to increase speed
        @tturbo for k in 1:length(inner_mz_values)
            # Loop over the MZ values (k) for this retention moment (i)
            for j in 1:length(unique_mz_values)
                # Check for all the unique masses (j)
                #pos != last_pos && break
                # If the current mass is equal to the current unique mass, safe the intensity
                test_outcome = (unique_mz_values[j] == inner_mz_values[k])
                inner_plot_matrix[j] += test_outcome * inner_mz_intensity[k]
                pos = j * test_outcome + !test_outcome
            end
        end
        plot_matrix[i, :] .= inner_plot_matrix
    end
    return unique_mz_values, plot_matrix
end

unique_mz_values, plot_matrix = mass_align(Rt, mz_vals, mz_int)

#Check the data visually
tic = sum(mz_int, dims=2)

p_tic = plot(Rt, tic, size=(1920, 1080), xlabel="Rt",
    ylabel="Intensity", left_margin=15Plots.mm, right_margin=7.5Plots.mm,
    bottom_margin=7.5Plots.mm, title="TIC of pest mix data", legend=false, c=:black)


p_surf = surface(unique_mz_values, Rt, plot_matrix, size=(1920, 1080), grid=true,
    camera=(45, 45),
    c=cgrad([:white, :navy, :indigo, :teal, :green, :yellow, :red], [0, 0.04, 1]),
    #color = :prism,
    #clims = (0,4000),
    title="sruface plot of pest mix data",
    ylabel="Rt",
    xlabel="m/z",
    zlabel="Intensity",
    legend=false
)

plot(p_tic, p_surf, layout=(2, 1))
#Function to calculate resolutions
function resolutions(SAFD_output)
    # Matrix with all retention time deltas
    Rt = SAFD_output[:, 4]
    D_Rt = zeros(length(SAFD_output[:, 1]), length(SAFD_output[:, 1]))
    for i = 1:length(SAFD_output[:, 1])
        for j = 1:length(SAFD_output[:, 1])
            D_Rt[i, j] = Rt[i] - Rt[j]
        end
    end
    D_Rt
    # Matrix with all mass deltas
    M = SAFD_output[:, 8]
    D_M = zeros(length(SAFD_output[:, 1]), length(SAFD_output[:, 1]))
    for i = 1:length(SAFD_output[:, 1])
        for j = 1:length(SAFD_output[:, 1])
            D_M[i, j] = M[i] - M[j]
        end
    end
    #Calculate the sigmas
    sigma = zeros(length(SAFD_output[:, 1]))
    for i = 1:length(SAFD_output[:, 1])
        sigma[i] = (SAFD_output[i, 5] - SAFD_output[i, 6]) / 4
    end
    sigma

    #Calculate the mass witdths
    mW = zeros(length(SAFD_output[:, 1]))
    for i = 1:length(SAFD_output[:, 1])
        mW[i] = (SAFD_output[i, 9] - SAFD_output[i, 10]) / 4
    end
    mW
    #Calculate resolutions
    res_Rt = zeros(length(SAFD_output[:, 1]), length(SAFD_output[:, 1]))
    res_M = zeros(length(SAFD_output[:, 1]), length(SAFD_output[:, 1]))
    for i = 1:length(SAFD_output[:, 1])
        @inbounds @simd for j = 1:length(SAFD_output[:, 1])
            res_Rt[i, j] = sqrt.((D_Rt[i, j] .^ 2 / (2 * (sigma[i] + sigma[j])) .^ 2))
            res_M[i, j] = sqrt.((D_M[i, j] .^ 2 / (2 * (mW[i] + mW[j])) .^ 2))
        end
    end

    return res_Rt, res_M
end

#Funtion to split windows for plotting

function window_split_Rt(Rt, wind_size)
    time_diff = (Rt[end] - Rt[1]) / wind_size
    Split = zeros(wind_size + 1)
    for i = 1:length(Split)
        Split[i] = time_diff * (i - 1)
    end
    Split[end] = Rt[end]
    Split
    return Split
end
Split = window_split_Rt(Rt, 10)

j = 6
i = 2
wind_size = 10
accepted_res = 1.5
SAFD_output = SAFD_output_100ppb
function Peaks_p_window(wind_size, Rt, SAFD_output)
    # This function splits the retention domain into n windows of equal length, then features that fall within
    # the retention time of each window are assigned into the different rows of Peaks_per_window
    Peaks_per_window = Vector{Vector{Int64}}(undef, wind_size)
    pos = 1
    #Divide the Rt domain in n equal parts 
    time_diff = (Rt[end] - Rt[1]) / wind_size
    #This loop checks in which of the defined windows the SAFD feature is by the center point of the feature
    for i = 1:wind_size
        inner_vector = Vector{Int64}(undef, 0)
        for j = pos:length(SAFD_output[:, 4])
            #if SAFD_output[j, 4] >= (time_diff * wind_size)
                #push!(inner_vector , j)
                #Peaks_per_window[end] = deepcopy(inner_vector)
            #end
            if SAFD_output[j, 4] >= time_diff * (i - 1) && SAFD_output[j, 4] <= time_diff * i
                push!(inner_vector, j)
            
                #When a given feature's RT exceeds the value of the next window, pos is updated to that feature and
                # the next loop starts to assign the next peaks
            elseif SAFD_output[j, 4] >= time_diff * i
                pos = j
                break
            end

            


        end

        Peaks_per_window[i] = deepcopy(inner_vector)
    end
    return Peaks_per_window, time_diff
end
P, T = Peaks_p_window(10, Rt, SAFD_output)
P

# Function that returns unresolved and resolved peaks
function unresolved_per_window_Rt_and_MS(Rt, SAFD_output, wind_size, accepted_res)
    #colors::Vector{Int32} = Vector{Int32}(undef, 0)
    #resolved_Rt::BitVector = BitVector(undef, 0)
    #Assign peaks to each of the windows
    Peaks_per_window, time_diff = Peaks_p_window(wind_size, Rt, SAFD_output)

    #Make some empty arrays to store data
    tot_unresolved_final = zeros(wind_size, 5)
    colors = ones(Int32, length(SAFD_output[:, 4]))
    #Keep track of peaks that have already been saved
    assigned_peaks = Set{Int64}()
    assigned_peaks_res = Set{Int64}()
    #This loop  will go calculate the resolution in Rt and MS for the features in each window
    #for i = 1:wind_size
    for i = 1:wind_size
        tot_unresolved_final[i, 1] = (time_diff * (i - 1))
        tot_unresolved_final[i, 2] = (time_diff * (i))
    
        if length(Peaks_per_window[i]) == 0
            tot_unresolved_final[i, 3] = 0
            
        elseif length(Peaks_per_window[i]) == 1
            for k = 1:length(Peaks_per_window[i])
                peak_index = Peaks_per_window[i][1]
                colors[peak_index] = 4
            end


        else
            #for handling the case when there are multiple peaks in the window

            #For counting unresolved peaks
            res_Rt, res_M = resolutions(SAFD_output[Peaks_per_window[i][1]:Peaks_per_window[i][end], :])
            res_Rt
            res_M
            result_Rt = res_Rt .>= accepted_res
            resolved_Rt = vec(sum(result_Rt, dims=1)) + vec(sum(result_Rt, dims=2))
            resolved_Rt_f = resolved_Rt .!= (size(res_M, 1)*2) - 2

            result_MS = res_M .>= accepted_res
            resolved_MS = vec(sum(result_MS, dims=1)) + vec(sum(result_MS, dims=2))
            resolved_MS_f = resolved_MS .!= (size(res_M, 1)*2) - 2
            
            result_f = (resolved_Rt_f.==1) .& (resolved_MS_f .==1)
            tot_unresolved_final[i, 3] = sum(result_f)
            

            #For coloring
            unresolved_Rt = vec(sum(result_Rt, dims=1)) + vec(sum(result_Rt, dims=2))
            unresolved_MS = vec(sum(result_MS, dims=1)) + vec(sum(result_MS, dims=2))
            unresolved_Rt_t = unresolved_Rt .!= size(res_M, 1)*2 - 2
            unresolved_Rt_f = sum(unresolved_Rt_t)
            unresolved_MS_t = unresolved_MS .!= size(res_M, 1)*2 - 2
            unresolved_MS_f = sum(unresolved_MS_t)
            
            # if 1, peak resolved in either RT or Mass domain
            #tot_resolved_final[i,3] = length(Peaks_per_window[i]) - sum(resolved)
            #tot_resolved_final[i, 3] = max(0, max(length(unique(Peaks_per_window[i]))) - sum(resolved))
            #max(length(Peaks_per_window[i]) - sum(resolved), 0)


            # This keeps track of which peaks are resolved or not for the coloring of peaks
            #It excludes peaks that have already appeared in other windows
            for k = 1:length(unresolved_Rt_t)
                peak_index = k + ((Peaks_per_window[i][1]) - 1)
                if unresolved_Rt_t[k] == 1 && unresolved_MS_t[k] == 1 && !in(peak_index, assigned_peaks)
                    colors[peak_index] = 1 #Red cross -> Resolved in neither
                    push!(assigned_peaks, peak_index)
                elseif unresolved_Rt_t[k] == 0 && unresolved_MS_t[k] == 1
                    !in(peak_index, assigned_peaks)
                    colors[peak_index] = 3 # Yellow cross -> Resolved in RT only
                    push!(assigned_peaks, peak_index)
                elseif unresolved_Rt_t[k] == 0 && unresolved_MS_t[k] == 0 && !in(peak_index, assigned_peaks)
                    colors[peak_index] = 4 #Green cross -> Resolved in both
                    push!(assigned_peaks, peak_index)
                elseif unresolved_Rt_t[k] == 1 && unresolved_MS_t[k] == 0 && !in(peak_index, assigned_peaks)
                    colors[peak_index] = 2 #Red circle -> Resolved in MS only
                    push!(assigned_peaks, peak_index)
                end
            end
        end
    end
    for i = 1:wind_size
        if length(Peaks_per_window[i]) > 1
            tot_unresolved_final[i,4] = (tot_unresolved_final[i,3]/length(Peaks_per_window[i]))*100
        end
        if length(Peaks_per_window[i]) > 1
            tot_unresolved_final[i,5] = (tot_unresolved_final[i,3]/length(SAFD_output[:,4]))*100
        end
    end
    

    final_df = DataFrame(Window_Start = tot_unresolved_final[:,1], Window_end = tot_unresolved_final[:,2],
                         Unresolved_peaks = tot_unresolved_final[:,3], Unresolved_compared_to_window = tot_unresolved_final[:,4],
                         Unresolved_compared_to_total = tot_unresolved_final[:,5])




    #The end result in the amount of peaks that have a Resolution in Rt and in MS lower than 1.5
    return tot_unresolved_final, colors, final_df

end

# Separated in both or either one of them
# Not separated Rt cause we can change it
# Not separated MS in red cause nothing can be done 
SAFD_output
@time result, colors, df = unresolved_per_window_Rt_and_MS(Rt, SAFD_output, 10, 1.5)
@code_warntype unresolved_per_window_Rt_and_MS(Rt, SAFD_output, 10, 1.5)
result
colors
df
sum(result[:,3])
#Final plotting

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
    bottom_margin=8.5Plots.mm
)
#color=[:red, :red, :yellow, :green]
# Create a scatter plot using the x and y coordinates and the colors and symbols vectors
paint = [:red, :orange, :yellow, :green]
colors_final = paint[colors]
mapping = Dict(1 => "Unresolved in RT and MS", 2 => "Resolved in MS only", 3 => "Resolved in RT only", 4 => "Fully resolved")
labels_final = map(x -> mapping[x], colors)
p2 = scatter!(SAFD_output_100ppb[:,4], SAFD_output_100ppb[:,8],
    #series_annotations = text.(1:length(SAFD_output[:,4]),size = 1),
    markershape = :xcross,
    color = colors_final,
    group = labels_final,
    legend = true,
    markersize=2.5,
    title="Pest mix, 2000 iterations, S/N = 3, r = 0.9, accepted_res = 1.5",
    left_margin=5Plots.mm, right_margin=7.5Plots.mm,
    bottom_margin=8.5Plots.mm,
    )

for i = 1:length(Split)
    p = plot!(ones(mz_thresh[end]) .* (Split[i]), collect(1:1:mz_thresh[end]), color=:red, legend=true, label = false)
    display(p)
end

savefig("/Users/tobias/Library/CloudStorage/OneDrive-UvA/For presentation/final algorithmheat map of pest mix.png")
savefig("/Users/tobias/Downloads/Elbow plot window nr Drug data.png")
#Visual resolution Testing
SAFD_output
peak_1 = 13
peak_2 = 12

m_1 = SAFD_output_100ppb[peak_1, 8]
m_2 = SAFD_output_100ppb[peak_2, 8]

s_1 = sqrt(((SAFD_output_100ppb[peak_1, 9] - SAFD_output_100ppb[peak_1, 10]) / 4)^2)
s_2 = sqrt(((SAFD_output_100ppb[peak_2, 9] - SAFD_output_100ppb[peak_2, 10]) / 4)^2)
D_1 = Normal(m_1, s_1)
D_2 = Normal(m_2, s_2)

M_R = sqrt((m_1 - m_2)^2 / (2 * (s_1 + s_2))^2)


plot(D_1)
p = plot!(D_2)

m_1 = SAFD_output_100ppb[peak_1, 4]
m_2 = SAFD_output_100ppb[peak_2, 4]

s_1 = sqrt(((SAFD_output_100ppb[peak_1, 5] - SAFD_output_100ppb[peak_1, 6]) / 4)^2)
s_2 = sqrt(((SAFD_output_100ppb[peak_2, 5] - SAFD_output_100ppb[peak_2, 6]) / 4)^2)
R_1 = Normal(m_1, s_1)
R_2 = Normal(m_2, s_2)

Rt_R = sqrt((m_1 - m_2)^2 / (2 * (s_1 + s_2))^2)



plot(R_1)
p_2 = plot!(R_2)

plot(p, p_2)

SAFD_output_100ppb
unresolved_peaks = zeros(150)
for i = 1:150
    result, colors, df = unresolved_per_window_Rt_and_MS(Rt, SAFD_output_Stef, i, 1.5)
    unresolved_peaks[i] = sum(result[:,3])
    @show i
end
windows = collect(0:5:150)
scatter(unresolved_peaks, size = (1920,1080), title = "Optimal nr of windows-Drug data", xlabel = "Nr of windows", ylabel = "Nr of unresolved peaks",left_margin=15Plots.mm, top_margin=7.5Plots.mm,
bottom_margin=8.5Plots.mm, grid = false, legend = false, xticks = windows)
plot!(unresolved_peaks)
