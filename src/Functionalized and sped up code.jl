
using Statistics, SAFD, CSV, DataFrames, LoopVectorization, StatsPlots, Distributions
using MS_Import, LinearAlgebra, CompCreate, JLD2

function mass_align(Rt::Vector{Float32}, Mz_values::Matrix{Float32}, Mz_intensity::Matrix{Float32})
    # Round the MS-values to the nearest integer
    Mz_values = round.(Mz_values .* 1) ./ 1
    # Get the unique masses (that were just rounded to the nearest integer)
    unique_mz_values::Vector{Float32} = sort(unique(Mz_values))
    # Define the intensity matrix
    plot_matrix::Matrix{Float32} = zeros(length(Rt), length(unique_mz_values))
    # Make variables to save memory allocations
    inner_mz_values = Vector{Float32}(undef, size(Mz_values, 2))
    inner_mz_intensity = Vector{Float32}(undef, size(Mz_values, 2))
    inner_plot_matrix = Vector{Float32}(undef, size(plot_matrix, 2))
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
function resolutions(SAFD_output::DataFrame)
    # Matrix with all retention time deltas

    Rt::Vector{Float32} = SAFD_output[:, 4]
    D_Rt = zeros(length(Rt), length(Rt))
    for i = 1:length(Rt)
        for j = i+1:length(Rt)
            D_Rt[i, j] = Rt[i] - Rt[j]
        end
    end

    # Matrix with all mass deltas
    M::Vector{Float32} = SAFD_output[:, 8]
    D_M = zeros(length(Rt), length(Rt))
    for i = 1:length(Rt)
        for j = i+1:length(Rt)
            D_M[i, j] = M[i] - M[j]
        end
    end
    #Calculate the sigmas
    Sigma_start::Vector{Float32} = SAFD_output[:,5]
    Sigma_end::Vector{Float32} = SAFD_output[:,6]
    sigma = zeros(length(Rt))
    for i = 1:length(Rt)
        sigma[i] = (Sigma_start[i] - Sigma_end[i]) / 4
    end

    #Calculate the mass witdths
    mW_start::Vector{Float32} = SAFD_output[:,9]
    mW_end::Vector{Float32} = SAFD_output[:,10]
    mW = zeros(length(Rt))
    for i = 1:length(Rt)
        mW[i] = (mW_start[i] - mW_end[i]) / 4
    end
    #Calculate resolutions
    res_Rt = zeros(length(Rt), length(Rt))
    res_M = zeros(length(Rt), length(Rt))
    for i = 1:length(Rt)
        @inbounds @simd for j = i+1:length(Rt)
            res_Rt[i, j] = sqrt.((D_Rt[i, j] .^ 2 / (2 * (sigma[i] + sigma[j])) .^ 2))
            res_M[i, j] = sqrt.((D_M[i, j] .^ 2 / (2 * (mW[i] + mW[j])) .^ 2))
        end
    end

    return res_Rt, res_M
end

function unresolved_per_window_Rt_and_MS(Rt::Vector{Float32}, SAFD_output::DataFrame, wind_size::Int64, accepted_res::Float64)

    #Assign peaks to each of the windows
    Peaks_per_window, time_diff = Peaks_p_window(wind_size, Rt, SAFD_output)
    #Make some empty arrays to store data
    #Store results
    tot_unresolved_final::Matrix{Float32} = zeros(wind_size, 6)
    #Assing colors for plotting
    colors = Vector{Int32}(undef, length(SAFD_output[:, 4]))
    #This loop  will go calculate the resolution in Rt and MS for the features in each window
    for i = 1:wind_size
        tot_unresolved_final[i, 1] = (time_diff * (i - 1))
        tot_unresolved_final[i, 2] = (time_diff * (i))
        #When there are no peaks in a window
        if length(Peaks_per_window[i]) == 0
            continue
        #When there is only one peak in a window
        elseif length(Peaks_per_window[i]) == 1
            peak_index = Peaks_per_window[i][1]
            colors[peak_index] = 4
        #for handling the case when there are multiple peaks in the window
        else 
            res_Rt, res_M = resolutions(SAFD_output[Peaks_per_window[i][1]:Peaks_per_window[i][end], :])
            #Initialize the bitmatrixes
            resolved_Rt_f = falses(size(res_M, 1))
            resolved_MS_f = falses(size(res_M, 1))
            result_f = falses(size(res_M, 1))
            #Retention RS time section
            result_Rt = res_Rt .>= accepted_res
            resolved_Rt = vec(sum(result_Rt, dims=1)) + vec(sum(result_Rt, dims=2))
            resolved_Rt_f = resolved_Rt .!= size(res_M, 1) -1
            #MS RS section
            result_MS = res_M .>= accepted_res
            resolved_MS = vec(sum(result_MS, dims=1)) + vec(sum(result_MS, dims=2))
            resolved_MS_f = resolved_MS .!= size(res_M, 1) -1
            #Both together
            result_f = (resolved_Rt_f) .& (resolved_MS_f)
            tot_unresolved_final[i, 4] = sum(result_f)

            #This code is to assign the colors in the plot
            for k = 1:length(resolved_Rt_f)
                peak_index = k + ((Peaks_per_window[i][1]) - 1)
                if resolved_Rt_f[k] && resolved_MS_f[k]
                    colors[peak_index] = 1 #Red cross -> Resolved in neither
                elseif !resolved_Rt_f[k] && resolved_MS_f[k]
                    colors[peak_index] = 2 #Orange cross ->Resolved only in Rt
                elseif !resolved_Rt_f[k] && !resolved_MS_f[k]
                    colors[peak_index] = 3 #Green cross -> Resolved in both
                else
                    colors[peak_index] = 4 #Yellow cross -> Resolved in MS only
                end
            end

        end
    end
    #This is just to calculate percentage unresolved based on 1)within window, 2) compared to all peaks
    for i = 1:wind_size
        if length(Peaks_per_window[i]) > 1
            tot_unresolved_final[i,5] = (tot_unresolved_final[i,4]/length(Peaks_per_window[i]))*100
        end
        if length(Peaks_per_window[i]) > 1
            tot_unresolved_final[i,6] = (tot_unresolved_final[i,4]/length(SAFD_output[:,4]))*100
        end
    end
    #To calculate the slope of the gradient in a specific window
    Split = window_split_Rt(Rt,wind_size)
    tot_unresolved_final[:,3] = Vector{Float64}(zeros((length(Split))-1))
    # first define pos at the first Split for the gradient 
    pos = findfirst(x->x>=Split[1], Rt)
    @time for i = 1:length(Split)-1
        # save the index of the i + 1 location in temp to not overwrite the pos variable
        tmp = findfirst(x->x>=Split[i+1], Rt)
        tot_unresolved_final[i,3] = (B_gradient_final[tmp] - B_gradient_final[pos])/(Split[i+1]-Split[i])
        # overwrite pos since the index in tmp is equal to the next position that we want to use
        pos = tmp
    end

    final_df::DataFrame = DataFrame(Window_Start = tot_unresolved_final[:,1], Window_end = tot_unresolved_final[:,2], Gradient_slope = 
                        tot_unresolved_final[:,3],Unresolved_peaks = tot_unresolved_final[:,4], Unresolved_compared_to_window = 
                        tot_unresolved_final[:,5],Unresolved_compared_to_total = tot_unresolved_final[:,6])

    #The end result in the amount of peaks that have a Resolution in Rt and in MS lower than 1.5
    return tot_unresolved_final, colors, final_df

end
function window_split_Rt(Rt::Vector{Float32}, wind_size::Int64)
    time_diff::Float32 = (Rt[end] - Rt[1]) / wind_size
    Split::Vector{Float64} = zeros(wind_size + 1)
    for i = 1:length(Split)
        Split[i] = time_diff * (i - 1)
    end
    Split[end] = Rt[end]
    Split
    return Split
end
function Peaks_p_window(wind_size::Int64, Rt::Vector{Float32}, SAFD_output::DataFrame)
    # This function splits the retention domain into n windows of equal length, then features that fall within
    # the retention time of each window are assigned into the different rows of Peaks_per_window
    Peaks_per_window = Vector{Vector{Int64}}(undef, wind_size)
    pos::Int64 = 1
    #Divide the Rt domain in n equal parts 
    time_diff = (Rt[end] - Rt[1]) / wind_size
    SAFD_Rt::Vector{Float32} = SAFD_output[:, 4]
    #This loop checks in which of the defined windows the SAFD feature is by the center point of the feature
    for i = 1:wind_size
        inner_vector = Vector{Int64}(undef, 0)
        for j = pos:length(SAFD_Rt)
            #if SAFD_output[j, 4] >= (time_diff * wind_size)
                #push!(inner_vector , j)
                #Peaks_per_window[end] = deepcopy(inner_vector)
            #end
            if SAFD_Rt[j] >= time_diff * (i - 1) && SAFD_Rt[j] <= time_diff * i
                push!(inner_vector, j)

                #When a given feature's RT exceeds the value of the next window, pos is updated to that feature and
                # the next loop starts to assign the next peaks
            elseif SAFD_Rt[j] >= time_diff * i
                pos = j
                break
            end


        end

        Peaks_per_window[i] = deepcopy(inner_vector)
    end
    return Peaks_per_window, time_diff
end
function plot_heatmap(SAFD::DataFrame, Rt::Vector{Float32}, unique_mz_values::Vector{Float32}, plot_matrix::Matrix{Float32}, wind_size::Int)
    split = window_split_Rt(Rt, wind_size)
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
    paint = [:Red, :Orange, :Green, :Yellow]
    colors_final::Vector{Symbol} = paint[colors]
    mapping = Dict(1 => "Unresolved in both", 2 => "Resolved only in RT", 3 => "Resolved in both MS and RT", 4 => "Resolved in MS only")
    labels_final::Vector{String} = map(x -> mapping[x], colors)
    p2 = scatter!(SAFD[:, 4], SAFD[:, 8],
        #series_annotations = text.(1:length(SAFD_output[:,4]),size = 1),
        markershape=:xcross,
        color=colors_final,
        group=labels_final,
        legend=true,
        markersize=2.5,
        title="Pest mix, 2000 iterations, S/N = 3, r = 0.9, accepted_res = 1.5",
        left_margin=5Plots.mm, right_margin=7.5Plots.mm,
        bottom_margin=8.5Plots.mm,
    )

    for i = 1:length(split)
        @show i
        p2 = plot!(ones(Int64(maximum(unique_mz_values))) .* (split[i]), collect(1:1:Int64(maximum(unique_mz_values))), color=:red, legend=true, label=false)
        display(p2)
    end
end
function calculate_percentage_coverage(intensities::Matrix{Float32}, threshold::Int64)
    # Apply threshold to intensities to remove noise
    intensities_masked::BitMatrix = intensities .>= threshold
    # Count the number of colored pixels in the heatmap plot
    total_pixels::Int32 = length(intensities_masked[:,1])*length(intensities_masked[1,:])
    #olred_pixels = sum(intensities_masked)
    colored_pixels::Int32 = length(findall(x->x==true, intensities_masked))
    # Calculate the percentage of the heatmap plot covered by the LC-MS data
    coverage_percentage::Float32 = colored_pixels / total_pixels * 100

    return coverage_percentage
end
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


