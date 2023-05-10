
using Statistics, SAFD, CSV, DataFrames, LoopVectorization, StatsPlots, Distributions
using MS_Import, LinearAlgebra, CompCreate, JLD2, VoronoiCells, GeometryBasics, Random

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
function gradient_curve(data::DataFrame, Rt::Vector{Float32})

    data = Matrix(data)
    b_modifiers = Vector{Vector{Float64}}(undef, (size(data, 1) - 1))

    for i in 1:length(b_modifiers)

        t_start = data[i, 1]
        t_end = data[i+1, 1]
        b_start = data[i, 2]
        b_end = data[i+1, 2]

        # find indices of retention times within current step
        if i == 1
            idx_start = findfirst(x -> x >= t_start, Rt)
        else
            idx_start = findfirst(x -> x >= t_start, Rt) + 1
        end
        idx_end = findfirst(x -> x >= t_end, Rt)
        if i != length(b_modifiers)

            # interpolate b_modifier values for retention times within current step
            if b_start == b_end # no change in %B
                b_modifiers[i] = Vector{Float32}(ones(length(idx_start:idx_end))) .* b_end
            else  # %B increases or decreases
                b_modifiers[i] = collect(b_start:((b_end-b_start)/(idx_end-idx_start)):b_end)

            end

        else
            #For the end part to make sure the length is correct
            if b_start == b_end  # no change in %B
                b_modifiers[i] = Vector{Float32}(ones(length(idx_start:length(Rt)))) .* b_end
            else   # %B increases or decreases
                b_modifiers[i] = collect(b_start:((b_end-b_start)/(length(Rt)-idx_start)):b_end)
            end
        end
    end
    b_modifiers_final = reduce(vcat, b_modifiers)
    return b_modifiers_final
end
function unresolved_per_window_Rt_and_MS(Rt::Vector{Float32}, SAFD_output::DataFrame, wind_size::Int64, accepted_res::Float64, gradient::DataFrame)

    #Assign peaks to each of the windows
    Peaks_per_window, time_diff = Peaks_p_window(wind_size, Rt, SAFD_output)
    #Calculate orthogonality score for weach window
    #Make some empty arrays to store data
    #Store results
    tot_unresolved_final::Matrix{Float32} = zeros(wind_size, 8)
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
            #Adding weights depending on the mass Distributions
            top_95 = zeros(wind_size)
            for i = 1:wind_size
                if length(Peaks_per_window[i]) > 1
                    Prior = fit_mle(Gamma, Vector{Float32}(SAFD_output[Peaks_per_window[i],8]))
                    top_95[i] = quantile(Prior, 0.95)
                end
            end
            Mass_weights = (top_95 .-minimum(top_95))./maximum(top_95).-minimum(top_95)
            #Calculating the orthogonality of the window with weights
            Rt_norm, MS_norm = normalize_lc_ms(SAFD_output[Peaks_per_window[i],:])
            tot_unresolved_final[i,7] = surface_voronoi(Rt_norm, MS_norm, 2.2) * Mass_weights[i]
            
            #############
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
            tot_unresolved_final[i,6] = (tot_unresolved_final[i,4]/length(SAFD_output[:,4]))*100
            #This sets the final score for a window
            tot_unresolved_final[i,8] = ((tot_unresolved_final[i,6]/100) + tot_unresolved_final[i,7])
        else
            tot_unresolved_final[i,8] = 1
        end
    end
    #To calculate the slope of the gradient in a specific window
    grad = gradient_curve(gradient,Rt)
    Split = window_split_Rt(Rt,wind_size)
    tot_unresolved_final[:,3] = Vector{Float64}(zeros((length(Split))-1))
    # first define pos at the first Split for the gradient 
    pos = findfirst(x->x>=Split[1], Rt)
    for i = 1:length(Split)-1
        # save the index of the i + 1 location in temp to not overwrite the pos variable
        tmp = findfirst(x->x>=Split[i+1], Rt)
        tot_unresolved_final[i,3] = (grad[tmp] - grad[pos])/(Split[i+1]-Split[i])
        # overwrite pos since the index in tmp is equal to the next position that we want to use
        pos = tmp
    end
    
    final_df::DataFrame = DataFrame(Window_Start = tot_unresolved_final[:,1], Window_end = tot_unresolved_final[:,2], Gradient_slope = 
                        tot_unresolved_final[:,3],Unresolved_peaks = tot_unresolved_final[:,4], Unresolved_compared_to_window = 
                        tot_unresolved_final[:,5],Unresolved_compared_to_total = tot_unresolved_final[:,6], Voronoi_surface_coverage = tot_unresolved_final[:,7], final_score = tot_unresolved_final[:,8])

    #The end result in the amount of peaks that have a Resolution in Rt and in MS lower than 1.5
    return colors, final_df, grad

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
function plot_heatmap(SAFD_output::DataFrame, Rt::Vector{Float32}, unique_mz_values::Vector{Float32}, plot_matrix::Matrix{Float32}, wind_size::Int, gradient::DataFrame, colors::Vector{Int32},filenames
                     ,pathin)
    # Extract the last part of the pathin variable
    pathin_parts = splitdir(pathin)
    pathin_last = pathin_parts[end]

    # Extract the filename without extension
    filename_parts = splitext(filenames[1])
    filename_no_ext = filename_parts[1]

    # Create the title
    num_features = length(SAFD_output[:,1])
    title_str = "$(pathin_last) -> $(filename_no_ext), $(num_features) features"
    title_str = replace(title_str, "_" => " ") # optional: replace underscores with spaces
    title_str = replace(title_str, "." => "") # optional: remove dots
    title_str = uppercase(title_str) # optional: convert to uppercase
    split = window_split_Rt(Rt, 12)
    heatmap(Rt, unique_mz_values, plot_matrix',
        #c = cgrad([:white,:navy,:indigo,:teal,:green,:yellow,:red],[0,0.04,1]),
        color=:plasma,
        clims=(20000, 80000),
        size=(1280, 720),
        xlabel="Rt (min)",
        ylabel="m/z",
        left_margin=5Plots.mm, right_margin=7.5Plots.mm,
        bottom_margin=8.5Plots.mm,
        colorbar = false,
        yticks = (0:(0.1*ceil(maximum(unique_mz_values))):ceil(maximum(unique_mz_values))),
        title=title_str

    )
    # Create a scatter plot using the x and y coordinates and the colors and symbols vectors
    paint = [:Red, :hotpink1, :Green, :Orange]
    colors_final = paint[colors]
    mapping = Dict(1 => "Unresolved in RT and MS", 2 => "Resolved in Rt only", 3 => "Resolved in both", 4 => "Resolved in MS only")
    labels_final = map(x -> mapping[x], colors)
    p2 = scatter!(SAFD_output[:, 4], SAFD_output[:, 8],
        #series_annotations = text.(1:length(SAFD_output[:,4]),size = 1),
        markershape=:xcross,
        color=colors_final,
        group=labels_final,
        legend=:topleft,
        markersize=2.5,
        #title="$(filenames[1]), $max_numb_iter iterations, S/N = $S2N, r = $r_thresh, accepted_res = 1.5, Componetization -> ($(length(SAFD_output[:,1])) features)",
        left_margin=5Plots.mm, right_margin=7.5Plots.mm,
        bottom_margin=8.5Plots.mm,
    )
    window_split_Rt(Rt,wind_size)
    for i = 1:length(split)
        @show i
        p2 = plot!(ones(Int32(ceil(maximum(unique_mz_values)))) .* (split[i]), collect(1:1:Int32(ceil(maximum(unique_mz_values)))), color=:red, label=false)
        display(p2)
    end
    gradient = gradient_curve(gradient, Rt)
    p2 = plot!(twinx(), Rt, gradient, yticks = (5:5:100), legend = false, ylabel = ("%B"), linewidth = 5, linestyle = :dot, xticks = (round.(split; digits = 1)))
    return p2
end
function surface_voronoi(x::Vector{Float64},y::Vector{Float64}, k)
    #Setting the maximum possible std using the model

    max_std = (0.5960069556116379 * length(x)^(-0.4724248759777625))/sqrt(length(x))

    #Real data
    rect = Rectangle(Point2(0, 0), Point2(1, 1))
    points = Point2.(x,y)
    tess = voronoicells(points, rect)
    a_data = voronoiarea(tess)
    std_real = std(a_data)/sqrt(length(x))

    ortho_score = tanh(k * (1 - (std_real / max_std))^7.5)/tanh(k)
    #ortho_score =  (1 - (std_real / max_std))^1/k
    return ortho_score
end
function load_and_prep_data(pathin, filenames, path2features)
    #Import MS data
    GC.gc()
    mz_thresh =[0,0]
    @time mz_val, mz_int, t0, t_end, m, path, msModel, msIonisation, msManufacturer,
    polarity, Rt = import_files_MS1(pathin, filenames, mz_thresh)

    #Componentization of features
    # Parameters for CompCreate
    
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

    return unique_mz_values,plot_matrix, SAFD_output, Rt
end
function normalize_lc_ms(SAFD_output::DataFrame)
    Rt_vals = SAFD_output[:,4]
    MS_vals = SAFD_output[:,8]
    Rt_norm = zeros(length(Rt_vals))
    MS_norm = zeros(length(MS_vals))
    #Test rand data
    for i = 1:length(Rt_vals)
        Rt_norm[i] = (Rt_vals[i]-minimum(Rt_vals))/(maximum(Rt_vals)-minimum(Rt_vals))
        MS_norm[i] = (MS_vals[i]-minimum(MS_vals))/(maximum(MS_vals) - minimum(MS_vals))
    end
    return Rt_norm, MS_norm
end
