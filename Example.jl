
using LC_MS_Resolutions
#This file is intended as a template for the worflow to go from an MZXML file to the results of the 
#LC_MS_Resolutions output 

#Code can be run either from scratch with an MZXML file or an already processed SAFD output file from a
#previous run:
#New MZXML file:
#Set filename
pathin = "path_to_your_LC_MS_file" 
filenames = ["name_of_your_file.mzXML"]
mz_thresh = [0, 0] #Sets threshold for mz values
int_thresh = 500 #Remove most of the noise
#Load the file for your gradient here
gradient_data = CSV.read("path_to_your_gradient/name_of_your_gradient.csv", DataFrame)

#Import MS data
GC.gc()

mz_val, mz_int, t0, t_end, m, path, msModel, msIonisation, msManufacturer,
polarity, Rt = LC_MS_Resolutions.MS_Import.import_files_MS1(pathin, filenames, mz_thresh, int_thresh)
FileName = m[1]
#Adjust the settings for SAFD here
max_numb_iter = 2000 #Roughly the number of features that will be found, if there are less than n_iter it will stop automatically
max_t_peak_w=300 # or 20
res=20000
min_ms_w=0.02
r_thresh=0.9 #How well the fitted gaussian overlaps the original peak
min_int=10000
sig_inc_thresh=5
S2N=3 #minimum signal to noise ratio for SAFD to take data as a feature
min_peak_w_s=3
GC.gc()
#Run SAFD (this may take some time depending on n iterations and min_int)
rep_table, SAFD_output = LC_MS_Resolutions.SAFD.safd_s3D(mz_val, mz_int, Rt, FileName, path, max_numb_iter,
max_t_peak_w, res, min_ms_w, r_thresh, min_int, sig_inc_thresh, S2N, min_peak_w_s)

#Componentization of features
#Add the file location of the features reported by SAFD to run CompCreate
#The report file should be located wherever your MZXML file is
#The following code fetches the name of the file automatically given that it is stored in the same place as your dataframe
#If you stored the report in a different location, create a variable called path2features leading to your report csv using CSV.read("yourpath.csv",DataFrame)
basename_pathin = basename(pathin)
filename_no_ext = splitext(filenames[1])[1]
path2features = joinpath(pathin*"/"*filename_no_ext*"_report.csv")
mass_win_per=0.8
ret_win_per=0.5
r_thresh=0.9
delta_mass=0.004
min_int = 750

chrom=LC_MS_Resolutions.MS_Import.import_files(pathin,filenames,mz_thresh,int_thresh)


## For only MS1
SAFD_output = LC_MS_Resolutions.CompCreate.comp_ms1(chrom,path2features,mass_win_per,ret_win_per,r_thresh,delta_mass, min_int)

wind_size = 12 #define the number of windows to split the Rt domain in
resolution = 1.5 #define the accepted resolution for two features to be considered resolved

#Align masses and run resolutions algorithm
unique_mz_values, plot_matrix = mass_align(Rt, mz_val, mz_int)
colors, results, gradient = unresolved_per_window_Rt_and_MS(Rt, SAFD_output, wind_size, resolution, gradient_data)
df   #This DataFrame contains all the results
#Save the dataframe as CSV
CSV.write("Path_to_where_you_want_your_results/name_of_your_file.csv" ,results)
#Plot the heatmap with features and windows
plot_heatmap(SAFD_output, Rt, unique_mz_values, plot_matrix, 12, gradient_data, colors, filenames, pathin)

savefig("Path_to_where_you_want_your_figure/name_of_your_figure.png") #save the figure as png, pdf or svg with the extension after .

#################
#END 


