
using LC_MS_Resolutions
#This file is intended as a template for the worflow to go from an MZXML file to the results of the 
#LC_MS_Resolutions output 

#Code can be run either from scratch with an MZXML file or an already processed SAFD output file from a
#previous run:
#New MZXML file:
#Set filename
pathin = "/Users/tobias/Documents/Tobi_data/20230324_pest_stds_SBC18_short_2" 
pathin = "/Users/tobias/Documents/Tobi_data/20230327_pest_stds_BonusRP" 
pathin = "/Users/tobias/Documents/Tobi_data/20230327_pest_stds_PFP" 
filenames = ["Std_04.mzXML"]
mz_thresh = [0, 0] #Sets threshold for mz values
int_thresh = 500 #Remove most of the noise
#Load the file for your gradient here
gradient_data = CSV.read("/Users/tobias/Documents/Tobi_data/Gradients/Gradient short_2.csv", DataFrame)
gradient_data = CSV.read("/Users/tobias/Documents/Tobi_data/Gradients/BONUS RP Gradient.csv", DataFrame)
gradient_data = CSV.read("/Users/tobias/Documents/Tobi_data/Gradients/POROSHELL PFP.csv", DataFrame)

#Import MS data
GC.gc()

@time mz_val, mz_int, t0, t_end, m, path, msModel, msIonisation, msManufacturer,
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
#Run SAFD (this may take some time depending on n iterations)
@time rep_table, SAFD_output = LC_MS_Resolutions.SAFD.safd_s3D(mz_val, mz_int, Rt, FileName, path, max_numb_iter,
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

#Align masses and run resolutions algorithm
unique_mz_values, plot_matrix = mass_align(Rt, mz_val, mz_int)
results, colors, df, gradient = @time unresolved_per_window_Rt_and_MS(Rt, SAFD_output, 12, 1.5, gradient_data)
df   #This DataFrame contains all the results
#Save the dataframe as CSV
CSV.write("/Users/tobias/Library/CloudStorage/OneDrive-UvA/Research project UvA/Figures and dataframes from US data/PFP, mix 4.csv" ,df)
#Plot the heatmap with features and windows
plot_heatmap(SAFD_output, Rt, unique_mz_values, plot_matrix, 12, gradient_data, colors, filenames, pathin)

savefig("/Users/tobias/Library/CloudStorage/OneDrive-UvA/Research project UvA/Figures and dataframes from US data/PFP, mix 4.png")

#################
#END 


