
using LC_MS_Resolutions
#This file is intended as a template for the worflow to go from an MZXML file to the results of the 
#LC_MS_Resolutions output 

#New MZXML file:
#Set filename
pathin = "/Users/tobias/Documents/Tobi_data/20230324_pest_stds_SBC18_short_2" 
filenames = ["Name of your file.mzXML"]
#Load the file for your gradient here
gradient_data = LC_MS_Resolutions.CSV.read("Path_to_your_gradientGradient/name of your gradient file.csv", LC_MS_Resolutions.DataFrame)

window_size = 12   #The number of windows the RT domain is split into (Default = 12)
accepted_resolution = 1.5  #The resolution accepted for two gaussian peaks to be considered separated (Baseline separation =1.5)
Resolved_peaks_algorithm(pathin, filenames,gradient_data,window_size, accepted_resolution)

#################
#END 