# LC-MS_Resolved_Peaks.jl

**LC_MS_Resolved_Peaks.jl** is a package developed for my research thesis. The package takes as input any mzXML file (LC-MS) along with a csv of the gradient used. The data is divided into the desired number of windows in the RT domain and, after feature detection, the number of overlapping peaks in each window is reported. The advantage of this algorithm is the fact that the entire XIC of the data is used to determine wheter a feature is successfully resolved. In addition, the algorithm takes into account how well the features are distributed in the LC-MS space using a novel Voronoi cell surface coverage calculation. A score is calculated using the % of unresolved peaks in each window compared to the total together with the surface coverage. Please keep in mind that the surface coverage calculation and the final score are still being optimized and do not represent the final version.

## Installation

In order to install the **LC_MS_Resolved_Peaks.jl** package in Julia, run the follwing: "]" to enter package manager and then "add https://github.com/tobihul/LC_MS_Resolved_Peaks"

Alternatively: 

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/tobihul/LC_MS_Resolved_Peaks"))

```
## Usage

First, import data using [**MS_Import.jl**](https://bitbucket.org/SSamanipour/ms_import.jl/src/master/) and import your gradient csv:

The gradient csv should be made in the following format:

Time | %B
------------ | -------------
0 | 5
10| 95
11| 95
11.1| 95
10| 5

### Data import, feature detection and componentization
```julia
using LC_MS_Resolutions

pathin = "path_to_your_mzXML_file"
filenames = "Name_of_your_mzXML_file.mzXML"
mz_thresh = [0, 0] #Sets threshold for mz values
int_thresh = 500 #Remove most of the noise
gradient_data = CSV.read("path_to_your_gradient.csv", DataFrame)

mz_val, mz_int, t0, t_end, m, path, msModel, msIonisation, msManufacturer,
polarity, Rt = LC_MS_Resolutions.MS_Import.import_files_MS1(pathin, filenames, mz_thresh, int_thresh)
```
Now we can run SAFD which has been described [elsewhere](https://pubs.acs.org/doi/full/10.1021/acs.analchem.9b02422):

Adjust settings for SAFD:

```julia
max_numb_iter = 2000 
max_t_peak_w=300 # or 20
res=20000
min_ms_w=0.02
r_thresh=0.9 #How well the fitted gaussian overlaps the original peak
min_int=10000
sig_inc_thresh=5
S2N=3 #minimum signal to noise ratio for SAFD to take data as a feature
min_peak_w_s=3
```

Run SAFD: 

```julia
rep_table, SAFD_output = LC_MS_Resolutions.SAFD.safd_s3D(mz_val, mz_int, Rt, FileName, path, max_numb_iter,
max_t_peak_w, res, min_ms_w, r_thresh, min_int, sig_inc_thresh, S2N, min_peak_w_s)
```

Run componentization, this can also be found [elsewhere](https://bitbucket.org/SSamanipour/compcreate.jl/src/master/):

```julia
basename_pathin = basename(pathin)
filename_no_ext = splitext(filenames[1])[1]
path2features = joinpath(pathin*"/"*filename_no_ext*"_report.csv")
mass_win_per=0.8
ret_win_per=0.5
r_thresh=0.9
delta_mass=0.004
min_int = 750
chrom=LC_MS_Resolutions.MS_Import.import_files(pathin,filenames,mz_thresh,int_thresh)

# For only MS1
SAFD_output = LC_MS_Resolutions.CompCreate.comp_ms1(chrom,path2features,mass_win_per,ret_win_per,r_thresh,delta_mass, min_int)
```
Once the Componentization has been done, the file can be saved to avoid having to rerun SAFD or **CompCreate** again

Now we can align the masses for plotting 

```julia
unique_mz_values, plot_matrix = mass_align(Rt, mz_val, mz_int)
```

## Determining resolved peaks per window and surface coverage

We can now run the main function of **LC_MS_Resolved_Peaks.jl**

###Inputs
* **Rt::Vector{Float32}** this is the list of retention times of features detected by SAFD
* **SAFD_output::DataFrame** this DataFrame contains all of the information about each feature detected in SAFD (Rt start, Rt end, MS start, MS end, peak area, etc.)
* **wind_size::Int** This determines the number of windows that the Rt domain will be split up in for resolution calculations
* **resolution::Int** This is the accepted resolution deciding if two peaks are resolved in Rt and/or MS domain. A resolution of 1.5 is generally accepted as baseline separation between two gaussian peaks
* **gradient_data::DataFrame** This is your csv containing the two columns with time and %B modifier describing the gradient used

```julia
wind_size = 12
resolution = 1.5
results, colors, df, gradient = unresolved_per_window_Rt_and_MS(Rt, SAFD_output, wind_size, resolution, gradient_data)
```
Finally, the results can be plotted to see where the unresolved peaks are in each window and what the gradient is doing

```julia
plot_heatmap(SAFD_output, Rt, unique_mz_values, plot_matrix, 12, gradient_data, colors, filenames, pathin)
```

This is an example of what the plot may look like:

![Alt Text](https://github.com/tobihul/LC_MS_Resolved_Peaks/blob/master/Short%202%2C%20mix%204.png?raw=true)

## Example

An example mzXML file and gradient can be found in the repository along with the **Example.jl** file to be used as a template

## Acknowledgements
I thank Saer Samanipour and Bob Pirok for supervising me during my thesis and providing insights throughout. I also thank Etienne Kant for providing the mass align and helping with bugs and code optimization. 
