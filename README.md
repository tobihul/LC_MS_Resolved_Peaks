# LC-MS_Resolved_Peaks.jl

** LC_MS_Resolved_Peaks.jl is a package developed for my research thesis. The package takes as input any mzXML file (LC-MS) along with a csv of the gradient used. The data is divided into the desired number of windows in the RT domain and, after feature detection, the number of overlapping peaks in each window is reported. The advantage of this algorithm is the fact that the entire XIC of the data is used to determine wheter a feature is successfully resolved. In addition, the algorithm takes into account how well the features are distributed in the LC-MS space using a novel Voronoi cell surface coverage calculation. A score is calculated using the % of unresolved peaks in each window compared to the total together with the surface coverage. Please keep in mind that the surface coverage calculation and the final score are still being optimized and do not represent the final version.

## Installation

In order to install the ** LC_MS_Resolved_Peaks.jl package in Julia, run the follwing: "]" to enter package manager and then "add https://github.com/tobihul/LC_MS_Resolved_Peaks"

Alternatively: 

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/tobihul/LC_MS_Resolved_Peaks"))

```



An example LC-MS dataset is included on this page along with its gradient CSV file. You can use your own data as long as it is in MZXML format
