# AdaptiveLASSO

AdaptiveLASSO.jl is a Julia package able to distinguish meaningful level shifts from typical signal fluctuations which, in a fiber monitoring context, is associated with the problem of identifying small losses within a noisy optical time-domain reflectometer (OTDR) profile. This is the Julia version of the original supplemental material for the paper ["Adaptive Filter for Automatic Identification of Mulitple Fiber Faults in a Noisy OTDR Profile"][linkartigo].

<!---
This package is a port of a Matlab code from ["adaLASSO Matlab Code"][linkcodematlab].


--->

## Installing

In order to use it, download the package in the julia prompt with the command:
```
Pkg.clone("git://github.com/FelipeCalliari/AdaptiveLASSO.jl.git")
```

Then simply include the command `using AdaptiveLASSO` to import the package.

## Example

```
# after initializing the package with
# Pkg.clone("git://github.com/FelipeCalliari/AdaptiveLASSO.jl.git")

# then you must include the library
using AdaptiveLASSO

# read the files
Y = readdlm("ExperimentalData.csv", ',', Float64, '\r')
Y = Y[:,2]

# run the algorithm
Selections, FilterResult, Beta_lasso = LSIAF(Y)
```

## Updating

There are two options to update the package:

### First

```
Pkg.checkout("AdaptiveLASSO")
```

### Second

The second option is to remove the package and add it again:

```
Pkg.rm("AdaptiveLASSO")
workspace()
Pkg.clone("git://github.com/FelipeCalliari/AdaptiveLASSO.jl.git")
```

[linkcodematlab]: https://github.com/Tungstenio/
[linkartigo]: http://ieeexplore.ieee.org/document/7471419/
