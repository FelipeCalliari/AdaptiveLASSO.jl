# AdaptiveLASSO

AdaptiveLASSO.jl is a Julia package for [...].

<!---
This package is a port of a Matlab code from ["adaLASSO Matlab Code"][linkcodematlab].

It is originally a supplemental material for our paper ["TITULO PAPER"][linkpaper].
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

## :blue_book: TODO

- [ ] Write something about this package
- [ ] Cite the url of the adaLASSO Matlab code
- [ ] Cite the paper

[linkcodematlab]: https://github.com/Tungstenio/
[linkpaper]: http://arxiv.org/pdf/
