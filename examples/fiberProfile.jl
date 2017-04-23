using Gadfly
set_default_plot_size(25cm, 16cm)

include("../src/utils.jl")
include("../src/ComposeMatrices.jl")
include("../src/CoordinateDescentTM.jl")
include("../src/CoordinateDescent3.jl")
include("../src/adaLASSO.jl")
include("../src/LSIAF.jl")

## Reading data
Y = readdlm("ExperimentalData.csv", ',', Float64, '\r')
Y = Y[:,2]

## Algorithm
Selections, FilterResult, Beta_lasso = LSIAF(Y)

## Plots
plot(
    layer(x=1:size(Y,1), y=Y, Geom.line, Theme(default_color=color("orange"))),
    layer(x=1:size(Y,1), y=FilterResult, Geom.line, Theme(default_color=color("blue"))),
    layer(x=1:size(Y,1), y=Beta_lasso, Geom.line, Theme(default_color=color("red"))) )
