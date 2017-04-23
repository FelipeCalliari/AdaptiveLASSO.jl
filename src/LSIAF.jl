"""
LSIAF{F<:AbstractFloat}(Y::Array{F})

The function `LSIAF` performs the filtering via coordinate descent and takes only following input:

* y : Signal or time-series to be filtered (a vector of reals).
"""
function LSIAF{F<:AbstractFloat}(Y::Array{F})
    ## File Input
    Y0 = Y
    Ym = mean(Y)
    Y  = Y0 - Ym
    n  = length(Y)

    ## Timming
    tic()

    ## Adaptive Least Absolute Shrinkage and Selection Algorithm
    β_lasso, ~, X = adaLASSO(Y)

    ## Post-Lasso
    # Determine Selected Variables
    selection = Int[]
    for i = 1:size(β_lasso,1)
        if abs(β_lasso[i,:])[1] > 1e-5
            push!(selection, i)
        end
    end

    #on error
    # if sum(selection == (1))==0
    #     selection = [selection 1]
    # end
    push!(selection, 1)

    # OLS to remove LASSO BIAS
    β_ols = inv(X[:,selection]' * X[:,selection]) * X[:,selection]' * Y

    ## Output
    OutputMatrix = Array{Number}(size(selection,1), 2)
    OutputMatrix[:, 1] = selection'
    OutputMatrix[:, 2] = β_ols

    Filter = X[:, selection] * β_ols + Ym

    ## Timming
    timeAdaLASSO = toq();
    error = sum( (Filter - Y0) .^ 2 )
    println("Finished LASSO with LMS approximation error of ", error, " and elapsed time of ", timeAdaLASSO, " seconds.")

    return OutputMatrix, Filter, β_lasso
end
