function CoordinateDescentTM{T<:Signed, F<:AbstractFloat}(Y::Array{F}, N::T; MaxIterλ::T = 100)
    ## Check variable type
    #readtype = typeof(Y[1])
    #if readtype == Float32
    #    N = Int32(N)
    #    MaxIterλ = Int32(MaxIterλ)
    #else
    #    N = Int64(N)
    #    MaxIterλ = Int64(MaxIterλ)
    #end

    ## Build LASSO-path

    ## Number of available lambdas: Granularity of the regulation parameter
    # MaxIterλ = 100

    ## Potential Penalization Parameters Grid
    passo = Float64(1 / MaxIterλ)

    ## Linear grid
    grid_linear = Array{Float64}(collect(1+0.001:passo*10:10))

    ## Logarithmic Grid
    LASSO_Path = Array{Float64}(1 - log10(grid_linear))

    ## Compose Matrices
    X, GMat, arrayInnerY, array_σ2, λ_max = ComposeMatrices(N, Y)

    ## LASSO Path
    println("Initializing LASSO path.")

    ## Initializing parameters
    β_tilde = zeros(N,1)
    β_ols = zeros(N,1)
    convergenceFlag = 0
    ## Active set
    A_actual = Int[]

    ## Preallocating vector β_lasso for speed
    β_lasso = zeros(N,length(LASSO_Path))

    ## First for runs through the LASSO path
    for cont = 1:length(LASSO_Path)-2
        ## cont
        ## Computing next lambda
        λ = LASSO_Path[cont] * λ_max
        lastchange = 0
        it_while = 1
        while (convergenceFlag != 2)
            ##  Reference active set
            A_reference = A_actual

            ##  Cycle
            if it_while == 1
                j_range = collect(1:N)'
            else
                j_range = A_actual
            end

            for j = j_range

                aux, λEff, β_ols = CoordinateDescentCore(A_actual, GMat, j, β_tilde, β_ols, N, arrayInnerY, λ, array_σ2)

                A_actual, β_tilde, lastchange, breakThis = updateCandidates(λEff, β_ols, j, A_actual, lastchange, β_tilde)

                if breakThis == 1
                    break
                end

            end
            ##  Check if the active set has been altered in the last cycle
            convergenceFlag = convCheck(A_reference, A_actual)
            it_while = it_while + 1

        end
        it_while
        convergenceFlag = 0

        β_lasso[:,cont] = β_tilde
    end

    println("LASSO path finished.")

    ## BIC
    bestB = calcBestB(N, Y, X, β_lasso, array_σ2)

    return β_lasso, LASSO_Path, bestB, X, array_σ2
end
