function ComposeMatrices{T<:Signed, F<:AbstractFloat}(N::T, Y::Array{F}; useCachedMatrices::Bool = false)

    readtype = typeof(Y[1])
    if readtype == Float32
        N = Int32(N)
    else
        N = Int64(N)
    end

    if useCachedMatrices
        check = isfile(string(N, "mu.csv")) && isfile(string(N, "sigma.csv")) && isfile(string(N, "U.csv")) && isfile(string(N, "L.csv")) && isfile(string(N, "X.csv")) &&  isfile(string(N, "GMat.csv")) && isfile(string(N, "sigma2.csv"))

        if check
            println("Reading pre-calculated matrices.")

            mu1         = readcsv(string(N, "mu.csv"), readtype)
            mu1         = mu1[1]
            sigma1      = readcsv(string(N, "sigma.csv"), readtype)
            sigma1      = sigma1[1]

            U           = readcsv(string(N, "U.csv"), readtype)
            L           = readcsv(string(N, "L.csv"), readtype)
            X           = readcsv(string(N, "X.csv"), readtype)
            GMat        = readcsv(string(N, "GMat.csv"), readtype)
            arraySigma2 = readcsv(string(N, "sigma2.csv"), readtype)
        else
            println("Pre-Calculating Matrices.. This might take a while.")

            ## Calculate Matrices
            U, L, mu1, sigma1, X, GMat, arraySigma2 = calculate_matrices(N)
            gc()

            writecsv(string(N, "U.csv"), U)
            writecsv(string(N, "L.csv"), L)
            writecsv(string(N, "mu.csv"), mu1)
            writecsv(string(N, "sigma.csv"), sigma1)
            writecsv(string(N, "X.csv"), X)
            writecsv(string(N, "GMat.csv"), GMat)
            writecsv(string(N, "sigma2.csv"), arraySigma2)
        end
    else
        println("Pre-Calculating Matrices.. This might take a while.")
        U, L, mu1, sigma1, X, GMat, arraySigma2 = calculate_matrices(N)
        gc()
    end

    println("Matrices finished. Ready to run LASSO path.");

    ## Compute maximum lambda
    arrayInnerY = smart_innerY(Y, N, U, L, sigma1, mu1)
    # TODO - verify lambda_max
    lambda_max, lambda_max_var = findmax(arraySigma2.*abs(arrayInnerY))
    lambda_max = lambda_max / arraySigma2[lambda_max_var] / N

    return X, GMat, arrayInnerY, arraySigma2, lambda_max
end
