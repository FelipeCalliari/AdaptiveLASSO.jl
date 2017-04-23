function adaLASSO{F<:AbstractFloat}(Y::Array{F})

    ## Initialize Parameters
    n      = length(Y)
    BICold = 10000000
    β      = Float64[]

    ## LASSO
    tic()
    B0TM, ~, bestB0TM, X, array_σ = CoordinateDescentTM(Y, n)
    TimeLASSO = toc()

    println("LASSO took ", TimeLASSO, "s")

    β_LASSO = B0TM[:, bestB0TM]

    ## Making real copies of data, so we can change and it'll be okay!
    brTM = copy(B0TM[:, bestB0TM])
    br   = copy(brTM)

    temp = 1:size(X,2)
    idx = abs(br ./ array_σ) .<= 1e-3
    temp = temp[!idx]

    if size(temp,1) >= 1
        if temp[1] != 1
            temp = [1; temp]
        end
    else
        temp = []
    end

    adaX = X[:, temp]
    adaSigma = array_σ[temp]

    ## OLS normalization
    adabr = br[temp] ./ adaSigma
    adabr = \(adaX'*adaX,adaX'*Y) ./ adaSigma

    ## Path
    λ_max, λ_max_var = findmax(adaSigma .* abs(adaX[:,1:end]' * Y) )
    λ_max = λ_max / adaSigma[λ_max_var] / n

    MaxIterλ = 100
    passo = Float64(1 / MaxIterλ)
    grid_linear = Array{Float64}(collect(1+0.001:passo*10:10))
    LASSOpath = Array{Float64}(1 - log10(grid_linear))

    println("Initializing adaLASSO path.")

    ## Gamma Grid
    BetaGamma = []
    for gamma = 1*sort(1-log10(1:1/4:10))+0.001
        wg = 1 ./ ( abs(adabr).^gamma )
        Blasso, bestB = CoordinateDescent3(adaX, Y, LASSOpath, λ_max, wg, adaSigma)

        Best = Blasso[:,bestB]
        BetaGamma = [BetaGamma; Best]
        BICnew = n * log(var( Y - adaX*Best )) + sum( abs(Best./adaSigma) .> 1e-3) * log(n)
        if BICnew < BICold
            BICold = BICnew
            β = Best
            lastgamma = gamma
        end
    end

    println("adaLASSO path finished.")

    ## TODO

    ## Output
    temp2 = 1:size(adaSigma,1)
    idx = abs(β ./ adaSigma) .> 1e-3
    temp2 = temp2[idx]

    if size(temp2,1) >= 1
        if temp2[1] != 1
            temp2 = [1; temp2]
        end
    else
        temp2 = []
    end

    selected = temp[temp2]
    β_adaLASSO = zeros(size(X,2),1)
    β_adaLASSO[selected] = β[temp2]

    return β_adaLASSO, β_LASSO, X, array_σ, TimeLASSO
end
