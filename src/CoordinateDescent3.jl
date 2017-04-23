function CoordinateDescent3{F<:AbstractFloat}(X::Array{F}, Y::Array{F}, LASSOpath::Array{F}, λ_max::F, w::Array{F}, σ::Array{F})
    ## Number of available lambdas: Granularity of the regulation parameter
    MaxIterλ = length(LASSOpath)

    n, p = size(X)

    ## Preallocating vector β_lasso for speed
    β_lasso = zeros(p, MaxIterλ)

    ## Initializing parameters
    β_tilde = zeros(p,1)

    ## Active set
    A= Int[]

    innerX = X' * X

    ## Calcula os produtos internos das colunas de X com Y
    innerY = sum(X.*repmat(Y,1,p),1)

    bestB = 1
    BICold = 100000000
    teste = 0

    β_ols = zeros(p,1)

    ## First for runs through the LASSO path
    for cont = 1:MaxIterλ
        ## Computing next lambda
        λ = LASSOpath[cont] * λ_max
        lastchange = 0
        it_while = 1

        while (teste != 2)
            ##  Reference active set
            A_anterior = A

            ##  Cycle
            if it_while == 1
                j_range = collect(1:p)'
            else
                #j_range = A_actual
                j_range = collect(1:p)'
            end

            for j = j_range

                aux, λEff, β_ols = CoordinateDescentCore(A, innerX, j, β_tilde, β_ols, n, innerY, λ, σ)



                #updateCandidates{T<:Signed, F<:AbstractFloat}(λEff::F, β_ols::Array{F},
                #                 j::T, A_actual::Array{T}, lastchange::T, β_tilde::Array{F}; w::Array{F} = [])

                A, β_tilde, lastchange, breakThis = updateCandidates(λEff, β_ols, j, A, lastchange, β_tilde, w = w)

                if breakThis == 1
                    break
                end

            end
            ##  Check if the active set has been altered in the last cycle
            teste = convCheck(A_anterior, A)
            it_while = it_while + 1

        end
        teste = 0

        BICnew = n * ( log(var( Y - X*β_tilde )) ) + sum( abs(β_tilde./σ) .> 1e-3 ) * log(n)
        if BICnew < BICold
            BICold = BICnew
            bestB = cont
        end

        β_lasso[:,cont] = β_tilde
    end


    return β_lasso, bestB
end
