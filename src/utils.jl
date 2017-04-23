# Meshgrid function (ndgrid)
function meshgrid{T}(vx::Array{T}, vy::Array{T})
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    (repmat(vx, m, 1), repmat(vy, 1, n))
end

function smart_innerY{T<:Signed, F<:AbstractFloat}(Y::Array{F}, N::T, U::Array{F}, L::Array{F}, σ1::F, μ1::F)
    ## Here, we assume that the vector Y is already normalized.

    ## The computation of lambda max involves taking all inner products of Y.
    ## The smart method allows for the computation of the i+1_th parcel based on
    ## the result of the i_th parcel.

    arraytype = typeof(Y[1])

    ## Preallocating the vector IP_vector for speed
    IP_vector = Array{arraytype}(zeros(N,1))

    index_array = collect(1:N)'
    IP_vector[1]= sum( Y.* (index_array'-μ1) / σ1 )
    sumY = 0

    for i = 2:N
        sumY = sumY + Y[i-1]
        IP_vector[i] = sumY * ( U[i-1] - L[i-1] )
    end

    return IP_vector
end

function calculate_matrices{T<:Int}(N::T)
    ## index
    i = collect(1:N-1)'

    ## calculate arrays mean and variance
    array_μ = (N-i')/N
    array_σ = sqrt( (1/(N-1)) *( (i').*array_μ .^2 + (N-i').*(1-array_μ).^2 ) )

    ## With these parameters, compose the Lower and Upper diagonal matrices
    U = - array_μ ./ array_σ
    L = (1 - array_μ) ./ array_σ

    ## Calculate μ1 and σ1
    μ1 = (N+1)/2
    σ1 = sqrt( (1/6)*(N+1)*(2*N+1)-μ1^2 )

    ## Concatenate
    array_σ2 = [σ1; array_σ]

    ## Build X
    repU, repLtranspose = meshgrid(U, L)
    #x1 = ( ((1:N)' - μ1) / σ1 )'
    #x2 = triu(repU) + tril(repLtranspose') - diagm(L[:])
    #x3 = [x1 [x2; L']]
    X = [( ((1:N)' - μ1) / σ1 )' [triu(repU) + tril(repLtranspose') - diagm(L[:]); L']]

    ## Calculate Gram Matrix (inner products matrix)
    GMat = triu( GramMatrixTM(N, U, L, σ1, μ1) );
    GMat = GMat + GMat' - N*eye(N);

    return U, L, μ1, σ1, X, GMat, array_σ2
end

function GramMatrixTM{T<:Int, F<:AbstractFloat}(N::T, U::Array{F}, L::Array{F}, σ1::F, μ1::F)
    ## builds the inner products matrix in the smart way
    ## (much faster than X'*X)
    I, J = meshgrid(collect(1:(N-1)), collect(1:(N-1)))
    Uaux, Laux = meshgrid(U, L)

    minimo = min(I, J);
    maximo = max(I, J);

    GMat1 = minimo .* (U*U') + (maximo-minimo) .* (Laux.*Uaux) + (N-maximo) .* (L*L')

    ## reduce the memory usage, since we don't need theses variables anymore...
    I, J, minimo, maximo = 0, 0, 0, 0
    gc()

    j = collect(1:N-1)';
    IP = (1/σ1) * ( j .* U' .* ((j+1)/2 - μ1) + (N-j) .* L' .* ((j+1+N)/2 - μ1) )

    GMat = [[N; IP'] [IP; GMat1]]
end

function calculate_μ{T<:Int}(N::T)
    i = 1:N-1
    (N-i')/(N)
end

function calculate_μ1{T<:Int}(N::T)
    (N+1)/2
end

function calculate_σ{T<:Int, F<:AbstractFloat}(N::T,array_μ::Array{F})
    i = 1:N-1
    sqrt( (1/(N- 1)) *( (i').*array_μ .^2 + (N-i').*(1-array_μ).^2 ) )
end

function calculate_σ1{T<:Int, F<:AbstractFloat}(N::T,μ1::F)
    sqrt( (1/6)*(N+1)*(2*N+1)-μ1^2 )
end

function CoordinateDescentCore{T<:Signed, F<:AbstractFloat}(A_actual::Array{T}, GMat::Array{F, 2}, j::T, β_tilde::Array{F}, β_ols::Array{F}, N::T, arrayInnerY::Array{F}, λ::F, array_σ2::Array{F})
    # aux = 0
    #
    # if size(A_actual)[1] >= 1
    #     aux = ( GMat[j,A_actual]' * β_tilde[A_actual] )[1]
    # end
    aux = sum(GMat[j,A_actual]' * β_tilde[A_actual])

    ##  Computation of the OLS estimator as simple regression
    β_ols[j] = 1 / N * (arrayInnerY[j] - aux) + β_tilde[j]

    ##  Soft-thresholding operator
    if j == 1
        λEff = λ
    else
        λEff = λ * array_σ2[j]
    end

    return aux, λEff, β_ols
end

function updateCandidates{T<:Signed, F<:AbstractFloat}(λEff::F, β_ols::Array{F}, j::T, A_actual::Array{T}, lastchange::T, β_tilde::Array{F}; w = [])
    breakThis = 0

    if size(w,1) == 0
        w = ones( size(β_ols) )
    end

    if λEff*w[j] >= abs(β_ols[j])
        if in(j, A_actual) # any(j == A_actual)
            splice!(A_actual,findin(A_actual,j)[1])
            lastchange = j
        elseif j == lastchange
            breakThis = 1
        end
        β_tilde[j] = 0
    else
        β_tilde[j] = sign(β_ols[j]) * (abs(β_ols[j]) - λEff*w[j])
        ## Includes an index
        if !in(j, A_actual) #not(any(j==A_actual))
            push!(A_actual, Int(j))
            sort!(A_actual)
            lastchange = j
        elseif j == lastchange
            breakThis = 1
        end
    end

    return A_actual, β_tilde, lastchange, breakThis
end

function convCheck{T<:Signed}(A_reference::Array{T}, A_actual::Array{T})
    if length(A_reference) != length(A_actual)
        ##  Did not converge
        convergenceFlag = 0
    elseif A_reference == A_actual && length(A_reference) == length(A_actual)
        ##  Converged!
        convergenceFlag = 2
    end

    return convergenceFlag
end

function calcBestB{T<:Signed, F<:AbstractFloat}(N::T, Y::Array{F}, X::Array{F}, β_lasso::Array{F}, array_σ2::Array{F})
    bestB = -1
    BICold = 1000000

    for i = 1:size(β_lasso, 2)
        BICnew = N * ( log(var( Y - X*β_lasso[:,i] )) ) + sum( abs(β_lasso[:,i]./array_σ2) .> 1e-3 ) * log(N)
        if BICnew < BICold
            BICold = BICnew
            bestB = i
        end
    end

    return bestB
end
