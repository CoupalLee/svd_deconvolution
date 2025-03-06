using LinearAlgebra
using StaticArrays
using DelimitedFiles
using BenchmarkTools
using Plots

#########
#= 
    This script contains the functions necessary to perform signal deconvolution using the singular
    value decomposition (SVD) approach.

    The convolution equation in its discretized form is given by:
    
    m = A * s,
    
    where A contains the transfer function and s is the signal to be deconvolved.
    Since A is a rectangular matrix, an inverse does not exist. The pseudo-inverse is used instead.
    
    The pseudo-inverse A^-1 is computed using the SVD of A.

    This script is a polished version of a homework for the course "OPTIMISATION NUMÉRIQUE ET APPLICATIONS"

    An effort is made to leverage Julia's in place features to optimize memory usage and speed.

=#
#########

function read_signal_data(filename::String)
    ## 2 columns data
    ## time points and signal values

    data = readdlm(filename)
    return data
end

function naive_gram_schmidt(V)

    ## Math background here
    ## U need not be orthogonal to decompose A 
    ## but if you want to leverage the fact that The
    ## inverse of U is U^T instead of inv(U) to deconvolute
    ## you need U to be orthogonal.

    Nrow = size(V)[1]
    Ncolumn = size(V)[2]    

    U = zeros(Nrow, Ncolumn)

    for column in 1:Ncolumn
        U[:,column] = V[:,column]
        for i in 1:column-1
            U[:,column] -= (U[:,i]' * V[:,column]) * U[:,i]
        end
        U[:,column] /= norm(U[:,column])
    end

    return U

end

function compute_orthogonal_complement(U)

    ## We have column u_1 to u_r that are orthonormal and span
    ## the subspace R^r 
    ## We compute the orthogonal complement of the subspace R^r
    ## to get the full orthonormal basis of R^m

    ## We build a matrix K whose rows are the columns u_1 to u_r
    ## and we compute the kernel of that matrix, i.e., the set
    ## of vectors that K maps to 0.
    
    #K = nullspace(U')
    #P = naive_gram_schmidt(K)
    #U[:,r+1:end] .= P

    #U = naive_gram_schmidt(U)

end

function naive_reduced_SVD(A)

    ## Potentially instable

    eigenvalues, V = eigen(A' * A)

    singular_values = sort(sqrt.(eigenvalues), rev=true)

    r = sum(singular_values .> 1e-3)

    Nm = size(A)[1]
    Ns = size(A)[2]

    U = zeros(Nm, r)
    D = diagm(singular_values[1:r])
    V = V[:,1:r]

    for (i, σ_i) in enumerate(singular_values[1:r])
        U[:,i] = A * V[:,i] / σ_i
    end

    U[:,1:r] = naive_gram_schmidt(U[:,1:r])


    return U, D, V'
end

function efficient_SVD(A)
    ## Golub and Kahan, 1965

end

function buid_transfer_matrix!(A, h, N)

    ## Transfer function : h

    ## The transfer matrix A is a Toeplitz matrix.
    ## and has the form
    ## A = [\vec{h} 0 ... 0  0]
    ##     [0 \vec{h} 0 ... 0] 
    ##     [0 0 \vec{h} 0 ... 0]
    ##     [0 0 0 ... 0 ... \vec{h}]
    
    Ns, Nm, Nh = N
    
    ## Column major algorithm
    counter = 0
    for column in 1:Ns
        for row in 1:Nm
            if row < Nh
                A[row + counter, column] = h[row]
            else
                break
            end
        end
        counter += 1
    end 
end

function deconvolution(U, Σ, V_T, m)

    s = U * Σ * V_T * m
    return m_test
    
end
function main()

    m = read_signal_data("")
    s = read_signal_data("")
    h = read_signal_data("")

    Nm = size(m)[1]
    Ns = size(s)[1]
    Nh = size(h)[1]

    A = zeros(Nm, Ns)
    

    buid_transfer_matrix!(A, h[:,2], (Ns, Nm, Nh))

    #############################
    ## SVD benchmark naive vs builtin svd
    ## @btime U, Σ, V_T = naive_SVD($A)
    ## @btime svd($A)
    ## 488.191 ms (28 allocations: 126.37 MiB)
    ## 869.665 ms (11 allocations: 102.55 MiB)
    ## Not too bad for a naive implementation
    #############################



    U, D, V_T = naive_reduced_SVD(A)

    A_test = U * D * V_T
    V = V_T'

    A_inv = V * inv(D) * U'

    s_test = A_inv * m[:,2]
    
    Plots.plot(s[:,1], s[:,2], label="s", c="red")
    Plots.plot!(s[:,1], s_test, label="s_test", c = "green")

    
end

main()

