# calculate eigenvalues of a matrix using the QR algorithm. This function is differntiable
function classical_gram_schmidt(A)
    n, m = size(A)
    Q = zeros(n, m)
    R = zeros(m, m)

    for j = 1:m
        v = A[:, j]
        for i = 1:j-1
            R[i, j] = Q[:, i]' * A[:, j]
            v -= R[i, j] * Q[:, i]
        end
        # don't use norm
        R[j, j] = sqrt(sum(v.^2))
        Q[:, j] = v / R[j, j]
    end

    return Q, R
end

function matrix_norm(A)
    return sqrt(sum(A .^ 2))
end
function diag1(A)
    n, m = size(A)
    if n != m
        error("Matrix must be square")
    end

    d = zeros(n)
    for i in 1:n
        d[i] = A[i, i]
    end

    return d
end

function qr_eigen(A; max_iters=10000, tolerance=1e-10)
    n, m = size(A)
    if n != m
        error("Matrix must be square")
    end

    Ak = copy(A)
    Qk = zeros(n, n)
    for i in 1:n
        Qk[i, i] = 1.0
    end

    for i in 1:max_iters
        Q, R = classical_gram_schmidt(Ak)
        Ak = R * Q
        Qk *= Q  

        # Check for convergence using the Frobenius norm of the off-diagonal elements
        off_diagonal_norm = sqrt(sum(Ak[i, j]^2 for i in 1:n for j in 1:n if i != j))
        if off_diagonal_norm < tolerance
            break
        end
    end

    eigenvalues = diag1(Ak)
    eigenvectors = Qk
    return eigenvalues, eigenvectors
end
