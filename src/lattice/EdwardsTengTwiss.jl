using LinearAlgebra # LinearAlgebra is not fully supported by Enzyme. E.g., det, eigen, inv, etc. are not supported.
abstract type AbstractTwiss end

"""
    lu_decomposition(A)

Performs an in-place LU decomposition of a square matrix A.
"""
function lu_decomposition(A::AbstractMatrix)
    n = size(A, 1)
    LU = copy(A)
    for k in 1:n
        for i in k+1:n
            LU[i, k] /= LU[k, k]
            for j in k+1:n
                LU[i, j] -= LU[i, k] * LU[k, j]
            end
        end
    end
    return LU
end

"""
    lu_solve(LU, b)

Solves the system Ax = b given the LU decomposition of A.
"""
function lu_solve(LU::AbstractMatrix, b::AbstractVector)
    n = size(LU, 1)
    x = zeros(eltype(b), n)
    y = zeros(eltype(b), n)

    # Forward substitution (solves Ly = b)
    for i in 1:n
        s = sum((LU[i, j] * y[j] for j in 1:i-1), init=zero(eltype(b)))
        y[i] = b[i] - s
    end

    # Backward substitution (solves Ux = y)
    for i in n:-1:1
        s = sum((LU[i, j] * x[j] for j in i+1:n), init=zero(eltype(b)))
        x[i] = (y[i] - s) / LU[i, i]
    end
    
    return x
end

function inv1(A::AbstractMatrix)
    n = size(A, 1)
    if size(A, 2) != n
        error("Matrix must be square.")
    end
    
    LU = lu_decomposition(A)
    I_mat = Matrix{eltype(A)}(I, n, n)
    A_inv = similar(A)
    
    for j in 1:n
        # Solve A * x = e_j for the j-th column of the inverse
        A_inv[:, j] = lu_solve(LU, I_mat[:, j])
    end
    
    return A_inv
end

function det1(A::Matrix)
    # determinant of a arbitrary matrix. Not used in the current implementation.
    # LU decomposition
    N = size(A, 1)
    U = zeros(Float64, N, N)
	for i = 1:N
		U[i, :] = A[i, :]
	end
    L = zeros(Float64, N, N)
	for i = 1:N
		L[i, i] = 1.0
	end
    swaps = 0

    # Manually create the identity matrix
    L = zeros(Float64, N, N)
    for i = 1:N
        L[i, i] = 1.0
    end

    for i = 1:N
        # Find the pivot
        pivot = i
        for j = i+1:N
            if abs(U[j, i]) > abs(U[pivot, i])
                pivot = j
            end
        end

        # Swap rows if necessary
        if pivot != i
            U[i, :], U[pivot, :] = U[pivot, :], U[i, :]
            L[i, 1:i-1], L[pivot, 1:i-1] = L[pivot, 1:i-1], L[i, 1:i-1]
            swaps += 1
        end

        # Check for singular matrix
        if U[i, i] == 0
            return 0.0
        end

        # Eliminate below
        for j = i+1:N
            L[j, i] = U[j, i] / U[i, i]
            U[j, i:end] -= L[j, i] * U[i, i:end]
        end
    end

    # Compute the determinant
    det_val = (-1.0)^swaps
    for i = 1:N
        det_val *= U[i, i]
    end

    return det_val
end

function det_small_matrix(A::Matrix)
	nrows, ncols = size(A)
	if nrows == 2
	    return A[1,1]*A[2,2] - A[1,2]*A[2,1]
	elseif nrows == 3
	    return A[1,1]*(A[2,2]*A[3,3] - A[2,3]*A[3,2]) -
           A[1,2]*(A[2,1]*A[3,3] - A[2,3]*A[3,1]) +
           A[1,3]*(A[2,1]*A[3,2] - A[2,2]*A[3,1])
	elseif nrows == 4
		m = A  # Alias for brevity

		# Compute the 4x4 determinant using the Laplace expansion
		det_val = m[1,1]*(
					 m[2,2]*(m[3,3]*m[4,4] - m[3,4]*m[4,3]) -
					 m[2,3]*(m[3,2]*m[4,4] - m[3,4]*m[4,2]) +
					 m[2,4]*(m[3,2]*m[4,3] - m[3,3]*m[4,2])
				 ) - m[1,2]*(
					 m[2,1]*(m[3,3]*m[4,4] - m[3,4]*m[4,3]) -
					 m[2,3]*(m[3,1]*m[4,4] - m[3,4]*m[4,1]) +
					 m[2,4]*(m[3,1]*m[4,3] - m[3,3]*m[4,1])
				 ) + m[1,3]*(
					 m[2,1]*(m[3,2]*m[4,4] - m[3,4]*m[4,2]) -
					 m[2,2]*(m[3,1]*m[4,4] - m[3,4]*m[4,1]) +
					 m[2,4]*(m[3,1]*m[4,2] - m[3,2]*m[4,1])
				 ) - m[1,4]*(
					 m[2,1]*(m[3,2]*m[4,3] - m[3,3]*m[4,2]) -
					 m[2,2]*(m[3,1]*m[4,3] - m[3,3]*m[4,1]) +
					 m[2,3]*(m[3,1]*m[4,2] - m[3,2]*m[4,1])
				 )
	
		return det_val
	end
	println(stderr, "Matrix size not supported.")
	return 0.0
end

struct EdwardsTengTwiss{T} <: AbstractTwiss 
	betax::T
	betay::T
	alphax::T
	alphay::T
	gammax::T
	gammay::T
	dx::T
	dpx::T
	dy::T
	dpy::T
	mux::T
	muy::T
	sinmux::T
	cosmux::T
	sinmuy::T
	cosmuy::T
	R::Matrix{T}
	mode::Int
end

"""
	EdwardsTengTwiss(betax::Float64, betay::Float64; alphax::Float64=0.0, alphay::Float64=0.0,
	dx::Float64=0.0, dy::Float64=0.0, dpx::Float64=0.0, dpy::Float64=0.0,
	mux::Float64=0.0, muy::Float64=0.0,
	R11::Float64=0.0, R12::Float64=0.0, R21::Float64=0.0, R22::Float64=0.0,
	mode::Int=1)

Construct a `EdwardsTengTwiss` object with betax and betay. All other parameters are optional.

# Arguments
- `betax::Float64`: Horizontal beta function.
- `betay::Float64`: Vertical beta function.
- `alphax::Float64=0.0`: Horizontal alpha function.
- `alphay::Float64=0.0`: Vertical alpha function.
- `dx::Float64=0.0`: Horizontal dispersion.
- `dy::Float64=0.0`: Vertical dispersion.
- `dpx::Float64=0.0`: derivative of horizontal dispersion.
- `dpy::Float64=0.0`: derivative of vertical dispersion.
- `mux::Float64=0.0`: Horizontal phase advance.
- `muy::Float64=0.0`: Vertical phase advance.
- `R11::Float64=0.0`: Matrix Element R11.
- `R12::Float64=0.0`: Matrix Element R12.
- `R21::Float64=0.0`: Matrix Element R21.
- `R22::Float64=0.0`: Matrix Element R22.
- `mode::Int=1`: mode for calculation.
"""											
function EdwardsTengTwiss(betax::Float64, betay::Float64;
						alphax::Float64 = 0.0, alphay::Float64 = 0.0,
						dx::Float64 = 0.0, dy::Float64 = 0.0,
						dpx::Float64 = 0.0, dpy::Float64 = 0.0,
						mux::Float64 = 0.0, muy::Float64 = 0.0,
						R11::Float64 = 0.0, R12::Float64 = 0.0,
						R21::Float64 = 0.0, R22::Float64 = 0.0,
						mode::Int = 1)
	gammax = (1.0 + alphax^2) / betax
	gammay = (1.0 + alphay^2) / betay
	sinmux = sin(mux)
	cosmux = cos(mux)
	sinmuy = sin(muy)
	cosmuy = cos(muy)
	dmux = mux  
	dmuy = muy  

	R = Matrix{Float64}(undef, 2, 2)
	R[1, 1] = R11
	R[1, 2] = R12
	R[2, 1] = R21
	R[2, 2] = R22

	# Return the struct instance
	return EdwardsTengTwiss{Float64}(betax, betay, alphax, alphay, gammax, gammay,
		dx, dpx, dy, dpy, dmux, dmuy,
		sinmux, cosmux, sinmuy, cosmuy, R, mode)
end

"""
	EdwardsTengTwiss(betax::DTPSAD{N,T}, betay::DTPSAD{N,T}; 
		alphax::DTPSAD{N,T}=zero(DTPSAD{N,T}),
		alphay::DTPSAD{N,T}=zero(DTPSAD{N,T}),
		dx::DTPSAD{N,T}=zero(DTPSAD{N,T}),
		dy::DTPSAD{N,T}=zero(DTPSAD{N,T}),
		dpx::DTPSAD{N,T}=zero(DTPSAD{N,T}),
		dpy::DTPSAD{N,T}=zero(DTPSAD{N,T}),
		mux::DTPSAD{N,T}=zero(DTPSAD{N,T}),
		muy::DTPSAD{N,T}=zero(DTPSAD{N,T}),
		R11::DTPSAD{N,T}=zero(DTPSAD{N,T}),
		R12::DTPSAD{N,T}=zero(DTPSAD{N,T}),
		R21::DTPSAD{N,T}=zero(DTPSAD{N,T}),
		R22::DTPSAD{N,T}=zero(DTPSAD{N,T}),
		mode::Int=1) where {N, T <: Number}
Construct a `EdwardsTengTwiss` object with betax and betay in TPSA (DTPSAD type) format. All other parameters are optional.
# Arguments
- `betax::DTPSAD{N,T}`: Horizontal beta function.
- `betay::DTPSAD{N,T}`: Vertical beta function.
- `alphax::DTPSAD{N,T}=zero(DTPSAD{N,T})`: Horizontal alpha function.
- `alphay::DTPSAD{N,T}=zero(DTPSAD{N,T})`: Vertical alpha function.
- `dx::DTPSAD{N,T}=zero(DTPSAD{N,T})`: Horizontal dispersion.
- `dy::DTPSAD{N,T}=zero(DTPSAD{N,T})`: Vertical dispersion.
- `dpx::DTPSAD{N,T}=zero(DTPSAD{N,T})`: derivative of horizontal dispersion.
- `dpy::DTPSAD{N,T}=zero(DTPSAD{N,T})`: derivative of vertical dispersion.
- `mux::DTPSAD{N,T}=zero(DTPSAD{N,T})`: Horizontal phase advance.
- `muy::DTPSAD{N,T}=zero(DTPSAD{N,T})`: Vertical phase advance.
- `R11::DTPSAD{N,T}=zero(DTPSAD{N,T})`: Matrix Element R11.
- `R12::DTPSAD{N,T}=zero(DTPSAD{N,T})`: Matrix Element R12.
- `R21::DTPSAD{N,T}=zero(DTPSAD{N,T})`: Matrix Element R21.
- `R22::DTPSAD{N,T}=zero(DTPSAD{N,T})`: Matrix Element R22.
- `mode::Int=1`: mode for calculation.
"""
function EdwardsTengTwiss(betax::DTPSAD{N,T}, betay::DTPSAD{N,T};
		alphax::DTPSAD{N,T} = zero(DTPSAD{N,T}),
		alphay::DTPSAD{N,T} = zero(DTPSAD{N,T}),
		dx::DTPSAD{N,T} = zero(DTPSAD{N,T}),
		dy::DTPSAD{N,T} = zero(DTPSAD{N,T}),
		dpx::DTPSAD{N,T} = zero(DTPSAD{N,T}),
		dpy::DTPSAD{N,T} = zero(DTPSAD{N,T}),
		mux::DTPSAD{N,T} = zero(DTPSAD{N,T}),
		muy::DTPSAD{N,T} = zero(DTPSAD{N,T}),
		R11::DTPSAD{N,T} = zero(DTPSAD{N,T}),
		R12::DTPSAD{N,T} = zero(DTPSAD{N,T}),
		R21::DTPSAD{N,T} = zero(DTPSAD{N,T}),
		R22::DTPSAD{N,T} = zero(DTPSAD{N,T}),
		mode::Int = 1) where {N, T <: Number}
		gammax = (1.0 + alphax^2) / betax
		gammay = (1.0 + alphay^2) / betay
		sinmux = sin(mux)
		cosmux = cos(mux)
		sinmuy = sin(muy)
		cosmuy = cos(muy)
		dmux = mux
		dmuy = muy
		R = Matrix{DTPSAD{N,T}}(undef, 2, 2)
		R[1, 1] = R11
		R[1, 2] = R12
		R[2, 1] = R21
		R[2, 2] = R22
		return EdwardsTengTwiss{DTPSAD{N,T}}(betax, betay, alphax, alphay, gammax, gammay,
			dx, dpx, dy, dpy, dmux, dmuy,
			sinmux, cosmux, sinmuy, cosmuy, R, mode)
end

"""
    find_closed_orbit(line::Vector{AbstractElement{Float64}}, dp::Float64=0.0; mass::Float64=m_e, energy::Float64=1e9,
	guess::Vector{Float64}=zeros(Float64, 6), max_iter::Int=20, tol::Float64=1e-8)
Calculates the closed orbit for a given lattice using Newton's method.
# Arguments
- `line::Vector{AbstractElement{Float64}}`: The lattice represented as a vector of elements.
- `dp::Float64=0.0`: Relative momentum deviation.
- `mass::Float64=m_e`: Mass of the particle.
- `energy::Float64=1e9`: Energy of the particle in eV.
- `guess::Vector{Float64}=zeros(Float64, 6)`: Initial guess for the closed orbit.
- `max_iter::Int=20`: Maximum number of iterations.
- `tol::Float64=1e-8`: Tolerance for convergence.
# Returns
- `x_closed::Vector{Float64}`: The closed orbit coordinates.
- `M::Matrix{Float64}`: The one-turn transfer matrix at the closed orbit
"""
function find_closed_orbit(line::Vector{AbstractElement{Float64}}, dp::Float64=0.0; mass::Float64=m_e, energy::Float64=1e9,
    guess::Vector{Float64}=zeros(Float64, 6), max_iter::Int=20, tol::Float64=1e-8)

    x = copy(guess)
    for i in 1:max_iter
        # Perform one-turn tracking
        M, x_out = track_a_turn_numeric(line, x, dp, mass=mass, energy=energy)
        
        # Check for convergence in 4D transverse space
        if norm(x_out[1:4] - x[1:4]) < tol
            return x_out[1:6], M
        end

        # Newton's method step in 4D using our pure Julia LU solver
        A = M[1:4, 1:4] - I
        b = -(x_out[1:4] - x[1:4])
        
        LU = lu_decomposition(A)
        delta_x = lu_solve(LU, b)
        
        x[1:4] = x[1:4] + delta_x
        
        # Keep longitudinal coordinates from the tracked output for the next guess
        x[5:6] = x_out[5:6]
    end
    
    @warn "Closed orbit did not converge after $max_iter iterations."
    M, x_out = track_a_turn_numeric(line, x, dp; mass=mass, energy=energy)
    return x_out[1:6], M
end

function track_particle(line::Vector{AbstractElement{Float64}}, x_in::Vector{Float64}, dp::Float64=0.0; 
		mass::Float64=m_e, energy::Float64=1e9)
    x = copy(x_in)
	x[6] = dp
    beam = Beam(reshape(x, 1, 6), energy=energy, mass=mass)
    linepass!(line, beam)
    return beam.r[1, :]
end
function track_a_turn_numeric(line::Vector{AbstractElement{Float64}}, x_in::Vector{Float64}, dp::Float64=0.0; 
	delta::Float64=1e-9, mass::Float64=m_e, energy::Float64=1e9)
    x_out_base = track_particle(line, x_in, dp; mass=mass, energy=energy)
    M = zeros(6, 6)
    for j in 1:6
        x_perturbed = copy(x_in)
        x_perturbed[j] += delta
        x_out_perturbed = track_particle(line, x_perturbed, dp; mass=mass, energy=energy)
        M[:, j] = (x_out_perturbed - x_out_base) / delta
    end
    return M, x_out_base
end

"""
	find_closed_orbit(line::Vector{AbstractElement{DTPSAD{N, T}}}, dp::Float64=0.0; mass::Float64=m_e, energy::Float64=1e9,
	guess::Vector{DTPSAD{N, T}}=zeros(DTPSAD{N, T}, 6), max_iter::Int=20, tol::Float64=1e-8) where {N, T <: Number}
Calculates the closed orbit for a given lattice using Newton's method in TPSA (DTPSAD type) format.
# Arguments
- `line::Vector{AbstractElement{DTPSAD{N, T}}}`: The lattice represented as a vector of elements.
- `dp::Float64=0.0`: Relative momentum deviation.
- `mass::Float64=m_e`: Mass of the particle.
- `energy::Float64=1e9`: Energy of the particle in eV.
- `guess::Vector{DTPSAD{N, T}}=zeros(DTPSAD{N, T}, 6)`: Initial guess for the closed orbit.
- `max_iter::Int=20`: Maximum number of iterations.
- `tol::Float64=1e-8`: Tolerance for convergence.
# Returns
- `x_closed::Vector{DTPSAD{N, T}}`: The closed orbit
- `M::Matrix{DTPSAD{N, T}}`: The one-turn transfer matrix at the closed orbit
"""
function find_closed_orbit(line::Vector{AbstractElement{DTPSAD{N, T}}}, dp::Float64=0.0; mass::Float64=m_e, energy::Float64=1e9,
    guess::Vector{DTPSAD{N, T}}=zeros(DTPSAD{N, T}, 6), max_iter::Int=20, tol::Float64=1e-8) where {N, T}

    x = copy(guess)
    for i in 1:max_iter
        # Perform one-turn tracking
        M, x_out = track_a_turn_numeric(line, x, dp, mass=mass, energy=energy)

        # Check for convergence in 4D transverse space
        norm_xout_x = DTPSAD(0.0)
        for j in 1:4
            norm_xout_x += (x_out[j] - x[j])^2
        end
        if sqrt(norm_xout_x) < tol
            return x_out[1:6], M
        end

        # Newton's method step in 4D using our pure Julia LU solver
        I_DPTSAD = zeros(DTPSAD{N, T}, 4, 4)
        for k in 1:4
            I_DPTSAD[k, k] = DTPSAD(1.0)
        end
        A = M[1:4, 1:4] - I_DPTSAD
        b = -(x_out[1:4] - x[1:4])
        
        LU = lu_decomposition(A)
        delta_x = lu_solve(LU, b)
        
        x[1:4] = x[1:4] + delta_x
        
        # Keep longitudinal coordinates from the tracked output for the next guess
        x[5:6] = x_out[5:6]
    end
    
    @warn "Closed orbit did not converge after $max_iter iterations."
    M, x_out = track_a_turn_numeric(line, x, dp, mass=mass, energy=energy)
    return x_out[1:6], M
end

function track_particle(line::Vector{AbstractElement{DTPSAD{N, T}}}, x_in::Vector{DTPSAD{N, T}}, dp::Float64=0.0; 
	mass::Float64=m_e, energy::Float64=1e9) where {N, T}
    x = copy(x_in)
    x[6] = dp
    beam = Beam(reshape(x, 1, 6), energy=energy, mass=mass)
    linepass!(line, beam)
    return beam.r[1, :]
end
function track_a_turn_numeric(line::Vector{AbstractElement{DTPSAD{N, T}}}, x_in::Vector{DTPSAD{N, T}}, dp::Float64=0.0; 
	mass::Float64=m_e, energy::Float64=1e9, delta::Float64=1e-9) where {N, T}
    x_out_base = track_particle(line, x_in, dp; mass=mass, energy=energy)
    M = zeros(DTPSAD{N, T}, 6, 6)
    for j in 1:6
        x_perturbed = copy(x_in)
        x_perturbed[j] += delta
        x_out_perturbed = track_particle(line, x_perturbed, dp; mass=mass, energy=energy)
        M[:, j] = (x_out_perturbed - x_out_base) / delta
    end
    return M, x_out_base
end

"""
	symplectic_conjugate_2by2(M::Matrix{T}) where T
Compute the symplectic conjugate of a 2x2 matrix M.
# Arguments
- `M::Matrix{T}`: The input 2x2 matrix.

# Returns
- `Matrix{T}`: The symplectic conjugate of the input matrix.
"""
function symplectic_conjugate_2by2(M::Matrix{T}) where T
	M_new = Matrix{T}(undef, 2, 2)
	M_new[1, 1] = M[2, 2]
	M_new[1, 2] = -M[1, 2]
	M_new[2, 1] = -M[2, 1]
	M_new[2, 2] = M[1, 1]
	return M_new
end

function matrixTransform_2by2(M::Matrix{Float64})
	m11 = M[1, 1]
	m21 = M[2, 1]
	m12 = M[1, 2]
	m22 = M[2, 2]
	M_new = zeros(Float64, 3, 3)
	M_new[1, 1] = m11^2
	M_new[1, 2] = -2 * m11 * m12
	M_new[1, 3] = m12^2
	M_new[2, 1] = -m11 * m21
	M_new[2, 2] = 1.0 + 2 * m12 * m21
	M_new[2, 3] = -m12 * m22
	M_new[3, 1] = m21^2
	M_new[3, 2] = -2 * m21 * m22
	M_new[3, 3] = m22^2
	# return [m11*m11 -2m11*m12 m12*m12
	# -m11*m21 1.0+2m12*m21 -m12*m22
	# m21*m21 -2m21*m22 m22*m22]
	return M_new
end

function matrixTransform_2by2(M::Matrix{DTPSAD{N,T}}) where {N,T}
	m11 = M[1, 1]
	m21 = M[2, 1]
	m12 = M[1, 2]
	m22 = M[2, 2]
	M_new = zeros(DTPSAD{N,T}, 3, 3)
	M_new[1, 1] = m11^2
	M_new[1, 2] = -2 * m11 * m12
	M_new[1, 3] = m12^2
	M_new[2, 1] = -m11 * m21
	M_new[2, 2] = 1.0 + 2 * m12 * m21
	M_new[2, 3] = -m12 * m22
	M_new[3, 1] = m21^2
	M_new[3, 2] = -2 * m21 * m22
	M_new[3, 3] = m22^2
	# return [m11*m11 -2m11*m12 m12*m12
	# -m11*m21 1.0+2m12*m21 -m12*m22
	# m21*m21 -2m21*m22 m22*m22]
	return M_new
end

"""
	twissPropagate(tin::EdwardsTengTwiss{Float64},M::Matrix{Float64})

Propagate the Twiss parameters through a matrix M.

# Arguments
- `tin::EdwardsTengTwiss{Float64}`: Input Twiss parameters.
- `M::Matrix{Float64}`: Transfer matrix.

# Returns
- `EdwardsTengTwiss{Float64}`: Output Twiss parameters.
"""
function twissPropagate(tin::EdwardsTengTwiss{Float64},M::Matrix{Float64})
	A= M[1:2,1:2]
	B= M[1:2,3:4]
	C= M[3:4,1:2]
	D= M[3:4,3:4]

	R1=tin.R
	_R1=symplectic_conjugate_2by2(R1)
	if tin.mode == 1
		X=A-B*R1
		t=det_small_matrix(X)
			if t>0.1
				R=(D*R1-C)*symplectic_conjugate_2by2(X)
				R/=t
				X/=sqrt(t)
				Y=D+C*_R1
				Y/=sqrt(det_small_matrix(Y))
				mode=1
			else
				X=C-D*R1
				X/=sqrt(det_small_matrix(X))
				Y=B+A*_R1
				t=det_small_matrix(Y)
				R=-(D+C*_R1)*symplectic_conjugate_2by2(Y)
				R/=t
				Y/=sqrt(t)
				mode=2
			end
	elseif tin.mode == Int(2) 
		X=B+A*_R1
		t=det_small_matrix(X)
			if t>0.1
				R=-(D+C*_R1)*symplectic_conjugate_2by2(X)
				R/=t
				X/=sqrt(t)
				Y=C-D*R1
				Y/=sqrt(det_small_matrix(Y))
				mode=1
			else
				X=D+C*_R1
				X/=sqrt(det_small_matrix(X))
				Y=A-B*R1
				t=det_small_matrix(Y)
				R=(D*R1-C)*symplectic_conjugate_2by2(Y)
				R/=t
				Y/=sqrt(t)
				mode=2
			end
	else
		#throw(AssertionError("Mode should be integer 1 or 2."))
		println(stderr,"Invalid mode for EdwardsTengTwiss.")
		# return EdwardsTengTwiss(;betax=1.0,betay=1.0,mode=0)
		# error("Invalid mode.")
	end

	Nx=matrixTransform_2by2(X)
	Ny=matrixTransform_2by2(Y)
	v1=Nx*[tin.betax;tin.alphax;tin.gammax]
	v2=Ny*[tin.betay;tin.alphay;tin.gammay]
	eta=( M[1:4,1:4])*[tin.dx,tin.dpx,tin.dy,tin.dpy]+( M[1:4,6])
	sin_dmux=X[1,2]/sqrt(v1[1]*tin.betax)
	cos_dmux=X[1,1]*sqrt(tin.betax/v1[1])-tin.alphax*sin_dmux
	sin_dmuy=Y[1,2]/sqrt(v2[1]*tin.betay)
	cos_dmuy=Y[1,1]*sqrt(tin.betay/v2[1])-tin.alphay*sin_dmuy

	smux0,cmux0,smuy0,cmuy0=tin.sinmux,tin.cosmux,tin.sinmuy,tin.cosmuy
	smux=sin_dmux*cmux0+cos_dmux*smux0
	cmux=cos_dmux*cmux0-sin_dmux*smux0
	smuy=sin_dmuy*cmuy0+cos_dmuy*smuy0
	cmuy=cos_dmuy*cmuy0-sin_dmuy*smuy0

	# Calculate the change in phase (Delta mux and Delta muy) in radians
	delta_mux = atan(sin_dmux, cos_dmux)
	delta_muy = atan(sin_dmuy, cos_dmuy)

	new_mux = tin.mux + delta_mux
	new_muy = tin.muy + delta_muy
	return EdwardsTengTwiss{Float64}(v1[1],v2[1],v1[2],v2[2],v1[3],v2[3],eta[1],eta[2],eta[3],eta[4],new_mux,new_muy,smux,cmux,smuy,cmuy,R,mode)
end

"""
	twissPropagate(tin::EdwardsTengTwiss{DTPSAD{N,T}},M::Matrix{DTPSAD{N,T}}) where {N,T}
Propagate the Twiss parameters through a matrix M in TPSA (DTPSAD type) format.
"""
function twissPropagate(tin::EdwardsTengTwiss{DTPSAD{N,T}},M::Matrix{DTPSAD{N,T}}) where {N,T}
	A= M[1:2,1:2]
	B= M[1:2,3:4]
	C= M[3:4,1:2]
	D= M[3:4,3:4]

	R1=tin.R
	_R1=symplectic_conjugate_2by2(R1)
	if tin.mode == 1
		X=A-B*R1
		t=det_small_matrix(X)
			if t>0.1
				R=(D*R1-C)*symplectic_conjugate_2by2(X)
				R/=t
				X/=sqrt(t)
				Y=D+C*_R1
				Y/=sqrt(det_small_matrix(Y))
				mode=1
			else
				X=C-D*R1
				X/=sqrt(det_small_matrix(X))
				Y=B+A*_R1
				t=det_small_matrix(Y)
				R=-(D+C*_R1)*symplectic_conjugate_2by2(Y)
				R/=t
				Y/=sqrt(t)
				mode=2
			end
	elseif tin.mode == Int(2) 
		X=B+A*_R1
		t=det_small_matrix(X)
			if t>0.1
				R=-(D+C*_R1)*symplectic_conjugate_2by2(X)
				R/=t
				X/=sqrt(t)
				Y=C-D*R1
				Y/=sqrt(det_small_matrix(Y))
				mode=1
			else
				X=D+C*_R1
				X/=sqrt(det_small_matrix(X))
				Y=A-B*R1
				t=det_small_matrix(Y)
				R=(D*R1-C)*symplectic_conjugate_2by2(Y)
				R/=t
				Y/=sqrt(t)
				mode=2
			end
	else
		error("Invalid mode for EdwardsTengTwiss.")
	end

	Nx=matrixTransform_2by2(X)
	Ny=matrixTransform_2by2(Y)
	v1=Nx*[tin.betax;tin.alphax;tin.gammax]
	v2=Ny*[tin.betay;tin.alphay;tin.gammay]
	eta=( M[1:4,1:4])*[tin.dx,tin.dpx,tin.dy,tin.dpy]+( M[1:4,6])
	sin_dmux=X[1,2]/sqrt(v1[1]*tin.betax)
	cos_dmux=X[1,1]*sqrt(tin.betax/v1[1])-tin.alphax*sin_dmux
	sin_dmuy=Y[1,2]/sqrt(v2[1]*tin.betay)
	cos_dmuy=Y[1,1]*sqrt(tin.betay/v2[1])-tin.alphay*sin_dmuy

	smux0,cmux0,smuy0,cmuy0=tin.sinmux,tin.cosmux,tin.sinmuy,tin.cosmuy
	smux=sin_dmux*cmux0+cos_dmux*smux0
	cmux=cos_dmux*cmux0-sin_dmux*smux0
	smuy=sin_dmuy*cmuy0+cos_dmuy*smuy0
	cmuy=cos_dmuy*cmuy0-sin_dmuy*smuy0

	# Calculate the change in phase (Delta mux and Delta muy) in radians
	delta_mux = atan(sin_dmux, cos_dmux)
	delta_muy = atan(sin_dmuy, cos_dmuy)

	new_mux = tin.mux + delta_mux
	new_muy = tin.muy + delta_muy
	return EdwardsTengTwiss{DTPSAD{N,T}}(v1[1],v2[1],v1[2],v2[2],v1[3],v2[3],eta[1],eta[2],eta[3],eta[4],new_mux,new_muy,smux,cmux,smuy,cmuy,R,mode)
end

"""
	findm66(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int; E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))

Find the 6x6 transfer matrix.

# Arguments
- `seq::Vector{<:AbstractElement{Float64}}`: Sequence of elements.
- `dp::Float64`: Relative momentum deviation.
- `order::Int`: Order of the TPSA. If `order == 0`, a fast numerical method is used.
- `E0::Float64=3e9`: Reference energy in eV.
- `m0::Float64=m_e`: Particle mass.
- `orb::Vector{Float64}=zeros(6)`: Initial orbit.		
# Returns
- `Matrix{Float64}`: 6x6 transfer matrix.
"""
function findm66(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int; E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
	map = zeros(Float64, 6, 6)
	if order == 0
		map .= fastfindm66(seq, dp, E0=E0, m0=m0, orb=orb)
		return map
	end
	x = CTPS(orb[1], 1, 6, order)
	px = CTPS(orb[2], 2, 6, order)
	y = CTPS(orb[3], 3, 6, order)
	py = CTPS(orb[4], 4, 6, order)
	z = CTPS(orb[5], 5, 6, order)
	delta = CTPS(dp, 6, 6, order)
	rin = [x, px, y, py, z, delta]
	# no radiation, cavity off
	linepass_TPSA!(seq, rin, E0=E0, m0=m0)

	for i in 1:6
		for j in 1:6
			map[i, j] = rin[i].map[j + 1]
		end
	end
	return map
end

"""
	findm66(seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, order::Int; E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6)) where {N, T}
Find the 6x6 transfer matrix for a DTPSAD lattice.
# Arguments
- `seq::Vector{<:AbstractElement{DTPSAD{N, T}}}`: Sequence of elements.
- `dp::Float64`: Relative momentum deviation.
- `order::Int`: Only `order == 0` is supported.
- `E0::Float64=3e9`: Reference energy in eV.
- `m0::Float64=m_e`: Particle mass.
- `orb::Vector{Float64}=zeros(6)`: Initial orbit.
# Returns
- `Matrix{Float64}`: 6x6 transfer matrix.
"""
function findm66(seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, order::Int; 
		E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6)) where {N, T}
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
	map = zeros(Float64, 6, 6)
	if order == 0
		map .= fastfindm66(seq, dp, E0=E0, m0=m0, orb=orb)
		return map
	else
		# only support order == 0
		println(stderr, "findm66: order > 0 is not supported for AbstractTPSAElement. Zero matrix will be returned.")
		return map
	end
end

function ADfindm66(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, changed_idx::Vector{Int}, changed_ele::Vector{<:AbstractElement{Float64}}; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
	map = zeros(Float64, 6, 6)
	if order == 0
		map .= ADfastfindm66(seq, dp, changed_idx, changed_ele, E0=E0, m0=m0)
		return map
	end
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
	x = CTPS(orb[1], 1, 6, order)
	px = CTPS(orb[2], 2, 6, order)
	y = CTPS(orb[3], 3, 6, order)
	py = CTPS(orb[4], 4, 6, order)
	z = CTPS(orb[5], 5, 6, order)
	delta = CTPS(dp, 6, 6, order)
	rin = [x, px, y, py, z, delta]
	# no radiation, cavity off
	ADlinepass_TPSA!(seq, rin, changed_idx, changed_ele, E0=E0, m0=m0)

	for i in 1:6
		for j in 1:6
			map[i, j] = rin[i].map[j + 1]
		end
	end
	return map
end

function ADfindm66_refpts(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, refpts::Vector{Int}, changed_idx::Vector{Int}, changed_ele::Vector{<:AbstractElement{Float64}}; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
	x = CTPS(orb[1], 1, 6, order)
	px = CTPS(orb[2], 2, 6, order)
	y = CTPS(orb[3], 3, 6, order)
	py = CTPS(orb[4], 4, 6, order)
	z = CTPS(orb[5], 5, 6, order)
	delta = CTPS(dp, 6, 6, order)
	rin = [x, px, y, py, z, delta]
	# no radiation, cavity off
	map_list = zeros(6, 6, length(refpts))
	for i in eachindex(refpts)
		if i == 1
			ADlinepass_TPSA!(seq[1:refpts[i]], rin, changed_idx, changed_ele, E0=E0, m0=m0)
		else
			ADlinepass_TPSA!(seq[refpts[i-1]+1:refpts[i]], rin, changed_idx, changed_ele, E0=E0, m0=m0)
		end
		map = zeros(Float64, 6, 6)
		for j in 1:6
			for k in 1:6
				map[j, k] = rin[j].map[k + 1]
			end
		end
		map_list[:,:,i] = map
		reassign!(rin[1], rin[1].map[1], 1)
		reassign!(rin[2], rin[2].map[1], 2)
		reassign!(rin[3], rin[3].map[1], 3)
		reassign!(rin[4], rin[4].map[1], 4)
		reassign!(rin[5], rin[5].map[1], 5)
		reassign!(rin[6], rin[6].map[1], 6)
	end

	return map_list
end

"""
	findm66_refpts(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, refpts::Vector{Int}; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
Find the 6x6 transfer matrix at specified reference points using high-order TPSA.
# Arguments
- `seq::Vector{<:AbstractElement{Float64}}`: Sequence of elements.
- `dp::Float64`: Relative momentum deviation.
- `order::Int`: Order of the TPSA. 
- `refpts::Vector{Int}`: Indices of reference points.
- `E0::Float64=3e9`: Reference energy in eV.
- `m0::Float64=m_e`: Particle mass.
- `orb::Vector{Float64}=zeros(6)`: Initial orbit.		
# Returns
- `Mat_list::Array{Float64,3}`: 6x6 transfer matrices at each reference point.
"""
function findm66_refpts(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, refpts::Vector{Int}; E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
	x = CTPS(orb[1], 1, 6, order)
	px = CTPS(orb[2], 2, 6, order)
	y = CTPS(orb[3], 3, 6, order)
	py = CTPS(orb[4], 4, 6, order)
	z = CTPS(orb[5], 5, 6, order)
	delta = CTPS(dp, 6, 6, order)
	rin = [x, px, y, py, z, delta]
	# no radiation, cavity off
	map_list = zeros(6, 6, length(refpts))
	num = 1
	for i in eachindex(refpts)
		if i == 1
			linepass_TPSA!(seq[1:refpts[i]], rin, E0=E0, m0=m0)
		else
			linepass_TPSA!(seq[refpts[i-1]+1:refpts[i]], rin, E0=E0, m0=m0)
		end
		map = zeros(Float64, 6, 6)
		for j in 1:6
			for k in 1:6
				map[j, k] = rin[j].map[k + 1]
			end
		end
		map_list[:,:,i] = map
		num += 1
		reassign!(rin[1], rin[1].map[1], 1)
		reassign!(rin[2], rin[2].map[1], 2)
		reassign!(rin[3], rin[3].map[1], 3)
		reassign!(rin[4], rin[4].map[1], 4)
		reassign!(rin[5], rin[5].map[1], 5)
		reassign!(rin[6], rin[6].map[1], 6)
	end

	return map_list
end


"""
	fastfindm66(LATTICE::Vector{<:AbstractElement{Float64}}, dp=0.0; E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
Find the 6x6 transfer matrix of a lattice using numerical differentiation.
# Arguments
- `LATTICE`: Beam line sequence.
- `dp::Float64=0.0`: Momentum deviation.
- `E0::Float64=3e9`: Reference energy in eV.
- `m0::Float64=m_e`: Particle mass.
- `orb::Vector{Float64}=zeros(6)`: Initial orbit.
# Returns
- `M66`: 6x6 transfer matrix.
"""
function fastfindm66(LATTICE::Vector{<:AbstractElement{Float64}}, dp=0.0; E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
	DPstep = 3e-8
    XYStep = 3e-8  # Default step size for numerical differentiation
    scaling =  [XYStep, XYStep, XYStep, XYStep, DPstep, DPstep]
    
    # find initial orbit
	orbitin = [orb[1] orb[2] orb[3] orb[4] orb[5] dp]

    # Build a diagonal matrix of initial conditions
	# manually build the matrix
    D6 = zeros(6, 6)
    for i in 1:6
        D6[i, i] = scaling[i] * 0.5
    end
    
    RIN = zeros(13, 6)
    for i in 1:6
        RIN[i, :] = orbitin[1,:] .+ D6[i, :]
    end
    for i in 1:6
        RIN[i+6, :] = orbitin[1,:] .- D6[i, :]
    end
    RIN[13, :] = orbitin

    beam = Beam(RIN, energy=E0, mass=m0)
    linepass!(LATTICE, beam)
    TMAT3 = transpose(beam.r)
    M66 = (TMAT3[:,1:6] - TMAT3[:,7:12]) ./ scaling

    return M66
end

"""
	fastfindm66(LATTICE::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64=0.0; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{DTPSAD{N, T}}=zeros(6)) where {N, T}
Find the 6x6 transfer matrix of a lattice using numerical differentiation for DTPSAD elements.
# Arguments
- `LATTICE`: Beam line sequence.
- `dp::Float64=0.0`: Momentum deviation.
- `E0::Float64=3e9`: Reference energy in eV.
- `m0::Float64=m_e`: Particle mass.
- `orb::Vector{DTPSAD{N, T}}=zeros(6)`: Initial orbit.
# Returns
- `M66`: 6x6 transfer matrix.
"""
function fastfindm66(LATTICE::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64=0.0; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{DTPSAD{N, T}}=zeros(6)) where {N, T}
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
	DPstep = 3e-8
    XYStep = 3e-8  # Default step size for numerical differentiation
    scaling =  [XYStep, XYStep, XYStep, XYStep, DPstep, DPstep]
    
    # find initial orbit
	orbitin = [orb[1] orb[2] orb[3] orb[4] orb[5] dp]
	orbitin = DTPSAD.(orbitin)

    # Build a diagonal matrix of initial conditions
	# manually build the matrix
    D6 = zeros(DTPSAD{NVAR(), Float64}, 6, 6)
    for i in 1:6
        D6[i, i] = scaling[i] * 0.5
    end
    
    RIN = zeros(DTPSAD{NVAR(), Float64}, 13, 6)
    for i in 1:6
        RIN[i, :] = orbitin[1,:] .+ D6[i, :]
    end
    for i in 1:6
        RIN[i+6, :] = orbitin[1,:] .- D6[i, :]
    end
    RIN[13, :] = orbitin

    beam = Beam(RIN, energy=E0, mass=m0)
    linepass!(LATTICE, beam)
    TMAT3 = transpose(beam.r)
    M66 = (TMAT3[:,1:6] - TMAT3[:,7:12]) ./ scaling

    return M66
end

"""
	fastfindm66_refpts(LATTICE::Vector{<:AbstractElement{Float64}}, dp::Float64, refpts::Vector{Int}; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
Find the 6x6 transfer matrix at specified reference points using numerical differentiation.
# Arguments
- `LATTICE`: Beam line sequence.
- `dp::Float64`: Momentum deviation.
- `refpts::Vector{Int}`: Indices of reference points.
- `E0::Float64=3e9`: Reference energy in eV.
- `m0::Float64=m_e`: Particle mass.
- `orb::Vector{Float64}=zeros(6)`: Initial orbit.
# Returns
- `M66_refpts`: 6x6 transfer matrices at each reference point.
"""
function fastfindm66_refpts(LATTICE::Vector{<:AbstractElement{Float64}}, dp::Float64, refpts::Vector{Int}; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
    # assume the closed orbit is zero
	DPstep = 3e-8
    XYStep = 3e-8  # Default step size for numerical differentiation
    scaling =  [XYStep, XYStep, XYStep, XYStep, DPstep, DPstep]
    
    # find initial orbit
	orbitin = [orb[1] orb[2] orb[3] orb[4] orb[5] dp]

    # Build a diagonal matrix of initial conditions
	# manually build the matrix
    D6 = zeros(6, 6)
    for i in 1:6
        D6[i, i] = scaling[i] * 0.5
    end
    
    RIN = zeros(13, 6)
    for i in 1:6
        RIN[i, :] = orbitin[1,:] .+ D6[i, :]
    end
    for i in 1:6
        RIN[i+6, :] = orbitin[1,:] .- D6[i, :]
    end
    RIN[13, :] = orbitin

    
    # rout = linepass!(LATTICE, beam, refpts)
    M66_refpts = zeros(6, 6, length(refpts))
    for i in 1:length(refpts)
		beam = Beam(copy(RIN), energy=E0, mass=m0)
		if i == 1
			linepass!(LATTICE[1:refpts[1]], beam)
		else
			linepass!(LATTICE[refpts[i-1]+1:refpts[i]], beam)
		end
        TMAT3 = transpose(beam.r)
        M66 = (TMAT3[:,1:6] - TMAT3[:,7:12]) ./ scaling
        M66_refpts[:, :, i] = M66
    end
    return M66_refpts
end

"""
	fastfindm66_refpts(LATTICE::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, refpts::Vector{Int}; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6)) where {N, T}
Find the 6x6 transfer matrix at specified reference points using numerical differentiation for DTPSAD elements.
# Arguments
- `LATTICE`: Beam line sequence.
- `dp::Float64`: Momentum deviation.
- `refpts::Vector{Int}`: Indices of reference points.
- `E0::Float64=3e9`: Reference energy in eV.
- `m0::Float64=m_e`: Particle mass.
- `orb::Vector{Float64}=zeros(6)`: Initial orbit.
# Returns
- `M66_refpts`: 6x6 transfer matrices at each reference point.
"""
function fastfindm66_refpts(LATTICE::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, refpts::Vector{Int}; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6)) where {N, T}
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
    # assume the closed orbit is zero
	DPstep = 3e-8
    XYStep = 3e-8  # Default step size for numerical differentiation
    scaling =  [XYStep, XYStep, XYStep, XYStep, DPstep, DPstep]
    
    # find initial orbit
	orbitin = [orb[1] orb[2] orb[3] orb[4] orb[5] dp]
	orbitin = DTPSAD.(orbitin)

    # Build a diagonal matrix of initial conditions
	# manually build the matrix
    D6 = zeros(DTPSAD{NVAR(), Float64}, 6, 6)
    for i in 1:6
        D6[i, i] = scaling[i] * 0.5
    end

    RIN = zeros(DTPSAD{NVAR(), Float64}, 13, 6)
    for i in 1:6
        RIN[i, :] = orbitin[1,:] .+ D6[i, :]
    end
    for i in 1:6
        RIN[i+6, :] = orbitin[1,:] .- D6[i, :]
    end
    RIN[13, :] = orbitin

    
    # rout = linepass!(LATTICE, beam, refpts)
    M66_refpts = zeros(DTPSAD{NVAR(), Float64}, 6, 6, length(refpts))
    for i in 1:length(refpts)
		beam = Beam(copy(RIN), energy=E0, mass=m0)
		if i == 1
			linepass!(LATTICE[1:refpts[1]], beam)
		else
			linepass!(LATTICE[refpts[i-1]+1:refpts[i]], beam)
		end
        TMAT3 = transpose(beam.r)
        M66 = (TMAT3[:,1:6] - TMAT3[:,7:12]) ./ scaling
        M66_refpts[:, :, i] = M66
    end
    return M66_refpts
end

function ADfastfindm66(LATTICE::Vector{<:AbstractElement{Float64}}, dp, changed_idx, changed_ele; E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
    # assume the closed orbit is zero
    NE = length(LATTICE)
	DPstep = 3e-8
    XYStep = 3e-8  # Default step size for numerical differentiation
    scaling =  [XYStep, XYStep, XYStep, XYStep, DPstep, DPstep]
    
    # find initial orbit
	orbitin = [orb[1] orb[2] orb[3] orb[4] orb[5] dp]

    # Build a diagonal matrix of initial conditions
	# manually build the matrix
    D6 = zeros(6, 6)
    for i in 1:6
        D6[i, i] = scaling[i] * 0.5
    end
    
    RIN = zeros(13, 6)
    for i in 1:6
        RIN[i, :] = orbitin[1,:] .+ D6[i, :]
    end
    for i in 1:6
        RIN[i+6, :] = orbitin[1,:] .- D6[i, :]
    end
    RIN[13, :] = orbitin

    beam = Beam(RIN, energy=E0, mass=m0)
    
    ADlinepass!(LATTICE, beam, changed_idx, changed_ele)
    TMAT3 = transpose(beam.r)
    M66 = (TMAT3[:,1:6] - TMAT3[:,7:12]) ./ scaling

    return M66
end

function ADfastfindm66_refpts(LATTICE, dp::Float64, refpts::Vector{Int}, changed_idx, changed_ele; E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
    # assume the closed orbit is zero
	DPstep = 3e-8
    XYStep = 3e-8  # Default step size for numerical differentiation
    scaling =  [XYStep, XYStep, XYStep, XYStep, DPstep, DPstep]
    
    # find initial orbit
	orbitin = [orb[1] orb[2] orb[3] orb[4] orb[5] dp]

    # Build a diagonal matrix of initial conditions
	# manually build the matrix
    D6 = zeros(6, 6)
    for i in 1:6
        D6[i, i] = scaling[i] * 0.5
    end
    
    RIN = zeros(13, 6)
    for i in 1:6
        RIN[i, :] = orbitin[1,:] .+ D6[i, :]
    end
    for i in 1:6
        RIN[i+6, :] = orbitin[1,:] .- D6[i, :]
    end
    RIN[13, :] = orbitin
    M66_refpts = zeros(6, 6, length(refpts))
    for i in 1:length(refpts)
		beam = Beam(RIN, energy=E0, mass=m0)
		if i == 1
            id_list = [j for j in 1:refpts[1]]
			ADlinepass!(LATTICE, id_list, beam, changed_idx, changed_ele)
		else
            id_list = [j for j in refpts[i-1]+1:refpts[i]]
			ADlinepass!(LATTICE, id_list, beam, changed_idx, changed_ele)
		end
        TMAT3 = transpose(beam.r)
        M66 = (TMAT3[:,1:6] - TMAT3[:,7:12]) ./ scaling
        M66_refpts[:, :, i] = M66
    end
    return M66_refpts
end

"""
	twissline(tin::EdwardsTengTwiss{Float64},seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, endindex::Int)
Propagate the Twiss parameters through a sequence of elements up to a specified index.
# Arguments
- `tin::EdwardsTengTwiss`: Input Twiss parameters.
- `seq::Vector`: Sequence of elements.
- `dp::Float64`: Momentum deviation.
- `order::Int`: Order of the map. 0 for finite difference, others for TPSA.
- `endindex::Int`: Index of the last element in the sequence.
# Returns
- `EdwardsTengTwiss`: Output Twiss parameters.
"""
function twissline(tin::EdwardsTengTwiss{Float64},seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, endindex::Int; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
	# obtain M through tracking
    ret = tin
    ss = 0.0
	used_seq = seq[1:endindex]
	# M = findm66(used_seq, dp, order)
	if order ==0
		M = fastfindm66(used_seq, dp, E0=E0, m0=m0, orb=orb)
	else
		M = findm66(used_seq, dp, order, E0=E0, m0=m0, orb=orb)
	end
	ret = twissPropagate(ret, M)
	return ret
end

"""	twissline(tin::EdwardsTengTwiss{DTPSAD{N,T}},seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, order::Int, endindex::Int) where {N, T}
Propagate the Twiss parameters through a sequence of DTPSAD elements up to a specified index
# Arguments
- `tin::EdwardsTengTwiss`: Input Twiss parameters.
- `seq::Vector`: Sequence of elements.
- `dp::Float64`: Momentum deviation.
- `order::Int`: Order of the map. Only 0 is supported for DTPSAD elements.
- `endindex::Int`: Index of the last element in the sequence.
# Returns
- `EdwardsTengTwiss`: Output Twiss parameters.
"""
function twissline(tin::EdwardsTengTwiss{DTPSAD{N,T}},seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, order::Int, endindex::Int; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6)) where {N, T}
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
	# obtain M through tracking
    ret = tin
    ss = 0.0
	used_seq = seq[1:endindex]
	# M = findm66(used_seq, dp, order)
	if order ==0
		M = fastfindm66(used_seq, dp, E0=E0, m0=m0, orb=orb)
	else
		error("Order > 0 is not supported for AbstractTPSAElement.")
	end
	ret = twissPropagate(ret, M)
	return ret
end

"""
	twissline(tin::EdwardsTengTwiss{Float64},seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, refpts::Vector{Int})
Propagate the Twiss parameters through a sequence of elements. Save the results at specified locations.
# Arguments
- `tin::EdwardsTengTwiss`: Input Twiss parameters.
- `seq::Vector`: Sequence of elements.
- `dp::Float64`: Momentum deviation.
- `order::Int`: Order of the map. 0 for finite difference, others for TPSA.
- `refpts::Vector{Int}`: Indices of elements where the Twiss parameters are calculated.

# Returns
- `Vector{EdwardsTengTwiss}`: Output Twiss parameters at specified locations.
"""
function twissline(tin::EdwardsTengTwiss{Float64},seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, refpts::Vector{Int}; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
	if refpts[end] > length(seq)
		error("Invalid reference point.")
	end
	# obtain M through tracking
	ret_vector = Vector{EdwardsTengTwiss{Float64}}(undef, length(refpts))
    ret = tin

	# use fastfindm66_refpts instead of findm66_refpts if order == 0
	if order == 0
		M_list = fastfindm66_refpts(seq, dp, refpts, E0=E0, m0=m0, orb=orb)
	else
		M_list = findm66_refpts(seq, dp, order, refpts, E0=E0, m0=m0, orb=orb)
	end
	for i in 1:length(refpts)
		ret = twissPropagate(ret, M_list[:, :, i])
		ret_vector[i] = ret
	end

	return ret_vector
end

"""
	twissline(tin::EdwardsTengTwiss{DTPSAD{N,T}},seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, order::Int, refpts::Vector{Int}) where {N, T}
Propagate the Twiss parameters through a sequence of DTPSAD elements up to a specified index

# Arguments
- `tin::EdwardsTengTwiss`: Input Twiss parameters.
- `seq::Vector`: Sequence of elements.
- `dp::Float64`: Momentum deviation.
- `order::Int`: Order of the map. Only 0 is supported for DTPSAD elements.
- `refpts::Vector{Int}`: Indices of elements where the Twiss parameters are calculated.

# Returns
- `EdwardsTengTwiss`: Output Twiss parameters.
"""
function twissline(tin::EdwardsTengTwiss{DTPSAD{N,T}},seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, order::Int, refpts::Vector{Int}; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6)) where {N, T}
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
	if refpts[end] > length(seq)
		error("Invalid reference point.")
	end
	# obtain M through tracking
	ret_vector = Vector{EdwardsTengTwiss{DTPSAD{N,T}}}(undef, length(refpts))
    ret = tin

	# use fastfindm66_refpts instead of findm66_refpts if order == 0
	if order == 0
		M_list = fastfindm66_refpts(seq, dp, refpts, E0=E0, m0=m0, orb=orb)
	else
		error("Order > 0 is not supported for AbstractTPSAElement.")
	end
	for i in 1:length(refpts)
		ret = twissPropagate(ret, M_list[:, :, i])
		ret_vector[i] = ret
	end

	return ret_vector
end

"""
	ADtwissline(tin::EdwardsTengTwiss,seq::Vector, dp::Float64, order::Int, refpts::Vector{Int}, changed_idx::Vector, changed_ele::Vector)

Propagate the Twiss parameters through a sequence of elements. Save the results at specified locations.
This function is used for automatic differentiation with Enzyme to avoid access issues.

# Arguments
- `tin::EdwardsTengTwiss`: Input Twiss parameters.
- `seq::Vector`: Sequence of elements.
- `dp::Float64`: Momentum deviation.
- `order::Int`: Order of the map. 0 for finite difference, others for TPSA.
- `refpts::Vector{Int}`: Indices of elements where the Twiss parameters are calculated.
- `changed_idx::Vector`: Indices of elements with changed parameters.
- `changed_ele::Vector`: Indices of changed parameters.

# Returns
- `Vector{EdwardsTengTwiss}`: Output Twiss parameters at specified locations.
"""
function ADtwissline(tin::EdwardsTengTwiss{Float64},seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, refpts::Vector{Int}, changed_idx::Vector{Int}, changed_ele::Vector{<:AbstractElement{Float64}}; 
	E0::Float64=3e9, m0::Float64=m_e, orb::Vector{Float64}=zeros(6))
	if dp == 0.0 && orb[6] != 0.0
		dp = orb[6]
	end
	ret_vector = Vector{EdwardsTengTwiss{Float64}}(undef, length(refpts))    
	ret = tin

	# use ADfastfindm66_refpts instead of ADfindm66_refpts if order == 0
	if order == 0
		M_list = ADfastfindm66_refpts(seq, dp, refpts, changed_idx, changed_ele, E0=E0, m0=m0, orb=orb)
	else
		M_list = ADfindm66_refpts(seq, dp, order, refpts, changed_idx, changed_ele, E0=E0, m0=m0, orb=orb)
	end
	for i in 1:length(refpts)
		ret = twissPropagate(ret, M_list[:, :, i])
		ret_vector[i] = ret
	end
	return ret_vector
end

function is_increasing(A::Array)
    # Check if the array has less than two elements
    if length(A) < 2
        return true
    end

    # Initialize the maximum with the first element
    max_val = A[1]

    # Iterate through the array starting from the second element
    for i in 2:length(A)
        # Check if the current element is not greater than the maximum so far
        if A[i] <= max_val
            return false
        end
        # Update the maximum value
        max_val = A[i]
    end

    return true
end


"""
	periodicEdwardsTengTwiss(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int)

Calculate the Twiss parameters for a periodic lattice.

# Arguments
- `seq::Vector{<:AbstractElement{Float64}}`: Sequence of elements.
- `dp::Float64`: Momentum deviation.
- `order::Int`: Order of the map. 0 for finite difference, others for TPSA.

# Returns
- `EdwardsTengTwiss`: Output Twiss parameters.
"""
function periodicEdwardsTengTwiss(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int; 
	E0::Float64=3e9, m0::Float64=m_e)
	orb, _ = find_closed_orbit(seq, dp, mass=m0, energy=E0)
	# M = findm66(seq, dp, order)
	if order == 0
		M = fastfindm66(seq, dp, E0=E0, m0=m0, orb=orb)
	else
		M = findm66(seq, dp, order, E0=E0, m0=m0, orb=orb)
	end
	A= M[1:2,1:2]
	B= M[1:2,3:4]
	C= M[3:4,1:2]
	D= M[3:4,3:4]
	invalid_ret=EdwardsTengTwiss(1.0,1.0,mode=0)

	Bbar_and_C=symplectic_conjugate_2by2(B)+C
	t1=0.5*(tr(A)-tr(D))
	Δ=t1*t1+det_small_matrix(Bbar_and_C)
	Δ<0.0 && (println(stderr,"Failed to decouple periodic transfer matrix. The linear matrix is unstable.");return invalid_ret)

	_sign= t1>0.0 ? Float64(-1) : 1.0

	t2=abs(t1)+sqrt(Δ)
	if t2==0.0
		R=Float64[0 0;0 0]
	else
		R=Bbar_and_C*(_sign/t2)
	end

	X=A-B*R
	Y=D+C*symplectic_conjugate_2by2(R)

	# It should be equal to 1
	(det_small_matrix(X)<Float64(0.9) || det_small_matrix(Y)<Float64(0.9))  && (println(stderr,"Failed to decouple the periodic transfer matrix with mode 1.");return invalid_ret)

	cmux=Float64(0.5)*(X[1,1]+X[2,2])
	cmuy=Float64(0.5)*(Y[1,1]+Y[2,2])
	(Float64(-1)<cmux<1.0 && Float64(-1)<cmuy<1.0) || (println(stderr,"Failed to get beta functions. The linear matrix is unstable.");return invalid_ret)

	smux=sqrt(1.0-cmux*cmux)*sign(X[1,2])
	smuy=sqrt(1.0-cmuy*cmuy)*sign(Y[1,2])
	betax=X[1,2]/smux
	gamx=-X[2,1]/smux
	betay=Y[1,2]/smuy
	gamy=-Y[2,1]/smuy

	alfx=Float64(0.5)*(X[1,1]-X[2,2])/smux
	alfy=Float64(0.5)*(Y[1,1]-Y[2,2])/smuy

	eta=inv1(Matrix{Float64}(I,(4,4))-( M[1:4,1:4]))*( M[1:4,6])
	return EdwardsTengTwiss{Float64}(betax,betay,alfx,alfy,gamx,gamy,eta[1],eta[2],eta[3],eta[4],0.0,0.0,smux,cmux,smuy,cmuy,R,1)
end

"""
	periodicEdwardsTengTwiss(seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, order::Int) where {N, T}
Calculate the Twiss parameters for a periodic lattice with DTPSAD elements.
# Arguments
- `seq::Vector{<:AbstractElement{DTPSAD{N, T}}}`: Sequence of elements.
- `dp::Float64`: Momentum deviation.
- `order::Int`: Order of the map. Only 0 is supported for DTPSAD elements.
# Returns
- `EdwardsTengTwiss`: Output Twiss parameters.
"""
function periodicEdwardsTengTwiss(seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, order::Int; 
	E0::Float64=3e9, m0::Float64=m_e) where {N, T}
	orb, _ = find_closed_orbit(seq, dp, mass=m0, energy=E0)
	# M = findm66(seq, dp, order)
	if order == 0
		M = fastfindm66(seq, dp, E0=E0, m0=m0, orb=orb)
	else
		error("only supports order 0 in periodicEdwardsTengTwiss for AbstractTPSAElement.")
	end
	A = M[1:2,1:2]
	B = M[1:2,3:4]
	C = M[3:4,1:2]
	D = M[3:4,3:4]
	invalid_ret = EdwardsTengTwiss(DTPSAD(1.0), DTPSAD(1.0), mode=0)

	Bbar_and_C=symplectic_conjugate_2by2(B)+C
	t1=0.5*(tr(A)-tr(D))
	Δ=t1*t1+det_small_matrix(Bbar_and_C)
	Δ<0.0 && (println(stderr,"Failed to decouple periodic transfer matrix. The linear matrix is unstable.");return invalid_ret)

	_sign= t1>0.0 ? Float64(-1) : 1.0

	t2=abs(t1)+sqrt(Δ)
	if t2==0.0
		R=Float64[0 0;0 0]
	else
		R=Bbar_and_C*(_sign/t2)
	end

	X=A-B*R
	Y=D+C*symplectic_conjugate_2by2(R)

	# It should be equal to 1
	(det_small_matrix(X)<Float64(0.9) || det_small_matrix(Y)<Float64(0.9))  && (println(stderr,"Failed to decouple the periodic transfer matrix with mode 1.");return invalid_ret)

	cmux=Float64(0.5)*(X[1,1]+X[2,2])
	cmuy=Float64(0.5)*(Y[1,1]+Y[2,2])
	(Float64(-1)<cmux<1.0 && Float64(-1)<cmuy<1.0) || (println(stderr,"Failed to get beta functions. The linear matrix is unstable.");return invalid_ret)

	smux=sqrt(1.0-cmux*cmux)*sign(X[1,2])
	smuy=sqrt(1.0-cmuy*cmuy)*sign(Y[1,2])
	betax=X[1,2]/smux
	gamx=-X[2,1]/smux
	betay=Y[1,2]/smuy
	gamy=-Y[2,1]/smuy

	alfx=Float64(0.5)*(X[1,1]-X[2,2])/smux
	alfy=Float64(0.5)*(Y[1,1]-Y[2,2])/smuy

	eta=inv1(Matrix{Float64}(I,(4,4))-M[1:4,1:4])* M[1:4,6]
	return EdwardsTengTwiss{DTPSAD{N,T}}(betax,betay,alfx,alfy,gamx,gamy,eta[1],eta[2],eta[3],eta[4],DTPSAD(0.0),DTPSAD(0.0),smux,cmux,smuy,cmuy,R,1)
end

function ADperiodicEdwardsTengTwiss(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, 
	changed_idx::Vector{Int}, changed_ele::Vector{<:AbstractElement{Float64}}; 
	E0::Float64=3e9, m0::Float64=m_e)
	orb, _ = find_closed_orbit(seq, dp, mass=m0, energy=E0)
	if order == 0
		M = ADfastfindm66(seq, dp, changed_idx, changed_ele, E0=E0, m0=m0, orb=orb)
	else
		M = ADfindm66(seq, dp, order, changed_idx, changed_ele, E0=E0, m0=m0, orb=orb)
	end
	A= M[1:2,1:2]
	B= M[1:2,3:4]
	C= M[3:4,1:2]
	D= M[3:4,3:4]
	invalid_ret=EdwardsTengTwiss(1.0,1.0,mode=0)

	Bbar_and_C=symplectic_conjugate_2by2(B)+C
	t1=0.5*(tr(A)-tr(D))
	Δ=t1*t1+det_small_matrix(Bbar_and_C)
	if Δ<0.0
		println("Failed to decouple periodic transfer matrix. The linear matrix is unstable.")
		return invalid_ret
	end

	_sign= t1>0.0 ? Float64(-1.0) : 1.0

	t2=abs(t1)+sqrt(Δ)
	if t2==0.0
		R=Float64[0 0;0 0]
	else
		R=Bbar_and_C*(_sign/t2)
	end

	X=A-B*R
	Y=D+C*symplectic_conjugate_2by2(R)

	# It should be equal to 1
	if (det_small_matrix(X)<Float64(0.9) || det_small_matrix(Y)<Float64(0.9))
		println("Failed to decouple the periodic transfer matrix with mode 1.")
		return invalid_ret
	end

	cmux=Float64(0.5)*(X[1,1]+X[2,2])
	cmuy=Float64(0.5)*(Y[1,1]+Y[2,2])

	if !(-1.0<cmux<1.0 && -1.0<cmuy<1.0)
		println("Failed to get beta functions. The linear matrix is unstable.")
		return invalid_ret
	end

	smux=sqrt(1.0-cmux*cmux)*sign(X[1,2])
	smuy=sqrt(1.0-cmuy*cmuy)*sign(Y[1,2])
	betax=X[1,2]/smux
	gamx=-X[2,1]/smux
	betay=Y[1,2]/smuy
	gamy=-Y[2,1]/smuy

	alfx=Float64(0.5)*(X[1,1]-X[2,2])/smux
	alfy=Float64(0.5)*(Y[1,1]-Y[2,2])/smuy
	Identity = zeros(Float64, 4, 4)
	for i in 1:4
		Identity[i, i] = 1.0
	end
	eta=inv1(Identity-( M[1:4,1:4]))*( M[1:4,6])
	return EdwardsTengTwiss{Float64}(betax,betay,alfx,alfy,gamx,gamy,eta[1],eta[2],eta[3],eta[4],0.0,0.0,smux,cmux,smuy,cmuy,R,1)
end

"""
	twissring(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int;
	E0::Float64=3e9, m0::Float64=m_e)

Calculate the periodic Twiss parameters along the ring.

# Arguments
- `seq::Vector`: Sequence of elements.
- `dp::Float64`: Momentum deviation.
- `order::Int`: Order of the map. 0 for finite difference, others for TPSA.
- `E0::Float64`: Beam energy in eV. Default is 3 GeV.
- `m0::Float64`: Particle rest mass in eV/c². Default is electron mass.

# Returns
- `twis`: Twiss parameters along the ring.
"""
function twissring(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int; E0::Float64=3e9, m0::Float64=m_e)
	twi0 = periodicEdwardsTengTwiss(seq, dp, order, E0=E0, m0=m0)
	nele = length(seq)
	refpts = [i for i in 1:nele]
	twis = twissline(twi0, seq, dp, order, refpts, E0=E0, m0=m0)
	return twis
end

"""
	twissring(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, refpts::Vector{Int};
	E0::Float64=3e9, m0::Float64=m_e)
Calculate the periodic Twiss parameters along the ring at specified reference points.
# Arguments
- `seq::Vector`: Sequence of elements.
- `dp::Float64`: Momentum deviation.
- `order::Int`: Order of the map. 0 for finite difference, others for TPSA.
- `refpts::Vector{Int}`: Reference points along the ring.
- `E0::Float64`: Beam energy in eV. Default is 3 GeV.
- `m0::Float64`: Particle rest mass in eV/c². Default is electron mass.

# Returns
- `twis`: Twiss parameters along the ring.
"""
function twissring(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, refpts::Vector{Int}; 
	E0::Float64=3e9, m0::Float64=m_e)
	twi0 = periodicEdwardsTengTwiss(seq, dp, order, E0=E0, m0=m0)
	nele = length(seq)
	twis = twissline(twi0, seq, dp, order, refpts, E0=E0, m0=m0)
	return twis
end

"""
	twissring(seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, order::Int;
	E0::Float64=3e9, m0::Float64=m_e) where {N, T}
Calculate the periodic Twiss parameters along the ring with DTPSAD elements.
# Arguments
- `seq::Vector`: Sequence of elements.
- `dp::Float64`: Momentum deviation.
- `order::Int`: Order of the map. 0 for finite difference, others for TPSA.
- `E0::Float64`: Beam energy in eV. Default is 3 GeV.
- `m0::Float64`: Particle rest mass in eV/c². Default is electron mass.

# Returns
- `twis`: Twiss parameters along the ring.
"""
function twissring(seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, order::Int; 
	E0::Float64=3e9, m0::Float64=m_e) where {N, T}
	twi0 = periodicEdwardsTengTwiss(seq, dp, order, E0=E0, m0=m0)
	nele = length(seq)
	refpts = [i for i in 1:nele]
	twis = twissline(twi0, seq, dp, order, refpts, E0=E0, m0=m0)
	return twis
end

"""
	twissring(seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, order::Int, refpts::Vector{Int};
	E0::Float64=3e9, m0::Float64=m_e) where {N, T}
Calculate the periodic Twiss parameters along the ring with DTPSAD elements at specified reference points.
# Arguments
- `seq::Vector`: Sequence of elements.
- `dp::Float64`: Momentum deviation.
- `order::Int`: Order of the map. 0 for finite difference, others for TPSA.
- `refpts::Vector{Int}`: Reference points along the ring.
- `E0::Float64`: Beam energy in eV. Default is 3 GeV.
- `m0::Float64`: Particle rest mass in eV/c². Default is electron mass.

# Returns
- `twis`: Twiss parameters along the ring.
"""
function twissring(seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64, order::Int, refpts::Vector{Int}; 
	E0::Float64=3e9, m0::Float64=m_e) where {N, T}
	twi0 = periodicEdwardsTengTwiss(seq, dp, order, E0=E0, m0=m0)
	twis = twissline(twi0, seq, dp, order, refpts, E0=E0, m0=m0)
	return twis
end

function ADtwissring(seq::Vector{<:AbstractElement{Float64}}, dp::Float64, order::Int, refpts::Vector{Int}, changed_idx::Vector{Int}, changed_ele::Vector{<:AbstractElement{Float64}}; 
	E0::Float64=3e9, m0::Float64=m_e)
	twi0 = ADperiodicEdwardsTengTwiss(seq, dp, order, changed_idx, changed_ele, E0=E0, m0=m0)
	nele = length(seq)
	# refpts = [i for i in 1:nele]
	twi = ADtwissline(twi0, seq, dp, order, refpts, changed_idx, changed_ele, E0=E0, m0=m0)
	return twi
end

function normalMatrix(tin::EdwardsTengTwiss{Float64})
	(tin.mode==1 || tin.mode==Int(2)) || begin
		println(stderr,"Warning: return identity matrix for unknown mode $(tin.mode) as the normal matrix (transformation matrix from normal space to physical space).")
		return 1.0*Matrix{Float64}(I,6,6)
	end
	D=[1.0 0.0 0.0 0.0 0.0 tin.dx
		0.0 1.0 0.0 0.0 0.0 tin.dpx
		0.0 0.0 1.0 0.0 0.0 tin.dy
		0.0 0.0 0.0 1.0 0.0 tin.dpy
		-tin.dpx tin.dx -tin.dpy tin.dy 1.0 0.0
		0.0 0.0 0.0 0.0 0.0 1.0]
	sbx = sqrt(tin.betax)
	sby = sqrt(tin.betay)
	B=[sbx 0.0 0.0 0.0 0.0 0.0
		-tin.alphax/sbx 1/sbx 0.0 0.0 0.0 0.0
		0.0 0.0 sby 0.0 0.0 0.0
		0.0 0.0 -tin.alphay/sby 1.0/sby 0.0 0.0
		0.0 0.0 0.0 0.0 1.0 0.0
		0.0 0.0 0.0 0.0 0.0 1.0]
	λ=1.0/sqrt(abs(1.0+det_small_matrix(tin.R)))
	R=λ*tin.R
	_R=symplectic_conjugate_2by2(R)
	O=[0.0 0.0;0.0 0.0]
	U=[λ 0.0;0.0 λ]
	if tin.mode==1
		V=[U _R O;-R U O;O O I]
	else
		V=[_R U O;U -R O;O O I]
	end
	return D*V*B
end

