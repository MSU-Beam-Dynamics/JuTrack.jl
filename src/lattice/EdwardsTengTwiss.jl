using LinearAlgebra # LinearAlgebra is not fully supported by Enzyme. E.g., det, eigen, inv, etc. are not supported.
abstract type AbstractTwiss end

function det1(A::Matrix{Float64})
    # Ensure the matrix is square
    if size(A, 1) != size(A, 2)
        error("Matrix must be square")
    end

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

function inv1(A::Matrix{Float64})
    nrows, ncols = size(A)
    if nrows != ncols
        error("Matrix must be square")
    end

    augmented_matrix = [A I]

    # Perform Gauss-Jordan elimination
    for i in 1:nrows
        # Pivot must be non-zero. If zero, try to swap with another row.
        if augmented_matrix[i, i] == 0
            for j in (i+1):nrows
                if augmented_matrix[j, i] != 0
                    augmented_matrix[i, :], augmented_matrix[j, :] = augmented_matrix[j, :], augmented_matrix[i, :]
                    break
                end
            end
        end

        # Make pivot 1 and eliminate all other entries in the column
        pivot = augmented_matrix[i, i]
        augmented_matrix[i, :] /= pivot
        for j in 1:nrows
            if i != j
                factor = augmented_matrix[j, i]
                augmented_matrix[j, :] -= factor * augmented_matrix[i, :]
            end
        end
    end

    invA = augmented_matrix[:, ncols+1:end]

    return invA
end


struct EdwardsTengTwiss <: AbstractTwiss
	betax::Float64
	betay::Float64
	alphax::Float64
	alphay::Float64
	gammax::Float64
	gammay::Float64
	dx::Float64
	dpx::Float64
	dy::Float64
	dpy::Float64
	dmux::Float64
	dmuy::Float64
	sinmux::Float64
	cosmux::Float64
	sinmuy::Float64
	cosmuy::Float64
	R::Matrix{Float64}
	mode::Int
end



EdwardsTengTwiss(;betax::Float64,betay::Float64,
			   alfx::Float64=0.0,alfy::Float64=0.0,
			   dx::Float64=0.0,dy::Float64=0.0,
			   dpx::Float64=0.0,dpy::Float64=0.0,
			   mux::Float64=0.0,muy::Float64=0.0,
			   R11::Float64=0.0,R12::Float64=0.0,
			   R21::Float64=0.0,R22::Float64=0.0,
			   mode::Int=1)=EdwardsTengTwiss(betax,betay,alfx,alfy,(1.0+alfx^2)/betax,(1.0+alfy^2)/betay,
														  dx,dpx,dy,dpy,mux,muy,sin(mux),cos(mux),sin(muy),cos(muy),[R11 R12;R21 R22],mode)
_symplectic_conjugate_2by2(M) = [M[2, 2] -M[1, 2]; -M[2, 1] M[1, 1]]
function _matrixTransform_2by2(M)
	m11 = M[1, 1]
	m21 = M[2, 1]
	m12 = M[1, 2]
	m22 = M[2, 2]
	return [m11*m11 -2m11*m12 m12*m12
	-m11*m21 1.0+2m12*m21 -m12*m22
	m21*m21 -2m21*m22 m22*m22]
end

function twissPropagate(tin::EdwardsTengTwiss,M::Matrix{Float64})
	A=@view M[1:2,1:2]
	B=@view M[1:2,3:4]
	C=@view M[3:4,1:2]
	D=@view M[3:4,3:4]

	R1=tin.R
	_R1=_symplectic_conjugate_2by2(R1)
	if tin.mode == 1
		X=A-B*R1
		t=det1(X)
			if t>0.1
				R=(D*R1-C)*_symplectic_conjugate_2by2(X)
				R/=t
				X/=sqrt(t)
				Y=D+C*_R1
				Y/=sqrt(det1(Y))
				mode=1
			else
				X=C-D*R1
				X/=sqrt(det1(X))
				Y=B+A*_R1
				t=det1(Y)
				R=-(D+C*_R1)*_symplectic_conjugate_2by2(Y)
				R/=t
				Y/=sqrt(t)
				mode=2
			end
	elseif tin.mode == Int(2) 
		X=B+A*_R1
		t=det1(X)
			if t>0.1
				R=-(D+C*_R1)*_symplectic_conjugate_2by2(X)
				R/=t
				X/=sqrt(t)
				Y=C-D*R1
				Y/=sqrt(det1(Y))
				mode=1
			else
				X=D+C*_R1
				X/=sqrt(det1(X))
				Y=A-B*R1
				t=det1(Y)
				R=(D*R1-C)*_symplectic_conjugate_2by2(Y)
				R/=t
				Y/=sqrt(t)
				mode=2
			end
	else
		#throw(AssertionError("Mode should be integer 1 or 2."))
		# println(stderr,"Invalid mode.")
		# return EdwardsTengTwiss(;betax=1.0,betay=1.0,mode=0)
		error("Invalid mode.")
	end

	Nx=_matrixTransform_2by2(X)
	Ny=_matrixTransform_2by2(Y)
	v1=Nx*[tin.betax;tin.alphax;tin.gammax]
	v2=Ny*[tin.betay;tin.alphay;tin.gammay]
	eta=(@view M[1:4,1:4])*[tin.dx,tin.dpx,tin.dy,tin.dpy]+(@view M[1:4,6])
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

	new_mux = tin.dmux + delta_mux
	new_muy = tin.dmuy + delta_muy
	return EdwardsTengTwiss(v1[1],v2[1],v1[2],v2[2],v1[3],v2[3],eta[1],eta[2],eta[3],eta[4],new_mux,new_muy,smux,cmux,smuy,cmuy,R,mode)
end

function findm66(seq::Vector{AbstractElement}, dp::Float64, order::Int)
	x = CTPS(0.0, 1, 6, order)
	px = CTPS(0.0, 2, 6, order)
	y = CTPS(0.0, 3, 6, order)
	py = CTPS(0.0, 4, 6, order)
	z = CTPS(dp, 5, 6, order)
	delta = CTPS(0.0, 6, 6, order)
	rin = [x, px, y, py, z, delta]
	# no radiation, cavity off
	linepass_TPSA!(seq, rin)

	map = zeros(Float64, 6, 6)
	for i in 1:6
		for j in 1:6
			map[i, j] = rin[i].map[j + 1]
		end
	end

	# Map66 = [x.map[2] x.map[3] x.map[4] x.map[5] x.map[6] x.map[7];
	# 		px.map[2] px.map[3] px.map[4] px.map[5] px.map[6] px.map[7];
	# 		y.map[2] y.map[3] y.map[4] y.map[5] y.map[6] y.map[7];
	# 		py.map[2] py.map[3] py.map[4] py.map[5] py.map[6] py.map[7];
	# 		delta.map[2] delta.map[3] delta.map[4] delta.map[5] delta.map[6] delta.map[7];
	# 		z.map[2] z.map[3] z.map[4] z.map[5] z.map[6] z.map[7]]
	return map
end
function ADfindm66(seq::Vector{AbstractElement}, dp::Float64, order::Int, changed_idx::Vector{Int}, changed_ele)
	x = CTPS(0.0, 1, 6, order)
	px = CTPS(0.0, 2, 6, order)
	y = CTPS(0.0, 3, 6, order)
	py = CTPS(0.0, 4, 6, order)
	z = CTPS(dp, 5, 6, order)
	delta = CTPS(0.0, 6, 6, order)
	rin = [x, px, y, py, z, delta]
	# no radiation, cavity off
	ADlinepass_TPSA!(seq, rin, changed_idx, changed_ele)

	map = zeros(Float64, 6, 6)
	for i in 1:6
		for j in 1:6
			map[i, j] = rin[i].map[j + 1]
		end
	end
	return map
end

function findm66_refpts(seq::Vector{AbstractElement}, dp::Float64, order::Int, refpts::Vector{Int})
	x = CTPS(0.0, 1, 6, order)
	px = CTPS(0.0, 2, 6, order)
	y = CTPS(0.0, 3, 6, order)
	py = CTPS(0.0, 4, 6, order)
	z = CTPS(dp, 5, 6, order)
	delta = CTPS(0.0, 6, 6, order)
	rin = [x, px, y, py, z, delta]
	# no radiation, cavity off
	map_list = Vector{Matrix{Float64}}(undef, length(refpts))
	num = 1
	for i in eachindex(refpts)
		if i == 1
			linepass_TPSA!(seq[1:refpts[i]], rin)
		else
			linepass_TPSA!(seq[refpts[i-1]+1:refpts[i]], rin)
		end
		map = zeros(Float64, 6, 6)
		for j in 1:6
			for k in 1:6
				map[j, k] = rin[j].map[k + 1]
			end
		end
		map_list[num] = map
		num += 1
		reassign!(rin[1], 0.0, 1)
		reassign!(rin[2], 0.0, 2)
		reassign!(rin[3], 0.0, 3)
		reassign!(rin[4], 0.0, 4)
		reassign!(rin[5], dp, 5)
		reassign!(rin[6], 0.0, 6)
	end

	return map_list
end

function Twissline(tin::EdwardsTengTwiss,seq::Vector{AbstractElement}, dp::Float64, order::Int, endindex::Int)
	# obtain M through tracking
    ret = tin
    ss = 0.0
	used_seq = seq[1:endindex]
	M = findm66(used_seq, dp, order)
	ret = twissPropagate(ret, M)
	ss = sum([mag.len for mag in used_seq])
	names = [mag.name for mag in used_seq]
	return ss,names,ret
end

function ADTwissline(tin::EdwardsTengTwiss,seq::Vector{AbstractElement}, dp::Float64, order::Int, endindex::Int, changed_idx::Vector{Int}, changed_ele)
	# obtain M through tracking
    ret = tin
    ss = 0.0
	used_seq = seq[1:endindex]
	M = ADfindm66(used_seq, dp, order, changed_idx, changed_ele)
	ret = twissPropagate(ret, M)
	ss = sum([mag.len for mag in used_seq])
	names = [mag.name for mag in used_seq]
	return ss,names,ret
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

function Twissline(tin::EdwardsTengTwiss,seq::Vector{AbstractElement}, dp::Float64, order::Int, refpts::Vector{Int})
	if !is_increasing(refpts)
		error("The reference points must be in increasing order.")
	end
	if refpts[end] > length(seq)
		error("Invalid reference point.")
	end
	# obtain M through tracking
	ret_vector = Vector{EdwardsTengTwiss}(undef, length(refpts))
    ret = tin

	M_list = findm66_refpts(seq, dp, order, refpts)
	for i in 1:length(refpts)
		ret = twissPropagate(ret, M_list[i])
		ret_vector[i] = ret
	end

	ss = zeros(Float64, length(refpts))
	names = Vector{String}(undef, length(refpts))
	for i in 1:length(refpts)
		ss[i] = sum([mag.len for mag in seq[1:refpts[i]]])
		names[i] = seq[refpts[i]].name
	end
	return ss,names,ret_vector
end

function periodicEdwardsTengTwiss(seq::Vector{AbstractElement}, dp, order::Int)
	M = findm66(seq, dp, order)
	A=@view M[1:2,1:2]
	B=@view M[1:2,3:4]
	C=@view M[3:4,1:2]
	D=@view M[3:4,3:4]
	invalid_ret=EdwardsTengTwiss(;betax=1.0,betay=1.0,mode=0)

	Bbar_and_C=_symplectic_conjugate_2by2(B)+C
	t1=0.5*(tr(A)-tr(D))
	Δ=t1*t1+det1(Bbar_and_C)
	Δ<0.0 && (println(stderr,"Failed to decouple periodic transfer matrix. The linear matrix is unstable.");return invalid_ret)

	_sign= t1>0.0 ? Float64(-1) : 1.0

	t2=abs(t1)+sqrt(Δ)
	if t2==0.0
		R=Float64[0 0;0 0]
	else
		R=Bbar_and_C*(_sign/t2)
	end

	X=A-B*R
	Y=D+C*_symplectic_conjugate_2by2(R)

	# It should be equal to 1
	(det1(X)<Float64(0.9) || det1(Y)<Float64(0.9))  && (println(stderr,"Failed to decouple the periodic transfer matrix with mode 1.");return invalid_ret)

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

	eta=inv1(Matrix{Float64}(I,(4,4))-(@view M[1:4,1:4]))*(@view M[1:4,6])
	return EdwardsTengTwiss(betax,betay,alfx,alfy,gamx,gamy,eta[1],eta[2],eta[3],eta[4],0.0,0.0,smux,cmux,smuy,cmuy,R,1)
end

function ADperiodicEdwardsTengTwiss(seq::Vector{AbstractElement}, dp, order::Int, changed_idx::Vector{Int}, changed_ele)
	M = ADfindm66(seq, dp, order, changed_idx, changed_ele)
	A=@view M[1:2,1:2]
	B=@view M[1:2,3:4]
	C=@view M[3:4,1:2]
	D=@view M[3:4,3:4]
	invalid_ret=EdwardsTengTwiss(;betax=1.0,betay=1.0,mode=0)

	Bbar_and_C=_symplectic_conjugate_2by2(B)+C
	t1=0.5*(tr(A)-tr(D))
	Δ=t1*t1+det1(Bbar_and_C)
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
	Y=D+C*_symplectic_conjugate_2by2(R)

	# It should be equal to 1
	if (det1(X)<Float64(0.9) || det1(Y)<Float64(0.9))
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
	eta=inv1(Identity-(@view M[1:4,1:4]))*(@view M[1:4,6])
	return EdwardsTengTwiss(betax,betay,alfx,alfy,gamx,gamy,eta[1],eta[2],eta[3],eta[4],0.0,0.0,smux,cmux,smuy,cmuy,R,1)
end

function twissring(seq::Vector{AbstractElement}, dp::Float64, order::Int)
	twi0 = periodicEdwardsTengTwiss(seq, dp, order)
	nele = length(seq)
	refpts = [i for i in 1:nele]
	ss, name, twi = Twissline(twi0, seq, dp, order, refpts)
	return twi, ss
end

function ADtwissring(seq::Vector{AbstractElement}, dp::Float64, order::Int, changed_idx::Vector{Int}, changed_ele)
	twi0 = ADperiodicEdwardsTengTwiss(seq, dp, order, changed_idx, changed_ele)
	nele = length(seq)
	refpts = [i for i in 1:nele]
	ss, name, twi = ADTwissline(twi0, seq, dp, order, nele, changed_idx, changed_ele)
	return twi, ss
end

function normalMatrix(tin::EdwardsTengTwiss)
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
	λ=1.0/sqrt(abs(1.0+det1(tin.R)))
	R=λ*tin.R
	_R=_symplectic_conjugate_2by2(R)
	O=[0.0 0.0;0.0 0.0]
	U=[λ 0.0;0.0 λ]
	if tin.mode==1
		V=[U _R O;-R U O;O O I]
	else
		V=[_R U O;U -R O;O O I]
	end
	return D*V*B
end

