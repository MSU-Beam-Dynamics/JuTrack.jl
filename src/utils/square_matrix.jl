# A code for extending the convergence map method based on square matrix method proposed by Li Hua Yu to 6-D space.
# Code created by Jinyu Wan, 03/24/2025.
using LinearAlgebra
using FFTW
using StaticArrays
# using TimerOutputs
# const to = TimerOutput()
# using Infiltrator
function uni_transform(ind_i::Int, ind_j::Int, ratio_left::ComplexF64, ratio_right::ComplexF64,
                       mat::Matrix{ComplexF64}, leftm::Matrix{ComplexF64}, rightm::Matrix{ComplexF64})
    rightm[:, ind_j] .+= rightm[:, ind_i] * ratio_right
    leftm[ind_i, :] .+= leftm[ind_j, :] * ratio_left
    temp1 = mat[ind_j, :] * ratio_left
    temp2 = mat[:, ind_i] * ratio_right
    temp3 = mat[ind_j, ind_i] * ratio_left * ratio_right
    mat[ind_i, :] .+= temp1
    mat[:, ind_j] .+= temp2
    mat[ind_i, ind_j] += temp3
    return nothing
end

function jordan_norm_form(sqrmat_ori::Matrix{ComplexF64}, ε::Float64=0.0)
    mat = copy(sqrmat_ori)
    dim = size(mat,1)
    size(mat,2) == dim || error("Matrix must be square")

    leftm = Matrix{ComplexF64}(I, dim, dim)
    rightm = Matrix{ComplexF64}(I, dim, dim)

    # First superdiagonal
    if abs(mat[1,2]) > ε
        siup = mat[1,2]
        uni_transform(1, 1, 1/siup - 1.0, siup - 1.0, mat, leftm, rightm)
    end

    for i in 3:dim
        nonzero = Int[]
        for j in 1:(i-1)
            if abs(mat[j,i]) > ε
                if j < i-1 && abs(mat[j,j+1] - 1) < 1e-10
                    uni_transform(j+1, i, mat[j,i], -mat[j,i], mat, leftm, rightm)
                else
                    push!(nonzero, j)
                end
            end
        end
        if isempty(nonzero)
            continue
        elseif length(nonzero) == 1
            j = nonzero[1]
            ratio = mat[j,i]
            uni_transform(i, i, ratio - 1, 1/ratio - 1, mat, leftm, rightm)
            js = j
        else
            function chain_len(k::Int)
                len = 1; kk = k
                while kk > 1 && abs(mat[kk-1,kk]) > ε
                    len += 1; kk -= 1
                end
                return len
            end
            je = nonzero[1]; js = -1
            for idx in 2:length(nonzero)
                js = nonzero[idx]
                le, ls = chain_len(je), chain_len(js)
                if le > ls
                    je, js = js, je; le, ls = ls, le
                end
                ratio = -mat[je,i]/mat[js,i]
                for k in 0:(le-1)
                    uni_transform(je-k, js-k, ratio, -ratio, mat, leftm, rightm)
                end
                je = js
            end
            ratio = mat[js,i]
            uni_transform(i, i, ratio - 1, 1/ratio - 1, mat, leftm, rightm)
        end
        if js < i-1
            perm = Matrix{ComplexF64}(I, dim, dim)
            perm[js+2:i, :] = perm[js+1:i-1, :]
            perm[js+1, :] = copy(perm[i, :])
            rightm .= rightm * transpose(perm)
            leftm .= perm * leftm
            mat .= perm * mat * transpose(perm)
        end
    end

    return leftm, mat, rightm
end

function jordan_chain_structure(jnf::Matrix{ComplexF64}, ε::Float64)
    dim = size(jnf,1)
    chains = Vector{Vector{Int}}()
    chain = [1]
    for i in 2:dim
        if abs(jnf[i-1,i]) <= ε
            push!(chains, chain)
            chain = [i]
        else
            push!(chain, i)
        end
    end
    push!(chains, chain)
    return chains
end

mutable struct TPSVar
    dim::Int
    order::Int
    vars::Vector
    sqrmat_dim::Int
    target::Int
    degenerate_list::Vector{Int}
    sqr_mat::Matrix{ComplexF64}
    sqr_mat_omegaI::Matrix{ComplexF64}
    left_mat::Matrix{ComplexF64}
    right_mat::Matrix{ComplexF64}
    reduced_mat::Matrix{ComplexF64}
    right_vector::Matrix{ComplexF64}
    left_vector::Matrix{ComplexF64}
    jnf_leftm::Matrix{ComplexF64}
    jnf_mat::Matrix{ComplexF64}
    jnf_rightm::Matrix{ComplexF64}
    epsilon::Float64
end

function TPSVar6D(order::Int)
    dim = 3
    init = 0.0 + 0.0im
    zx = CTPS(init, 1, 6, order)
    zxs  = CTPS(init, 2, 6, order)
    zy = CTPS(init, 3, 6, order)
    zys  = CTPS(init, 4, 6, order)
    zz = CTPS(init, 5, 6, order)
    zzs  = CTPS(init, 6, 6, order)

    vars = [zx, zxs, zy, zys, zz, zzs]
    d = length(zx.map)
    I_mat = Matrix{ComplexF64}(I, d, d)
    return TPSVar(dim, order, vars, d, 0, Int[], zeros(ComplexF64,d,d), 
    zeros(ComplexF64,d,d), I_mat, I_mat, zeros(ComplexF64,0,0), zeros(ComplexF64,0,0), 
    zeros(ComplexF64,0,0), zeros(ComplexF64,0,0), zeros(ComplexF64,0,0), zeros(ComplexF64,0,0), 1e-14)
end

function TPSVar4D(order::Int)
    dim = 2
    init = 0.0 + 0.0im
    zx = CTPS(init, 1, 4, order)
    zxs  = CTPS(init, 2, 4, order)
    zy = CTPS(init, 3, 4, order)
    zys  = CTPS(init, 4, 4, order)

    vars = [zx, zxs, zy, zys]
    d = length(zx.map)
    I_mat = Matrix{ComplexF64}(I, d, d)
    return TPSVar(dim, order, vars, d, 0, Int[], zeros(ComplexF64,d,d), 
    zeros(ComplexF64,d,d), I_mat, I_mat, zeros(ComplexF64,0,0), zeros(ComplexF64,0,0), 
    zeros(ComplexF64,0,0), zeros(ComplexF64,0,0), zeros(ComplexF64,0,0), zeros(ComplexF64,0,0), 1e-14)
end

function get_degenerate_list!(tpsvar::TPSVar, target::Int; resonance::Union{Nothing, Vector{Int}}=nothing)
    tpsvar.target = target
    if target > tpsvar.dim || target <= 0
        println("Wrong target specified, should be one of the $(tpsvar.dim) dimensions")
        return
    end

    powerlist = [getindexmap(tpsvar.vars[1].polymap[], i) for i in 1:tpsvar.sqrmat_dim]
    powerlist = hcat(powerlist...)'
    fns = zeros(Int, tpsvar.sqrmat_dim, tpsvar.dim)

    for i in 1:tpsvar.dim
        fns[:, i] .= powerlist[:, 2 * i] .- powerlist[:, 2 * i + 1]
    end
    fns[:, target] .-= 1

    if resonance !== nothing
        res = Int.(resonance)
        if length(resonance) != tpsvar.dim
            println("Resonance dimension $(length(resonance)) is not compatible with the variable dimension $(tpsvar.dim)")
            return
        end
        if any(res .!= 0)
            absres2 = sum(res .^ 2)
            absfn2 = sum(fns .^ 2, dims=2)
            fndotres = fns * res
            sign = sign.(fndotres)
            mask1 = absfn2 .> 0
            mask2 = (fndotres .^ 2) .== absres2 .* absfn2
            mask = mask1 .& mask2

            ratio = sign[mask] .* floor.(sqrt.(absfn2[mask] ./ absres2))
            fns[mask, :] .-= ratio .* res'
        end
    end

    mask = .!any(fns .!= 0, dims=2)
    mask_id = [i for i in 1:tpsvar.sqrmat_dim if mask[i]]
    tpsvar.degenerate_list = mask_id
end

function get_variables(tpsvar::TPSVar; betax=1.0, betay=1.0, alphax=0.0, alphay=0.0, betaz=1.0, alphaz=0.0)
    vars_real = []
    for i in 1:tpsvar.dim
        if i == 1
            push!(vars_real, (tpsvar.vars[2*i-1] + tpsvar.vars[2*i]) / 2.0 * sqrt(betax))
            push!(vars_real, (tpsvar.vars[2*i-1] - tpsvar.vars[2*i]) / (-2.0im) /sqrt(betax) -
            (tpsvar.vars[2*i-1] + tpsvar.vars[2*i]) / 2.0 * alphax / sqrt(betax))
        elseif i == 2
            push!(vars_real, (tpsvar.vars[2*i-1] + tpsvar.vars[2*i]) / 2.0 * sqrt(betay))
            push!(vars_real, (tpsvar.vars[2*i-1] - tpsvar.vars[2*i]) / (-2.0im) /sqrt(betay) -
            (tpsvar.vars[2*i-1] + tpsvar.vars[2*i]) / 2.0 * alphay / sqrt(betay))
        elseif i == 3
            push!(vars_real, (tpsvar.vars[2*i-1] + tpsvar.vars[2*i]) / 2.0 * sqrt(betaz))
            push!(vars_real, (tpsvar.vars[2*i-1] - tpsvar.vars[2*i]) / (-2.0im) /sqrt(betaz) -
            (tpsvar.vars[2*i-1] + tpsvar.vars[2*i]) / 2.0 * alphaz / sqrt(betaz)) 
        else
            error("The dimension is either 2 or 3")  
        end
    end
    return vars_real
end

function construct_sqr_matrix(tpsvar::TPSVar, periodic_map)
    if length(periodic_map) != 2 * tpsvar.dim
        println("The length of the periodic map is not compatible with the variable dimension")
        return
    end
    sqrmat = zeros(ComplexF64, tpsvar.sqrmat_dim, tpsvar.sqrmat_dim)
    for i in 1:tpsvar.sqrmat_dim
        pind = getindexmap(tpsvar.vars[1].polymap[], i)
        newvar = CTPS(1.0 + 0.0im, 2 * tpsvar.dim, tpsvar.order)
        for j in 1:length(periodic_map)
            for k in 1:(pind[j + 1])
                newvar *= periodic_map[j]
            end
        end
        indices = newvar.map
        sqrmat[i, 1:length(indices)] .= indices
    end
    tpsvar.sqr_mat = sqrmat
    return nothing
end

function jordan_form!(tpsvar::TPSVar)
    if length(tpsvar.degenerate_list) == 1
        tpsvar.right_vector = reshape(tpsvar.right_mat[:, tpsvar.degenerate_list[1]],1, :)
        tpsvar.left_vector  = reshape(tpsvar.left_mat[tpsvar.degenerate_list[1], :],1, :)
        return nothing
    end

    tpsvar.jnf_leftm, tpsvar.jnf_mat, tpsvar.jnf_rightm = jordan_norm_form(tpsvar.reduced_mat, tpsvar.epsilon)
    if tpsvar.epsilon > 0
        digits = Int(-log10(tpsvar.epsilon))
        tpsvar.jnf_leftm .= round.(tpsvar.jnf_leftm, digits=digits)
        tpsvar.jnf_rightm .= round.(tpsvar.jnf_rightm, digits=digits)
        tpsvar.jnf_mat .= round.(tpsvar.jnf_mat, digits=digits)
    end

    temp_mat = zeros(ComplexF64, tpsvar.sqrmat_dim, length(tpsvar.degenerate_list))
    for i in 1:length(tpsvar.degenerate_list)
        temp_mat[tpsvar.degenerate_list[i], :] .= tpsvar.jnf_rightm[i, :]
    end
    tpsvar.right_vector = tpsvar.right_mat * temp_mat

    temp_mat = zeros(ComplexF64, length(tpsvar.degenerate_list), tpsvar.sqrmat_dim)
    for i in 1:length(tpsvar.degenerate_list)
        temp_mat[:, tpsvar.degenerate_list[i]] .= tpsvar.jnf_leftm[:, i]
    end
    tpsvar.left_vector = temp_mat * tpsvar.left_mat
    return nothing
end

function sqrmat_reduction!(tpsvar)
    dim = tpsvar.sqrmat_dim
    # subtract ωI
    tpsvar.right_mat = Matrix{ComplexF64}(I, dim, dim)
    tpsvar.left_mat  = Matrix{ComplexF64}(I, dim, dim)
    deg_idxs = tpsvar.degenerate_list
    idx = tpsvar.target * 2
    tpsvar.sqr_mat_omegaI = tpsvar.sqr_mat .- tpsvar.sqr_mat[idx, idx] * Matrix{ComplexF64}(I, dim, dim)
    # Pass 1: eliminate above‑diagonal entries
    for ind_i in reverse(deg_idxs)
        for ind_j in (ind_i+1):dim
            val = tpsvar.sqr_mat_omegaI[ind_i, ind_j]
            if abs(val) ≤ tpsvar.epsilon
                tpsvar.sqr_mat_omegaI[ind_i, ind_j] = 0
                continue
            elseif ind_j in deg_idxs
                continue
            end

            ratio = -val / tpsvar.sqr_mat_omegaI[ind_j, ind_j]
            uni_transform(ind_i, ind_j, ratio, -ratio,
                           tpsvar.sqr_mat_omegaI, tpsvar.left_mat, tpsvar.right_mat)
            if abs(ratio) > 1/tpsvar.epsilon
                @warn "Too large ratio observed in reducing function pass 1, ratio=$ratio"
            end
        end
    end

    # Pass 2: eliminate below‑diagonal entries
    for ind_j in deg_idxs
        for ind_i in (ind_j-1):-1:1
            val = tpsvar.sqr_mat_omegaI[ind_i, ind_j]
            if abs(val) ≤ tpsvar.epsilon
                tpsvar.sqr_mat_omegaI[ind_i, ind_j] = 0
                continue
            elseif ind_i in deg_idxs
                continue
            end

            ratio = val / tpsvar.sqr_mat_omegaI[ind_i, ind_i]
            uni_transform(ind_i, ind_j, ratio, -ratio,
                           tpsvar.sqr_mat_omegaI, tpsvar.left_mat, tpsvar.right_mat)
            if abs(ratio) > 1/tpsvar.epsilon
                @warn "Too large ratio observed in reducing function pass 2, ratio=$ratio"
            end
        end
    end

    # Extract the reduced submatrix
    tpsvar.reduced_mat = tpsvar.sqr_mat_omegaI[deg_idxs, deg_idxs]

    # Round everything to avoid floating‑point noise
    if tpsvar.epsilon > 0
        digits = Int(-log10(tpsvar.epsilon))
        tpsvar.reduced_mat .= round.(tpsvar.reduced_mat, digits=digits)
        tpsvar.left_mat    .= round.(tpsvar.left_mat,    digits=digits)
        tpsvar.right_mat   .= round.(tpsvar.right_mat,   digits=digits)
    end

    return nothing
end

function evaluate(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    return function U(inivalue)
        if length(inivalue) != TPS_Dim
            error("Inconsistent dimension to evaluate CTPS")
        end
        sum = ctps.map[1]
        for i in 2:ctps.terms
            temp = getindexmap(ctps.polymap[], i)
            product = 1.0 + 0.0im
            for j in 1:TPS_Dim
                dimpower = inivalue[j]^temp[j+1]
                product *= dimpower
            end
            sum += product * ctps.map[i]
        end
        return sum
    end
end

struct PolyEvaluator{T,D}
    coeffs::Vector{T}
    exps::Vector{NTuple{D,Int}}
end

function PolyEvaluator(ctps::CTPS{T,D,M}) where {T,D,M}
    nterms = ctps.terms
    coeffs = Vector{T}(undef, nterms)
    exps   = Vector{NTuple{D,Int}}(undef, nterms)

    for i in 1:nterms
        coeffs[i] = ctps.map[i]
        idx = getindexmap(ctps.polymap[], i)
        exps[i] = tuple(idx[2:end]...)  # drop total‑degree slot
    end

    return PolyEvaluator{T,D}(coeffs, exps)
end

@inline function evaluate(pe::PolyEvaluator{T,D}, x) where {T,D}
    s = pe.coeffs[1]
    @inbounds for i in 2:length(pe.coeffs)
        prod = one(T)
        e = pe.exps[i]
        @inbounds for j in 1:D
            prod *= x[j]^e[j]
        end
        s += prod * pe.coeffs[i]
    end
    return s
end

function compute_action_angle_polynomials(tpsvar::TPSVar, dim::Int, order::Int)
    # Perform the Jordan decomposition and compute the polynomial forms
    # for the approximate action-angle variables.
    vectors = []
    for i in 1:dim
        get_degenerate_list!(tpsvar, i)
        println("Degenerate indices (set $i):", tpsvar.degenerate_list)
        sqrmat_reduction!(tpsvar)
        jordan_form!(tpsvar)
        leftvect = tpsvar.left_vector[1,:]
        jordan_chain_structure(tpsvar.jnf_mat, 1e-15)

        push!(vectors, leftvect)
    end
    return vectors
end

function linear(f::CTPS{T,D,M}) where {T,D,M}
    lin = CTPS(f)                        # copy
    for i in eachindex(lin.map)
        idx = getindexmap(lin.polymap[], i)
        lin.map[i] = sum(idx[2:end]) == 1 || i == 1 ? lin.map[i] : zero(T)
    end
    return lin
end
function compose(f::CTPS{T,D,M}, xs::Vector{CTPS{T,D,M}}) where {T,D,M}
    result = CTPS(zero(T), D, M)
    for i in eachindex(f.map)
        c = f.map[i]
        c == zero(T) && continue
        idx = getindexmap(f.polymap[], i)
        term = CTPS(c, D, M)
        for j in 1:D
            if idx[j+1] != 0
                term *= xs[j]^idx[j+1]
            end
        end
        result += term
    end
    return result
end
function inverse_map(fs::Vector{CTPS{T,D,M}}, order,
    initial::Union{Nothing,Vector{CTPS{T,D,M}}}=nothing,
    iteration_limit::Int=0) where {T,D,M}

    n = length(fs)
    # dimension check
    # if get(TPS_Dim, fs[1]) != n
    # error("Dimension mismatch: length(fs) = $n, TPS dim = $(get(TPS_Dim,fs[1]))")
    # end

    # Build linear Jacobian and nonlinear remainder
    A = zeros(T,n,n)
    nl = Vector{CTPS{T,D,M}}(undef,n)
    for i in 1:n
        li = linear(fs[i])
        A[i,:] .= li.map[2:n+1]
        nl[i] = fs[i] - li
    end

    Ainv = zeros(T,n,n)
    try 
        Ainv = inv(A)
    catch e
        println("Matrix is singular")
        return nothing
    end
    Iv   = [ CTPS(zero(T), i, D, M) for i in 1:n ]
    res  = isnothing(initial) ? copy(Iv) : initial

    max_iter = iteration_limit > 0 ? iteration_limit : order
    for _ in 1:max_iter
        temp = [ Iv[i] - compose(nl[i], res) for i in 1:n ]
        res  = [ sum(Ainv[i,j]*temp[j] for j in 1:n) for i in 1:n ]
    end
    return res
end

function conjugate(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, mode::Int=1) where {T, TPS_Dim, Max_TPS_Degree}
    # Must be even dimension
    if TPS_Dim % 2 != 0
        return CTPS(ctps)
    end

    out = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)
    out.map[1] = ctps.map[1]  # copy constant term

    for i in 2:length(ctps.map)
        idx = getindexmap(ctps.polymap[], i)
        if mode == 1
            for p in 1:(TPS_Dim ÷ 2)
                pos1 = 2*p
                pos2 = 2*p + 1
                idx[pos1], idx[pos2] = idx[pos2], idx[pos1]
            end
        end
        newi = findindex(ctps, idx)
        out.map[newi] = conj(ctps.map[i])
    end

    # recompute_degree!(out)
    return out
end

function compute_inverse_maps(vectors::Vector, dim::Int, order::Int)
    w0z = [CTPS(0.0im, dim*2, order) for i in 1:dim*2]
    # w0z and their conjugate
    for i in 1:dim
        w0z[2*i-1].map .= vectors[i]
        w0z[2*i] = conjugate(w0z[2*i-1])
    end 

    # w0zfunc = [evaluate(w0z[i]) for i in 1:dim*2]
    w0zfunc = [PolyEvaluator(w0z[i]) for i in 1:dim*2]
    function ztow(zs)
        return [evaluate(w0zfunc[i], zs) for i in 1:dim*2]
    end

    invw0z = inverse_map(w0z, order)
    # invw0zfunc = [evaluate(invw0z[i]) for i in 1:dim*2]
    invw0zfunc = [PolyEvaluator(invw0z[i]) for i in 1:dim*2]
    function wtoz(ws)
        # return [invw0zfunc[i](ws) for i in 1:dim*2]
        return [evaluate(invw0zfunc[i], ws) for i in 1:dim*2]
    end
    return ztow, wtoz, w0z
end

function derivative(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, ndim::Int, order::Int) where {T, TPS_Dim, Max_TPS_Degree}
    if order ≤ 0
        return CTPS(ctps)
    end

    degree = 0
    for i in 1:length(ctps.map)
        poly = getindexmap(ctps.polymap[], i)
        current_degree = sum(poly[2:end])
        if ctps.map[i] != zero(T) && current_degree > degree
            degree = current_degree
        end
    end
    # If the current degree is lower than the derivative order, return the zero TPS.
    if degree < order
        return CTPS(zero(T), TPS_Dim, Max_TPS_Degree)
    end

    if !(1 ≤ ndim ≤ TPS_Dim)
        error("Inconsistent dimension to take derivative")
    end

    temp = CTPS(ctps) * zero(T)

    for i in 2:length(ctps.map)
        # Get the index map for term i. 
        # (indexlist[1] is total degree; indexlist[j+1] is the exponent for variable j.)
        indexlist = getindexmap(ctps.polymap[], i)
        if indexlist[ndim+1] ≥ order
            thisdim = indexlist[ndim+1]
            indexlist[ndim+1] -= order
            indexlist[1]       -= order
            new_i = findindex(ctps, indexlist)
            temp.map[new_i] = ctps.map[i] * binomial(thisdim, order)
        end
    end
    # temp.degree = ctps.degree - order
    return temp
end

function numerical_inverse_6D_fast(funcs, tpsinv, jacinv, wvals; max_iter=3, tol=1e-11)
    n = length(wvals[1])
    zapprox = tpsinv.(SVector.(wvals...))
    # Pre‑allocate
    zdiff = Matrix{ComplexF64}(undef, n, 6)
    wdiff = Matrix{ComplexF64}(undef, n, 6)
    for it in 1:max_iter
        # compute wdiff in-place
        Threads.@threads for i in 1:n
            vals = funcs(zapprox[i])   # returns a 6‑element tuple/array
            @inbounds for j in 1:6
                wdiff[i, j] = vals[j] - wvals[j][i]
            end
        end
        # fill & solve jacobian in-place
        Threads.@threads for i in 1:n
            J = MMatrix{6,6,ComplexF64}(undef)
            # jacinv returns 36-element vector row-major
            vals = [evaluate(jacinv[j], (zapprox[i])) for j in 1:36]
            J .= transpose(reshape(vals,6,6))
            invJ = inv(J)
            zdiff[i, :] = -invJ * @view wdiff[i, :]
        end

        # update zapprox
        Threads.@threads for i in 1:n
            zapprox[i] += SVector(zdiff[i, :]...)
        end

        if sum(abs, zdiff) < tol
            break
        end
    end
    return transpose(reduce(hcat, zapprox))
end

function itearation_freq(dtheta, theta0, freq, nmap, z_to_w, w_to_z; jacobian=nothing, dist_avoid_res=0.01)
    n1, n2, n3 = size(dtheta[1,:,:,:])
    t = theta0 .+ dtheta

    w_x = exp.(1.0im .* t[1,:,:,:])
    w_xc = conj.(w_x)
    w_y = exp.(1.0im .* t[2,:,:,:])
    w_yc = conj.(w_y)
    w_z = exp.(1.0im .* t[3,:,:,:])
    w_zc = conj.(w_z)

    if isnothing(jacobian)
        xcur, pxcur, ycur, pycur, zcur, pzcur = w_to_z([w_x, w_xc, w_y, w_yc, w_z, w_zc])
    else
        flat_wvals = [vec(permutedims(arr, (3, 2, 1))) for arr in (w_x, w_xc, w_y, w_yc, w_z, w_zc)]
        zapprox = numerical_inverse_6D_fast(z_to_w, w_to_z, jacobian, flat_wvals)
        xcur = zapprox[:, 1]
        pxcur = zapprox[:, 2]
        ycur = zapprox[:, 3]
        pycur = zapprox[:, 4]
        zcur = zapprox[:, 5]
        pzcur = zapprox[:, 6]
        xcur = permutedims(reshape(xcur, n1, n2, n3), (3, 2, 1))  
        pxcur = permutedims(reshape(pxcur, n1, n2, n3), (3, 2, 1))
        ycur = permutedims(reshape(ycur, n1, n2, n3), (3, 2, 1))  
        pycur = permutedims(reshape(pycur, n1, n2, n3), (3, 2, 1))
        zcur = permutedims(reshape(zcur, n1, n2, n3), (3, 2, 1))
        pzcur = permutedims(reshape(pzcur, n1, n2, n3), (3, 2, 1))
    end
    # one-turn map
    # nx, npx, ny, npy, nz, npz = nmap([xcur, pxcur, ycur, pycur, zcur, pzcur]) # not work in Julia
    Ys = broadcast((a,b,c,d,e,f) -> nmap([a,b,c,d,e,f]), xcur, pxcur, ycur, pycur, zcur, pzcur)
    nx = getindex.(Ys, 1)
    npx = getindex.(Ys, 2)
    ny = getindex.(Ys, 3)
    npy = getindex.(Ys, 4)
    nz = getindex.(Ys, 5)
    npz = getindex.(Ys, 6)

    # w-space
    # n_wx, n_wxc, n_wy, n_wyc, n_wz, n_wzc = z_to_w([nx, npx, ny, npy, nz, npz]) # not work in Julia
    Ws = broadcast((a,b,c,d,e,f) -> z_to_w([a,b,c,d,e,f]), nx, npx, ny, npy, nz, npz)
    n_wx = getindex.(Ws, 1)
    # n_wxc = getindex.(Ws, 2)
    n_wy = getindex.(Ws, 3)
    # n_wyc = getindex.(Ws, 4)
    n_wz = getindex.(Ws, 5)
    # n_wzc = getindex.(Ws, 6)

    phix = (-1.0im).*log.(n_wx./w_x) .- freq[1]
    phiy = (-1.0im).*log.(n_wy./w_y) .- freq[2]
    phiz = (-1.0im).*log.(n_wz./w_z) .- freq[3]

    phix_fft = fft(phix)
    phiy_fft = fft(phiy)
    phiz_fft = fft(phiz)

    # Define a helper function similar to numpy.fft.fftfreq.
    function fftfreq_np(n::Integer, d::Real=1.0)
        if iseven(n)
            freqs = vcat(
                0:(n÷2-1),
                -n÷2,
                -(n÷2-1):-1
            )
        else
            freqs = vcat(
                0:((n-1)÷2),
                -((n-1)÷2):-1
            )
        end
        return freqs ./ (n*d)
    end

    fl1 = fftfreq_np(n1, 1.0/n1)
    fl2 = fftfreq_np(n2, 1.0/n2)
    fl3 = fftfreq_np(n3, 1.0/n3)

    freq1list = [a for a in fl1, b in fl2, c in fl3]
    freq2list = [b for a in fl1, b in fl2, c in fl3]
    freq3list = [c for a in fl1, b in fl2, c in fl3]

    freq1list[1, 1, 1] = 1.0
    freq2list[1, 1, 1] = 1.0
    freq3list[1, 1, 1] = 1.0

    newfreq = [0.0 for i in 1:3]
    newfreq[1] = real(freq[1] + phix_fft[1,1,1] / (n1*n2*n3))
    newfreq[2] = real(freq[2] + phiy_fft[1,1,1] / (n1*n2*n3))
    newfreq[3] = real(freq[3] + phiz_fft[1,1,1] / (n1*n2*n3))
    
    res_term = exp.(im .* freq1list .* newfreq[1] .+ im .* freq2list .* 
                newfreq[2] .+ im .* freq3list .* newfreq[3]) .- 1.0
    min_abs_res_term = minimum(abs.(res_term))

    trial = 0
    newfreq_trials = []
    resmin_trials = []

    while min_abs_res_term < 1e-3
        trial += 1
        newfreq[1] = real(freq[1] + (1.0 + dist_avoid_res - 2 * dist_avoid_res * rand()) * phix_fft[1,1,1] / (n1*n2*n3))
        newfreq[2] = real(freq[2] + (1.0 + dist_avoid_res - 2 * dist_avoid_res * rand()) * phiy_fft[1,1,1] / (n1*n2*n3))
        newfreq[3] = real(freq[3] + (1.0 + dist_avoid_res - 2 * dist_avoid_res * rand()) * phiz_fft[1,1,1] / (n1*n2*n3))
        
        push!(newfreq_trials, copy(newfreq))
        
        res_term = exp.(im * freq1list * newfreq[1] + im * freq2list * newfreq[2] + im * freq3list * newfreq[3]) .- 1
        min_abs_res_term = minimum(abs.(res_term))
        
        push!(resmin_trials, min_abs_res_term)
        
        println("Avoiding resonance, trial ", trial)
        
        if trial > 10
            max_value = maximum(resmin_trials)
            max_index = argmax(resmin_trials)
            newfreq = newfreq_trials[max_index]
            println("Selected trial index ", max_index, " with min value ", resmin_trials[max_index])
            break
        end
    end
    
    theta_mx = phix_fft ./ res_term
    theta_my = phiy_fft ./ res_term
    theta_mz = phiz_fft ./ res_term
    theta_mx[1,1,1] = 0.0+0.0im
    theta_my[1,1,1] = 0.0+0.0im
    theta_mz[1,1,1] = 0.0+0.0im

    dt = zeros(ComplexF64, 3, n1, n2, n3)
    dt[1, :, :, :] = ifft(theta_mx)
    dt[2, :, :, :] = ifft(theta_my)
    dt[3, :, :, :] = ifft(theta_mz)
    return dt, newfreq, theta_mx, theta_my, theta_mz
end

function CMscan(transfer_map, dim, order, tunes, 
            x_min, x_max, y_min, y_max, z_min, z_max, 
            n_x, n_y, n_z, pxini, pyini, pzini)
    x_vals = range(x_min, x_max, length=n_x)
    y_vals = range(y_min, y_max, length=n_y)
    z_vals = range(z_min, z_max, length=n_z)

    conv_metrix = zeros(n_x, n_y, n_z)

    tunex, tuney, tunez = tunes
    freq_init = [2*pi*tunex, 2*pi*tuney, 2*pi*tunez]
    hp, tpsMap = transfer_map(dim, order, tunes)

    leftvects = compute_action_angle_polynomials(hp, dim, order)
    ztow, wtoz, w0z = compute_inverse_maps(leftvects, dim, order)

    N = 4 * dim * dim
    wlist = Vector{PolyEvaluator{ComplexF64,2*dim}}(undef, N)
    idx = 1
    
    for i in 1:dim
        tps = CTPS(0.0im, 2*dim, order)
        tps.map .= leftvects[i]
    
        for ct in (tps, conjugate(tps))
            for j in 1:(2*dim)
                wlist[idx] = PolyEvaluator(derivative(ct, j, 1))
                idx += 1
            end
        end
    end

    nalpha, nbeta, ngamma = 16, 16, 16

    # Create ranges that mimic np.linspace(0, 2π, num, endpoint=false)
    alphas = range(0, step=2π/nalpha, length=nalpha)
    betas  = range(0, step=2π/nbeta, length=nbeta)
    gammas = range(0, step=2π/ngamma, length=ngamma)
    
    # Build 3D arrays using comprehensions
    aa = [a for a in alphas, b in betas, g in gammas]
    bb = [b for a in alphas, b in betas, g in gammas]
    gg = [g for a in alphas, b in betas, g in gammas]

    wx0z = CTPS(0.0im, 2*dim, order)
    wy0z = CTPS(0.0im, 2*dim, order)
    wz0z = CTPS(0.0im, 2*dim, order)
    wx0z.map .= leftvects[1]
    wy0z.map .= leftvects[2]
    wz0z.map .= leftvects[3]
    wx0z = PolyEvaluator(wx0z)
    wy0z = PolyEvaluator(wy0z)
    wz0z = PolyEvaluator(wz0z)
    itertimes = 10
    for ix in 1:n_x
        xini = x_vals[ix]
        for iy in 1:n_y
            yini = y_vals[iy]
            for iz in 1:n_z
                start_time = time()
                zini = z_vals[iz]
                zxini = xini - 1.0im * pxini
                zxcini = xini + 1.0im * pxini
                zyini = yini - 1.0im * pyini
                zycini = yini + 1.0im * pyini
                zzini = zini - 1.0im * pzini
                zzcini = zini + 1.0im * pzini

                wxini = evaluate(wx0z, [zxini, zxcini, zyini, zycini, zzini, zzcini])
                wyini = evaluate(wy0z, [zxini, zxcini, zyini, zycini, zzini, zzcini])
                wzini = evaluate(wz0z, [zxini, zxcini, zyini, zycini, zzini, zzcini])

                theta1ini = -1.0im * log(wxini)
                theta2ini = -1.0im * log(wyini)
                theta3ini = -1.0im * log(wzini)

                thetas = zeros(ComplexF64, 3, nalpha, nbeta, ngamma)
                thetas[1, :, :, :] .= theta1ini .+ aa
                thetas[2, :, :, :] .= theta2ini .+ bb
                thetas[3, :, :, :] .= theta3ini .+ gg

                dt = zeros(ComplexF64, 3, nalpha, nbeta, ngamma)
                freqs = copy(freq_init)

                dtdiff_final = nothing
                for _ in 1:itertimes
                    dtold = copy(dt)
                    dt, freqs, _, _, _ = itearation_freq(dt, thetas, freqs, 
                                                        tpsMap, ztow, wtoz, jacobian=wlist)
                    thetas[1, :, :, :] = aa .+ theta1ini .- dt[1, 1, 1, 1]
                    thetas[2, :, :, :] = bb .+ theta2ini .- dt[2, 1, 1, 1]
                    thetas[3, :, :, :] = gg .+ theta3ini .- dt[3, 1, 1, 1]
                    dtdiff = sqrt(sum(abs.(dt .- dtold).^2) / (nalpha * nbeta * ngamma))
                    dtdiff_final = dtdiff
                    if isnan(dtdiff_final) || dtdiff_final > 1e4
                        println("NaN or large value encountered, breaking out of loop.")
                        break
                    end
                end
                if dtdiff_final <= 0.0
                    conv_metrix[ix, iy, iz] = -30
                else
                    conv_metrix[ix, iy, iz] = log10(dtdiff_final)
                end
                end_time = time()
                println("Time taken for (x, y, z) = ($xini, $yini, $zini): $(end_time - start_time) seconds")
            end
        end
    end

    # save the convergence matrix to a file
    open("conv_metrix.txt", "w") do f
        for ix in 1:n_x
            for iy in 1:n_y
                for iz in 1:n_z
                    write(f, "$(conv_metrix[ix, iy, iz]) ")
                end
                write(f, "\n")
            end
            write(f, "\n")
        end
    end
    return conv_metrix
end

function schur_transform!(tpsvar::TPSVar)
    S = schur(tpsvar.sqr_mat)
    Q = S.Z
    T = S.T   # T is upper triangular
    tpsvar.sqr_mat .= Q' * tpsvar.sqr_mat * Q
    return Q
end

function twiss_from_6x6(M::AbstractMatrix{<:Real})
    function twiss2(M2)
        cosμ = (M2[1,1] + M2[2,2]) / 2
        sinμ = sign(M2[1,2]) * sqrt(max(0, 1 - cosμ^2))
        μ    = atan(sinμ, cosμ)      
        β    = M2[1,2] / sinμ
        α    = (M2[1,1] - M2[2,2]) / (2*sinμ)
        γ    = (1 + α^2) / β
        return (α, β, γ, μ)
    end

    Mx = M[1:2, 1:2]
    My = M[3:4, 3:4]
    Mz = M[5:6, 5:6]

    return (
        horizontal = twiss2(Mx),
        vertical   = twiss2(My),
        longitudinal = twiss2(Mz)
    )
end

# function Heono_test(dim, order, tunes)
#     hp = TPSVar4D(order)
#     x, px, y, py = get_variables(hp)
#     tunex = tunes[1]
#     tuney = 1.01 - 2.0*tunex
#     xmu = 2.0*pi*tunex
#     ymu = 2.0*pi*tuney
#     cmx = cos(xmu)
#     smx = sin(xmu)
#     cmy = cos(ymu)
#     smy = sin(ymu)

#     pxm = px - x*x + y*y
#     pym = py + 2.0*x*y

#     xmap = x * cmx + pxm * smx
#     pxmap = -x * smx + pxm * cmx
#     ymap = y * cmy + pym * smy
#     pymap = -y * smy + pym * cmy

#     z1map = xmap - 1.0im * pxmap
#     z1cmap = xmap + 1.0im * pxmap
#     z2map = ymap - 1.0im * pymap
#     z2cmap = ymap + 1.0im * pymap

#     # z1func = numpify(z1map, dim * 2)
#     # z1cfunc = numpify(z1cmap, dim * 2)
#     # z2func = numpify(z2map, dim * 2)
#     # z2cfunc = numpify(z2cmap, dim * 2)

#     z1func = evaluate(z1map)
#     z1cfunc = evaluate(z1cmap)
#     z2func = evaluate(z2map)
#     z2cfunc = evaluate(z2cmap)

#     construct_sqr_matrix(hp, [z1map, z1cmap, z2map, z2cmap])
#     return hp, function(zs)
#         return [
#             z1func(zs),
#             z1cfunc(zs),
#             z2func(zs),
#             z2cfunc(zs)
#         ]
#     end
# end

# function crab_cavity_map(dim, order, tunes)
#     theta_c = 25e-3 
#     f_cc = 197e6
#     beta_cc = 1300.0
#     beta_IP = 0.9
#     b2 = 0.0
#     b3 = 50000.0
#     tune_x = tunes[1]
#     tune_y = tunes[2]
#     tune_z = tunes[3]
#     alpha_c = 1.5e-2
#     V_rf = 15.8e6
#     f_rf = 591e6
#     E = 275e9
#     h = 7560.0
#     phi_s = 0.0

#     c = 299792458.0
#     k_c = 2 * pi * f_cc / c
#     k_rf = 2 * pi * f_rf / c
#     e = 1.602176634e-19
#     gamma = E / (938.272046e6)
#     beta = sqrt(1.0 - 1.0 / gamma^2)
#     eta = alpha_c - 1.0 / gamma^2

#     map = TPSVar6D(order)
#     x, px, y, py, z, pz = get_variables(map)

#     # crabbing kicks
#     delta_px = -tan(theta_c) * sin(k_c * z) / (k_c * sqrt(beta_cc*beta_IP))
#     # delta_px += -b2 * x * sin(k_c * z) # b2 = 0
#     delta_px += b3 * (x^2 - y^2) * sin(k_c * z)
#     p_x_cc = px + delta_px

#     delta_py = -2.0 * b3 * x * y * sin(k_c * z)
#     p_y_cc = py + delta_py

#     delta_pz = -x * tan(theta_c) * cos(k_c * z) / (sqrt(beta_cc*beta_IP))
#     # delta_pz += (b2*k_c/2.0) * (x^2 + y^2) * sin(k_c * z) # b2 = 0
#     delta_pz += (b3*k_c/3.0) * (x^3 - 3.0*x*y^2) * cos(k_c * z)
#     p_z_cc = pz + delta_pz

#     # transverse map
#     mu_x = 2.0 * pi * tune_x
#     x_rot = x*cos(mu_x) + p_x_cc*sin(mu_x)
#     p_x_rot = -x*sin(mu_x) + p_x_cc*cos(mu_x)

#     mu_y = 2.0 * pi * tune_y
#     y_rot = y*cos(mu_y) + p_y_cc*sin(mu_y)
#     p_y_rot = -y*sin(mu_y) + p_y_cc*cos(mu_y)

#     # longitudinal map
#     mu_z = 2.0 * pi * tune_z
#     z_drift = z*cos(mu_z) + p_z_cc*sin(mu_z)
#     p_z_drift = -z*sin(mu_z) + p_z_cc*cos(mu_z)

#     x_map = x_rot
#     px_map = p_x_rot
#     y_map = y_rot
#     py_map = p_y_rot
#     z_map = z_drift
#     pz_map = p_z_drift

#     z1 = x_map - 1.0im * px_map
#     z1c = x_map + 1.0im * px_map
#     z2 = y_map - 1.0im * py_map
#     z2c = y_map + 1.0im * py_map
#     z3 = z_map - 1.0im * pz_map
#     z3c = z_map + 1.0im * pz_map

#     z1mapfunc = evaluate(z1)
#     z1cmapfunc = evaluate(z1c)
#     z2mapfunc = evaluate(z2)
#     z2cmapfunc = evaluate(z2c)
#     z3mapfunc = evaluate(z3)
#     z3cmapfunc = evaluate(z3c)
#     construct_sqr_matrix(map, [z1, z1c, z2, z2c, z3, z3c])
#     return map, function(zs)
#         return [
#             z1mapfunc(zs),
#             z1cmapfunc(zs),
#             z2mapfunc(zs),
#             z2cmapfunc(zs),
#             z3mapfunc(zs),
#             z3cmapfunc(zs)
#         ]
#     end
# end
# CMscan(crab_cavity_map, 3, 3, [0.26, 0.23, 0.005], 
#     -0.0002, -0.0002, 0.0001, 0.0001, -0.0001, -0.0001,
#     1, 1, 1,
#     0.0, 0.0, 0.0)