# using LinearAlgebra
# using Distributions

# linear algebra functions
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

function diagm1(v)
    n = length(v)
    M = zeros(n, n)  
    for i in 1:n
        M[i, i] = v[i] 
    end
    return M
end

function get_centroid!(beam::Beam)
    beam.centroid[1] = sum(beam.r[:, 1]) / beam.nmacro
    beam.centroid[3] = sum(beam.r[:, 3]) / beam.nmacro
    beam.centroid[2] = sum(beam.r[:, 2]) / beam.nmacro
    beam.centroid[4] = sum(beam.r[:, 4]) / beam.nmacro
    
    beam.centroid[5] = sum(beam.r[:, 5]) / beam.nmacro
    beam.centroid[6] = sum(beam.r[:, 6]) / beam.nmacro
    return nothing
end


function get_2nd_moment!(beam::Beam)
    sums = zeros(6,6)
    
    for c in 1:beam.nmacro
        r6 = @view beam.r[c, :]
        for i in 1:6
            for j in 1:6
                sums[i, j] += r6[i] * r6[j]
            end
        end
    end

    for i in 1:6
        for j in 1:6
            beam.moment2nd[i, j] = sums[i, j] / beam.nmacro
        end
    end
    return nothing
end


function get_emittance!(beam::Beam)
    get_centroid!(beam)
    get_2nd_moment!(beam)
    @inbounds for i in 1:6
        @inbounds for j in 1:6
            beam.moment2nd[i, j] -= beam.centroid[i] * beam.centroid[j]
        end
    end
    @inbounds for i in 1:3
        beam.emittance[i] = sqrt(det1(beam.moment2nd[2*i-1:2*i,2*i-1:2*i]))
    end
    @inbounds for i in 1:6
        beam.beamsize[i] = sqrt(beam.moment2nd[i, i])
    end
    beam.emittance[3] = beam.emittance[3] * beam.beta * beam.energy / 2.99792458e8
    return nothing
end

function initilize_6DGaussiandist!(beam::Beam, optics::AbstractOptics4D, lmap::AbstractLongitudinalMap, cutoff::Float64=5.0)
    # 6D Gaussian distribution

    # dist=Truncated(Normal(0.0,1.0),-cutoff,cutoff)

    # beam.r = rand(dist, beam.nmacro, 6)
    beam.r = randn(beam.nmacro, 6)

    get_centroid!(beam)

    for c in 1:beam.nmacro
        beam.r[c, 1] -= beam.centroid[1]
        beam.r[c, 2] -= beam.centroid[2]
        beam.r[c, 3] -= beam.centroid[3]
        beam.r[c, 4] -= beam.centroid[4]
        beam.r[c, 5] -= beam.centroid[5]
        beam.r[c, 6] -= beam.centroid[6]
    end

    get_2nd_moment!(beam)
    
    eigval,eigvec=qr_eigen(beam.moment2nd)
    
    mscale=eigvec * diagm1(1.0 ./ sqrt.(eigval)) * eigvec'
    # # #beam.dist=mscale*beam.dist
    for c in 1:beam.nmacro
        beam.temp1[c]  = mscale[1,1] * beam.r[c, 1] + mscale[1,2] * beam.r[c, 2] + mscale[1,3] * beam.r[c, 3] + mscale[1,4] * beam.r[c, 4] + mscale[1,5] * beam.r[c, 5] + mscale[1,6] * beam.r[c, 6]
        beam.temp2[c]  = mscale[2,1] * beam.r[c, 1] + mscale[2,2] * beam.r[c, 2] + mscale[2,3] * beam.r[c, 3] + mscale[2,4] * beam.r[c, 4] + mscale[2,5] * beam.r[c, 5] + mscale[2,6] * beam.r[c, 6]
        beam.temp3[c]  = mscale[3,1] * beam.r[c, 1] + mscale[3,2] * beam.r[c, 2] + mscale[3,3] * beam.r[c, 3] + mscale[3,4] * beam.r[c, 4] + mscale[3,5] * beam.r[c, 5] + mscale[3,6] * beam.r[c, 6]
        beam.temp4[c]  = mscale[4,1] * beam.r[c, 1] + mscale[4,2] * beam.r[c, 2] + mscale[4,3] * beam.r[c, 3] + mscale[4,4] * beam.r[c, 4] + mscale[4,5] * beam.r[c, 5] + mscale[4,6] * beam.r[c, 6]
        beam.temp5[c]  = mscale[5,1] * beam.r[c, 1] + mscale[5,2] * beam.r[c, 2] + mscale[5,3] * beam.r[c, 3] + mscale[5,4] * beam.r[c, 4] + mscale[5,5] * beam.r[c, 5] + mscale[5,6] * beam.r[c, 6]
        beam.r[c, 6] = mscale[6,1] * beam.r[c, 1] + mscale[6,2] * beam.r[c, 2] + mscale[6,3] * beam.r[c, 3] + mscale[6,4] * beam.r[c, 4] + mscale[6,5] * beam.r[c, 5] + mscale[6,6] * beam.r[c, 6]
    end
    for c in 1:beam.nmacro
        beam.r[c, 1] = beam.temp1[c] * sqrt(beam.emittance[1]*optics.optics_x.beta)
        beam.r[c, 2] = beam.temp2[c] * sqrt(beam.emittance[1]/optics.optics_x.beta)
        beam.r[c, 3] = beam.temp3[c] * sqrt(beam.emittance[2]*optics.optics_y.beta)
        beam.r[c, 4] = beam.temp4[c] * sqrt(beam.emittance[2]/optics.optics_y.beta)
        beam.r[c, 2] += beam.r[c, 1] * (optics.optics_x.alpha/optics.optics_x.beta)
        beam.r[c, 4] += beam.r[c, 3] * (optics.optics_y.alpha/optics.optics_y.beta)
    end

    
    # #generate longtiudinal distribution based on small amplitude approximation
    eta_p=lmap.alphac-1.0/beam.gamma^2
    
    Qs=sqrt(lmap.RF.volt*lmap.RF.h*abs(eta_p*cos(lmap.RF.phis))/2/π/beam.beta^2/beam.energy)

    emit_deltap_z=beam.emittance[3]*2.99792458e8/beam.beta/beam.energy
    invbeta_deltap_z=Qs*lmap.RF.k/lmap.RF.h/abs(eta_p)
    for c in 1:beam.nmacro
        beam.r[c, 5] = beam.temp5[c] * sqrt(emit_deltap_z*invbeta_deltap_z)
        beam.r[c, 6] = beam.r[c, 6] * sqrt(emit_deltap_z*invbeta_deltap_z)
    end
    
    return nothing
end

function histogram1DinZ!(beam::Beam, nbins::Int64, inzindex, zhist, zhist_edges) 
    # histogram in z
    num_macro=beam.nmacro
    zmax=maximum(abs.(beam.r[:, 5]))
    zmin=-zmax
    zhist .= 0.0

    total_range_start = zmin - (zmax - zmin) / nbins
    total_range_end = zmax + (zmax - zmin) / nbins
    zstep = (total_range_end - total_range_start) / nbins
    @inbounds for i in 0:nbins
        zhist_edges[i+1] = total_range_start + i * zstep
    end
    zsep=(zhist_edges[end]-zhist_edges[1])/nbins
    @inbounds for i in 1:num_macro
        ibin = (beam.r[:, 5][i] - zhist_edges[1])/zsep  # number of bin from 0
        dx = round(ibin) - ibin
        binnum=Int64(floor(ibin)+1)
        inzindex[i] = binnum 
        neighbor = binnum + Int64(sign(dx))
        ratio = (0.5 - abs(dx))/0.5
        weight_neighbor = 0.5 * ratio^2
        zhist[binnum] += 1.0 - weight_neighbor
        zhist[neighbor] += weight_neighbor
    end
    return nothing
end

function histogram1DinZ!(beam::Beam)
    # histogram in z
    histogram1DinZ!(beam, beam.znbin, beam.inzindex, beam.zhist, beam.zhist_edges)
    return nothing
end