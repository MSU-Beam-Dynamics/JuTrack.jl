include("../src/JuTrack.jl")
using .JuTrack
using Plots  
using LaTeXStrings
using Enzyme
Enzyme.API.runtimeActivity!(true)


function f(x)
L12_ID = 0.2584003731494474
# k1_1 = -4.611667042711146e-5
# k1_2 = 6.4800348490202825e-6
# k1_3 = -6.480034849020284e-6
L01_ID = 2.7260771047243093
L23_ID = 0.891400429082984
a2 = 0.010394537621095908
a1 = (0.02418798691622576 - a2) / 2.0

S = KSEXT(len=1.0, k2=x)

OSB12_ID = DRIFT(name="OSB12", len=L12_ID)
# EDGE1_ID = thinMULTIPOLE(name="EDGE1", PolynomB=[0.0, k1_1, 0.0, 0.0])
# EDGE2_ID = thinMULTIPOLE(name="EDGE2", PolynomB=[0.0, k1_2, 0.0, 0.0])
# EDGE3_ID = thinMULTIPOLE(name="EDGE3", PolynomB=[0.0, k1_3, 0.0, 0.0])

D01A_ID = SBEND(name="D01A", len=L01_ID, angle=a1)
D23_ID = SBEND(name="D23", len=L23_ID, angle=a2)
D01B_ID = SBEND(name="D01B", len=L01_ID, angle=a1)
println("build line")
line = [S, D01A_ID, OSB12_ID, D23_ID, OSB12_ID, D01B_ID]
println("build beam")
ebeam = Beam(randn(10000,6)./1e4)
# opIP=optics4DUC(4.0,0.0,1.0,0.0)
# ebeam = Beam(zeros(10000,6), emittance=[20e-9/0.8, 1.3e-9, 1.36e-4])
# mainRF=AccelCavity(197e6, 1e6, 2520.0, 0.0)
# αc=1e-3
# lmap=LongitudinalRFMap(αc, mainRF)
# initilize_6DGaussiandist!(ebeam, opIP, lmap)
# get_emittance!(ebeam)

println("pass line")
linepass!(line, ebeam)
println("get emittance")
get_emittance!(ebeam)
return ebeam.emittance
end

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
function initilize_6DGaussiandist1!(beam::Beam, betx, bety, alphax, alphay)
    # 6D Gaussian distribution

    # dist=Truncated(Normal(0.0,1.0),-cutoff,cutoff)

    # beam.r = rand(dist, beam.nmacro, 6)
    # beam.r = randn(beam.nmacro, 6) # do it in this function will be problematic

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
    # beam.r[:, 1] *= x
    # #beam.dist=mscale*beam.dist
    for c in 1:beam.nmacro
        beam.temp1[c]  = mscale[1,1] * beam.r[c, 1] + mscale[1,2] * beam.r[c, 2] + mscale[1,3] * beam.r[c, 3] + mscale[1,4] * beam.r[c, 4] + mscale[1,5] * beam.r[c, 5] + mscale[1,6] * beam.r[c, 6]
        beam.temp2[c]  = mscale[2,1] * beam.r[c, 1] + mscale[2,2] * beam.r[c, 2] + mscale[2,3] * beam.r[c, 3] + mscale[2,4] * beam.r[c, 4] + mscale[2,5] * beam.r[c, 5] + mscale[2,6] * beam.r[c, 6]
        beam.temp3[c]  = mscale[3,1] * beam.r[c, 1] + mscale[3,2] * beam.r[c, 2] + mscale[3,3] * beam.r[c, 3] + mscale[3,4] * beam.r[c, 4] + mscale[3,5] * beam.r[c, 5] + mscale[3,6] * beam.r[c, 6]
        beam.temp4[c]  = mscale[4,1] * beam.r[c, 1] + mscale[4,2] * beam.r[c, 2] + mscale[4,3] * beam.r[c, 3] + mscale[4,4] * beam.r[c, 4] + mscale[4,5] * beam.r[c, 5] + mscale[4,6] * beam.r[c, 6]
        beam.temp5[c]  = mscale[5,1] * beam.r[c, 1] + mscale[5,2] * beam.r[c, 2] + mscale[5,3] * beam.r[c, 3] + mscale[5,4] * beam.r[c, 4] + mscale[5,5] * beam.r[c, 5] + mscale[5,6] * beam.r[c, 6]
        beam.r[c, 6] = mscale[6,1] * beam.r[c, 1] + mscale[6,2] * beam.r[c, 2] + mscale[6,3] * beam.r[c, 3] + mscale[6,4] * beam.r[c, 4] + mscale[6,5] * beam.r[c, 5] + mscale[6,6] * beam.r[c, 6]
    end
    for c in 1:beam.nmacro
        beam.r[c, 1] = beam.temp1[c] * sqrt(beam.emittance[1]*betx)
        beam.r[c, 2] = beam.temp2[c] * sqrt(beam.emittance[1]/betx)
        beam.r[c, 3] = beam.temp3[c] * sqrt(beam.emittance[2]*bety)
        beam.r[c, 4] = beam.temp4[c] * sqrt(beam.emittance[2]/bety)
        beam.r[c, 2] += beam.r[c, 1] * (alphax/betx)
        beam.r[c, 4] += beam.r[c, 3] * (alphay/bety)
    end

    
    # #generate longtiudinal distribution based on small amplitude approximation
    # eta_p=lmap.alphac-1.0/beam.gamma^2
    
    # Qs=sqrt(lmap.RF.volt*lmap.RF.h*abs(eta_p*cos(lmap.RF.phis))/2/π/beam.beta^2/beam.energy)

    # emit_deltap_z=beam.emittance[3]*2.99792458e8/beam.beta/beam.energy
    # invbeta_deltap_z=Qs*lmap.RF.k/lmap.RF.h/abs(eta_p)
    # for c in 1:beam.nmacro
    #     beam.r[c, 5] = beam.temp5[c] * sqrt(emit_deltap_z*invbeta_deltap_z)
    #     beam.r[c, 6] = beam.r[c, 6] * sqrt(emit_deltap_z*invbeta_deltap_z)
    # end
    
    return nothing
end

function f1(ex)

    particles = randn(10000,6)/1e4
    # particles[:, 1] *= ex # works
    ebeam = Beam(particles, emittance=[ex, 1.3e-9, 1.36e-4])
    # ebeam.r[:, 1] *= ex # works
    initilize_6DGaussiandist1!(ebeam, 0.45*0.8, 0.056, 0.0, 0.0)
    get_emittance!(ebeam)
    return ebeam.emittance
end

println(f1(2.5e-9))
grad = autodiff(Forward, f1, Duplicated, Duplicated(2.5e-9, 1.0)) 
println(grad)