using JuTrack
function fastfindm66(LATTICE, dp=0.0)
    # assume the closed orbit is zero
    NE = length(LATTICE)
	DPstep = 3e-6
    XYStep = 3e-8  # Default step size for numerical differentiation
    scaling =  [XYStep, XYStep, XYStep, XYStep, DPstep, DPstep]
    
    # find initial orbit
	orbitin = [0.0 0.0 0.0 0.0 dp 0.0]

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

    beam = Beam(RIN)
    
    linepass!(LATTICE, beam)
    TMAT3 = transpose(beam.r)
    M66 = (TMAT3[:,1:6] - TMAT3[:,7:12]) ./ scaling

    return M66
end

function fastfindm66_refpts(LATTICE, dp::Float64, refpts::Vector{Int})
    # assume the closed orbit is zero
    NE = length(LATTICE)
	DPstep = 3e-6
    XYStep = 3e-8  # Default step size for numerical differentiation
    scaling =  [XYStep, XYStep, XYStep, XYStep, DPstep, DPstep]
    
    # find initial orbit
	orbitin = [0.0 0.0 0.0 0.0 dp 0.0]

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

    beam = Beam(RIN)
    
    rout = linepass!(LATTICE, beam, refpts)
    M66_refpts = zeros(6, 6, length(refpts))
    for i in 1:length(refpts)
        TMAT3 = transpose(rout[i])
        M66 = (TMAT3[:,1:6] - TMAT3[:,7:12]) ./ scaling
        M66_refpts[:, :, i] = M66
    end
    return M66_refpts
end

function ADfastfindm66(LATTICE, dp, changed_idx, changed_ele)
    # assume the closed orbit is zero
    NE = length(LATTICE)
	DPstep = 3e-6
    XYStep = 3e-8  # Default step size for numerical differentiation
    scaling =  [XYStep, XYStep, XYStep, XYStep, DPstep, DPstep]
    
    # find initial orbit
	orbitin = [0.0 0.0 0.0 0.0 dp 0.0]

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

    beam = Beam(RIN)
    
    ADlinepass!(LATTICE, beam, changed_idx, changed_ele)
    TMAT3 = transpose(beam.r)
    M66 = (TMAT3[:,1:6] - TMAT3[:,7:12]) ./ scaling

    return M66
end

function ADfastfindm66_refpts(LATTICE, dp::Float64, refpts::Vector{Int}, changed_idx, changed_ele)
    # assume the closed orbit is zero
    NE = length(LATTICE)
	DPstep = 3e-6
    XYStep = 3e-8  # Default step size for numerical differentiation
    scaling =  [XYStep, XYStep, XYStep, XYStep, DPstep, DPstep]
    
    # find initial orbit
	orbitin = [0.0 0.0 0.0 0.0 dp 0.0]

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

    beam = Beam(RIN)
    
    rout = ADlinepass!(LATTICE, beam, refpts, changed_idx, changed_ele)
    M66_refpts = zeros(6, 6, length(refpts))
    for i in 1:length(refpts)
        TMAT3 = transpose(rout[i])
        M66 = (TMAT3[:,1:6] - TMAT3[:,7:12]) ./ scaling
        M66_refpts[:, :, i] = M66
    end
    return M66_refpts
end
# function findorbit(ring, dp)
#     # Default orbit initial condition 
#     orbitin = xorbit_dp(ring, dp)
#     return orbitin
# end


# function xorbit_dp(ring, dp=0.0)
#     XYStep = 3e-8   # Step size for numerical differentiation
#     dps = 1e-12    # Convergence threshold
#     max_iterations = 20  # Max. iterations
#     Ri = [0 0 0 0 dp 0]  # Default initial guess for the orbit

#     scaling = XYStep .* [1, 1, 1, 1]
#     D = zeros(5, 6)
#     for i in 1:4
#         D[i, i] = scaling[i]
#     end

#     # change = Inf
#     I = zeros(4, 4)
#     for i in 1:4
#         I[i, i] = 1
#     end
#     for itercount in 1:max_iterations
#         RMATi = ones(5,1) * Ri + D
#         beam = Beam(RMATi)

#         linepass!(ring, beam)
        
#         RMATf = transpose(beam.r)
#         Rf = RMATf[:,end]
#         # Compute the transverse part of the Jacobian
#         J4 = (RMATf[1:4,1:4] .- RMATf[1:4,5]) ./ scaling
#         Ri_next = Ri + transpose([ (I - J4) \ reshape((Rf[1:4] - Ri[1:4]), 4, 1) ;0 ;0])
#         change = sqrt(sum((Ri_next - Ri).^2))
#         Ri = Ri_next
#         if change < dps
#             break
#         end
#     end

#     return Ri
# end

D1 = DRIFT(name="D1", len=0.5)
D2 = DRIFT(name="D2", len=0.5)
Q1 = QUAD(name="Q1", len=0.5, k1=1.0)
Kick = HKICKER(name="Kick", len=1.0, xkick=0.2)
line = [D1, Q1,  Kick]
M66_refpts = fastfindm66_refpts(line, 0.0, [1, 2, 3])
function f(k)
    global line
    changed_idx = findelem(line, :name, "Q1")
    changed_ele = [QUAD(len=0.5, k1=k) for i in 1:length(changed_idx)]
    m66 = ADfastfindm66_refpts(line, 0.0, [1,2,3], changed_idx, changed_ele)
    return m66[1, 1, end] + m66[2, 2, end]
end

function f_TPSA(k)
    global line
    changed_idx = findelem(line, :name, "Q1")
    changed_ele = [QUAD(len=0.5, k1=k) for i in 1:length(changed_idx)]
    m66 = ADfindm66(line, 0.0, 1, changed_idx, changed_ele)
    return m66[1, 1] + m66[2, 2]
end

@time println(f(1.2))

using BenchmarkTools
@btime grad = autodiff(Forward, f, DuplicatedNoNeed, Duplicated(1.2, 1.0))
@btime grad = autodiff(Forward, f_TPSA, DuplicatedNoNeed, Duplicated(1.2, 1.0))
println(autodiff(Forward, f, DuplicatedNoNeed, Duplicated(1.2, 1.0)))
println(autodiff(Forward, f_TPSA, DuplicatedNoNeed, Duplicated(1.2, 1.0)))