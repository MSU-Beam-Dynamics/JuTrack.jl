using StatsBase
using StatsBase: Weights
using PyCall
np = pyimport("numpy")
RegularGridInterpolator = pyimport("scipy.interpolate").RegularGridInterpolator
function space_charge_ML!(r_in, K, num_particles, le, lost_flags, model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)
    x = r_in[1:6:end]
    y = r_in[3:6:end]
    x_survi = x[iszero.(lost_flags)]
    y_survi = y[iszero.(lost_flags)]
    H = fit(Histogram, (x_survi, y_survi), (xedges, yedges), closed=:left)
    H_counts = float(H.weights)
    H_avg = H_counts ./ num_particles
    input_tensor = zeros(1, 1, 128, 128)
    input_tensor[1, 1, :, :] .= (H_avg .- x_mean) ./ x_std
    input_tensor = convert(Array{Float32}, input_tensor)
    input = Dict("input" => input_tensor)
    output = model(input)
    u_output = (output["output"] .* y_std) .+ y_mean

    dx = zeros(1, 128, 128)
    dy = zeros(1, 128, 128)
    smoothed_H = GAN_model.smooth_prediction(u_output)

    dx[1,:,:] .= np.gradient(smoothed_H[1,1,:,:], delta, axis=0)
    dy[1,:,:] .= np.gradient(smoothed_H[1,1,:,:], delta, axis=1)

    dh_dx, dh_dy = derivative_histogram(smoothed_H[1,1,:,:], delta, delta, xaxis, xaxis, x, y)

    r_in[2:6:end] .-= (le * K).* dh_dx
    r_in[4:6:end] .-= (le * K).* dh_dy
    return nothing
end
mutable struct DRIFT_SC_ML <: AbstractElement
    name::String
    len::Float64
    T1::Array{Float64,1}
    T2::Array{Float64,1}
    R1::Array{Float64,2}
    R2::Array{Float64,2}        
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    a::Float64
    b::Float64
    Nl::Int64
    Nm::Int64
    Nsteps::Int64 # Number of steps for space charge calculation. One step represents a half-kick-half.
    eletype::String

    DRIFT_SC_ML(;name::String = "DRIFT_SC", len::Float64 = 0.0, T1::Array{Float64,1} = zeros(6), 
        T2::Array{Float64,1} = zeros(6), R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6), 
        RApertures::Array{Float64,1} = zeros(6), EApertures::Array{Float64,1} = zeros(6), a::Float64 = 1.0, b::Float64 = 1.0,
        Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1) = new(name, len, T1, T2, R1, R2, RApertures, EApertures, a, b, Nl, Nm, Nsteps, "DRIFT_SC")
end
mutable struct KQUAD_SC_ML <: AbstractElement
    name::String
    len::Float64
    k1::Float64
    PolynomA::Array{Float64,1}
    PolynomB::Array{Float64,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    FringeIntM0::Array{Float64,1}
    FringeIntP0::Array{Float64,1}
    T1::Array{Float64,1}
    T2::Array{Float64,1}
    R1::Array{Float64,2}
    R2::Array{Float64,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{Float64,1}
    a::Float64
    b::Float64
    Nl::Int64
    Nm::Int64
    Nsteps::Int64
    eletype::String

    function KQUAD_SC_ML(;name::String = "Quad", len::Float64 = 0.0, k1::Float64 = 0.0, 
                    PolynomA::Array{Float64,1} = zeros(Float64, 4), 
                    PolynomB::Array{Float64,1} = zeros(Float64, 4), MaxOrder::Int64=1, 
                    NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                    FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
                    FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
                    T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
                    R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                    EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2),
                    a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)
        if k1 != 0.0 && PolynomB[2] == 0.0
            PolynomB[2] = k1
        end
        new(name, len, k1, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit, 
            FringeIntM0, FringeIntP0, T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, "KQUAD_SC")
    end
end


function pass!(ele::DRIFT_SC_ML, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam, model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    lstep = ele.len / ele.Nsteps
    for i in 1:ele.Nsteps
        DriftPass_SC_ML!(r_in, lstep, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags,
            ele.a, ele.b, ele.Nl, ele.Nm, K, model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)
    end
    return nothing
end

function DriftPass_SC_ML!(r_in::Array{Float64,1}, le::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64, 2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1}, a, b, Nl, Nm, K, model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)

    # half-kick-half
    Threads.@threads for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            drift6!(r6, le/2.0)
            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end

    space_charge_ML!(r_in, K, num_particles, le, lost_flags, model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)

    Threads.@threads for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            drift6!(r6, le/2.0)
            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end

    return nothing
end

function derivative_histogram(H, dx, dy, xaxis, yaxis, xvec, yvec)
    df_dx = np.gradient(H, dx, axis=0)
    df_dy = np.gradient(H, dy, axis=1)
    
    xy = np.transpose(np.array([xvec, yvec]))

    interpolator_dx = RegularGridInterpolator((xaxis, yaxis), df_dx)
    interpolator_dy = RegularGridInterpolator((xaxis, yaxis), df_dy)

    df_dx_v = interpolator_dx(xy)
    df_dy_v = interpolator_dy(xy)
    return df_dx_v, df_dy_v
end

function StrMPoleSymplectic4Pass_SC!(r::Array{Float64,1}, le::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    fringeIntM0::Array{Float64,1},  # I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 
    fringeIntP0::Array{Float64,1},  # I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1}, a::Float64, b::Float64, Nl::Int, Nm::Int, K::Float64, Nsteps::Int,
    model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)
    # Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # and [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    lstep = le / Nsteps
    SL = lstep / 2.0 /num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if FringeQuadEntrance==2 #&& !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleEntrance = 1
    else
        useLinFrEleEntrance = 0
    end
    if FringeQuadExit==2 #&& !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleExit = 1
    else
        useLinFrEleExit = 0
    end

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le 
    end

    for step in 1:Nsteps
        Threads.@threads for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[(c-1)*6+1:c*6]
            NormL1 = L1 / (1.0 + r6[6])
            NormL2 = L2 / (1.0 + r6[6])
            
            if step == 1
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end

                if FringeQuadEntrance != 0 && B[2] != 0
                    if useLinFrEleEntrance == 1
                        linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
                    else
                        QuadFringePassP!(r6, B[2])
                    end
                end
            end

            # Integrator
            for m in 1:num_int_step
                fastdrift!(r6, NormL1, L1)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkick!(r6, A, B, K2, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL1, L1)
            end

            if check_lost(r6)
                lost_flags[c] = 1
            end
        end

        space_charge_ML!(r, K, num_particles, lstep, lost_flags, model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)

        Threads.@threads for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[(c-1)*6+1:c*6]

            NormL1 = L1 / (1.0 + r6[6])
            NormL2 = L2 / (1.0 + r6[6])
    
            # Integrator
            for m in 1:num_int_step
                fastdrift!(r6, NormL1, L1)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkick!(r6, A, B, K2, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL1, L1)
            end
            
            if step == Nsteps
                if FringeQuadExit != 0 && B[2] != 0
                    if useLinFrEleExit == 1
                        linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
                    else
                        QuadFringePassN!(r6, B[2])
                    end
                end
    
                # Misalignment at exit
                if !iszero(R2)
                    multmv!(r6, R2)
                end
                if !iszero(T2)
                    addvv!(r6, T2)
                end
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end
    if le > 0
        B[1] += sin(KickAngle[1]) / le
        A[1] -= sin(KickAngle[2]) / le 
    end

    return nothing
end

function pass!(ele::KQUAD_SC_ML, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam, model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)
    # ele: KQUAD_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    PolynomB = zeros(4)
    E0 = particles.energy
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[2] = ele.k1
        StrMPoleSymplectic4Pass_SC!(r_in,  ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps, model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        StrMPoleSymplectic4Pass_SC!(r_in,  ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps, model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)
    end
    return nothing
end

function linepass_ML!(line::Vector, particles::Beam, model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    if length(particles6) != np*6
        error("The number of particles does not match the length of the particle array")
    end
    for i in eachindex(line)
        pass!(line[i], particles6, np, particles, model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)        
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end