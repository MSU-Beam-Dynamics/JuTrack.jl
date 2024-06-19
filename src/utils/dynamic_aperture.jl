function dynamic_aperture(RING, nturns, amp_max, amp_step, angle_steps, E=3.0e9)
    # DA calculation
    # RING: lattice
    # nturns: number of turns
    # amp_max: maximum amplitude [m]
    # amp_step: step size of amplitude scan [m]
    # angle_steps: number of lines. 11, 13, 19, 21, 37 etc.
    # E: beam energy [eV]
    # return: boundary points of dynamic aperture, survived particles

    println("number of threads using for parallel comptuing: ", Threads.nthreads())
    angle_list = [pi/(angle_steps-1) * i for i in 0:angle_steps-1]
    amp_list = [amp_step * i for i in 1:Int(amp_max/amp_step)]

    particles = zeros(length(angle_list) * length(amp_list), 6)
    for i in 1:length(angle_list)
        for j in 1:length(amp_list)
            particles[(i-1)*length(amp_list) + j, 1] = amp_list[j] * cos(angle_list[i])
            particles[(i-1)*length(amp_list) + j, 3] = amp_list[j] * sin(angle_list[i])
        end
    end

    beam = Beam(copy(particles), energy=E)
    pringpass!(RING, beam, nturns)
    survived = findall(x -> x == 0, beam.lost_flag)

    N_survived = length(survived)
    survived_particles = zeros(N_survived, 6)
    for i in 1:N_survived
        survived_particles[i, :] = particles[survived[i], :]
    end

    DA = zeros(length(angle_list), 2)
    for i in 1:length(angle_list)
        for j in 1:length(amp_list)
            c = survived_particles[:, 1] ./ sqrt.(survived_particles[:, 1].^2 + survived_particles[:, 3].^2)
            idx = findall(x -> x < 0.00001, abs.(c .- cos(angle_list[i])))
            for j in 1:length(idx) -1
                DA[i, 1] = survived_particles[idx[j], 1]
                DA[i, 2] = survived_particles[idx[j], 3]
                if abs(sqrt(survived_particles[idx[j], 1]^2 + survived_particles[idx[j], 3]^2) - 
                    sqrt(survived_particles[idx[j+1], 1]^2 + survived_particles[idx[j+1], 3]^2)) > amp_step * 1.1
                    break
                end
            end
        end
    end
    return DA, survived_particles
end
