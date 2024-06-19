include("../src/JuTrack.jl")
include("../src/demo/ssrf_ring.jl")
using .JuTrack
using Serialization
using Plots  
using LinearAlgebra
using LaTeXStrings
using DelimitedFiles

# RING = deserialize("src/demo/spear3.jls")
RING = ssrf(-1.063770, 0)
# x: -0.04 to 0.04
angle_list = [ 0 + pi/36 * i for i in 0:36]
amp_list = [0.001 + 0.001 * i for i in 0:40]

N = length(angle_list) * length(amp_list)
particles = zeros(length(angle_list) * length(amp_list), 6)
for i in 1:length(angle_list)
    for j in 1:length(amp_list)
        particles[(i-1)*length(amp_list) + j, 1] = amp_list[j] * cos(angle_list[i])
        particles[(i-1)*length(amp_list) + j, 3] = amp_list[j] * sin(angle_list[i])
    end
end

beam = Beam(copy(particles), energy=3.0e9)
ringpass!(RING, beam, 1)
survived = findall(x -> x == 0, beam.lost_flag)

N_survived = length(survived)
survived_particles = zeros(N_survived, 6)
for i in 1:N_survived
    survived_particles[i, :] = particles[survived[i], :]
end

function TPSA_track_jacobian(amp, angle)
    x = amp * cos(angle)
    y = amp * sin(angle)
    X = CTPS(x, 1, 6, 1)
    PX = CTPS(0.0, 2, 6, 1)
    Y = CTPS(y, 3, 6, 1)
    PY = CTPS(0.0, 4, 6, 1)
    Z = CTPS(0.0, 5, 6, 1)
    DELTA = CTPS(0.0, 6, 6, 1)
    rin = [X, PX, Y, PY, Z, DELTA]
    ringpass_TPSA!(RING, rin, 1)
    jaco = zeros(6, 6)
    for i in 1:6
        jaco[i, :] = rin[i].map[2:7]
    end
    # eigenvalues, eigenvectors = qr_eigen(jaco)
    eigenvalues, eigenvectors = eigen(jaco)
    return eigenvalues, eigenvectors
end
function TPSA_track_jacobian(amp, angle)
    x = amp * cos(angle)
    y = amp * sin(angle)
    X = CTPS(x, 1, 6, 10)
    PX = CTPS(0.0, 2, 6, 10)
    Y = CTPS(y, 3, 6, 10)
    PY = CTPS(0.0, 4, 6, 10)
    Z = CTPS(0.0, 5, 6, 10)
    DELTA = CTPS(0.0, 6, 6, 10)
    rin = [X, PX, Y, PY, Z, DELTA]
    E1 = KQUAD(name="Q1", k1=2.0, len=0.34)
    ringpass_TPSA!([E1], rin, 1)
    jaco = zeros(6, 6)
    for i in 1:6
        jaco[i, :] = rin[i].map[2:7]
    end
    eigenvalues, eigenvectors = eigen(jaco)
    deteminant = det(jaco)
    return jaco, deteminant, eigenvalues, eigenvectors
end
for i in 1:N_survived
    angle = atan(survived_particles[i, 3], survived_particles[i, 1])
    amplitude = sqrt(survived_particles[i, 1]^2 + survived_particles[i, 3]^2)
    jaco, det, eigenvalues, eigenvectors = TPSA_track_jacobian(amplitude, angle)
    println("Progress: ", i, "/", N_survived, " Determinant: ", det, "eigen: ", abs.(eigenvalues))
end


function normalize(v)
    return v / norm(v)
end
function cosine_similarity(v1, v2)
    return abs(dot(v1, v2)) / (norm(v1) * norm(v2))
end

function track_eigen(Jacobian_previous, Jacobian_current)
    eigenvalues_prev, eigenvectors_prev = eigen(Jacobian_previous)
    eigenvalues_curr, eigenvectors_curr = eigen(Jacobian_current)
    
    # Normalize eigenvectors
    eigenvectors_prev = [normalize(vec) for vec in eachcol(eigenvectors_prev)]
    eigenvectors_curr = [normalize(vec) for vec in eachcol(eigenvectors_curr)]
    
    matched_eigenvalues = Complex{Float64}[]
    matched_eigenvectors = zeros(Complex{Float64}, 6, 6)
    
    i = 1
    for vec_prev in eigenvectors_prev
        best_match_idx = -1
        best_similarity = -Inf
        
        # Compare with current eigenvectors
        for (j, vec_curr) in enumerate(eigenvectors_curr)
            similarity = cosine_similarity(vec_prev, vec_curr)
            if similarity > best_similarity
                best_similarity = similarity
                best_match_idx = j
            end
        end
        
        # Match eigenvalues based on best eigenvector alignment
        push!(matched_eigenvalues, eigenvalues_curr[best_match_idx])
        matched_eigenvectors[:, i] = eigenvectors_curr[best_match_idx]
        i += 1
    end
    
    return matched_eigenvalues, matched_eigenvectors
end

function track_eigenvalues(eigenvectors_prev, Jacobian_current)
    eigenvalues_curr, eigenvectors_curr = eigen(Jacobian_current)
    
    eigenvectors_prev = [normalize(vec) for vec in eachcol(eigenvectors_prev)]
    eigenvectors_curr = [normalize(vec) for vec in eachcol(eigenvectors_curr)]
    
    matched_eigenvalues = Complex{Float64}[]
    matched_eigenvectors = zeros(Complex{Float64}, 6, 6)
    
    for vec_prev in eigenvectors_prev
        best_match_idx = -1
        best_similarity = -Inf
        
        # Compare with current eigenvectors
        for (j, vec_curr) in enumerate(eigenvectors_curr)
            similarity = cosine_similarity(vec_prev, vec_curr)
            if similarity > best_similarity
                best_similarity = similarity
                best_match_idx = j
            end
        end
        
        # Match eigenvalues based on best eigenvector alignment
        push!(matched_eigenvalues, eigenvalues_curr[best_match_idx])
    end
    
    return matched_eigenvalues
end

x = amp_list[1] * cos(angle_list[1])
y = amp_list[1] * sin(angle_list[1])
rin = [CTPS(x, 1, 6, 1), CTPS(0.0, 2, 6, 1), CTPS(y, 3, 6, 1), CTPS(0.0, 4, 6, 1), CTPS(0.0, 5, 6, 1), CTPS(0.0, 6, 6, 1)]
ringpass_TPSA!(RING, rin, 1)
jaco_0 = zeros(6, 6)
for i in 1:6
    jaco_0[i, :] = rin[i].map[2:7]
end

eigens = []
for i in 1:length(angle_list)
    if i > 0
        x = amp_list[1] * cos(angle_list[i])
        y = amp_list[1] * sin(angle_list[i])
        rin = [CTPS(x, 1, 6, 1), CTPS(0.0, 2, 6, 1), CTPS(y, 3, 6, 1), CTPS(0.0, 4, 6, 1), CTPS(0.0, 5, 6, 1), CTPS(0.0, 6, 6, 1)]
        ringpass_TPSA!(RING, rin, 1)
        jaco_1 = zeros(6, 6)
        for i in 1:6
            jaco_1[i, :] = rin[i].map[2:7]
        end
        e_1, v_1 = track_eigen(jaco_0, jaco_1)
    end
    for j in 1:length(amp_list)
        x = amp_list[j] * cos(angle_list[i])
        y = amp_list[j] * sin(angle_list[i])
        beam = Beam([x 0.0 y 0.0 0.0 0.0], energy=3.0e9)
        ringpass!(RING, beam, 1)
        if beam.lost_flag[1] == 1
            continue
        end
        rin = [CTPS(x, 1, 6, 1), CTPS(0.0, 2, 6, 1), CTPS(y, 3, 6, 1), CTPS(0.0, 4, 6, 1), CTPS(0.0, 5, 6, 1), CTPS(0.0, 6, 6, 1)]
        ringpass_TPSA!(RING, rin, 1)
        jaco = zeros(6, 6)
        for i in 1:6
            jaco[i, :] = rin[i].map[2:7]
        end
        push!(eigens, track_eigenvalues(v_1, jaco))
    end
end

eigens_array = zeros(Complex{Float64}, length(eigens), 6)
for i in 1:length(eigens)
    eigens_array[i, :] = [eigens[i][1], eigens[i][2], eigens[i][3], eigens[i][4], eigens[i][5], eigens[i][6]]
end

survived_particles = hcat(survived_particles, eigens_array)
# writedlm("Spear3_eigenvalues.csv", eigens, ',')
