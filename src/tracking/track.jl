function matrix_to_array(matrix::Matrix{Float64})
    particles = zeros(Float64, size(matrix, 1)*6)
    for i in 1:size(matrix, 1)
        for j in 1:6
            particles[(i-1)*6+j] = matrix[i, j]
        end
    end

    return particles
end

function array_to_matrix(array::Vector{Float64}, n::Int)
    particles = zeros(Float64, n, 6)
    for i in 1:n
        for j in 1:6
            particles[i, j] = array[(i-1)*6+j]
        end
    end
    return particles
end

function linepass!(line, particles::Beam)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    if length(particles6) != np*6
        error("The number of particles does not match the length of the particle array")
    end
    for i in eachindex(line)
        pass!(line[i], particles6, np, particles)        
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end

function linepass!(line, particles::Beam, refpts::Vector)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    if length(particles6) != np*6
        error("The number of particles does not match the length of the particle array")
    end
    saved_particles = []
    for i in eachindex(line)
        # ele = line[i]
        pass!(line[i], particles6, np, particles)        
        if i in refpts
            push!(saved_particles, copy(array_to_matrix(particles6, np)))
        end
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return saved_particles
end


function ADlinepass!(line, particles::Beam, changed_idx::Vector, changed_ele::Vector)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    if length(particles6) != np*6
        error("The number of particles does not match the length of the particle array")
    end
    count = 1
    for i in eachindex(line)
        # ele = line[i]
        if i in changed_idx
            pass!(changed_ele[count], particles6, np, particles)
            count += 1
        else
            pass!(line[i], particles6, np, particles)        
        end
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end

function ADlinepass!(line, particles::Beam, refpts::Vector, changed_idx::Vector, changed_ele::Vector)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    if length(particles6) != np*6
        error("The number of particles does not match the length of the particle array")
    end
    count = 1
    saved_particles = []
    for i in eachindex(line)
        # ele = line[i]
        if i in changed_idx
            pass!(changed_ele[count], particles6, np, particles)
            count += 1
        else
            pass!(line[i], particles6, np, particles)        
        end
        if i in refpts
            push!(saved_particles, copy(array_to_matrix(particles6, np)))
        end

    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return saved_particles
end

function AD_ringpass!(line, particles::Beam, nturn::Int, changed_idx::Vector, changed_ele::Vector)
    for i in 1:nturn
        ADlinepass!(line, particles, changed_idx, changed_ele)    
    end
    return nothing
end
function AD_ringpass!(line, particles::Beam, nturn::Int, changed_idx::Vector, changed_ele::Vector, save::Bool)
    save_beam = []
    for i in 1:nturn
        ADlinepass!(line, particles, changed_idx, changed_ele)    
        if save
            push!(save_beam, copy(particles.r))
        end
    end
    return save_beam
end

function ringpass!(line, particles::Beam, nturn::Int)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    for i in 1:nturn
        linepass!(line, particles)    
    end
    return nothing
end

function ringpass!(line, particles::Beam, nturn::Int, save::Bool)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    save_beam = []
    for i in 1:nturn
        linepass!(line, particles)    
        if save
            push!(save_beam, copy(particles.r))
        end
    end
    return save_beam
end

function linepass_TPSA!(line, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    if length(rin) != 6
        error("The length of TPSA must be 6")
    end

    for i in eachindex(line)
        # ele = line[i]
        pass_TPSA!(line[i], rin)        
    end
    return nothing
end
function ADlinepass_TPSA!(line, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, changed_idx::Vector, changed_ele::Vector) where {T, TPS_Dim, Max_TPS_Degree}
    if length(rin) != 6
        error("The length of TPSA must be 6")
    end
    count = 1
    for i in eachindex(line)
        # ele = line[i]
        if i in changed_idx
            pass_TPSA!(changed_ele[count], rin)
            count += 1
        else
            pass_TPSA!(line[i], rin)        
        end
    end
    return nothing
end

function ringpass_TPSA!(line, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, nturn::Int) where {T, TPS_Dim, Max_TPS_Degree}
    if length(rin) != 6
        error("The length of TPSA must be 6")
    end
    for i in 1:nturn
        linepass_TPSA!(line, rin)    
    end
    return nothing
end

function check_lost(r6)
    if isnan(r6[1]) || isinf(r6[1])
        return true
    end
    if maximum(abs.(r6[1:4])) > CoordLimit || abs(r6[6]) > CoordLimit
        return true
    end
    return false
end