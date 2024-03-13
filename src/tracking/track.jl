function matrix_to_array(matrix::Matrix{Float64})
    particles = zeros(Float64, size(matrix, 1)*size(matrix, 2))
    for i in 1:size(matrix, 1)
        for j in 1:size(matrix, 2)
            particles[(i-1)*size(matrix, 2)+j] = matrix[i, j]
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
        # ele = line[i]
        pass!(line[i], particles6, np, particles)        
        if isnan(particles6[1]) || isinf(particles6[1])
            println("The particle is lost at element ", i, "element name is ", line[i].name)
            rout = array_to_matrix(particles6, np)
            particles.r = rout
            return nothing
        end
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end

function ADlinepass!(line, particles::Beam, changed_idx::Vector{Int}, changed_ele)
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
        if isnan(particles6[1]) || isinf(particles6[1])
            println("The particle is lost at element ", i, "element name is ", line[i].name)
            rout = array_to_matrix(particles6, np)
            particles.r = rout
            return nothing
        end
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
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
function ADlinepass_TPSA!(line, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, changed_idx::Vector{Int}, changed_ele) where {T, TPS_Dim, Max_TPS_Degree}
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

