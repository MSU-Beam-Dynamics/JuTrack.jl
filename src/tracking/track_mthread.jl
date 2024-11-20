"""
    plinepass!(line::Vector, particles::Beam)

Pass particles through the line element by element by implementing multi-threading. The number of threads is determined by the environment variable `JULIA_NUM_THREADS`.

# Arguments
- line::Vector: a vector of beam line elements
- particles::Beam: a beam object
"""
function plinepass!(line::Vector, particles::Beam)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    if length(particles6) != np*6
        error("The number of particles does not match the length of the particle array")
    end
    for i in eachindex(line)
        pass_P!(line[i], particles6, np, particles)        
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end

function ADplinepass!(line::Vector, particles::Beam, changed_idx::Vector, changed_ele::Vector)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    # if length(particles6) != np*6
    #     error("The number of particles does not match the length of the particle array")
    # end
    count = 1
    for i in eachindex(line)
        if i in changed_idx
            pass_P!(changed_ele[count], particles6, np, particles)
            count += 1
        else
            pass_P!(line[i], particles6, np, particles)        
        end
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end
function ADpringpass!(line::Vector, particles::Beam, nturn::Int, changed_idx::Vector, changed_ele::Vector)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    for i in 1:nturn
        ADplinepass!(line, particles, changed_idx, changed_ele)    
    end
    return nothing
end

"""
    pringpass!(line::Vector, particles::Beam, nturn::Int)

Pass particles through the ring by implementing multi-threading. The number of threads is determined by the environment variable `JULIA_NUM_THREADS`.

# Arguments
- line::Vector: a vector of beam line elements
- particles::Beam: a beam object
- nturn::Int: number of turns
"""
function pringpass!(line::Vector, particles::Beam, nturn::Int)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    for i in 1:nturn
        plinepass!(line, particles)    
    end
    return nothing
end

function pringpass!(line::Vector, particles::Beam, nturn::Int, save::Bool)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    save_beam = []
    for i in 1:nturn
        plinepass!(line, particles)    
        if save
            push!(save_beam, copy(particles.r))
        end
    end
    return save_beam
end
