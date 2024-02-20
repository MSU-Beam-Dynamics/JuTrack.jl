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

function plinepass!(line, particles::Beam)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    if length(particles6) != np*6
        error("The number of particles does not match the length of the particle array")
    end
    for i in eachindex(line)
        # ele = line[i]
        pass_P!(line[i], particles6, np, particles)        
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end

# function plinepass!(line, particles::Beam)
#     # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
#     # Check if the particle is lost by checking the lost_flag
#     np = particles.nmacro
#     particles6 = matrix_to_array(particles.r)
#     if length(particles6) != np*6
#         error("The number of particles does not match the length of the particle array")
#     end

#     for i in eachindex(line)
#         # ele = line[i]
#         if line[i] isa AbstractElement
#             pass_P!(line[i], particles6, np, particles)      
#         else 
#             for j in eachindex(line[i])
#                 if line[i][j] isa AbstractElement
#                     pass_P!(line[i][j], particles6, np, particles)      
#                 else
#                     for k in eachindex(line[i][j])
#                         if line[i][j][k] isa AbstractElement
#                             pass_P!(line[i][j][k], particles6, np, particles)      
#                         else
#                            error("The element is not an AbstractElement")
#                         end
#                     end
#                 end     
#             end  
#         end     
#     end
#     rout = array_to_matrix(particles6, np)
#     particles.r = rout
#     return nothing
# end

function pringpass!(line::Vector{AbstractElement}, particles::Beam, nturn::Int)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    for i in 1:nturn
        plinepass!(line, particles)    
    end
    return nothing
end


