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

# function linepass!(line, particles::Beam)
#     # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
#     # Check if the particle is lost by checking the lost_flag
#     np = particles.nmacro
#     particles6 = matrix_to_array(particles.r)
#     if length(particles6) != np*6
#         error("The number of particles does not match the length of the particle array")
#     end
#     noTarray = zeros(6)
#     noRmatrix = zeros(6, 6)
#     for i in eachindex(line)
#         # ele = line[i]
#         pass!(line[i], particles6, np, particles, noTarray, noRmatrix)        
#     end
#     rout = array_to_matrix(particles6, np)
#     particles.r = rout
#     return nothing
# end

function linepass!(line, particles::Beam)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    if length(particles6) != np*6
        error("The number of particles does not match the length of the particle array")
    end
    noTarray = zeros(6)
    noRmatrix = zeros(6, 6)
    for i in eachindex(line)
        # ele = line[i]
        if line[i] isa AbstractElement
            pass!(line[i], particles6, np, particles, noTarray, noRmatrix)      
        else 
            for j in eachindex(line[i])
                if line[i][j] isa AbstractElement
                    pass!(line[i][j], particles6, np, particles, noTarray, noRmatrix)      
                else
                    for k in eachindex(line[i][j])
                        if line[i][j][k] isa AbstractElement
                            pass!(line[i][j][k], particles6, np, particles, noTarray, noRmatrix)      
                        else
                           error("The element is not an AbstractElement")
                        end
                    end
                end     
            end  
        end     
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end


function ringpass!(line::Vector{AbstractElement}, particles::Beam, nturn::Int)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    for i in 1:nturn
        linepass!(line, particles)    
    end
    return nothing
end

function linepass_TPSA!(line::Vector{AbstractElement}, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    if length(rin) != 6
        error("The length of TPSA must be 6")
    end
    noTarray = zeros(6)
    noRmatrix = zeros(6, 6)
    for i in eachindex(line)
        # ele = line[i]
        pass_TPSA!(line[i], rin, noTarray, noRmatrix)        
    end
    return nothing
end

function ringpass_TPSA!(line::Vector{AbstractElement}, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, nturn::Int) where {T, TPS_Dim, Max_TPS_Degree}
    if length(rin) != 6
        error("The length of TPSA must be 6")
    end
    for i in 1:nturn
        linepass_TPSA!(line, rin)    
    end
    return nothing
end

