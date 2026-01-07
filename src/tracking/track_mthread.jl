"""
    plinepass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64})

Pass particles through the line element by element by implementing multi-threading. 
The number of threads is determined by the environment variable `JULIA_NUM_THREADS`.

# Arguments
- line::Vector{<:AbstractElement{Float64}}: a vector of beam line elements
- particles::Beam{Float64}: a beam object
"""
function plinepass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64})
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    for i in eachindex(line)
        pass_P!(line[i], particles.r, np, particles)        
    end
    return nothing
end

function ADplinepass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, 
    changed_idx::Vector{Int}, changed_ele::Vector{<:AbstractElement{Float64}})
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    count = 1
    for i in eachindex(line)
        if i in changed_idx
            pass_P!(changed_ele[count], particles.r, np, particles)
            count += 1
        else
            pass_P!(line[i], particles.r, np, particles)        
        end
    end
    return nothing
end
function ADpringpass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, nturn::Int, 
    changed_idx::Vector{Int}, changed_ele::Vector{<:AbstractElement{Float64}})
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    for i in 1:nturn
        ADplinepass!(line, particles, changed_idx, changed_ele)    
    end
    return nothing
end

"""
    pringpass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, nturn::Int)

Pass particles through the ring by implementing multi-threading. 
The number of threads is determined by the environment variable `JULIA_NUM_THREADS`.

# Arguments
- line::Vector{<:AbstractElement{Float64}}: a vector of beam line elements
- particles::Beam{Float64}: a beam object
- nturn::Int: number of turns
"""
function pringpass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, nturn::Int)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares
    # Check if the particle is lost by checking the lost_flag
    for i in 1:nturn
        plinepass!(line, particles)    
    end
    return nothing
end

function pringpass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, nturn::Int, save::Bool)
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
