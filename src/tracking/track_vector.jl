"""
    linepass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64})

Pass the beam through the line element by element.

# Arguments
- line::Vector{<:AbstractElement{Float64}}: a vector of beam line elements
- particles::Beam{Float64}: a beam object
"""
function linepass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64})
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    for i in eachindex(line)
        pass!(line[i], particles.r, np, particles)        
    end
    return nothing
end

"""
    linepass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, refpts::Vector{Int})

Pass particles through the line element by element. Save the particles at the reference points.

# Arguments
- line::Vector{<:AbstractElement{Float64}}: a vector of beam line elements
- particles::Beam{Float64}: a beam object
- refpts::Vector{Int}: a vector of reference points

# Returns
- saved_particles::Vector: a vector of saved particles at the reference points
"""
function linepass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, refpts::Vector{Int})
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    saved_particles = []
    for i in eachindex(line)
        pass!(line[i], particles.r, np, particles)        
        if i in refpts
            push!(saved_particles, copy(particles.r))
        end
    end
    return saved_particles
end

"""
    ADlinepass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, 
    changed_idx::Vector{Int}, changed_ele::Vector{<:AbstractElement{Float64}})

Pass particles through the line element by element. The elements in the `changed_idx` will be replaced by the elements in `changed_ele`.
This is a convinent function to implement automatic differentiation that avoid directly changing parameters in `line``.

# Arguments
- line::Vector{<:AbstractElement{Float64}}: a vector of beam line elements
- particles::Beam{Float64}: a beam object
- changed_idx::Vector{Int}: a vector of indices of the elements to be changed
- changed_ele::Vector{<:AbstractElement{Float64}}: a vector of elements to replace the elements in `changed_idx`
"""
function ADlinepass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, 
    changed_idx::Vector{Int}, changed_ele::Vector{<:AbstractElement{Float64}})
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    count = 1
    for i in eachindex(line)
        if i in changed_idx
            pass!(changed_ele[count], particles.r, np, particles)
            count += 1
        else
            pass!(line[i], particles.r, np, particles)        
        end
    end
    return nothing
end

function ADlinepass!(line::Vector{<:AbstractElement{Float64}}, id_list::Vector{Int}, particles::Beam{Float64}, 
    changed_idx::Vector{Int}, changed_ele::Vector{<:AbstractElement{Float64}})
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    count = 1
    for i in eachindex(line)
        if i in id_list
            if i in changed_idx
                pass!(changed_ele[count], particles.r, np, particles)
                count += 1
            else
                pass!(line[i], particles.r, np, particles)        
            end
        end
    end
    return nothing
end

function ADlinepass!(line::Vector, particles::Beam{Float64}, refpts::Vector, changed_idx::Vector, changed_ele::Vector)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    count = 1
    saved_particles = []
    for i in eachindex(line)
        if i in changed_idx
            pass!(changed_ele[count], particles.r, np, particles)
            count += 1
        else
            pass!(line[i], particles.r, np, particles)        
        end
        if i in refpts
            push!(saved_particles, copy(particles.r))
        end
    end
    return saved_particles
end


"""
    ringpass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, nturn::Int)

Pass particles through the ring for `nturn` turns.

# Arguments
- line::Vector{<:AbstractElement{Float64}}: a vector of beam line elements
- particles::Beam{Float64}: a beam object
- nturn::Int: number of turns
"""
function ringpass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, nturn::Int)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    for i in 1:nturn
        linepass!(line, particles)    
    end
    return nothing
end

"""
    ringpass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, nturn::Int, save::Bool)

Pass particles through the ring for `nturn` turns. Save the particles at each turn.

# Arguments
- line::Vector{<:AbstractElement{Float64}}: a vector of beam line elements
- particles::Beam{Float64}: a beam object
- nturn::Int: number of turns
- save::Bool: Flag
"""
function ringpass!(line::Vector{<:AbstractElement{Float64}}, particles::Beam{Float64}, nturn::Int, save::Bool)
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

"""
    ADringpass!(line::Vector, particles::Beam{Float64}, nturn::Int, changed_idx::Vector, changed_ele::Vector)

Pass particles through the ring for `nturn` turns. The elements in the `changed_idx` will be replaced by the elements in `changed_ele`.
This is a convinent function to implement automatic differentiation that avoid directly changing parameters in `line``.

# Arguments
- line::Vector: a vector of a ring
- particles::Beam{Float64}: a beam object
- nturn::Int: number of turns
- changed_idx::Vector: a vector of indices of the elements to be changed
- changed_ele::Vector: a vector of elements to replace the elements in `changed_idx`
"""
function ADringpass!(line::Vector, particles::Beam{Float64}, nturn::Int, changed_idx::Vector, changed_ele::Vector)
    for i in 1:nturn
        ADlinepass!(line, particles, changed_idx, changed_ele)    
    end
    return nothing
end
function ADringpass!(line::Vector, particles::Beam{Float64}, nturn::Int, changed_idx::Vector, changed_ele::Vector, save::Bool)
    save_beam = []
    for i in 1:nturn
        ADlinepass!(line, particles, changed_idx, changed_ele)    
        if save
            push!(save_beam, copy(particles.r))
        end
    end
    return save_beam
end

"""
    linepass_TPSA!(line::Vector{<:AbstractElement{Float64}}, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}};
    E0::Float64=3e9, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}

Pass 6-D high-order TPSA coordinates through the line element by element.

# Arguments
- line::Vector{<:AbstractElement{Float64}}: a vector of beam line elements
- rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}: a vector of 6-D high-order TPSA coordinates
- E0::Float64=3e9: reference energy in eV
- m0::Float64=m_e: rest mass in eV/c^2
"""
function linepass_TPSA!(line::Vector{<:AbstractElement{Float64}}, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; 
    E0::Float64=3e9, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    if length(rin) != 6
        error("The length of TPSA must be 6")
    end

    for i in eachindex(line)
        # ele = line[i]
        pass_TPSA!(line[i], rin, E0=E0, m0=m0)       
        if isnan(rin[1].map[2])
            println("The particle is lost at element $(i)")
        end 
    end
    return nothing
end

"""
    ADlinepass_TPSA!(line::Vector{<:AbstractElement{Float64}}, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, 
    changed_idx::Vector, changed_ele::Vector; E0::Float64=3e9, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}

Pass 6-D high-order TPSA coordinates through the line element by element. The elements in the `changed_idx` will be replaced by the elements in `changed_ele`.
This is a convinent function to implement automatic differentiation that avoid directly changing parameters in `line``.

# Arguments
- line::Vector: a vector of beam line elements
- rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}: a vector of 6-D high-order TPSA coordinates
- changed_idx::Vector: a vector of indices of the elements to be changed
- changed_ele::Vector: a vector of elements to replace the elements in `changed_idx`
- E0::Float64=3e9: reference energy in eV
- m0::Float64=m_e: rest mass in eV/c^2
"""
function ADlinepass_TPSA!(line::Vector{<:AbstractElement{Float64}}, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, 
    changed_idx::Vector, changed_ele::Vector; E0::Float64=3e9, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    if length(rin) != 6
        error("The length of TPSA must be 6")
    end
    count = 1
    for i in eachindex(line)
        # ele = line[i]
        if i in changed_idx
            pass_TPSA!(changed_ele[count], rin, E0=E0, m0=m0)
            count += 1
        else
            pass_TPSA!(line[i], rin, E0=E0, m0=m0)        
        end
    end
    return nothing
end

"""
    ringpass_TPSA!(line::Vector{<:AbstractElement{Float64}}, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, nturn::Int;
    E0::Float64=3e9, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}

Pass 6-D high-order TPSA coordinates through the ring for `nturn` turns.

# Arguments
- line::Vector: a vector of beam line elements
- rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}: a vector of 6-D high-order TPSA coordinates
- nturn::Int: number of turns
- E0::Float64=3e9: reference energy in eV
- m0::Float64=m_e: rest mass in eV/c^2
"""
function ringpass_TPSA!(line::Vector{<:AbstractElement{Float64}}, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, nturn::Int; 
    E0::Float64=3e9, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    if length(rin) != 6
        error("The length of TPSA must be 6")
    end
    for i in 1:nturn
        linepass_TPSA!(line, rin, E0=E0, m0=m0)    
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
    # sqrt(1.0 + 2.0*r[6]*beti + r[6]^2 - r[2]^2 - r[4]^2) must be real
    if r6[2]^2 + r6[4]^2 > 1.0  + r6[6]^2
        return true
    end
    return false
end