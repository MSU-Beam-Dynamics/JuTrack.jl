mutable struct Beam
    r::Matrix{Float64} 
    np::Int
    energy::Float64
    lost_flag::Vector{Int}
end
function Beam(r::Matrix{Float64}, energy::Float64)
    np = size(r, 1)
    lost_flag = zeros(Int, np)
    return Beam(r, np, energy, lost_flag)
end
function Beam(r::Matrix{Float64})
    np = size(r, 1)
    energy = 1e9
    lost_flag = zeros(Int, np)
    return Beam(r, np, energy, lost_flag)
end