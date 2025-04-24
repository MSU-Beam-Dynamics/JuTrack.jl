function get_len(ele::AbstractElement)
    return get_len_value(ele.len)
end
function get_len_value(L::Float64)
    return L[1]
end

function total_length(ring::Vector)
    leng = 0.0
    for i in eachindex(ring)
        leng += get_len(ring[i])
    end
    return leng
end

"""
    spos(ring::Vector)

Calculate the s position of each element in the lattice.

# Arguments
- ring::Vector: a vector of beam line elements

# Return
- pos::Vector{Float64}: a vector of s positions
"""
function spos(ring::Vector)
    pos = zeros(length(ring))
    len = 0.0
    for i in eachindex(ring)
        len += get_len(ring[i])
        pos[i] = len
    end
    return pos
end

"""
    spos(ring::Vector, idx::Vector)

Calculate the s position of the specified elements in the lattice.

# Arguments
- ring::Vector: a vector of beam line elements
- idx::Vector: a vector of indices of the elements

# Return
- pos::Vector{Float64}: a vector of s positions of the specified elements
"""
function spos(ring::Vector, idx::Vector)
    pos_all = spos(ring)
    pos = zeros(length(idx))
    for i in 1:length(idx)
        pos[i] = pos_all[idx[i]]
    end
    return pos
end
# function findelem(ring::Vector, field::Symbol, value)
#     # Example: findelem(ring, :name, "QF")
#     # warning: may not compatible with autodifferentiation
#     ele_index = []
#     for i in eachindex(ring)
#         if field in fieldnames(typeof(ring[i])) && getfield(ring[i], field) == value
#             push!(ele_index, i)
#         end
#     end
#     return ele_index
# end
""" 
    findelem(ring::Vector, field::Symbol, value)

Find the index of elements with the specified field value in the ring.

# Arguments
- ring::Vector: a vector of beam line elements
- field::Symbol: the field name
- value: the value of the field

# Return
- ele_index::Vector{Int}: a vector of indices of the elements with the specified field value

# Example
```julia
ele_index = findelem(ring, :name, "QF")
```
"""
function findelem(ring::Vector, field::Symbol, value)
    c = 0
    for i in eachindex(ring)
        if field in fieldnames(typeof(ring[i])) && getfield(ring[i], field) == value
            c += 1
        end
    end
    ele_index = zeros(Int, c)
    c = 0
    for i in eachindex(ring)
        if field in fieldnames(typeof(ring[i])) && getfield(ring[i], field) == value
            c += 1
            ele_index[c] = i
        end
    end
    return ele_index
end

# function findelem(ring::Vector, type::Type)
#     # Example: findelem(ring, DRFIT)    
#     # warning: may not compatible with autodifferentiation
#     ele_index = []
#     for i in eachindex(ring)
#         if typeof(ring[i]) == type
#             push!(ele_index, i)
#         end
#     end
#     return ele_index
# end
"""
    findelem(ring::Vector, type::Type)

Find the index of elements with the specified type in the ring.

# Arguments
- ring::Vector: a vector of beam line elements
- type::Type: the type of the element

# Return
- ele_index::Vector{Int}: a vector of indices of the elements with the specified type

# Example
```julia
ele_index = findelem(ring, DRIFT)
```
"""
function findelem(ring::Vector, type::Type)
    c = 0
    for i in eachindex(ring)
        if typeof(ring[i]) == type
            c += 1
        end
    end
    ele_index = zeros(Int, c)
    c = 0
    for i in eachindex(ring)
        if typeof(ring[i]) == type
            c += 1
            ele_index[c] = i
        end
    end
    return ele_index
end

function use_exact_drift(flag)
    if flag == 1
        global use_exact_Hamiltonian = 1
    else
        global use_exact_Hamiltonian = 0
    end
end

function array_optics(Twi)
    beta = zeros(length(Twi), 2)
    beta[:, 1] = [Twi[i].betax for i in eachindex(Twi)]
    beta[:, 2] = [Twi[i].betay for i in eachindex(Twi)]
    alpha = zeros(length(Twi), 2)
    alpha[:, 1] = [Twi[i].alphax for i in eachindex(Twi)]
    alpha[:, 2] = [Twi[i].alphay for i in eachindex(Twi)]
    gamma = zeros(length(Twi), 2)
    gamma[:, 1] = [Twi[i].gammax for i in eachindex(Twi)]
    gamma[:, 2] = [Twi[i].gammay for i in eachindex(Twi)]
    mu = zeros(length(Twi), 2)
    mu[:, 1] = [Twi[i].dmux for i in eachindex(Twi)]
    mu[:, 2] = [Twi[i].dmuy for i in eachindex(Twi)]
    dp = zeros(length(Twi), 4)
    dp[:, 1] = [Twi[i].dx for i in eachindex(Twi)]
    dp[:, 2] = [Twi[i].dpx for i in eachindex(Twi)]
    dp[:, 3] = [Twi[i].dy for i in eachindex(Twi)]
    dp[:, 4] = [Twi[i].dpy for i in eachindex(Twi)]
    return beta, alpha, gamma, mu, dp
end


function symplectic(M66::Array{Float64,2})
    # check if a transfer map is symplectic
    # the canonical coordinates are (x, px/p0, y, py/p0, tau, -dE/beta0/c/p0)
    # z = tau, dp/p0 ~ dE/beta0/c/p0. Therefore, the last block of J is -J2
    J2 = [0.0 1.0; -1.0 0.0]
    J = [J2 zeros((2,2)) zeros((2,2)); zeros((2,2)) J2 zeros((2,2)); zeros((2,2)) zeros((2,2)) -J2]
    delta = transpose(M66) * J * M66 .- J
    println("max deviation: ", maximum(abs.(delta)))
    return maximum(abs.(delta))
end