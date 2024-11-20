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

function spos(ring::Vector)
    pos = zeros(length(ring))
    len = 0.0
    for i in eachindex(ring)
        len += get_len(ring[i])
        pos[i] = len
    end
    return pos
end

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

function insert_space_charge(lattice, dphi, a, b, Nl, Nm)
    # Insert space charge elements into the lattice every dphi phase advance
    twi = twissring(lattice, 0.0, 1)
    phi0 = 0.0
    s0 = 0.0
    s = spos(lattice)
    insert_idx = [] 
    insert_len = []

    if twi[end].dmux < dphi
        println("The phase advance of the whole ring is less than ", dphi, ", only one space charge element is inserted.")
        sc = SPACECHARGE(effective_len=s[end], a=a, b=b, Nl=Nl, Nm=Nm)
        new_lattice = [lattice..., sc]
        return new_lattice
    end
    for i in 1:length(lattice)
        if twi[i].dmux - phi0 > dphi
            len = s[i] - s0
            push!(insert_idx, i)
            push!(insert_len, len)
            s0 = s[i]
            phi0 = twi[i].dmux
        end
    end

    new_lattice = []
    for i in 1:length(lattice)
        push!(new_lattice, lattice[i])
        if i in insert_idx
            num = findfirst(x -> x == i, insert_idx)
            sc = SPACECHARGE(effective_len=insert_len[num], a=a, b=b, Nl=Nl, Nm=Nm)
            push!(new_lattice, sc)
        end
    end
    println("Space charge elements are inserted at: ", insert_idx)
    return new_lattice
end
