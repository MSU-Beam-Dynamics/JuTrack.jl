function total_length(ring)
    length = 0.0
    for i in eachindex(ring)
        length += ring[i].len
    end
    return length
end

function spos(ring)
    pos = zeros(length(ring))
    for i in eachindex(ring)
        pos[i] = total_length(ring[1:i])
    end
    return pos
end

function findelem(ring, field::Symbol, value)
    # Example: findelem(ring, :name, "QF")
    # warning: may not compatible with autodifferentiation
    ele_index = []
    for i in eachindex(ring)
        if field in fieldnames(typeof(ring[i])) && getfield(ring[i], field) == value
            push!(ele_index, i)
        end
    end
    return ele_index
end

function findelem(ring, type::Type)
    # Example: findelem(ring, DRFIT)    
    # warning: may not compatible with autodifferentiation
    ele_index = []
    for i in eachindex(ring)
        if typeof(ring[i]) == type
            push!(ele_index, i)
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

function optics(Twi)
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
    dp[:, 2] = [Twi[i].dy for i in eachindex(Twi)]
    dp[:, 3] = [Twi[i].dpx for i in eachindex(Twi)]
    dp[:, 4] = [Twi[i].dpy for i in eachindex(Twi)]
    return beta, alpha, gamma, mu, dp
end

# function plot_optics(RING)
#     Twi,s = twissring(RING, 0.0, 1)
#     using Plots
#     beta, alpha, gamma, mu, dp = optics(Twi)
#     p1=plot(s, beta[:, 1], label="betax", xlabel="s", ylabel="beta")
#     plot!(s, beta[:, 2], label="betay", xlabel="s", ylabel="beta",)
#     p2=plot(s, dp[:, 1], label="dispersion", xlabel="s", ylabel="dispersion")
#     plot(p1, p2, layout=(2,1))
# end