include("../src/JuTrack.jl")
using .JuTrack
using Distributions, Plots
using ProgressMeter
include("../src/demo/ssrf_ring.jl")
# RING = ssrf(-1.063770, 0)

function insert_space_charge(lattice, dphi, a, b, Nl, Nm)
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
    # new list of AbstractElement
    new_lattice = Vector{AbstractElement}()
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


# new_RING = insert_space_charge(RING, pi/4, 10e-3, 10e-3, 15, 15)

function f(x1)
    D1 = DRIFT(len=0.2)
    D2 = DRIFT(len=0.4)
    D3 = DRIFT(len=0.2)
    Q1 = KQUAD(len=0.1, k1=x1)
    Q2 = KQUAD(len=0.1, k1=-x1)
    
    a = 13e-3
    b = 13e-3
    nl = 12
    nm = 12
    SC_D1 = SPACECHARGE(effective_len=0.2, a=a, b=b, Nl=nl, Nm=nm)
    SC_D2 = SPACECHARGE(effective_len=0.4, a=a, b=b, Nl=nl, Nm=nm)
    SC_D3 = SPACECHARGE(effective_len=0.2, a=a, b=b, Nl=nl, Nm=nm)
    SC_Q1 = SPACECHARGE(effective_len=0.1, a=a, b=b, Nl=nl, Nm=nm)
    SC_Q2 = SPACECHARGE(effective_len=0.1, a=a, b=b, Nl=nl, Nm=nm)
    
    
    line_sc = [D1, SC_D1, Q1, SC_Q1, D2, SC_D2, Q2, SC_Q2, D3, SC_D3]
    
    beam = Beam(zeros(10000, 6), energy=1.0e5, current=200.0, mass=m_p, charge=1.0, emittance=[1e-6, 1e-6, 0.0])
    beta = beam.beta
    gamma = beam.gamma
    emit_norm = 1e-6
    emit_phys = emit_norm # / (beta * gamma)
    beam.emittance = [emit_phys, emit_phys, 0.0]
    
    vbase=3.42*8.5e6
    ϕs=10.0
    vact=vbase/cos(ϕs*π/180.0)
    freq = 591e6
    mainRFe=AccelCavity(freq, vact, 7560.0, π-ϕs*π/180.0)
    tunex, tuney=50.08, 44.14
    αc=3.42/tunex/tunex
    lmap=LongitudinalRFMap(αc, mainRFe)
    opt=optics4DUC(1.0, 0.0, 1.0, 0.0)
    initilize_6DGaussiandist!(beam, opt, lmap)
    
    N = 5000
    for i in 1:N
        linepass!(line_sc, beam)
        if i % 100 == 0
            println("i: ", i)
        end
    end
    get_emittance!(beam)
    emit = beam.emittance
    return emit
end
    
k1 = 29.6
grad = autodiff(Forward, f, Duplicated, Duplicated(k1, 1.0))
k2 = 31.6
grad2 = autodiff(Forward, f, Duplicated, Duplicated(k2, 1.0))