include("../lattice/canonical_elements.jl")
include("../tracking/EDrift.jl")
include("../tracking/csbend.jl")
include("../tracking/multipole.jl")

function linepass(line, particles, Po, p_error, sigmaDelta2)
    nline = length(line)
    np = length(particles)
    for i in 1:nline
        ele = line[i]
        if ele isa EDRIFT
            particles = EDrift(particles, np, ele.len)
        elseif ele isa KQUAD || ele isa KSEXT || ele isa KOCT
            particles = multipole_tracking2(particles, np, ele, Po)
        elseif ele isa CSBEND
            particles = track_through_csbend(particles, np, ele, p_error, Po, sigmaDelta2)
        else
            error("Unknown element type")
        end
    end
    return particles
end

parts = [Float64[0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0], Float64[0.001, 0.0, 0.0, 0.0, 0.0, 0.0]]
D1 = EDRIFT("D1", 0.5)
D2 = EDRIFT("D2", 0.5)
Q1 = KQUAD("Q1", 0.5, 2)
Q2 = KQUAD("Q2", 0.5, -3)
B1 = CSBEND("B1", 0.5, 0.1)
line = [D1, Q1, D2, Q2, B1]
Po = 1000.0
p_error = 0.0
sigmaDelta2 = 0.0
pout = linepass(line, parts, Po, p_error, sigmaDelta2)
println(pout)