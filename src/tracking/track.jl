include("../lattice/canonical_elements.jl")
include("../tracking/EDrift.jl")
include("../tracking/csbend.jl")
include("../tracking/multipole.jl")
include("../tracking/rfca.jl")
include("../tracking/multipole_TPSA.jl")
include("../tracking/csbend_TPSA.jl")
include("../tracking/rfca_TPSA.jl")

function linepass(line, particles, Po, p_error, sigmaDelta2)
    nline = length(line)
    np = length(particles)
    z = 0.0
    see_rf = 0
    sigmaDelta2 = nothing
    for i in 1:nline
        ele = line[i]
        
        if ele isa EDRIFT
            particles = EDrift(particles, np, ele.len)
        elseif ele isa KQUAD || ele isa KSEXT || ele isa KOCT
            particles = multipole_tracking2(particles, np, ele, Po, sigmaDelta2)
        elseif ele isa CSBEND
            particles = track_through_csbend(particles, np, ele, p_error, Po, sigmaDelta2)
        elseif ele isa RFCA
            if see_rf == 0
                fiducial_seen = 0
                n_references = 0.0
                reference_ref_number = 0
                reference_phase = 0.0
                reference_flags = 0
                n_references = 0
                rfca_phase_reference = ele.phase_reference
                see_rf = 1
            end
            particles, fiducial_seen, Po, n_references, rfca_phase_reference, reference_ref_number, reference_phase, reference_flags = 
                simple_rf_cavity(particles, np, ele, Po, z, fiducial_seen, rfca_phase_reference, reference_ref_number, 
                reference_phase, reference_flags, n_references)
        else
            error("Unknown element type")
        end
        z = z + ele.len
    end
    return particles
end

function linepass(line, x, xp, y, yp, z, delta, Po, p_error, sigmaDelta2)
    nline = length(line)
    z_total = 0.0
    see_rf = 0
    for i in 1:nline
        ele = line[i]
        
        if ele isa EDRIFT
            x, xp, y, yp, z, delta = EDrift(x, xp, y, yp, z, delta, ele.len)
        elseif ele isa KQUAD || ele isa KSEXT || ele isa KOCT
            x, xp, y, yp, z, delta = multipole_tracking2(x, xp, y, yp, z, delta, 1, ele, Po)
        elseif ele isa CSBEND
            x, xp, y, yp, z, delta = track_through_csbend(x, xp, y, yp, z, delta, 1, ele, p_error, Po, sigmaDelta2)
        elseif ele isa RFCA
            if see_rf == 0
                fiducial_seen = 0
                n_references = 0.0
                reference_ref_number = 0
                reference_phase = 0.0
                reference_flags = 0
                n_references = 0
                rfca_phase_reference = ele.phase_reference
                see_rf = 1
            end
            x, xp, y, yp, z, delta, fiducial_seen, Po, n_references, rfca_phase_reference, reference_ref_number, reference_phase, reference_flags = 
                simple_rf_cavity(x, xp, y, yp, z, delta, 
                ele, Po, z_total, fiducial_seen, rfca_phase_reference, reference_ref_number, 
                reference_phase, reference_flags, n_references)
        else
            error("Unknown element type")
        end
        # println("No. ", i, " element: ", ele.name)
        z_total = z_total + ele.len
    end
    return x, xp, y, yp, z, delta
end

# function f(volt, freq)
#     parts = [Float64[0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0], Float64[0.001, 0.0, 0.0, 0.0, 0.0, 0.0]]
#     D1 = EDRIFT("D1", 0.5)
#     D2 = EDRIFT("D2", 0.5)
#     Q1 = KQUAD("Q1", 0.5, 2, 0.0, 0.0, 0.0, 0, 0, 10)
#     Q2 = KQUAD("Q2", 0.5, -3, 0.0, 0.0, 0.0, 0, 0, 10)
#     B1 = CSBEND("B1", 0.5, 0.1, 0.0, 0.0, 10)
#     RFC = RFCA("RFCA", 0.1, freq, volt, 90, 10)
#     line = [D1, Q1, D2, Q2, B1, RFC]
#     Po = 19569.507622969009
#     p_error = 0.0
#     sigmaDelta2 = 0.0
#     pout = linepass(line, parts, Po, p_error, sigmaDelta2)
#     return pout[1]
# end
# println("Forward pass: ", f(1e6, 500e6))
# using Zygote
# grad = Zygote.jacobian(f, 1e6, 500e6)
# println("Gradient: ", grad)
# freq = 500e6
# volt = 1e6
# function f(freq, volt)
#     parts = [Float64[0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0], Float64[0.001, 0.0, 0.0, 0.0, 0.0, 0.0]]
#     D1 = EDRIFT("D1", 0.5)
#     D2 = EDRIFT("D2", 0.5)
#     Q1 = KQUAD("Q1", 0.5, 2, 0.0, 0.0, 0.0, 0, 0, 10)
#     Q2 = KQUAD("Q2", 0.5, -3, 0.0, 0.0, 0.0, 0, 0, 10)
#     B1 = CSBEND("B1", 0.5, 0.1, 0.0, 0.0, 10)
#     RFC = RFCA("RFCA", 0.1, freq, volt, 90, 10)
#     line = [D1, Q1, D2, Q2, B1, RFC]
# x = CTPS(0.0, 1, 6, 2)
# xp = CTPS(0.0, 2, 6, 2)
# y = CTPS(0.0, 3, 6, 2)
# yp = CTPS(0.0, 4, 6, 2)
# delta = CTPS(0.0, 5, 6, 2)
# z = CTPS(0.0, 6, 6, 2)
# Po = 19569.507622969009
# x, xp, y, yp, z, delta = linepass(line, x, xp, y, yp, z, delta, Po, 0.0, 0.0)
# println(x)
# println(xp)
# println(y)
# println(yp)
# println(z)
# println(delta)

# xvalue = evaluate(x, [0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0])
# yvalue = evaluate(y, [0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0])
# println(xvalue)
# println(yvalue)
# Map66 = [x.map[2] x.map[3] x.map[4] x.map[5] x.map[6] x.map[7];
# xp.map[2] xp.map[3] xp.map[4] xp.map[5] xp.map[6] xp.map[7];
# 			y.map[2] y.map[3] y.map[4] y.map[5] y.map[6] y.map[7];
# 			yp.map[2] yp.map[3] yp.map[4] yp.map[5] yp.map[6] yp.map[7];
# 			delta.map[2] delta.map[3] delta.map[4] delta.map[5] delta.map[6] delta.map[7];
# 			z.map[2] z.map[3] z.map[4] z.map[5] z.map[6] z.map[7]]
# println(Map66)
# return x.map
# end

# grad = Zygote.jacobian(f, 500e6, 1e6)
# println("Gradient: ", grad)