include("../lattice/canonical_elements.jl")
include("EDrift_Enzyme.jl")
include("multipole_Enzyme.jl")
include("multipole_TPSA_Enzyme.jl")
include("csbend_Enzyme.jl")
include("csbend_TPSA_Enzyme.jl")
include("rfca_Enzyme.jl")
include("rfca_TPSA_Enzyme.jl")
include("../TPSA_Enzyme/TPSA_fixedmap.jl")


function linepass!(line, particles, Po, p_error)
    nline = length(line)
    np = length(particles)
    z = 0.0
    see_rf = 0
    sigmaDelta2 = 0.0
    dzLoss = 0.0
    for i in 1:nline
        # ele = line[i]
        pass!(particles, line[i], Po, sigmaDelta2)        
        # if ele isa EDRIFT
        #     EDrift(particles, np, ele.len)
        # elseif ele isa KQUAD || ele isa KSEXT || ele isa KOCT
        #     multipole_tracking2(particles, np, ele, Po, sigmaDelta2)
        # elseif ele isa CSBEND
        #     particles = track_through_csbend(particles, np, ele, p_error, Po, sigmaDelta2)
        # elseif ele isa RFCA
        #     if see_rf == 0
        #         fiducial_seen = 0
        #         n_references = 0.0
        #         reference_ref_number = 0
        #         reference_phase = 0.0
        #         reference_flags = 0
        #         n_references = 0
        #         rfca_phase_reference = ele.phase_reference
        #         see_rf = 1
        #     end
        #     particles, fiducial_seen, Po, n_references, rfca_phase_reference, reference_ref_number, reference_phase, reference_flags = 
        #         simple_rf_cavity(particles, np, ele, Po, z, fiducial_seen, rfca_phase_reference, reference_ref_number, 
        #         reference_phase, reference_flags, n_references)
        # else
        #     error("Unknown element type")
        # end
        # z += line[i].len
    end
    # return particles
end

function linepass_TPSA(line, x, xp, y, yp, z, delta, Po, p_error)
    nline = length(line)
    # z = 0.0
    see_rf = 0
    sigmaDelta2 = 0.0
    dzLoss = 0.0
    for i in 1:nline
        x, xp, y, yp, z, delta = pass_TPSA(x, xp, y, yp, z, delta, line[i], Po, sigmaDelta2)
    end
    return x, xp, y, yp, z, delta
end

# function test_f(x)
#     D1 = EDRIFT(len=1.0)
#     D2 = EDRIFT(len=1.0)
#     QF = KQUAD(name="QF",len=1.0,k1=x[1])
#     QD = KQUAD(name="QD",len=1.0,k1=-x[1])
#     S = KSEXT(len=1.0,k2=x[2])
#     freq = 500e6
#     volt = 1000000.0
#     rfca1 = RFCA(name="rfca", freq=freq, volt=volt, phase=90.0, len=0.1, nSlices=10)
#     rfca2 = RFCA(name="rfca", freq=freq, volt=volt, phase=90.0, len=0.1, nSlices=10)
#     line = [rfca1, rfca2] # it works if we only use 1 element, or we use two same elements. It fails if we use two different elements.
#     # particles = [Float64[0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0], Float64[0.001, 0.0, 0.0, 0.0, 1.0, 0.0]]
#     particles = [0.001 0.0001 0.0005 0.0002 0.0 0.0; 0.001 0.0 0.0 0.0 0.0 0.0]
#     Po = 19569.50762296901
#     linepass!(line, particles, Po, 0.0)
#     return particles[1,:]
# end
# println(test_f([0.1, 0.1]))
# println((test_f([0.1, 0.1])[6]+1.0)*19569.50762296901)


# function test_TPSA(xx)
#     D1 = EDRIFT(len=1.0)
#     D2 = EDRIFT(len=1.0)
#     QF = KQUAD(name="QF",len=1.0,k1=xx[1],synch_rad=0)
#     QD = KQUAD(name="QD",len=1.0,k1=xx[2],synch_rad=0)
#     S = KSEXT(len=1.0,k2=xx[3],synch_rad=0)
#     O = KOCT(len=1.0,k3=xx[4],synch_rad=0)
#     B = CSBEND(len=1.0, angle=xx[5], synch_rad=1)
#     rfca1 = RFCA(name="rfca", freq=xx[6], volt=xx[7], phase=90.0, len=0.1, nSlices=10)
#     line = [QF,D1,QD,D2,S,O,B,rfca1] # it works if we only use 1 element, or we use two same elements. It fails if we use two different elements.
#     x = CTPS(0.0, 1, 6, 3)
#     xp = CTPS(0.0, 2, 6, 3)
#     y = CTPS(0.0, 3, 6, 3)
#     yp = CTPS(0.0, 4, 6, 3)
#     delta = CTPS(0.0, 5, 6, 3)
#     z = CTPS(0.0, 6, 6, 3) 
#     Po = 19569.50762296901
#     x, xp, y, yp, z, delta = linepass_TPSA(line, x, xp, y, yp, z, delta, Po, 0.0)
#     println(x.map)
#     return x.map[3]
# end
# using Enzyme
# using BenchmarkTools
# Enzyme.API.runtimeActivity!(true)
# Enzyme.API.strictAliasing!(false)
x = [2.0, -2.0, 2.0, -1.0, 0.1571, 500e6, 1e6]
# println(test_TPSA(x))
# @btime begin
#     grad = gradient(Forward, test_TPSA, x)
# end
# println(grad)
# include("ssrf_ring.jl")

# function test_TPSA_ssrf(xx)
#     # line = ssrf(K1)
#     # D1 = EDRIFT(name="D1", len=xx[2])
#     # DL = EDRIFT(name="DL", len=6.0)    
#     # QL1 = KQUAD(name="QL1", len=0.32, k1=xx[1],  synch_rad=1)
#     line = ssrf(xx[1], xx[2])
#     x = CTPS(0.0, 1, 6, 3)
#     xp = CTPS(0.0, 2, 6, 3)
#     y = CTPS(0.0, 3, 6, 3)
#     yp = CTPS(0.0, 4, 6, 3)
#     delta = CTPS(0.0, 5, 6, 3)
#     z = CTPS(0.0, 6, 6, 3) 
#     Po = 19569.50762296901
#     x, xp, y, yp, z, delta = linepass_TPSA(line, x, xp, y, yp, z, delta, Po, 0.0)
#     # println(x.map)
#     return x.map[3]
# end
# x = [-1.063770, 0.34]
# println(test_TPSA_ssrf(x))
# grad = gradient(Forward, test_TPSA_ssrf, x)
# println(grad)
# @btime begin
#     grad = gradient(Forward, test_TPSA_ssrf, [-1.063770])
# end