include("../lattice/canonical_elements.jl")
include("../tracking/EDrift.jl")
include("../tracking/csbend.jl")
include("../tracking/multipole.jl")
include("../tracking/rfca.jl")

function linepass(line, particles, Po, p_error, sigmaDelta2)
    nline = length(line)
    np = length(particles)
    z = 0.0
    for i in 1:nline
        ele = line[i]
        see_rf = 0
        if ele isa EDRIFT
            particles = EDrift(particles, np, ele.len)
        elseif ele isa KQUAD || ele isa KSEXT || ele isa KOCT
            particles = multipole_tracking2(particles, np, ele, Po)
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
                see_rf = 1
            end
            particles, fiducial_seen, Po, n_references, reference_ref_number, reference_phase, reference_flags = 
                simple_rf_cavity(particles, np, ele, Po, z, fiducial_seen, n_references, reference_ref_number, reference_phase, 
                                reference_flags, n_references)
        else
            error("Unknown element type")
        end
        z = z + ele.len
    end
    return particles
end

