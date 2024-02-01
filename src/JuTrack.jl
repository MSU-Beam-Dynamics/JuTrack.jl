module JuTrack
CoordLimit = 1.0
AngleLimit = 1.0

include("TPSA_Enzyme/TPSA_fixedmap.jl")
include("lattice/beam.jl")
include("lattice/canonical_elements.jl")
include("tracking_AT/bend_AT.jl")
include("tracking_AT/drift_AT.jl")
include("tracking_AT/multipole_AT.jl")
include("tracking_AT/rfcavity_AT.jl")
include("tracking_AT/track.jl")
include("tracking_AT/bend_TPSA_struct.jl")
include("tracking_AT/drift_TPSA_struct.jl")
include("tracking_AT/multipole_TPSA_struct.jl")
include("tracking_AT/rfcavity_TPSA_struct.jl")
include("tracking_AT/track.jl")
include("lattice/EdwardsTengTwiss.jl")


export CTPS, cst, findindex, PolyMap, getindexmap, tadd, tminus, tmult, tdiv, tpow, tsqrt, tsin, tcos, ttan, tcosh, tsinh, reassign!
export DRIFT, KQUAD, KSEXT, KOCT, SBEND, RFCA, AbstractElement
export Beam
export EdwardsTengTwiss, AbstractTwiss, twissPropagate, findm66, periodicEdwardsTengTwiss, Twissline
export linepass!, pass!, ringpass!, linepass_TPSA!, pass_TPSA!, ringpass_TPSA!

end
