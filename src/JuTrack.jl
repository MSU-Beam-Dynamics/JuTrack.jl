module JuTrack
CoordLimit = 1.0
AngleLimit = 1.0

include("TPSA_Enzyme/TPSA_fixedmap.jl")
include("lattice/beam.jl")
include("lattice/canonical_elements.jl")
include("tracking/bend.jl")
include("tracking/drift.jl")
include("tracking/multipole.jl")
include("tracking/rfcavity.jl")
include("tracking/track.jl")
include("tracking/bend_TPSA.jl")
include("tracking/drift_TPSA.jl")
include("tracking/multipole_TPSA.jl")
include("tracking/rfcavity_TPSA.jl")
include("tracking/track.jl")
include("lattice/EdwardsTengTwiss.jl")


export CTPS, cst, findindex, PolyMap, getindexmap, tadd, tminus, tmult, tdiv, tpow, tsqrt, tsin, tcos, ttan, tcosh, tsinh, reassign!
export DRIFT, KQUAD, KSEXT, KOCT, SBEND, RFCA, AbstractElement
export Beam
export EdwardsTengTwiss, AbstractTwiss, twissPropagate, findm66, periodicEdwardsTengTwiss, Twissline
export linepass!, pass!, ringpass!, linepass_TPSA!, pass_TPSA!, ringpass_TPSA!

end
