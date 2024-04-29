module JuTrack
const CoordLimit = 1.0
const AngleLimit = 1.0
const m_e = 0.51099895e6
const m_p = 938.27208816e6
const m_goldion = 931.49410242e6 # charge 79, atomic number 197
const CGAMMA =	8.846056192e-05
use_exact_Hamiltonian = 0

include("TPSA_Enzyme/TPSA.jl")
include("lattice/beam.jl")
include("lattice/canonical_elements.jl")
include("tracking/bend.jl")
include("tracking/drift.jl")
include("tracking/multipole.jl")
include("tracking/rfcavity.jl")
include("tracking/thinmultipole.jl")
include("tracking/corrector.jl")
include("tracking/wakefield.jl")

include("tracking/bend_TPSA.jl")
include("tracking/drift_TPSA.jl")
include("tracking/multipole_TPSA.jl")
include("tracking/rfcavity_TPSA.jl")
include("tracking/solenoid.jl")
include("tracking/solenoid_TPSA.jl")
include("tracking/crabcavity_TPSA.jl")
include("tracking/thinmultipole_TPSA.jl")
include("tracking/corrector_TPSA.jl")
include("tracking/track.jl")
include("lattice/EdwardsTengTwiss.jl")

# multi-threading
include("tracking/track_mthread.jl")

include("tracking/lorentz.jl")
include("tracking/crabcavity.jl")
include("tracking/strongbb.jl")
include("lattice/optics.jl")
include("lattice/bunchedbeam.jl")

include("utils/lattice_utils.jl")
include("utils/matrix.jl")
export Beam
export m_e, m_p, m_goldion, CGAMMA, CoordLimit, AngleLimit, use_exact_Hamiltonian, use_exact_drift
export qr_eigen, diag1
export CRABCAVITY, easyCRABCAVITY, AccelCavity, LorentzBoost, InvLorentzBoost, StrongGaussianBeam, StrongThinGaussianBeam, AbstractStrongBeamBeam
export LongitudinalRFMap, AbstractLongitudinalRFMap, AbstractTransferMap, AbstractTransverseMap
export LongitudinalRLCWake, LongitudinalWake
export AbstractOptics, AbstractOptics2D, AbstractOptics4D, optics2D, optics4DUC
export initilize_6DGaussiandist!, get_emittance!, get_2nd_moment!, get_centroid!, histogram1DinZ!
export initilize_zslice!

export CTPS, cst, findindex, PolyMap, getindexmap, tadd, tminus, tmult, tdiv, tpow, tsqrt, tsin, tcos, ttan, tcosh, tsinh, reassign!
export AbstractElement, DRIFT, KQUAD, KSEXT, KOCT, SBEND, RBEND, RFCA, SOLENOID, MARKER, CORRECTOR, HKICKER, VKICKER, thinMULTIPOLE
export EdwardsTengTwiss, AbstractTwiss, twissPropagate, findm66, periodicEdwardsTengTwiss, Twissline, ADTwissline, twissring, ADfindm66, ADtwissring, ADperiodicEdwardsTengTwiss
export linepass!, pass!, ringpass!, linepass_TPSA!, pass_TPSA!, ringpass_TPSA!
export plinepass!, pringpass!, pass_P!, ADlinepass!, ADlinepass_TPSA!
export matrix_to_array, array_to_matrix
export total_length, spos, findelem, plot_optics
export ADfindm66_refpts
end
