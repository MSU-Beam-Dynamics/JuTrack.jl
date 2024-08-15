module JuTrack
using Enzyme
Enzyme.API.runtimeActivity!(true) # this is temporarily used
const CoordLimit = 1.0
const AngleLimit = 1.0
const m_e = 0.51099895e6
const m_p = 938.27208816e6
const m_goldion = 931.49410242e6 # charge 79, atomic number 197
const CGAMMA =	8.846056192e-05
const epsilon_0 = 8.854187817e-12
const speed_of_light = 2.99792458e8 # m/s
const charge_e = 1.602176634e-19 # C
use_exact_Hamiltonian = 0

include("TPSA/TPSA.jl")
include("lattice/beam.jl")
include("lattice/canonical_elements.jl")
include("tracking/bend.jl")
include("tracking/drift.jl")
include("tracking/multipole.jl")
include("tracking/rfcavity.jl")
include("tracking/thinmultipole.jl")
include("tracking/corrector.jl")
include("tracking/wakefield.jl")
include("tracking/quad.jl")
include("tracking/space_charge.jl")

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
include("lattice/ResonanceDrivingTerms.jl")

include("tracking/fringe.jl")
include("tracking/fringe_TPSA.jl")

# multi-threading
include("tracking/track_mthread.jl")

include("tracking/lorentz.jl")
include("tracking/crabcavity.jl")
include("tracking/strongbb.jl")
include("lattice/bunchedbeam.jl")

include("utils/lattice_utils.jl")
include("utils/matrix.jl")
include("utils/dynamic_aperture.jl")
include("utils/fma.jl")

export Beam
export m_e, m_p, m_goldion, speed_of_light, epsilon_0, CGAMMA, CoordLimit, AngleLimit, use_exact_Hamiltonian, use_exact_drift
export qr_eigen, diag1
export CRABCAVITY, easyCRABCAVITY, AccelCavity, LorentzBoost, InvLorentzBoost, StrongGaussianBeam, 
    StrongThinGaussianBeam, AbstractStrongBeamBeam, crab_crossing_setup!, pass_lumi!, pass_lumi_P!, Bassetti_Erskine!
export LongitudinalRFMap, AbstractLongitudinalRFMap, AbstractTransferMap, AbstractTransverseMap
export LongitudinalRLCWake, LongitudinalWake, wakefieldfunc_RLCWake
export AbstractOptics, AbstractOptics2D, AbstractOptics4D, optics2D, optics4DUC
export initilize_6DGaussiandist!, get_emittance!, get_2nd_moment!, get_centroid!, histogram1DinZ!
export initilize_zslice!, twiss_2d, twiss_beam

export CTPS, cst, findindex, PolyMap, getindexmap, reassign!
export AbstractElement, DRIFT, KQUAD, KSEXT, KOCT, SBEND, RBEND, RFCA, SOLENOID, MARKER, CORRECTOR, HKICKER, VKICKER, thinMULTIPOLE
export QUAD, buildlatt
export SPACECHARGE
export EdwardsTengTwiss, AbstractTwiss, twissPropagate, findm66, periodicEdwardsTengTwiss, twissline, ADtwissline, twissring, ADfindm66, ADtwissring, ADperiodicEdwardsTengTwiss
export fastfindm66, fastfindm66_refpts, ADfastfindm66_refpts
export linepass!, pass!, ringpass!, linepass_TPSA!, pass_TPSA!, ringpass_TPSA!, check_lost
export plinepass!, pringpass!, pass_P!, ADlinepass!, ADlinepass_TPSA!
export matrix_to_array, array_to_matrix
export total_length, spos, findelem, insert_space_charge, array_optics, get_len
export ADfindm66_refpts

function Duplicated(x::Float64, dx::Base.RefValue{Float64})
    return Duplicated(x, dx[])
end
export autodiff, Forward, gradient, jacobian, Duplicated, DuplicatedNoNeed, Const, Val, Enzyme, BatchDuplicated

export dynamic_aperture, naff, FMA, computeRDT, ADcomputeRDT
end
