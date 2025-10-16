module JuTrack
using PyCall
using Enzyme
# Enzyme.API.runtimeActivity!(true) # this is temporarily used
const CoordLimit = 1.0
const AngleLimit = 1.0
const m_e = 0.51099895069e6
const m_p = 938.27208816e6
const m_goldion = 931.49410242e6 # charge 79, atomic number 197
const CGAMMA =	8.846273768691263e-5
const __RE = 2.8179403205e-15 # classical electron radius
const RAD_CONST_E = 2.0/3.0*__RE
const RAD_CONST_P = RAD_CONST_E * m_e/m_p 
const epsilon_0 = 8.854187817e-12
const speed_of_light = 2.99792458e8 # m/s
const charge_e = 1.602176634e-19 # C
use_exact_Hamiltonian = 1 # use exact pz
use_exact_beti = 0 # use delta p/p0 as the sixth coordinate. Change it to 1 to use delta E/p0 
const _jlplotlib_available = let ok
    try
        pyimport("matplotlib.pyplot")
        ok = true
    catch
        ok = false
    end
    ok
end
# include("TPSA/TPSA.jl")
include("TPSA/fast_TPSA_module.jl")
include("TPSA/TPSA.jl")
using .HighOrderTPS
using .TPSAadStatic
export DTPSAD, NVAR, set_tps_dim, Gradient, Jacobian
export Number2TPSAD, TPSAD2Number, to_TPSAD, to_Number
include("utils/erfc.jl")
using .PureErrorFunctions
export erfcx, erf, erfinv
include("lattice/beam.jl")
include("lattice/canonical_elements.jl")
include("TPSA/fast_TPSA_tracking.jl")
include("tracking/bend.jl")
include("tracking/drift.jl")
include("tracking/multipole.jl")
include("tracking/rfcavity.jl")
include("tracking/thinmultipole.jl")
include("tracking/corrector.jl")
include("tracking/wakefield.jl")
include("tracking/quad.jl")
include("tracking/space_charge.jl")
include("tracking/refobt.jl")
include("tracking/bendlinear.jl")

include("tracking/solenoid.jl")
include("tracking/track_vector.jl")
include("lattice/EdwardsTengTwiss.jl")
include("lattice/ResonanceDrivingTerms.jl")

include("tracking/fringe.jl")

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
if _jlplotlib_available
    include("utils/lattice_plot.jl")
else
    @warn "Matplotlib is not available. Lattice plotting functions will not work."
end
# include("utils/lattice_plot.jl")

export Beam
export m_e, m_p, m_goldion, charge_e, speed_of_light, epsilon_0, CGAMMA, CoordLimit, AngleLimit, use_exact_Hamiltonian, use_exact_drift, use_exact_beti
export qr_eigen, diag1, randn_approx
# export Lattice, add!, buildlattice
export CRABCAVITY,CRABCAVITY_K2, easyCRABCAVITY, AccelCavity, LorentzBoost, InvLorentzBoost, StrongGaussianBeam, 
    StrongThinGaussianBeam, AbstractStrongBeamBeam, crab_crossing_setup!, pass_lumi!, pass_lumi_P!, Bassetti_Erskine!
export LongitudinalRFMap, AbstractLongitudinalRFMap, AbstractTransferMap, AbstractTransverseMap
export LongitudinalRLCWake, LongitudinalWake, wakefieldfunc_RLCWake
export AbstractOptics, AbstractOptics2D, AbstractOptics4D, optics2D, optics4DUC
export initilize_6DGaussiandist!, get_emittance!, get_2nd_moment!, get_centroid!, histogram1DinZ!
export initilize_zslice!, twiss_2d, twiss_beam, Gauss3_Dist

export CTPS, cst, findindex, PolyMap, getindexmap, reassign!, assign!
export AbstractElement, DRIFT, KQUAD, KSEXT, KOCT, SBEND, RBEND, RFCA, SOLENOID, MARKER, CORRECTOR, HKICKER, VKICKER, thinMULTIPOLE
export QUAD, LBEND, ESBEND, ERBEND, buildlatt
export TRANSLATION, YROTATION
export SPACECHARGE, QUAD_SC, DRIFT_SC, KQUAD_SC, KSEXT_SC, KOCT_SC, SBEND_SC, RBEND_SC, calculate_K
export EdwardsTengTwiss, AbstractTwiss, twissPropagate, findm66, periodicEdwardsTengTwiss, twissline, ADtwissline, twissring, ADfindm66, ADtwissring, ADperiodicEdwardsTengTwiss
export find_closed_orbit, fastfindm66, fastfindm66_refpts, ADfastfindm66_refpts, findm66_refpts
export linepass!, pass!, ringpass!, linepass_TPSA!, pass_TPSA!, ringpass_TPSA!, check_lost
export plinepass!, pringpass!, pass_P!, ADlinepass!, ADlinepass_TPSA!, ADringpass!, ADpringpass!, ADplinepass!
export matrix_to_array, array_to_matrix
export plot_lattice
export total_length, spos, findelem, insert_space_charge, array_optics, get_len, symplectic
export find_closed_orbit_6d, find_closed_orbit_4d, tracking_U0, integral_U0, rad_on!, rad_off!, fast_closed_orbit_4d, fast_closed_orbit_6d
export gettune, getchrom
# export CMscan, TPSVar6D, TPSVar4D, twiss_from_6x6, get_variables, evaluate, construct_sqr_matrix
export ADfindm66_refpts


function Duplicated(x::Float64, dx::Base.RefValue{Float64})
    return Duplicated(x, dx[])
end
export autodiff, Forward, ForwardWithPrimal, gradient, jacobian, Duplicated, set_runtime_activity, Const, Val, Enzyme, BatchDuplicated

export dynamic_aperture, computeRDT, ADcomputeRDT, FMA, plot_fma, fma_map_from_segments
export drift6!, multmv!, addvv!, linearQuadFringeElegantEntrance!, QuadFringePassP!, fastdrift!, strthinkick!
end
