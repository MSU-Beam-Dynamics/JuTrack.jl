module JuTrack
# include("TPSA_Enzyme/arrayTPSA_fixedmap.jl")
include("lattice/canonical_elements.jl")
# include("lattice/EdwardsTengTwiss.jl") # TPSA still not stable with Enzyme
include("tracking_Enzyme/EDrift_Enzyme.jl")
include("tracking_Enzyme/multipole_Enzyme.jl")
# include("tracking_Enzyme/multipole_TPSA_Enzyme.jl")
include("tracking_Enzyme/csbend_Enzyme.jl")
# include("tracking_Enzyme/csbend_TPSA_Enzyme.jl")
include("tracking_Enzyme/rfca_Enzyme.jl")
# include("tracking_Enzyme/rfca_TPSA_Enzyme.jl")
include("tracking_Enzyme/track.jl")


export CTPS, cst, findindex, PolyMap, getindexmap, tadd, tminus, tmult, tdiv, tpow, tsqrt, tsin, tcos, ttan, tcosh, tsinh
export EDRIFT, KQUAD, KSEXT, KOCT, CSBEND, RFCA, AbstractElement
# export EdwardsTengTwiss, AbstractTwiss, twissPropagate, findm66, periodicEdwardsTengTwiss
export linepass!, pass!, ringpass!#, linepass_TPSA

end
