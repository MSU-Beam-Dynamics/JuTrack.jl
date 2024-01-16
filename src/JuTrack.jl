module JuTrack
include("TPSA_Enzyme/TPSA_fixedmap.jl")
include("lattice/canonical_elements.jl")
# include("lattice/EdwardsTengTwiss.jl") # TPSA still not stable with Enzyme
include("tracking_Enzyme/track.jl")


export CTPS, cst, findindex, findpower, redegree, assign!,  pow, PolyMap, getindexmap # element, evaluate, derivative, integrate,
export EDRIFT, KQUAD, KSEXT, KOCT, CSBEND, RFCA, AbstractElement
# export EdwardsTengTwiss, AbstractTwiss, twissPropagate, findm66, periodicEdwardsTengTwiss
export linepass!, linepass_TPSA

end
