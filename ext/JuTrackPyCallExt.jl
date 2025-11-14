module JuTrackPyCallExt

using JuTrack
import JuTrack: _pycall_forced_off, _register_pycall_features!

using PyCall

if _pycall_forced_off
    @info "PYJUTRACK_DISABLE_PYCALL is set; skipping PyCall-backed JuTrack features."
else
    JuTrack._register_pycall_features!(pyimport)
end

end # module
