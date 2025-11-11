"""
Wakefield Demo

This demonstrates wakefield simulation with automatic differentiation.
We compute the effect of longitudinal RLC wakefield on beam emittance
and calculate derivatives w.r.t. the wakefield frequency.
"""

import numpy as np
import matplotlib.pyplot as plt
import pyJuTrack as jt

jt.set_tps_dim(1)

def simulate_wakefield(freq):
    # Handle both direct DTPSAD input and list input from Jacobian
    if isinstance(freq, list):
        freq_val = freq[0]
    else:
        freq_val = freq
    
    # Convert to DTPSAD if needed (Jacobian will pass DTPSAD values)
    if not hasattr(freq_val, '_jl_callmethod'):
        freq_dtpsad = jt.convert_to_DTPSAD(freq_val)
    else:
        freq_dtpsad = freq_val
    
    # Create RLC wakefield element
    RLCwake = jt.LongitudinalRLCWake(
        "RLCWake",
        freq=freq_dtpsad,
        Rshunt=jt.DTPSAD(5.5e3),
        Q0=jt.DTPSAD(3.0)
    )
    
    D1 = jt.DRIFT("D1", jt.DTPSAD(0.1))
    D2 = jt.DRIFT("D2", jt.DTPSAD(0.1))
    
    np.random.seed(42)  
    initial_coords = jt.to_dtpsad_matrix(np.random.randn(500, 6) * 1e-4)
    ebeam = jt.Beam(initial_coords, energy=10e9)
    
    jt.histogram1DinZ(ebeam)
    
    line = [D1, RLCwake]
    
    jt.linepass(line, ebeam)
    emitx, emity, emitz = jt.get_emittance(ebeam)
    return emitx, emity, emitz


freq0 = 180e9 
freq0_dtpsad = jt.DTPSAD(freq0, 1)
emitx, emity, emitz = simulate_wakefield([freq0_dtpsad])


jac, results = jt.Jacobian(simulate_wakefield, [freq0], return_value=True)

print(f"  Emittances: {results}")

print("\nJacobian (d[emitx, emity, emitz] / d[frequency]):")
print(jac)


