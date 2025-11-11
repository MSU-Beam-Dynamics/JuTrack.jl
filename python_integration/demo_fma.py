"""
Frequency Map Analysis (FMA) Demo

This script demonstrates how to perform Frequency Map Analysis using pyJuTrack.
"""

import numpy as np
import pyJuTrack as jt


def create_test_ring(rad=0):
    """Create a test storage ring for FMA analysis"""
    
    # Markers
    markers = [jt.MARKER(f"CENOFSTR{i:02d}") for i in range(1, 13)]
    
    # Drifts
    D11 = jt.DRIFT("D11", 0.5)
    D11A = jt.DRIFT("D11A", 0.535)
    DX = jt.DRIFT("DX", 0.0375)
    D12 = jt.DRIFT("D12", 0.075)
    D15 = jt.DRIFT("D15", 0.225)
    
    # Quadrupoles
    QF1 = jt.KQUAD("QF1", 0.18, 13.84045189)
    QD1 = jt.KQUAD("QD1", 0.14, -13.77444068)
    QF2 = jt.KQUAD("QF2", 0.19, 10.19367)
    QF3 = jt.KQUAD("QF3", 0.115, 10.94517)
    QF4 = jt.KQUAD("QF4", 0.305, 15.32064)
    QF5 = jt.KQUAD("QF5", 0.305, 15.8)
    QF6 = jt.KQUAD("QF6", 0.305, 15.68564)
    
    # Sextupoles
    SHH = jt.KSEXT("SHH", 0.075, 3.514648)
    SHH2 = jt.KSEXT("SHH2", 0.075, -929.5772)
    SD = jt.KSEXT("SD", 0.28, -1367.68410852)
    SF = jt.KSEXT("SF", 0.28, 1610.76982434)

    # Bends
    BEND1 = jt.SBEND("BEND1", 0.34, 0.05817765, PolynomB=[0.0,-2.827967,0.0,0.0])
    BEND2 = jt.SBEND("BEND2", 0.5, 0.05817765, PolynomB=[0.0,-7.057813,0.0,0.0])
    BEND3 = jt.SBEND("BEND3", 0.5, 0.05817765, PolynomB=[0.0,-7.057813,0.0,0.0])
    
    # Straight sections
    STR_A = [D11A, D11, D11, D11, D11]
    STR_B = [D11, D11, D11, D11, D11A]
    
    # Arc section
    ARC = [
        DX, SHH, D12, QF1, DX, DX, QD1, DX, DX, SHH2, D12, BEND1,
        DX, DX, SD, D12, QF2, D12, SF, DX, DX, QF3, D15, BEND2,
        DX, DX, QF4, DX, DX, BEND3, DX, DX, QF5, D12, BEND3,
        DX, DX, QF6, DX, DX, BEND3, DX, DX, QF6, DX, DX, BEND3,
        DX, DX, QF5, D12, BEND3, D12, QF4, DX, DX, BEND2,
        D15, QF3, D12, SF, DX, DX, QF2, D12, SD, DX, DX, BEND1,
        DX, DX, SHH2, DX, DX, QD1, DX, DX, QF1, D12, SHH, DX
    ]
    
    # Build complete ring (12 superperiods)
    ring_elements = []
    for i, marker in enumerate(markers):
        ring_elements.append(marker)
        ring_elements.extend(STR_B)
        ring_elements.extend(ARC)
        ring_elements.extend(STR_A)
    
    return jt.Lattice(ring_elements)

RING = create_test_ring()
print(f"  Total length: {RING.total_length():.3f} m")

# Calculate tunes
try:
    nux, nuy = jt.gettune(RING, dp=0.0)
    print(f"\nBetatron tunes:")
    print(f"  νx = {nux:.6f}")
    print(f"  νy = {nuy:.6f}")
except Exception as e:
    print(f"\nCould not calculate tunes: {e}")

# Perform FMA
print("\nPerforming Frequency Map Analysis...")
try:
    num_turns = 256    
    rows = jt.FMA(RING, num_turns)

    
    # Use JuTrack's plot_fma function to visualize results
    try:
        jt.plot_fma(rows, 
                    figsize=(10, 4),
                    s=10,
                    x_min=-3e-3, x_max=3e-3,
                    y_min=0.0, y_max=3e-3,
                    resonance_lines=True,
                    resonance_orders=[1, 2, 3, 4],
                    filepath="fma_result.png")
        print("Saved FMA plot to 'fma_result.png'")
    except Exception as e:
        print(f"\n Could not create plot: {e}")
        import traceback
        traceback.print_exc()

except Exception as e:
    print(f"\nFMA failed: {e}")
    import traceback
    traceback.print_exc()
