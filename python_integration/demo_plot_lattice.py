"""
Lattice Plotting Demo

This demonstrates how to visualize an accelerator lattice using pyJuTrack.
"""

import numpy as np
import pyJuTrack as jt

# Create lattice elements
D1 = jt.DRIFT("D1", 0.2)
Q1 = jt.QUAD("Q1", length=0.1, k1=29.6)
Q2 = jt.QUAD("Q2", length=0.3, k1=-29.6)
D2 = jt.DRIFT("D2", 0.4)
D3 = jt.DRIFT("D3", 0.2)
S1 = jt.KSEXT("S1", length=0.2, k2=0.0)
O1 = jt.KOCT("O1", length=0.1, k3=0.0)
B1 = jt.SBEND("B1", length=0.6, angle=np.pi/6)
B2 = jt.RBEND("B2", length=0.4, angle=-np.pi/6)
B3 = jt.RBEND("B3", length=0.2, angle=np.pi/6)
B4 = jt.SBEND("B4", length=0.2, angle=-np.pi/6)
RF = jt.RFCA("RF", length=0.2, volt=3.42*8.5e6, freq=591e6)

# Build beamline
line = [RF, D1, Q1, D2, Q2, D3, B1, D1, Q1, D2, Q2, D3,
        B2, D1, Q1, D2, Q2, D3, B3, D1, Q1, D2, Q2, D3,
        B4, D1, Q1, D2, Q2, D3, S1, D2, O1, D2, S1, D2]

lattice = jt.Lattice(line)

# Count element types
element_types = {}
for i in range(len(lattice)):
    elem = lattice[i]
    elem_type = type(elem).__name__
    element_types[elem_type] = element_types.get(elem_type, 0) + 1

for elem_type, count in sorted(element_types.items()):
    print(f"  {elem_type}: {count}")


jt.plot_lattice(line, 0.3, True, savepath="lattice_plot.png")

