"""
This demonstrates RDT optimization using Python version of JuTrack.jl.
Objective: Minimize sum of RDT terms (h21000, h10110, h30000, h10200, h10020)
by tuning sextupole strengths (SDM and SFM) in the SPEAR3 lattice.
"""

import numpy as np
import matplotlib.pyplot as plt
import pyJuTrack as jt

jt.set_tps_dim(2)

# Load SPEAR3 lattice
jt.include("src/demo/SPEAR3/spear3.jl")
ring = jt.Lattice(jt.call_julia("spear3"))

sdm_idx = jt.findelem(ring, name="SDM")
sfm_idx = jt.findelem(ring, name="SFM")
print(f"Found {len(sdm_idx)} SDM and {len(sfm_idx)} SFM sextupoles")

idx_marker = jt.findelem(ring, element_type=jt.element_types.MARKER)
print(f"Found {len(idx_marker)} MARKER elements")

ring_tpsa = jt.Number2TPSAD(ring)


def objective(x):
    """
    Pure Python objective function.
    
    Args:
        x: DTPSAD Vector [sdm_k2, sfm_k2]
    
    Returns:
        Sum of RDT terms (h21000, h10110, h30000, h10200, h10020)
    """
    x1 = x[0]  
    x2 = x[1]  
    
    for i in sdm_idx:
        elem = ring_tpsa[i-1]
        elem.k2 = x1
        elem.PolynomB[2] = x1
    
    for i in sfm_idx:
        elem = ring_tpsa[i-1]
        elem.k2 = x2
        elem.PolynomB[2] = x2
    
    dlist, tune = jt.computeRDT(ring_tpsa, idx_marker, E0=3e9, m0=jt.m_e)
    
    tot_21000 = jt.DTPSAD(0.0)
    tot_10110 = jt.DTPSAD(0.0)
    tot_30000 = jt.DTPSAD(0.0)
    tot_10200 = jt.DTPSAD(0.0)
    tot_10020 = jt.DTPSAD(0.0)
    
    for i in range(len(dlist)):
        tot_21000 += dlist[i].h21000[0] 
        tot_10110 += dlist[i].h10110[0]
        tot_30000 += dlist[i].h30000[0]
        tot_10200 += dlist[i].h10200[0]
        tot_10020 += dlist[i].h10020[0]
    
    return tot_21000 + tot_10110 + tot_30000 + tot_10200 + tot_10020

# Initial values
x = np.array([-17.0, 15.0])
print(f"\nInitial strengths: SDM k2 = {x[0]:.2f}, SFM k2 = {x[1]:.2f}")
print("\nStarting optimization...")

learning_rate = 0.001
num_iterations = 1000

x_history = [x.copy()]
obj_history = []
grad_history = []

for i in range(num_iterations):
    grad, obj_val = jt.Gradient(objective, x, return_value=True)
    
    obj_history.append(obj_val)
    grad_history.append(grad.copy())
    
    x = x - learning_rate * grad
    x_history.append(x.copy())
    
    if i % 10 == 0 or i < 3:
        print(f"Iter {i:3d}: obj = {obj_val:.6e}, x = [{x[0]:+7.2f}, {x[1]:+7.2f}], "
              f"grad = [{grad[0]:+.2e}, {grad[1]:+.2e}]")

print(f"  Initial objective: {obj_history[0]:.6e}")
print(f"  Final objective:   {obj_history[-1]:.6e}")
print(f"  Reduction:         {(1 - obj_history[-1]/obj_history[0])*100:.1f}%")
print(f"  Final SDM k2:      {x[0]:.2f}")
print(f"  Final SFM k2:      {x[1]:.2f}")


fig, axes = plt.subplots(1, 3, figsize=(13, 4))
ax = axes[0]
ax.semilogy(obj_history, 'b-', linewidth=2)
ax.set_xlabel('Iteration')
ax.set_ylabel('Objective (sum of RDTs)')
ax.grid(True, alpha=0.3)

ax = axes[1]
x_arr = np.array(x_history)
ax.plot(x_arr[:, 0], x_arr[:, 1], 'g-o', markersize=3, alpha=0.6)
ax.plot(x_arr[0, 0], x_arr[0, 1], 'ro', markersize=5, label='Start')
ax.plot(x_arr[-1, 0], x_arr[-1, 1], 'bs', markersize=5, label='End')
ax.set_xlabel(r'SDM $k_2$ [$m^{-3}$]')
ax.set_ylabel(r'SFM $k_2$ [$m^{-3}$]')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[2]
grad_arr = np.array(grad_history)
ax.plot(grad_arr[:, 0], 'b-', label='∂f/∂(SDM k2)', linewidth=2)
ax.plot(grad_arr[:, 1], 'r-', label='∂f/∂(SFM k2)', linewidth=2)
ax.set_xlabel('Iteration')
ax.set_ylabel('Gradient')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('rdt_optimization.png', dpi=150, bbox_inches='tight')
print("Saved plot to 'rdt_optimization.png'")

plt.show()