"""
This script demonstrates using pyJuTrack to tune the ESR lattice for a specific 
phase advance between two crab cavities using auto-differentiation.

Warning: This code is for demonstration purposes. Simply tuning the quadrupoles 
may result in unstable solutions.
"""

import numpy as np
import pyJuTrack as jt
import sys

E0 = 17.846262619763e9  # Beam energy in eV
jt.set_tps_dim(7) # 7 variables for the 7 quadrupoles we'll tune

# Load the ESR main lattice
# ESR_crab = jt.load_lattice("src/demo/ESR/esr_main.jls")
from juliacall import Main as jl
# RING is defined in esr_main.jl
jl.seval('include("src/demo/ESR/esr-main.jl")')
ESR_crab = jt.Lattice(jl.RING)
print(f"Loaded ESR lattice with {len(ESR_crab)} elements")
print(f"Total length: {ESR_crab.total_length():.3f} m")

ESR_crab = jt.Number2TPSAD(ESR_crab)

# Function to perturb quadrupoles to create some errors
def Q_perturb(lattice):
    """Add random perturbations to quadrupole strengths"""
    for i in range(len(lattice)):
        elem = lattice[i]
        elem_type = type(elem).__name__
        if 'KQUAD' in str(type(elem)):
            k1_val = elem.k1
            if hasattr(k1_val, 'val'):
                k1 = float(k1_val.val)
            else:
                k1 = float(k1_val)
            
            k1_perturbed = k1 * (1 + 0.001 * np.random.randn())
            
            if hasattr(k1_val, 'val'):
                new_KQ = jt.KQUAD(name=str(elem.name), length=elem.len, k1=jt.DTPSAD(k1_perturbed))
            else:
                new_KQ = jt.KQUAD(name=str(elem.name), length=float(elem.len), k1=k1_perturbed)
            lattice[i] = new_KQ
    return lattice

ESR_perturb = Q_perturb(ESR_crab)

# Find crab cavity indices
crab_idx = jt.findelem(ESR_crab, element_type=jt.element_types.CRABCAVITY)

# Indices of quadrupoles to be tuned.
# !!! Note: These are 1-based indices in Julia. Subtract 1 for 0-based Python indexing.
changed_idx = [9, 13, 17, 23, 27, 31, 5537]
print(f"\nQuadrupoles to be tuned: {changed_idx}")

# objective function
def get_phase14_zero(X):
    working_lattice = jt.Lattice(ESR_perturb.julia_object)
    
    for i, idx in enumerate(changed_idx):
        working_lattice[idx-1].k1 = X[i]
    
    refpts = list(range(1, len(working_lattice) + 1))
    twi = jt.twissring(working_lattice, dp=0.0, refpts=refpts, E0=E0, m0=jt.m_e)
    
    # Phase advance between crab cavities
    phase41 = twi[35].mux + twi[-1].mux - twi[5533].mux
    
    # Subtract 2*pi
    return phase41 - jt.DTPSAD(2.0 * np.pi)

def multi_val_optimization(x0, niter, step):
    """
    Perform gradient descent optimization to minimize phase advance error.
    
    Args:
        x0: Initial quadrupole strengths (list of 7 values)
        niter: Maximum number of iterations
        step: Step size for gradient descent
    
    Returns:
        x0_vals: History of quadrupole strengths
        goal_vals: History of objective function values
        grad_vals: History of gradients
    """
    target = 0.015  
    x0_vals = np.zeros((7, niter))
    goal_vals = []
    grad_vals = np.zeros((7, niter))
    
    grads_init, g0 = jt.Gradient(get_phase14_zero, x0, return_value=True)
    print(f"\nInitial phase advance error: {g0:.6f} rad")
    print(f"Target: {target:.6f} rad")

    for i in range(niter):
        grads, new_phase = jt.Gradient(get_phase14_zero, x0, return_value=True)
        
        for j in range(7):
            x0[j] -= step * grads[j]
        
        x0_vals[:, i] = x0
        goal_vals.append(new_phase)
        grad_vals[:, i] = grads
        
        if (i + 1) % 10 == 0 or i < 5:  # Print first 5 and every 10th iteration
            print(f"Step {i+1:3d}: Phase error = {new_phase:+.6f} rad")
        
        if abs(new_phase) < target:
            x0_vals = x0_vals[:, :i+1]
            grad_vals = grad_vals[:, :i+1]
            break
    else:
        print(f"Maximum iterations ({niter}) reached without convergence")
        print(f"Final phase error: {new_phase:.6f} rad (target: {target:.6f} rad)")
    
    return x0_vals, goal_vals, grad_vals

# Initial guess for quadrupole strengths
xinit = [-1e-6, 1e-6, -1e-6, 1e-6, -1e-6, 1e-6, -1e-6]
print(f"\nInitial quadrupole strengths: {xinit}")

# Run optimization
x0_vals, goal_vals, grad_vals = multi_val_optimization(xinit.copy(), niter=100, step=1e-4)

# Print final results
print("\nFinal quadrupole strengths:")
for i, k in enumerate(x0_vals[:, -1]):
    print(f"  k{i+1} (index {changed_idx[i]}): {k:+.6e} m^-2")
print(f"Improvement: {(goal_vals[0] - goal_vals[-1]):.6f} rad")

try:
    import matplotlib.pyplot as plt
    
    plot_steps = len(goal_vals)
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10))
    
    # Plot 1: Evolution of quadrupole strengths
    for i in range(7):
        ax1.plot(range(1, plot_steps+1), x0_vals[i, :plot_steps], 
                marker='o', linestyle='--', label=f'k{i+1}', markersize=4)
    ax1.set_xlabel('Iterations')
    ax1.set_ylabel('Strength (m$^{-2}$)')
    ax1.set_title('Evolution of Quadrupole Strengths')
    ax1.legend(ncol=2)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Evolution of phase advance error
    ax2.plot(range(1, plot_steps+1), goal_vals, 
            marker='o', linestyle='--', color='blue', linewidth=2, markersize=4)
    ax2.axhline(y=0.015, color='r', linestyle=':', linewidth=2, label='Target')
    ax2.set_xlabel('Iterations')
    ax2.set_ylabel('Phase advance error (rad)')
    ax2.set_title('Evolution of Objective Function')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Evolution of gradients
    for i in range(7):
        ax3.plot(range(1, plot_steps+1), grad_vals[i, :plot_steps], 
                marker='o', linestyle='--', label=f'∂Δφ/∂k{i+1}', markersize=4)
    ax3.set_xlabel('Iterations')
    ax3.set_ylabel('Gradient')
    ax3.set_title('Evolution of Gradients')
    ax3.legend(ncol=2)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig('phase_advance_tuning_results.png', dpi=150, bbox_inches='tight')
    print(f"Saved plot to 'phase_advance_tuning_results.png'")
    
    # plt.show()
    
except ImportError:
    print("\nmatplotlib not available - skipping plots")
except Exception as e:
    print(f"\nError creating plots: {e}")


