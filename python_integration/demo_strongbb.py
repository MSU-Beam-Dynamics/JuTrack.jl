"""
Strong Beam-Beam Interaction Demo

This demonstrates strong beam-beam interaction simulation with automatic differentiation.
We calculate the final emittance of the electron beam after interaction with a proton beam,
and compute derivatives w.r.t. the proton beam emittance at the interaction point.

The related Python functions are not fully developed yet. Here we use jl = jt.get_julia_module() 
to call the necessary Julia functions directly. 
This method must be used with caution and will be replaced by pure Python APIs in the future.
"""

import numpy as np
import matplotlib.pyplot as plt
import pyJuTrack as jt

# Setup for 3-var optimization
jt.set_tps_dim(3)

def simulate_beambeam(emitx_p, emity_p, sigz_p):
    vbase = 3.42 * 8.5e6
    phis = 10.0  
    vact = vbase / np.cos(phis * np.pi / 180.0)
    
    # Call Julia functions directly
    jl = jt.get_julia_module()
    mainRFe = jl.AccelCavity(freq=jt.DTPSAD(591e6), volt=vact, h=7560.0, 
                             phis=np.pi - phis * np.pi / 180.0)
    
    tunex, tuney = 50.08, 44.14
    alphac = 3.42 / tunex / tunex
    
    lmap = jl.LongitudinalRFMap(jt.DTPSAD(alphac), mainRFe)
    
    # Optics at IP
    opIPp = jl.optics4DUC(0.8, 0.0, 0.072, 0.0)  
    opIPe = jl.optics4DUC(0.45, 0.0, 0.056, 0.0) 
    
    pstrong = jl.StrongGaussianBeam(
        jt.DTPSAD(1.0),              # charge
        jt.DTPSAD(jt.m_p),           # mass (proton)
        jt.DTPSAD(1.0),              # beta
        int(0.688e11),               # number of particles
        jt.DTPSAD(275e9),            # energy (275 GeV)
        opIPp,                       # optics
        jt._to_julia_dtpsad_vector([emitx_p, emity_p, sigz_p]),  # emittances
        9                            # number of slices
    )
    
    jl.initilize_zslice_b(pstrong, jl.Symbol("gaussian"), jl.Symbol("evennpar"), 7.0)
    
    crab_ratio = 0.33
    overcrab = 1.0
    pcrab1st = jt.easyCRABCAVITY("cc1", freq=197.0e6, 
                                  halfthetac=overcrab * 12.5e-3 * (1 + crab_ratio))
    pcrab2nd = jt.easyCRABCAVITY("cc2", freq=197.0e6 * 2.0,
                                  halfthetac=-overcrab * 12.5e-3 * crab_ratio)
    
    jl.crab_crossing_setup_b(pstrong, 12.5e-3, pcrab1st, pcrab2nd)
    
    particles_dtpsad = jt.zeros(jt.DTPSAD, 5000, 6)
    ebeam = jt.Beam(
        particles_dtpsad,
        np=int(1.72e11 * 3),
        energy=10e9,
        emittance=[20e-9, 1.3e-9, 1.36e-4]
    )
    
    jl.initilize_6DGaussiandist_b(ebeam.julia_object, opIPe, lmap)
    
    jt.linepass([pstrong], ebeam)
    jl.get_emittance_b(ebeam.julia_object)
    
    emit = ebeam.julia_object.emittance
    return emit[0], emit[1], emit[2]

emitx_p0 = 95e-6  
emity_p0 = 8.5e-6  
sigz_p0 = 0.06     

x1 = jt.DTPSAD(emitx_p0, 1)
x2 = jt.DTPSAD(emity_p0, 2)
x3 = jt.DTPSAD(sigz_p0, 3)

emitx_e, emity_e, emitz_e = simulate_beambeam(x1, x2, x3)

print("\nElectron beam emittances after interaction:")
print(f"  Horizontal: {emitx_e} m")
print(f"  Vertical:   {emity_e} m")
print(f"  Longitudinal: {emitz_e} m")

# Compute Jacobian: d(emit_e) / d(emit_p)
jac = jt.Jacobian(simulate_beambeam, [emitx_p0, emity_p0, sigz_p0])

print("\nJacobian matrix (d[emitx_e, emity_e, emitz_e] / d[emitx_p, emity_p, sigz_p]):")
print(jac)

print("\nInterpretation:")
print(f"  ∂emitx_e/∂emitx_p = {jac[0,0]:.6e}")
print(f"  ∂emity_e/∂emity_p = {jac[1,1]:.6e}")
print(f"  ∂emitz_e/∂sigz_p  = {jac[2,2]:.6e}")

# Visualize Jacobian
fig, ax = plt.subplots(figsize=(10, 8))

im = ax.imshow(jac, cmap='RdBu_r', aspect='auto')
ax.set_xticks([0, 1, 2])
ax.set_yticks([0, 1, 2])
ax.set_xticklabels(['εx(p)', 'εy(p)', 'σz(p)'])
ax.set_yticklabels(['εx(e)', 'εy(e)', 'εz(e)'])
ax.set_xlabel('Proton beam parameters', fontsize=12)
ax.set_ylabel('Electron beam emittances', fontsize=12)

# Add colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Sensitivity', fontsize=10)

# Add text annotations
for i in range(3):
    for j in range(3):
        text = ax.text(j, i, f'{jac[i, j]:.2e}',
                      ha="center", va="center", color="black", fontsize=10)

plt.tight_layout()
plt.savefig('strongbb_jacobian.png', dpi=150, bbox_inches='tight')
print("\nSaved Jacobian plot to 'strongbb_jacobian.png'")

plt.show()
