# pyJuTrack - Python Interface to JuTrack.jl

A Python wrapper for JuTrack.jl.

## Installation

```bash
# Clone JuTrack.jl repository
git clone https://github.com/MSU-Beam-Dynamics/JuTrack.jl.git
cd JuTrack.jl/python_integration

# Install as Python package
pip install -e .

# Verify installation
python test_installation.py
```


## Quick Start

```python
import numpy as np
import pyJuTrack as jt

# Create lattice elements
d1 = jt.DRIFT("D1", 1.0)
qf = jt.KQUAD("QF", 0.5, k1=1.2)
d2 = jt.DRIFT("D2", 2.0)
qd = jt.KQUAD("QD", 0.5, k1=-1.2)

# Build lattice (accepts Python lists!)
fodo = jt.Lattice([d1, qf, d2, qd, d1])
print(f"Total length: {fodo.total_length():.2f} m")

# Create beam (accepts NumPy arrays!)
particles = np.random.randn(100, 6) * 1e-4
beam = jt.Beam(particles, energy=1e9)

# Track particles
jt.linepass(fodo, beam)

# Calculate Twiss parameters
twiss = jt.twissring(fodo, dp=0.0, E0=1e9, m0=jt.m_e)

# Get tunes
nux, nuy = jt.gettune(fodo, dp=0.0, E0=1e9, m0=jt.m_e)
print(f"Tunes: νx={nux:.4f}, νy={nuy:.4f}")
```


### Lattice Elements

**Basic:**
- `DRIFT(name, length)`
- `MARKER(name)`

**Magnets:**
- `QUAD(name, length, k1)` - Thick lens quadrupole
- `KQUAD(name, length, k1)` - Kick-drift-kick quadrupole
- `KSEXT(name, length, k2)` - Sextupole
- `KOCT(name, length, k3)` - Octupole
- `SBEND(name, length, angle)` - Sector bend
- `RBEND(name, length, angle)` - Rectangular bend
- `LBEND(name, length, angle)` - Linear bend
- `ESBEND(name, length, angle)` - Exact sector bend
- `ERBEND(name, length, angle)` - Exact rectangular bend
- `SOLENOID(name, length, ks)` - Solenoid
- `thinMULTIPOLE(name, ...)` - Thin multipole
- `CORRECTOR(name, length)`
- `HKICKER(name, length)` - Horizontal kicker
- `VKICKER(name, length)` - Vertical kicker

**RF and Special:**
- `RFCA(name, length, voltage, frequency)` - RF cavity
- `CRABCAVITY(name, ...)` - Crab cavity
- `CRABCAVITY_K2(name, ...)` - Crab cavity (K2 model)

**Space Charge Elements:**
- `SPACECHARGE`
- `SPACECHARGE2P5D`
- `insert_space_charge_2p5d(lattice, sc_path_length_min, ...)`
- `QUAD_SC`, `QUAD_SC2P5D`
- `DRIFT_SC`, `DRIFT_SC2P5D`
- `KQUAD_SC`, `KQUAD_SC2P5D`
- `KSEXT_SC`, `KSEXT_SC2P5D`
- `KOCT_SC`, `KOCT_SC2P5D`
- `SBEND_SC`, `SBEND_SC2P5D`
- `RBEND_SC`, `RBEND_SC2P5D`

**Transformations:**
- `TRANSLATION(name, dx, dy, dz)`
- `YROTATION(name, angle)`

### Tracking Functions

```python
# Single-pass tracking
jt.linepass(lattice, beam)

# Multi-turn tracking
jt.ringpass(lattice, beam, num_turns=1000)

# Parallel tracking (Using Julia's multi-threaded)
jt.plinepass(lattice, beam, nthreads=8)
jt.pringpass(lattice, beam, num_turns=1000, nthreads=8)
```

### 2.5-D Space Charge

JuTrack's 2.5-D model follows the coordinate convention
`[x, px, y, py, z, dp/p0]`. In Python you can either create
`SPACECHARGE2P5D` nodes explicitly, or insert them automatically into an
existing numeric lattice.

```python
import numpy as np
import pyJuTrack as jt

d1 = jt.DRIFT("D1", 0.2)
q1 = jt.KQUAD("Q1", 0.1, k1=29.6)
d2 = jt.DRIFT("D2", 0.4)
q2 = jt.KQUAD("Q2", 0.1, k1=-29.6)
d3 = jt.DRIFT("D3", 0.2)

cell = jt.Lattice([d1, q1, d2, q2, d3])
cell_sc = jt.insert_space_charge_2p5d(
    cell,
    0.05,
    periodic=True,
    xsize=32,
    ysize=32,
    zsize=32,
    pipe_radius=13e-3,
    xy_ratio=1.0,
    long_avg_n=3,
)

particles = np.zeros((1000, 6))
particles[:, 0] = np.random.randn(1000) * 2e-4
particles[:, 2] = np.random.randn(1000) * 2e-4
particles[:, 5] = np.random.randn(1000) * 1e-5
particles[:, 4] = np.random.randn(1000) * 3e-2

beam = jt.Beam(
    particles,
    energy=1.0e9,
    np=2_000_000_000,
    current=20.0,
    mass=jt.m_p,
    charge=1.0,
)

jt.linepass(cell_sc, beam)
```

Notes:
- `effective_length` is the path length represented by each thin kick.
- For thick native 2.5-D elements, use `DRIFT_SC2P5D`, `QUAD_SC2P5D`, `KQUAD_SC2P5D`, `KSEXT_SC2P5D`, `KOCT_SC2P5D`, `SBEND_SC2P5D`, or `RBEND_SC2P5D`.
- `insert_space_charge_2p5d(...)` currently inserts nodes into numeric lattices. If you need DTPSAD/AD, insert the nodes first and then call `jt.Number2TPSAD(...)`.
- The 2.5-D model is inactive if the bunch has zero longitudinal extent, so `z` must not be identical for every particle.

### Twiss Parameters and Optics

```python
# Calculate Twiss parameters for a ring
twiss = jt.twissring(ring, dp=0.0, refpts=None, E0=1e9, m0=jt.m_e)

# Access Twiss parameters
for tw in twiss:
    print(tw.betax, tw.betay, tw.mux, tw.muy, tw.etax)

# Calculate for beamline
twiss = jt.twissline(ring, dp=0.0, refpts=[0, 10, 20], E0=1e9, m0=jt.m_e)

# AD-compatible versions (for optimization)
twiss = jt.ADtwissring(ring, dp=0.0, changed_idx=[1,2], changed_ele=[q1, q2])

# Periodic Twiss
periodic_twiss = jt.periodicEdwardsTengTwiss(ring, E0=1e9, m0=jt.m_e)
```

### Transfer Matrices

```python
# Calculate 6x6 one-turn matrix
m66 = jt.findm66(ring, dp=0.0, E0=1e9)

# Fast calculation
m66 = jt.fastfindm66(ring, dp=0.0, E0=1e9)

# Matrices at specific points
matrices = jt.findm66_refpts(ring, refpts=[0, 10, 20], dp=0.0)
```

### Closed Orbit

```python
# Find closed orbit (6D)
orbit = jt.find_closed_orbit(ring, dp=0.0, E0=1e9, m0=jt.m_e)
```

### Lattice Utilities

```python
# Total length
length = jt.total_length(ring)

# S-positions
spos_all = jt.spos(ring)
spos_sel = jt.spos(ring, indices=[0, 5, 10])

# Find elements by type
jl = jt.get_julia_module()
quad_indices = jt.findelem(ring, element_type=jl.QUAD)

# Find by name
marker_indices = jt.findelem(ring, name="M1")

# Build lattice
ring = jt.Lattice([d1, q1, d2, q2])
```

### Serialization

```python
# Load lattice from Julia .jls file
ring = jt.load_lattice("path/to/lattice.jls")

# Save lattice
jt.save_lattice(ring, "path/to/output.jls")
```

### DTPSAD and Setting Element Fields. Currently, we cannot directly assign a DTPSAD variable in Python to Julia lattice element.

```python
# Set TPSA dimension (number of variables)
jt.set_tps_dim(6)

# Convert lattice to TPSA format for automatic differentiation
ring_tpsa = jt.Number2TPSAD(elements)

# Create DTPSAD variables
len_var = jt.DTPSAD(5.27, 1)   # value=5.27, derivative w.r.t. variable 1
k1_var  = jt.DTPSAD(1.2,  2)   # value=1.2, derivative w.r.t. variable 2

# Set fields on TPSA elements using set_field()
# (Direct assignment like `elem.len = val` does not work for DTPSAD values)
jt.set_field(ring_tpsa[1], 'len', len_var)
jt.set_field(ring_tpsa[3], 'k1',  k1_var)

# Convert back to numeric lattice
ring_num = jt.TPSAD2Number(ring_tpsa)

# Compute gradient and Jacobian
grad = jt.Gradient(func, x)
jac  = jt.Jacobian(func, x)
```

### Advanced: Direct Julia Access

```python
# Get Julia module for advanced operations
jl = jt.get_julia_module()

# Evaluate Julia code
result = jt.seval("2 + 2")

# Access Julia types
quad_type = jl.QUAD
crab_type = jl.CRABCAVITY
```

## Contributing

To add new JuTrack.jl functionality:
1. Add wrapper function in `pyJuTrack.py`
2. Handle type conversions (Python ↔ Julia)
3. Add to `__all__` export list

## License

Same as JuTrack.jl
