"""
pyJuTrack - Comprehensive Python Wrapper for JuTrack.jl

This module provides a Python interface to JuTrack.jl particle tracking library.
All type conversions between Python and Julia are handled automatically.

Author: Jinyu Wan
Version: 1.0.0
"""

import os
import sys
import warnings
from typing import List, Union, Optional, Tuple, Any

import numpy as np

os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")
os.environ.setdefault("PYJUTRACK_DISABLE_PYCALL", "1")

import juliacall


_julia_initialized = False
_jl = None
_python_functions = {}


def _pycall_disabled() -> bool:
    return os.environ.get("PYJUTRACK_DISABLE_PYCALL", "1").strip().lower() in ("1", "true", "yes", "on")


def _ensure_conda_pycall_ready(julia_handle):
    """Ensure Conda.jl and PyCall.jl are built for the active Python interpreter."""
    if _pycall_disabled():
        return

    # Make PyCall reuse the currently running Python interpreter instead of
    # trying to bootstrap its own Conda-based Python (which often fails on HPC).
    julia_handle.ENV["PYTHON"] = sys.executable

    try:
        julia_handle.seval(
            """
            import Pkg

            function _pyjutrack_build(pkg::String)
                path = Base.find_package(pkg)
                path === nothing && return

                deps_file = joinpath(dirname(dirname(path)), "deps", "deps.jl")
                if !isfile(deps_file)
                    @info "pyJuTrack: building $(pkg) for the current environment"
                    Pkg.build(pkg)
                end
            end

            _pyjutrack_build("Conda")
            _pyjutrack_build("PyCall")
            nothing
            """
        )
    except Exception as exc:  # pragma: no cover - best-effort setup step
        warnings.warn(
            "pyJuTrack could not auto-configure Conda/PyCall. "
            "Please run `Pkg.build(\"Conda\")` and `Pkg.build(\"PyCall\")` manually.",
            RuntimeWarning,
        )
        raise

def _initialize_julia():
    """Initialize the Julia runtime and load JuTrack.jl"""
    global _julia_initialized, _jl
    
    if _julia_initialized:
        return _jl
    
    import io
    import contextlib
    
    print("Initializing JuTrack.jl...")
    
    # Activate the JuTrack.jl project
    project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    juliacall.Pkg.activate(project_path)
    
    # Import Julia Main module
    from juliacall import Main as jl
    _jl = jl

    # Ensure Julia-side PyCall can see the same Python interpreter and has its
    # Conda deps built before JuTrack (which depends on PyCall) is loaded.
    _ensure_conda_pycall_ready(jl)
    
    # Automatically fix version mismatches by updating the manifest
    # This resolves the StyledStrings issue when switching between Julia versions
    try:
        jl.seval('import Pkg; Pkg.resolve()')
        jl.seval('import Pkg; Pkg.instantiate()')
    except:
        # If resolve/instantiate fail, try more aggressive approach
        try:
            # Update packages to be compatible with current Julia version
            jl.seval('import Pkg; Pkg.update()')
        except:
            pass  # Continue even if this fails
    
    # Load JuTrack (suppress all warnings and output)
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        jl.seval('''
        import Logging
        old_logger = Logging.global_logger()
        Logging.global_logger(Logging.NullLogger())
        ''')
        
        try:
            jl.seval('using JuTrack')
        except Exception as first_error:
            # If JuTrack fails to load, try to fix automatically
            jl.seval('Logging.global_logger(old_logger)')
            
            # Attempt automatic fix: update all packages for current Julia version
            try:
                jl.seval('import Pkg; Pkg.update(); Pkg.resolve(); Pkg.instantiate()')
                # Try loading JuTrack again after fix
                jl.seval('Logging.global_logger(Logging.NullLogger())')
                jl.seval('using JuTrack')
            except Exception as second_error:
                # If still fails, it's a more serious issue
                jl.seval('Logging.global_logger(old_logger)')
                print("\n" + "="*70, file=sys.stderr)
                print("ERROR: Failed to load JuTrack after attempting automatic fix", file=sys.stderr)
                print("="*70, file=sys.stderr)
                print("\nOriginal error:", str(first_error)[:200], file=sys.stderr)
                print("\nAutomatic package update was attempted but failed.", file=sys.stderr)
                print("Please check your Julia installation and JuTrack dependencies.", file=sys.stderr)
                print("="*70 + "\n", file=sys.stderr)
                raise
        except:
            # Catch any other exception and re-raise
            jl.seval('Logging.global_logger(old_logger)')
            raise
        finally:
            jl.seval('Logging.global_logger(old_logger)')
    
    _julia_initialized = True
    
    # Try to import PyCall if available (optional - only needed for some plotting features)
    # If it fails, pyJuTrack will still work for tracking and optimization
    try:
        jl.seval('import PyCall')
        jl.seval('global _python_callbacks = Dict{String, Any}()')
    except:
        # PyCall not available - this is OK, core functionality will still work
        # Users can still use pyJuTrack for tracking, optimization, etc.
        warnings.warn(
            "PyCall could not be loaded. Some advanced plotting features may not work. "
            "Core pyJuTrack functionality (tracking, optimization, etc.) is unaffected.",
            UserWarning
        )
        # Create empty callback dict anyway for compatibility
        jl.seval('global _python_callbacks = Dict{String, Any}()')
    
    return jl

_jl = _initialize_julia()

def _convert_kwargs_to_julia(kwargs):
    """
    Convert Python types in kwargs to Julia types automatically.
    
    This handles common type conversions:
    - Python lists -> Julia Vectors (properly typed)
    - NumPy arrays -> Julia Arrays
    - Python strings -> Julia Strings
    
    Args:
        kwargs: Dictionary of keyword arguments
    
    Returns:
        Dictionary with converted values
    """
    converted = {}
    for key, value in kwargs.items():
        if isinstance(value, list):
            # Convert Python list to properly typed Julia Vector
            # Check if list contains DTPSAD objects
            if value and hasattr(value[0], '__class__') and 'DTPSAD' in str(type(value[0])):
                # List of DTPSAD objects
                converted[key] = _to_julia_dtpsad_vector(value)
            elif all(isinstance(x, (int, float, np.integer, np.floating)) for x in value):
                # Create Float64 vector
                jl_vec = _jl.seval('Float64[]')
                for v in value:
                    _jl.push_b(jl_vec, float(v))
                converted[key] = jl_vec
            else:
                # Fall back to generic Vector
                converted[key] = _jl.Vector(value)
        elif isinstance(value, np.ndarray):
            # Convert NumPy array to Julia array
            if value.ndim == 1:
                # 1D array -> Vector{Float64}
                jl_vec = _jl.seval('Float64[]')
                for v in value:
                    _jl.push_b(jl_vec, float(v))
                converted[key] = jl_vec
            else:
                converted[key] = _to_julia_matrix(value)
        else:
            # Keep as-is
            converted[key] = value
    return converted

def _to_julia_vector(python_list, element_type=None):
    """Convert Python list to properly typed Julia vector"""
    if element_type is not None:
        jl_vec = _jl.seval(f"{element_type}[]")
        for elem in python_list:
            _jl.push_b(jl_vec, elem)
        return jl_vec
    else:
        return _jl.Vector(python_list)

def _to_julia_int_vector(python_list):
    """Convert Python list of integers to Julia Vector{Int64}"""
    # Create the vector directly in Julia and return it
    jl_vec = _jl.seval('Int64[]')
    for val in python_list:
        _jl.push_b(jl_vec, int(val))
    return jl_vec

def _to_julia_float_vector(python_list):
    """Convert Python list of floats to Julia Vector{Float64}"""
    # Create the vector directly in Julia and return it
    jl_vec = _jl.seval('Float64[]')
    for val in python_list:
        _jl.push_b(jl_vec, float(val))
    return jl_vec

def _to_julia_dtpsad_vector(python_list):
    """
    Convert Python list of DTPSAD objects to Julia Vector{DTPSAD{N,T}}.
    
    Args:
        python_list: List of DTPSAD objects (from jt.DTPSAD)
    
    Returns:
        Julia Vector{DTPSAD{N,T}} that can be used in Julia functions
    """
    if not python_list:
        # Return empty DTPSAD vector
        return _jl.seval('DTPSAD[]')
    
    # Get the type of the first element to create properly typed vector
    first_elem = python_list[0]
    elem_type = _jl.typeof(first_elem)
    
    # Create empty vector with the specific DTPSAD type
    jl_vec = _jl.Vector[elem_type]()
    
    # Push each DTPSAD element
    for val in python_list:
        _jl.push_b(jl_vec, val)
    
    return jl_vec

def _to_julia_matrix(numpy_array):
    """Convert NumPy array to Julia Matrix"""
    return _jl.Matrix(np.asarray(numpy_array, dtype=np.float64))

def _to_numpy_array(julia_array):
    """Convert Julia array to NumPy array"""
    return np.array(julia_array)

class Lattice:
    """
    Python wrapper for a Julia lattice vector.
    Handles automatic type conversions for lattice elements.
    Automatically detects Float64 vs DTPSAD parametric types.
    """
    def __init__(self, elements: Union[List, Any] = None):
        """
        Create a lattice from a list of elements or wrap existing Julia lattice.
        
        Args:
            elements: List of lattice elements or Julia lattice object
        """
        if elements is None:
            self._jl_lattice = _jl.seval("AbstractElement{Float64}[]")
        elif isinstance(elements, list):
            if not elements:
                # Empty list - default to Float64
                self._jl_lattice = _jl.seval("AbstractElement{Float64}[]")
            else:
                # Check the first element to determine lattice type
                first_elem = elements[0]
                first_type_str = str(_jl.typeof(first_elem))
                
                # Determine if DTPSAD-parametric or Float64
                is_dtpsad = 'DTPSAD' in first_type_str
                
                # Validate that ALL elements have consistent parametric type
                for i, elem in enumerate(elements):
                    elem_type_str = str(_jl.typeof(elem))
                    elem_is_dtpsad = 'DTPSAD' in elem_type_str
                    if elem_is_dtpsad != is_dtpsad:
                        raise TypeError(
                            f"Inconsistent element types in lattice. "
                            f"First element is {'DTPSAD' if is_dtpsad else 'Float64'}-parametric, "
                            f"but element {i} ({elem.name if hasattr(elem, 'name') else 'unnamed'}) "
                            f"is {'DTPSAD' if elem_is_dtpsad else 'Float64'}-parametric. "
                            f"All elements in a lattice must have the same parametric type."
                        )
                
                if is_dtpsad:
                    # DTPSAD-parametric lattice
                    self._jl_lattice = _jl.seval(f"AbstractElement{{DTPSAD{{NVAR(), Float64}}}}[]")
                else:
                    # Float64 lattice
                    self._jl_lattice = _jl.seval("AbstractElement{Float64}[]")
                
                # Push all elements
                for elem in elements:
                    _jl.push_b(self._jl_lattice, elem)
        else:
            # Assume it's already a Julia lattice
            self._jl_lattice = elements
    
    @property
    def julia_object(self):
        """Get the underlying Julia lattice object"""
        return self._jl_lattice
    
    def __len__(self):
        return len(self._jl_lattice)
    
    def __getitem__(self, index):
        """Access elements by index (0-based Python indexing)"""
        if isinstance(index, slice):
            start, stop, step = index.indices(len(self))
            return [self._jl_lattice[i] for i in range(start, stop, step)]
        else:
            # Julia uses 1-based indexing
            return self._jl_lattice[index]
    
    def __setitem__(self, index, value):
        """Set element by index (0-based Python indexing)"""
        self._jl_lattice[index] = value
    
    def append(self, element):
        """Append an element to the lattice"""
        _jl.push_b(self._jl_lattice, element)
    
    def insert(self, index, element):
        """Insert element at index"""
        _jl.insert_b(self._jl_lattice, index + 1, element)  # Julia 1-based
    
    def total_length(self) -> float:
        """Calculate total length of the lattice"""
        return float(_jl.total_length(self._jl_lattice))
    
    def spos(self, indices: Optional[List[int]] = None):
        """
        Get s-positions of elements.
        
        Args:
            indices: Optional list of element indices (0-based)
        
        Returns:
            NumPy array of s-positions
        """
        if indices is None:
            result = _jl.spos(self._jl_lattice)
        else:
            # Convert to 1-based indexing and create properly typed Julia vector
            indices_1based = [i + 1 for i in indices]
            julia_indices = _to_julia_int_vector(indices_1based)
            result = _jl.spos(self._jl_lattice, julia_indices)
        return _to_numpy_array(result)
    
    def __repr__(self):
        return f"Lattice(n_elements={len(self)}, length={self.total_length():.3f}m)"


class Beam:
    """
    Python wrapper for a Julia Beam object.
    
    Matches Julia constructor:
    Beam(r::Matrix; energy=1e9, np=size(r,1), charge=-1.0, mass=m_e, atn=1.0,
         emittance=zeros(3), centroid=zeros(6), T0=0.0, znbin=99, current=0.0)
    """
    def __init__(self, r, energy: float = 1e9, **kwargs):
        """
        Create a beam from particle coordinates matrix.
        
        Args:
            r: Particle coordinate matrix (n_particles * 6)
               Can be NumPy array, or Julia matrix (Float64 or DTPSAD from jt.zeros)
            energy: Beam energy in eV (default: 1e9)
            **kwargs: Additional keyword arguments:
                - np: Number of real particles (default: size(r,1))
                - charge: Particle charge in units of e (default: -1.0)
                - mass: Particle mass in eV (default: m_e)
                - atn: Atomic number (default: 1.0)
                - emittance: List/array of 3 emittances [ex, ey, ez] (default: [0,0,0])
                - centroid: List/array of 6 centroid values (default: [0,0,0,0,0,0])
                - T0: Revolution period (default: 0.0)
                - znbin: Number of z bins (default: 99)
                - current: Beam current in A (default: 0.0)
        
        Examples:
            >>> # Float64 beam
            >>> particles = np.zeros((5000, 6))
            >>> beam = jt.Beam(particles, energy=10e9, emittance=[20e-9, 1.3e-9, 1.36e-4])
            
            >>> # DTPSAD beam for automatic differentiation
            >>> jt.set_tps_dim(3)
            >>> particles = jt.zeros(jt.DTPSAD, 5000, 6)
            >>> beam = jt.Beam(particles, energy=10e9, np=int(1.72e11*3))
        """
        kwargs = _convert_kwargs_to_julia(kwargs)
        
        if hasattr(r, '_jl_callmethod'):
            r_julia = r
        else:
            r_julia = _to_julia_matrix(r)
        
        self._jl_beam = _jl.Beam(r_julia, energy=energy, **kwargs)
    
    @property
    def julia_object(self):
        """Get the underlying Julia Beam object"""
        return self._jl_beam
    
    def get_particles(self) -> np.ndarray:
        """Get particle coordinates as NumPy array"""
        return _to_numpy_array(self._jl_beam.r)
    
    @property
    def r(self) -> np.ndarray:
        """Particle coordinates as NumPy array"""
        return self.get_particles()
    
    @property
    def energy(self) -> float:
        """Beam energy"""
        return float(self._jl_beam.energy)
    
    @property
    def nmacro(self) -> int:
        """Number of macro particles"""
        return int(self._jl_beam.nmacro)
    
    @property
    def lost_flag(self) -> np.ndarray:
        """
        Lost particle flags as boolean NumPy array.
        
        Returns:
            Boolean array where True indicates a lost particle
        
        Example:
            >>> beam = jt.Beam(particles, energy=3e9)
            >>> jt.linepass(ring, beam)
            >>> lost = beam.lost_flag
            >>> n_lost = np.sum(lost)
            >>> print(f"Lost {n_lost} particles")
        """
        return _to_numpy_array(self._jl_beam.lost_flag)
    
    def __repr__(self):
        particles = self.get_particles()
        return f"Beam(n_particles={particles.shape[0]}, energy={self.energy:.2e}eV)"


def DRIFT(name: str, length: float, **kwargs):
    """Create a drift space"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.DRIFT(name=name, len=length, **kwargs)

def MARKER(name: str, **kwargs):
    """Create a marker element"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.MARKER(name=name, **kwargs)

def QUAD(name: str, length: float, k1: float, **kwargs):
    """Create a quadrupole (thick lens)"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.QUAD(name=name, len=length, k1=k1, **kwargs)

def KQUAD(name: str, length: float, k1: float, **kwargs):
    """Create a quadrupole (kick-drift-kick model)"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.KQUAD(name=name, len=length, k1=k1, **kwargs)

def KSEXT(name: str, length: float, k2: float, **kwargs):
    """Create a sextupole"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.KSEXT(name=name, len=length, k2=k2, **kwargs)

def KOCT(name: str, length: float, k3: float, **kwargs):
    """Create an octupole"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.KOCT(name=name, len=length, k3=k3, **kwargs)

def SBEND(name: str, length: float, angle: float, **kwargs):
    """Create a sector bending magnet"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.SBEND(name=name, len=length, angle=angle, **kwargs)

def RBEND(name: str, length: float, angle: float, **kwargs):
    """Create a rectangular bending magnet"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.RBEND(name=name, len=length, angle=angle, **kwargs)

def LBEND(name: str, length: float, angle: float, **kwargs):
    """Create a linear bending magnet"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.LBEND(name=name, len=length, angle=angle, **kwargs)

def ESBEND(name: str, length: float, angle: float, **kwargs):
    """Create an exact sector bend"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.ESBEND(name=name, len=length, angle=angle, **kwargs)

def ERBEND(name: str, length: float, angle: float, **kwargs):
    """Create an exact rectangular bend"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.ERBEND(name=name, len=length, angle=angle, **kwargs)

def SOLENOID(name: str, length: float, ks: float, **kwargs):
    """Create a solenoid"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.SOLENOID(name=name, len=length, ks=ks, **kwargs)

def thinMULTIPOLE(name: str, **kwargs):
    """Create a thin multipole"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.thinMULTIPOLE(name=name, **kwargs)

def CORRECTOR(name: str, length: float = 0.0, **kwargs):
    """Create a corrector element"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.CORRECTOR(name=name, len=length, **kwargs)

def HKICKER(name: str, length: float = 0.0, **kwargs):
    """Create a horizontal kicker"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.HKICKER(name=name, len=length, **kwargs)

def VKICKER(name: str, length: float = 0.0, **kwargs):
    """Create a vertical kicker"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.VKICKER(name=name, len=length, **kwargs)

def RFCA(name: str, length: float, volt: float, freq: float, **kwargs):
    """Create an RF cavity"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.RFCA(name=name, len=length, volt=volt, freq=freq, **kwargs)

def CRABCAVITY(name: str, **kwargs):
    """Create a crab cavity"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.CRABCAVITY(name=name, **kwargs)

def CRABCAVITY_K2(name: str, **kwargs):
    """Create a crab cavity with K2 model"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.CRABCAVITY_K2(name=name, **kwargs)

def easyCRABCAVITY(name: str, **kwargs):
    """Create a simplified crab cavity"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.easyCRABCAVITY(name=name, **kwargs)

def LongitudinalRLCWake(name: str = "RLCWake", **kwargs):
    """
    Create a longitudinal RLC wakefield element.
    
    Args:
        name: Element name (default: "RLCWake")
        **kwargs: Wakefield parameters
            - freq: Resonant frequency (Hz)
            - Rshunt: Shunt impedance (Ohms)
            - Q0: Quality factor
    
    Example:
        >>> wake = jt.LongitudinalRLCWake("wake1", freq=180e9, Rshunt=5.5e3, Q0=3.0)
    """
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.LongitudinalRLCWake(name=name, **kwargs)

def SPACECHARGE(name: str, **kwargs):
    """Create a space charge element"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.SPACECHARGE(name=name, **kwargs)

def QUAD_SC(name: str, length: float, k1: float, **kwargs):
    """Create a quadrupole with space charge"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.QUAD_SC(name=name, len=length, k1=k1, **kwargs)

def DRIFT_SC(name: str, length: float, **kwargs):
    """Create a drift with space charge"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.DRIFT_SC(name=name, len=length, **kwargs)

def KQUAD_SC(name: str, length: float, k1: float, **kwargs):
    """Create a KQUAD with space charge"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.KQUAD_SC(name=name, len=length, k1=k1, **kwargs)

def KSEXT_SC(name: str, length: float, k2: float, **kwargs):
    """Create a sextupole with space charge"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.KSEXT_SC(name=name, len=length, k2=k2, **kwargs)

def KOCT_SC(name: str, length: float, k3: float, **kwargs):
    """Create an octupole with space charge"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.KOCT_SC(name=name, len=length, k3=k3, **kwargs)

def SBEND_SC(name: str, length: float, angle: float, **kwargs):
    """Create a sector bend with space charge"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.SBEND_SC(name=name, len=length, angle=angle, **kwargs)

def RBEND_SC(name: str, length: float, angle: float, **kwargs):
    """Create a rectangular bend with space charge"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.RBEND_SC(name=name, len=length, angle=angle, **kwargs)

class ElementTypes:
    """Julia element type constants for use with findelem"""
    DRIFT = _jl.DRIFT
    MARKER = _jl.MARKER
    QUAD = _jl.QUAD
    KQUAD = _jl.KQUAD
    KSEXT = _jl.KSEXT
    KOCT = _jl.KOCT
    SBEND = _jl.SBEND
    RBEND = _jl.RBEND
    LBEND = _jl.LBEND
    ESBEND = _jl.ESBEND
    ERBEND = _jl.ERBEND
    SOLENOID = _jl.SOLENOID
    CORRECTOR = _jl.CORRECTOR
    HKICKER = _jl.HKICKER
    VKICKER = _jl.VKICKER
    RFCA = _jl.RFCA
    CRABCAVITY = _jl.CRABCAVITY
    thinMULTIPOLE = _jl.thinMULTIPOLE
    SPACECHARGE = _jl.SPACECHARGE

element_types = ElementTypes()

def TRANSLATION(name: str, dx: float = 0.0, dy: float = 0.0, dz: float = 0.0, **kwargs):
    """Create a translation element"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.TRANSLATION(name=name, dx=dx, dy=dy, dz=dz, **kwargs)

def YROTATION(name: str, angle: float, **kwargs):
    """Create a Y-rotation element"""
    kwargs = _convert_kwargs_to_julia(kwargs)
    return _jl.YROTATION(name=name, angle=angle, **kwargs)

def linepass(lattice: Union[Lattice, List], beam: Beam, refpts=None):
    """
    Track particles through a beamline once.
    
    Args:
        lattice: Lattice object or list of elements
        beam: Beam object
        refpts: Optional list of element indices (0-based) where to save particle coordinates.
                If None, no coordinates are saved (normal linepass).
                If provided, returns list of particle coordinate arrays at those positions.
    
    Returns:
        None if refpts is None, otherwise List[np.ndarray] of particle coordinates at each refpt
    
    Example:
        >>> # Normal tracking without saving
        >>> jt.linepass(ring, beam)
        
        >>> # Track and save at specific positions
        >>> saved = jt.linepass(ring, beam, refpts=[0, 5, 10])
        >>> print(saved[0].shape)  # (n_particles, 6)
    """
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    elif isinstance(lattice, list):
        jl_lattice = Lattice(lattice).julia_object
    else:
        jl_lattice = lattice
    
    if refpts is not None:
        # Convert Python 0-based indices to Julia 1-based indices
        refpts_1based = [i + 1 for i in refpts]
        julia_refpts = _to_julia_int_vector(refpts_1based)
        
        # Call Julia linepass! with refpts - returns saved particle arrays
        saved_particles = _jl.linepass_b(jl_lattice, beam.julia_object, julia_refpts)
        
        # Convert each saved particle matrix to NumPy array
        result = [_to_numpy_array(particles) for particles in saved_particles]
        return result
    else:
        # Normal linepass without saving
        _jl.linepass_b(jl_lattice, beam.julia_object)

def ringpass(lattice: Union[Lattice, List], beam: Beam, num_turns: int):
    """
    Track particles through a ring for multiple turns.
    
    Args:
        lattice: Lattice object or list of elements
        beam: Beam object
        num_turns: Number of turns
    """
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    elif isinstance(lattice, list):
        jl_lattice = Lattice(lattice).julia_object
    else:
        jl_lattice = lattice
    
    _jl.ringpass_b(jl_lattice, beam.julia_object, num_turns)

def plinepass(lattice: Union[Lattice, List], beam: Beam, nthreads: int = None):
    """Parallel linepass tracking"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    elif isinstance(lattice, list):
        jl_lattice = Lattice(lattice).julia_object
    else:
        jl_lattice = lattice
    
    if nthreads is not None:
        _jl.plinepass_b(jl_lattice, beam.julia_object, nthreads)
    else:
        _jl.plinepass_b(jl_lattice, beam.julia_object)

def pringpass(lattice: Union[Lattice, List], beam: Beam, num_turns: int, nthreads: int = None):
    """Parallel ringpass tracking"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    elif isinstance(lattice, list):
        jl_lattice = Lattice(lattice).julia_object
    else:
        jl_lattice = lattice
    
    if nthreads is not None:
        _jl.pringpass_b(jl_lattice, beam.julia_object, num_turns, nthreads)
    else:
        _jl.pringpass_b(jl_lattice, beam.julia_object, num_turns)

def check_lost(beam: Beam) -> np.ndarray:
    """Check which particles are lost"""
    return _to_numpy_array(_jl.check_lost(beam.julia_object))

def set_tps_dim(dim: int):
    """Set TPSA dimension"""
    _jl.set_tps_dim(dim)

def DTPSAD(value: float, var_idx: int = None):
    """
    Create a DTPSAD (dual number) with specified value.
    
    Julia's multiple dispatch allows two forms:
    - DTPSAD(value): Creates a constant (no derivatives)
    - DTPSAD(value, var_idx): Creates a variable with unit derivative at var_idx
    
    Args:
        value: Initial value
        var_idx: Variable index (1-based, Julia convention). 
                 If None, creates a constant.
                 If provided, creates a variable with unit derivative.
    
    Returns:
        DTPSAD object that can be used in arithmetic operations
    
    Examples:
        >>> # Create constant (for accumulators)
        >>> tot = jt.DTPSAD(0.0)
        >>> tot += dlist[0].h21000[0]  # Add DTPSAD values
        
        >>> # Create variable (for optimization parameters)
        >>> x1 = jt.DTPSAD(5.0, 1)  # First variable
        >>> x2 = jt.DTPSAD(3.0, 2)  # Second variable
    """
    if var_idx is None:
        # DTPSAD(value) - create constant
        return _jl.DTPSAD(value)
    else:
        # DTPSAD(value, var_idx) - create variable
        return _jl.DTPSAD(value, var_idx)

def pow(x, n):
    """
    Power function that works with both Python numbers and Julia DTPSAD objects.
    
    This allows Python syntax like pow(x, 2) to work with DTPSAD values.
    Alternative to x**2 which doesn't work with Julia objects in Python.
    
    Args:
        x: Base (Python number or Julia DTPSAD)
        n: Exponent (number)
    
    Returns:
        x^n (same type as input)
    
    Examples:
        >>> x1 = jt.DTPSAD(2.0, 1)
        >>> y = jt.pow(x1, 2)  # Works!
        >>> # Instead of: y = x1**2  # Would fail
    """
    if hasattr(x, '_jl_callmethod'):
        # It's a Julia object, use Julia's ^ operator
        return _jl.seval("^")(x, n)
    else:
        # Python number, use Python's **
        return x ** n

def to_dtpsad_matrix(arr):
    """
    Convert a NumPy array or Python list to a DTPSAD matrix.
    
    Each element will be converted to DTPSAD(value) - a constant with no derivatives.
    Useful for initializing beam distributions from NumPy arrays when using automatic
    differentiation.
    
    Args:
        arr: NumPy array or nested Python list (n*m shape)
    
    Returns:
        Julia Matrix{DTPSAD{NVAR(),Float64}}
    
    Examples:
        >>> import numpy as np
        >>> jt.set_tps_dim(1)
        >>> # Create random distribution in NumPy
        >>> coords = np.random.randn(500, 6) * 1e-4
        >>> # Convert to DTPSAD matrix for AD tracking
        >>> dtpsad_coords = jt.to_dtpsad_matrix(coords)
        >>> ebeam = jt.Beam(dtpsad_coords, energy=10e9)
    """
    # Convert to NumPy array if needed
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr)
    
    n, m = arr.shape
    
    # First convert NumPy array to Julia Float64 matrix
    float_matrix = _jl.Matrix(arr)
    
    # Then convert each element to DTPSAD
    julia_code = f"""
    begin
        mat = Matrix{{DTPSAD{{NVAR(), Float64}}}}(undef, {n}, {m})
        for i in 1:{n}
            for j in 1:{m}
                mat[i, j] = DTPSAD(float_matrix[i, j])
            end
        end
        mat
    end
    """
    # Store the float_matrix in Julia namespace temporarily
    _jl.seval("float_matrix = nothing")  # Initialize
    _jl.float_matrix = float_matrix
    result = _jl.seval(julia_code)
    
    return result

def zeros(dtype, n: int, m: int):
    """
    Create a matrix of zeros with specified type.
    
    Matches Julia's zeros(Type, n, m) syntax.
    
    Args:
        dtype: Data type - use jt.DTPSAD for DTPSAD{NVAR,Float64}, or 'Float64' for Float64
        n: Number of rows
        m: Number of columns
    
    Returns:
        Julia matrix of zeros with specified type
    
    Examples:
        >>> # DTPSAD matrix for automatic differentiation
        >>> jt.set_tps_dim(3)
        >>> particles = jt.zeros(jt.DTPSAD, 5000, 6)
        >>> ebeam = jt.Beam(particles, np=int(1.72e11*3), energy=10e9)
        
        >>> # Float64 matrix
        >>> particles = jt.zeros('Float64', 5000, 6)
        >>> ebeam = jt.Beam(particles, energy=10e9)
    """
    # Check if dtype is DTPSAD type
    if dtype is DTPSAD or (hasattr(dtype, '__name__') and 'DTPSAD' in str(dtype)):
        # Call Julia: zeros(DTPSAD{NVAR(), Float64}, n, m)
        return _jl.seval(f'zeros(DTPSAD{{NVAR(), Float64}}, {n}, {m})')
    elif dtype == 'Float64' or dtype == float:
        # Call Julia: zeros(Float64, n, m)
        return _jl.seval(f'zeros(Float64, {n}, {m})')
    else:
        # Generic case - pass dtype as string
        return _jl.seval(f'zeros({dtype}, {n}, {m})')
    
def conj(x):
    """
    Conjugate function that works with both Python numbers and Julia DTPSAD objects.
    
    This allows Python syntax like conj(x) to work with DTPSAD values.
    
    Args:
        x: Input (Python number or Julia DTPSAD)
    
    Returns:
        Conjugate of x (same type as input)
    
    Examples:
        >>> x1 = jt.DTPSAD(2.0, 1)
        >>> y = jt.conj(x1)  # Works!
    """
    if hasattr(x, '_jl_callmethod'):
        # It's a Julia object, use Julia's conj function
        return _jl.conj(x)
    else:
        # Python number, use Python's complex conjugate
        return np.conj(x)
    
def NVAR():
    """
    Get current TPSA dimension (number of variables).
    
    Returns:
        Integer number of TPSA variables
    
    Example:
        >>> jt.set_tps_dim(3)
        >>> print(jt.NVAR())  # Output: 3
    """
    return int(_jl.NVAR())

def Number2TPSAD(obj):
    """Convert numbers to TPSA format"""
    if isinstance(obj, Lattice):
        tpsa_lattice = _jl.Number2TPSAD(obj.julia_object)
        result = Lattice.__new__(Lattice)
        result._jl_lattice = tpsa_lattice
        return result
    else:
        return _jl.Number2TPSAD(obj)

def convert_to_DTPSAD(value):
    """
    Convert a value to DTPSAD if it's a plain number.
    If already DTPSAD, return as-is.
    
    This helps with automatic type conversion when assigning to DTPSAD lattice elements.
    
    Args:
        value: Number (float, int, or already DTPSAD)
    
    Returns:
        DTPSAD version of the value
    """
    if hasattr(value, '_jl_callmethod'):
        return value
    else:
        return _jl.DTPSAD(float(value))

def TPSAD2Number(obj):
    """Convert from TPSA format to numbers"""
    if isinstance(obj, Lattice):
        num_lattice = _jl.TPSAD2Number(obj.julia_object)
        result = Lattice.__new__(Lattice)
        result._jl_lattice = num_lattice
        return result
    else:
        return _jl.TPSAD2Number(obj)

# Dictionary to store Python callback functions
_python_objective_functions = {}
_function_counter = 0

def Gradient(func, x, return_value=False):
    """
    Calculate gradient using TPSA.
    
    Args:
        func: Python callable OR Julia function
        x: Point to evaluate gradient (list or array)
        return_value: If True, also return function value
    
    Returns:
        gradient (and optionally function value)
    
    Note:
        For best results with TPSA, define objective functions in Julia.
        Python callbacks may not preserve TPSA derivative tracking correctly.
    """
    if isinstance(x, (list, np.ndarray)):
        x_julia = _to_julia_float_vector(x)
    else:
        x_julia = x
    
    if hasattr(func, '_jl_callmethod'):
        result = _jl.Gradient(func, x_julia, return_value)
    else:
        # Detect if function expects a single array or individual arguments
        import inspect
        sig = inspect.signature(func)
        num_params = len(sig.parameters)
        
        global _function_counter
        func_id = f"_pyfunc_{_function_counter}"
        _function_counter += 1
        
        # Store the Python function in Julia's global dict
        _jl.seval(f'global _python_callbacks["{func_id}"] = nothing')
        _jl._python_callbacks[func_id] = func
        
        if num_params == 1:
            # Function expects single array argument: objective(x)
            julia_code = f"""
function {func_id}(args...)
    py_func = _python_callbacks["{func_id}"]
    # Pass as array
    result = py_func([args...])
    return result
end
"""
        else:
            # Function expects individual arguments: objective(x1, x2, ...)
            julia_code = f"""
function {func_id}(args...)
    py_func = _python_callbacks["{func_id}"]
    # Pass DTPSAD objects directly as individual arguments
    result = py_func(args...)
    return result
end
"""
        _jl.seval(julia_code)
        julia_func = getattr(_jl, func_id)
        result = _jl.Gradient(julia_func, x_julia, return_value)
    
    if return_value:
        grad, val = result
        return _to_numpy_array(grad), val
    else:
        return _to_numpy_array(result)

def Jacobian(func, x, return_value=False):
    """
    Calculate Jacobian matrix using TPSA.
    
    Args:
        func: Python callable function(x) -> tuple OR function(x1, x2, ...) -> tuple
              Returns tuple/list of results
        x: Point to evaluate Jacobian (list or array)
        return_value: If True, return (jacobian, function_values), else just jacobian
    
    Returns:
        If return_value=False: Jacobian matrix (numpy array)
        If return_value=True: (Jacobian matrix, list of function values)
    
    Examples:
        >>> # Style 1: Individual arguments
        >>> def my_func(x1, x2):
        ...     return (x1*x1 + x2, x2*x1)
        >>> jac = jt.Jacobian(my_func, [1.0, 2.0])
        
        >>> # Style 2: Array argument
        >>> def my_func(x):
        ...     return (x[0]*x[0] + x[1], x[1]*x[0])
        >>> jac = jt.Jacobian(my_func, [1.0, 2.0])
    """
    global _function_counter
    
    if isinstance(x, (list, np.ndarray)):
        x_julia = _to_julia_float_vector(x)
    else:
        x_julia = x
    
    # Check if func is already a Julia function
    if hasattr(func, '_jl_callmethod'):
        # It's a Julia function, use directly
        result = _jl.Jacobian(func, x_julia, return_value)
    else:
        # It's a Python function - detect signature
        import inspect
        sig = inspect.signature(func)
        num_params = len(sig.parameters)
        
        # Store Python function in Julia's global dictionary
        func_id = f"_pyfunc_{_function_counter}"
        _function_counter += 1
        
        # Store the function in Julia's _python_callbacks dictionary
        _jl._python_callbacks[func_id] = func
        
        if num_params == 1:
            # Function expects single array: func(x)
            julia_code = f"""
function {func_id}(args...)
    py_func = _python_callbacks["{func_id}"]
    result_py = py_func([args...])
    n = pylen(result_py)
    result_vec = [pyconvert(Any, result_py[i-1]) for i in 1:n]
    return tuple(result_vec...)
end
"""
        else:
            # Function expects individual arguments: func(x1, x2, ...)
            julia_code = f"""
function {func_id}(args...)
    py_func = _python_callbacks["{func_id}"]
    result_py = py_func(args...)
    n = pylen(result_py)
    result_vec = [pyconvert(Any, result_py[i-1]) for i in 1:n]
    return tuple(result_vec...)
end
"""
        _jl.seval(julia_code)
        julia_func = getattr(_jl, func_id)
        result = _jl.Jacobian(julia_func, x_julia, return_value)
    
    if return_value:
        # Julia returns (Jacobian_matrix, [val1, val2, ...])
        jac, vals = result
        return _to_numpy_array(jac), vals
    else:
        return _to_numpy_array(result)

# Helper function accessible from Julia
def get_python_function(func_id):
    """Internal: Get stored Python function by ID"""
    return _python_objective_functions[func_id]

def twissring(lattice, dp=0.0, refpts=None, **kwargs):
    """
    Calculate Twiss parameters around a ring.
    
    Args:
        lattice: Lattice object
        dp: Momentum deviation
        refpts: Reference points (list of indices)
        **kwargs: E0, m0, etc.
    
    Returns:
        List of Twiss objects
    """
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    if refpts is None:
        refpts = list(range(1, len(jl_lattice) + 1))
        julia_refpts = _to_julia_int_vector(refpts)
    elif isinstance(refpts, (list, np.ndarray)):
        julia_refpts = _to_julia_int_vector([int(r) for r in refpts])
    else:
        julia_refpts = refpts
    
    return _jl.twissring(jl_lattice, dp, 0, julia_refpts, **kwargs)

def twissline(lattice, dp=0.0, refpts=None, **kwargs):
    """Calculate Twiss parameters along a beamline"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
        
    if refpts is None:
        refpts = list(range(1, len(jl_lattice) + 1))
        julia_refpts = _to_julia_int_vector(refpts)
    elif isinstance(refpts, (list, np.ndarray)):
        julia_refpts = _to_julia_int_vector([int(r) for r in refpts])
    else:
        julia_refpts = refpts
    
    return _jl.twissline(jl_lattice, dp, 0, julia_refpts, **kwargs)

def ADtwissring(lattice, dp=0.0, refpts=None, changed_idx=None, changed_ele=None, **kwargs):
    """AD-compatible twissring"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    if refpts is None:
        refpts = list(range(1, len(jl_lattice) + 1))
        julia_refpts = _to_julia_int_vector(refpts)
    elif isinstance(refpts, (list, np.ndarray)):
        julia_refpts = _to_julia_int_vector([int(r) for r in refpts])
    else:
        julia_refpts = refpts
    
    if changed_idx is not None and changed_ele is not None:
        return _jl.ADtwissring(jl_lattice, dp, 0, julia_refpts, changed_idx, changed_ele, **kwargs)
    else:
        return _jl.ADtwissring(jl_lattice, dp, 0, julia_refpts, **kwargs)

def periodicEdwardsTengTwiss(lattice, dp, order, **kwargs):
    """Calculate periodic Twiss parameters"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    return _jl.periodicEdwardsTengTwiss(jl_lattice, dp, order, **kwargs)

def ADperiodicEdwardsTengTwiss(lattice, dp, order, changed_idx=None, changed_ele=None, **kwargs):
    """AD-compatible periodic Twiss"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    if changed_idx is not None and changed_ele is not None:
        return _jl.ADperiodicEdwardsTengTwiss(jl_lattice, dp, order, changed_idx, changed_ele, **kwargs)
    else:
        return _jl.ADperiodicEdwardsTengTwiss(jl_lattice, dp, order, **kwargs)

def findm66(lattice, dp=0.0, order=1, **kwargs):
    """Find 6x6 transfer matrix
    
    Args:
        lattice: Lattice object or Julia lattice
        dp: Momentum deviation
        order: TPSA order (default: 1)
        **kwargs: Additional arguments (E0, m0, orb)
    """
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    result = _jl.findm66(jl_lattice, dp, order, **kwargs)
    return _to_numpy_array(result)

def fastfindm66(lattice, dp=0.0, **kwargs):
    """Fast calculation of 6x6 transfer matrix"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    result = _jl.fastfindm66(jl_lattice, dp, **kwargs)
    return _to_numpy_array(result)

def findm66_refpts(lattice, refpts, dp=0.0, order=1, **kwargs):
    """Find transfer matrices at reference points
    
    Args:
        lattice: Lattice object or Julia lattice
        refpts: Reference point indices
        dp: Momentum deviation
        order: TPSA order (default: 1)
        **kwargs: Additional arguments (E0, m0, orb)
    """
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    if isinstance(refpts, (list, np.ndarray)):
        julia_refpts = _to_julia_int_vector([int(r) for r in refpts])
    else:
        julia_refpts = refpts
    
    return _jl.findm66_refpts(jl_lattice, dp, order, julia_refpts, **kwargs)

def find_closed_orbit(lattice, dp=0.0, **kwargs):
    """Find closed orbit"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    result = _jl.find_closed_orbit(jl_lattice, dp, **kwargs)
    return _to_numpy_array(result)

def find_closed_orbit_4d(lattice, dp=0.0, **kwargs):
    """Find 4D closed orbit"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice

    result = _jl.find_closed_orbit_4d(jl_lattice, dp=dp, **kwargs)
    return _to_numpy_array(result)

def find_closed_orbit_6d(lattice, **kwargs):
    """Find 6D closed orbit"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    result = _jl.find_closed_orbit_6d(jl_lattice, **kwargs)
    return _to_numpy_array(result)

def total_length(lattice) -> float:
    """Calculate total length of lattice"""
    if isinstance(lattice, Lattice):
        return lattice.total_length()
    else:
        return float(_jl.total_length(lattice))

def spos(lattice, indices: Optional[List[int]] = None):
    """Get s-positions"""
    if isinstance(lattice, Lattice):
        return lattice.spos(indices)
    else:
        if indices is None:
            return _to_numpy_array(_jl.spos(lattice))
        else:
            indices_1based = [i + 1 for i in indices]
            julia_indices = _to_julia_int_vector(indices_1based)
            result = _jl.spos(lattice, julia_indices)
            return _to_numpy_array(result)

def findelem(lattice, element_type=None, name=None, field=None, value=None):
    """
    Find elements in lattice.
    
    Args:
        lattice: Lattice object
        element_type: Type of element (e.g., _jl.QUAD)
        name: Element name
        field: Field name
        value: Field value
    
    Returns:
        List of indices (1-based Julia convention)
    """
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    if element_type is not None:
        return _jl.findelem(jl_lattice, element_type)
    elif name is not None:
        return _jl.findelem(jl_lattice, _jl.Symbol("name"), name)
    elif field is not None and value is not None:
        return _jl.findelem(jl_lattice, _jl.Symbol(field), value)
    else:
        raise ValueError("Must provide element_type, name, or field+value")

def buildlatt(elements: List) -> Lattice:
    """Build a lattice from list of elements"""
    return Lattice(elements)

def load_lattice(filename: str):
    """Load lattice from .jls file"""
    filename = filename.replace('\\', '/')
    jl_lattice = _jl.seval(f'using Serialization; deserialize("{filename}")')
    
    if hasattr(jl_lattice, '__len__'):
        return Lattice(jl_lattice)
    else:
        return jl_lattice

def save_lattice(lattice, filename: str):
    """Save lattice to .jls file"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    filename = filename.replace('\\', '/')
    
    # Declare global variable and assign in one step
    _jl.seval('global _temp_lattice_tosave = nothing')
    _jl._temp_lattice_tosave = jl_lattice
    _jl.seval(f'using Serialization; serialize("{filename}", Main._temp_lattice_tosave)')
    # Clean up
    _jl.seval('global _temp_lattice_tosave = nothing')

def initilize_6DGaussiandist(beam: Beam, **kwargs):
    """Initialize 6D Gaussian distribution"""
    _jl.initilize_6DGaussiandist_b(beam.julia_object, **kwargs)

def initilize_zslice(beam: Beam, **kwargs):
    """Initialize z-slices for beam"""
    _jl.initilize_zslice_b(beam.julia_object, **kwargs)

def histogram1DinZ(beam: Beam):
    """
    Calculate 1D histogram in longitudinal (z) direction.
    
    Creates histogram bins for the beam distribution along z-axis.
    Required before using wakefield elements.
    
    Args:
        beam: Beam object
    
    Example:
        >>> beam = jt.Beam(particles, energy=10e9)
        >>> jt.histogram1DinZ(beam)
    """
    _jl.histogram1DinZ_b(beam.julia_object)

def get_emittance(beam: Beam) -> Tuple[float, float, float]:
    """Get beam emittances (ex, ey, ez)"""
    _jl.get_emittance_b(beam.julia_object)
    # Function modifies beam in place, return the emittance values
    emittance = beam.julia_object.emittance
    return (emittance[0], emittance[1], emittance[2])

def get_centroid(beam: Beam) -> np.ndarray:
    """Get beam centroid"""
    _jl.get_centroid_b(beam.julia_object)
    # Function modifies beam in place, return the centroid
    return _to_numpy_array(beam.julia_object.centroid)

def get_2nd_moment(beam: Beam) -> np.ndarray:
    """Get 2nd moments of beam"""
    _jl.get_2nd_moment_b(beam.julia_object)
    # Function modifies beam in place, return the moment2nd matrix
    return _to_numpy_array(beam.julia_object.moment2nd)

def dynamic_aperture(lattice, **kwargs):
    """Calculate dynamic aperture"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    return _jl.dynamic_aperture(jl_lattice, **kwargs)

def FMA(lattice, num_turns: int, **kwargs):
    """Frequency Map Analysis"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice

    return _jl.FMA(jl_lattice, num_turns, **kwargs)

def plot_fma(rows, **kwargs):
    """
    Plot Frequency Map Analysis results.
    
    Creates a two-panel plot:
    - Left panel: Initial particle positions (x, y) colored by diffusion
    - Right panel: Tune diagram (nux, nuy) with resonance lines
    
    Args:
        rows: FMA results from FMA() function (list of NamedTuples)
        **kwargs: Optional plotting parameters
            - figsize: Tuple (width, height) in inches, default (10, 4)
            - s: Marker size, default 10
            - x_min, x_max: X-axis limits for initial conditions plot
            - y_min, y_max: Y-axis limits for initial conditions plot
            - resonance_lines: Bool, whether to plot resonance lines, default True
            - resonance_orders: List of resonance orders to plot, default [1,2,3,4]
            - filepath: String, output file path, default "fma_plot.png"
    
    Example:
        >>> rows = jt.FMA(ring, 256)
        >>> jt.plot_fma(rows, figsize=(12, 5), resonance_orders=[1, 2, 3])
    
    Note:
        Requires PyCall and matplotlib installed in Julia environment.
        If PyCall is not available, this function will raise an error.
    """
    try:
        return _jl.plot_fma(rows, **kwargs)
    except Exception as e:
        if 'PyCall' in str(e) or 'matplotlib' in str(e):
            raise RuntimeError(
                "plot_fma requires PyCall and matplotlib in Julia environment. "
                "These dependencies are optional and not available on this system. "
                "Consider exporting FMA data and plotting with Python's matplotlib directly."
            ) from e
        else:
            raise

def plot_lattice(lattice, scale=0.25, axis=True, savepath=None):
    """
    Plot lattice layout showing element positions and types.
    
    Visualizes the beamline layout with different colors for different element types.
    
    Args:
        lattice: Lattice object or list of elements
        scale: Float, vertical scale factor for element heights, default 0.25
        axis: Bool, whether to show axis, default True
    
    Example:
        >>> line = jt.Lattice([drift1, quad1, drift2, bend1])
        >>> jt.plot_lattice(line, scale=0.3, axis=True)
    
    Note:
        Requires PyCall and matplotlib installed in Julia environment.
        If PyCall is not available, this function will raise an error.
    """
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    try:
        return _jl.plot_lattice(jl_lattice, scale, axis, savepath)
    except Exception as e:
        if 'PyCall' in str(e) or 'matplotlib' in str(e):
            raise RuntimeError(
                "plot_lattice requires PyCall and matplotlib in Julia environment. "
                "These dependencies are optional and not available on this system. "
                "Consider using other lattice visualization tools."
            ) from e
        else:
            raise


def computeRDT(lattice, indices=None, **kwargs):
    """
    Compute Resonance Driving Terms.
    
    This function automatically uses TPSA-based AD if lattice is DTPSAD type.
    
    Args:
        lattice: Lattice object (Float64 or DTPSAD type)
        indices: List of element indices to compute RDTs at (optional)
                 If None, uses all MARKER elements
        **kwargs: Additional keyword arguments
            - E0: Energy (default 3e9)
            - m0: Mass (default m_e)
            - chromatic, coupling, geometric1, geometric2, tuneshifts: bool flags
    
    Returns:
        Tuple (dlist, tune) where:
            - dlist: List of DrivingTerms objects (one per index)
            - tune: Tuple of tunes (nux, nuy)
        
    Example (TPSA-based AD for optimization):
        >>> # Convert to TPSA format
        >>> ring_tpsa = jt.Number2TPSAD(ring)
        >>> # Compute RDTs at first marker
        >>> dlist, tune = jt.computeRDT(ring_tpsa, [1], E0=3e9, m0=jt.m_e)
        >>> h11001 = dlist[0].h11001[0]  # Access first RDT's h11001 value
    """
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    # If no indices provided, find all markers
    if indices is None:
        indices = _jl.findelem(jl_lattice, _jl.MARKER)
        if len(indices) == 0:
            indices = [1]  # Default to first element
    
    if isinstance(indices, (list, np.ndarray)):
        indices = _to_julia_int_vector(indices)
    
    return _jl.computeRDT(jl_lattice, indices, **kwargs)

def ADcomputeRDT(lattice, index, changed_ids, changed_elems, **kwargs):
    """
    Enzyme AD-compatible RDT calculation.
    
    This is for Enzyme-based automatic differentiation (Julia-side AD).
    For Python TPSA optimization, use computeRDT() with DTPSAD lattice instead.
    
    Args:
        lattice: Lattice object
        index: Indices to compute RDTs at
        changed_ids: IDs of changed elements
        changed_elems: Changed element objects
        **kwargs: Additional keyword arguments (E0, m0, chromatic, etc.)
    
    Returns:
        RDT results
        
    Note:
        This function is primarily for Julia-side Enzyme AD.
        Python users should use the TPSA approach with computeRDT().
    """
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    return _jl.ADcomputeRDT(jl_lattice, index, changed_ids, changed_elems, **kwargs)

def gettune(lattice, dp=0.0, **kwargs) -> Tuple[float, float]:
    """Get betatron tunes (nux, nuy)"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice

    result = _jl.gettune(jl_lattice, dp=dp, **kwargs)
    return (float(result[0]), float(result[1]))

def getchrom(lattice, dp=0.0, **kwargs) -> Tuple[float, float]:
    """Get chromaticities (xi_x, xi_y)"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice

    result = _jl.getchrom(jl_lattice, dp=dp, **kwargs)
    return (float(result[0]), float(result[1]))

def rad_on():
    """Turn on radiation"""
    _jl.rad_on_b()

def rad_off():
    """Turn off radiation"""
    _jl.rad_off_b()

def tracking_U0(lattice, **kwargs) -> float:
    """Calculate U0 via tracking"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    return float(_jl.tracking_U0(jl_lattice, **kwargs))

def integral_U0(lattice, **kwargs) -> float:
    """Calculate U0 via integration"""
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    return float(_jl.integral_U0(jl_lattice, **kwargs))

def ringpara(lattice, energy: float = 3e9, Vrf: float = 0.0, harm: int = 1,
             freq_rf: float = 0.0, dp: float = 0.0, print_summary: bool = False):
    """
    Calculate and optionally print ring parameters including equilibrium emittance,
    damping times, energy spread, and RF-dependent parameters.
    
    This function computes radiation integrals from dipoles, quadrupoles, and wigglers,
    then calculates equilibrium beam parameters based on synchrotron radiation theory.
    
    Args:
        lattice: Lattice object or list of elements
        energy: Beam energy in eV (default: 3e9)
        Vrf: RF voltage per cell [V] (default: 0.0)
        harm: Harmonic number (default: 1)
        freq_rf: RF frequency in Hz (default: 0.0)
        dp: Momentum deviation (default: 0.0)
        print_summary: If True, print parameter summary (default: True)
    
    Returns:
        Dictionary with all calculated parameters including:
        - Basic parameters: energy, gamma, beta, circumference, etc.
        - Tunes: nux, nuy, chromx, chromy
        - Radiation integrals: I1, I2, I3, I4, I5, I6, Iv
        - Damping: Jx, Jy, Je, damping times and rates
        - Equilibrium: emittance, energy spread, U0
        - RF parameters: phi_s, nus, delta_max, bunchlength (if Vrf > 0)
    
    Examples:
        >>> # Print summary
        >>> jt.ringpara(ring, energy=3e9, Vrf=3.2e6, harm=372, freq_rf=476e6)
        
        >>> # Get parameters as dictionary
        >>> params = jt.ringpara(ring, energy=3e9, Vrf=3e6, harm=372,
        ...                      freq_rf=476e6, print_summary=False)
        >>> print(f"Natural emittance: {params['emittx']*1e9} nm⋅rad")
    
    Reference:
        Based on AT's ringpara function. See also:
        - H. Wiedemann, "Particle Accelerator Physics"
        - M. Sands, "The Physics of Electron Storage Rings"
    """
    if isinstance(lattice, Lattice):
        jl_lattice = lattice.julia_object
    else:
        jl_lattice = lattice
    
    # Call Julia function and get results as NamedTuple
    jl_results = _jl.ringpara(jl_lattice, energy=energy, Vrf=Vrf, harm=harm,
                              freq_rf=freq_rf, dp=dp, print_summary=print_summary)
    
    # Convert Julia NamedTuple to Python dictionary
    # Access NamedTuple fields using property access
    results = {
        'E0': jl_results.E0,
        'E0_GeV': jl_results.E0_GeV,
        'gamma': jl_results.gamma,
        'beta': jl_results.beta,
        'R': jl_results.R,
        'circumference': jl_results.circumference,
        'T0': jl_results.T0,
        'frev': jl_results.frev,
        'alphac': jl_results.alphac,
        'etac': jl_results.etac,
        'nux': jl_results.nux,
        'nuy': jl_results.nuy,
        'chromx': jl_results.chromx,
        'chromy': jl_results.chromy,
        'I1': jl_results.I1,
        'I2': jl_results.I2,
        'I3': jl_results.I3,
        'I4': jl_results.I4,
        'I5': jl_results.I5,
        'I6': jl_results.I6,
        'Iv': jl_results.Iv,
        'Jx': jl_results.Jx,
        'Jy': jl_results.Jy,
        'Je': jl_results.Je,
        'U0': jl_results.U0,
        'sigma_E': jl_results.sigma_E,
        'emittx': jl_results.emittx,
        'emitty_d': jl_results.emitty_d,
        'emitty_lim': jl_results.emitty_lim,
        'dampingtime_x': jl_results.dampingtime_x,
        'dampingtime_y': jl_results.dampingtime_y,
        'dampingtime_E': jl_results.dampingtime_E,
        'dampingalpha_x': jl_results.dampingalpha_x,
        'dampingalpha_y': jl_results.dampingalpha_y,
        'dampingalpha_E': jl_results.dampingalpha_E,
        'Vrf': jl_results.Vrf,
        'harm': jl_results.harm,
        'freq_rf': jl_results.freq_rf,
    }
    
    # Add RF-dependent parameters (may be NaN)
    try:
        results['phi_s'] = float(jl_results.phi_s)
        results['nus'] = float(jl_results.nus)
        results['delta_max'] = float(jl_results.delta_max)
        results['bunchlength'] = float(jl_results.bunchlength)
    except:
        results['phi_s'] = float('nan')
        results['nus'] = float('nan')
        results['delta_max'] = float('nan')
        results['bunchlength'] = float('nan')
    
    return results

# Physical constants
m_e = float(_jl.m_e)  # Electron mass (eV)
m_p = float(_jl.m_p)  # Proton mass (eV)
m_goldion = float(_jl.m_goldion)  # Gold ion mass (eV)
charge_e = float(_jl.charge_e)  # Elementary charge (C)
speed_of_light = float(_jl.speed_of_light)  # Speed of light (m/s)
epsilon_0 = float(_jl.epsilon_0)  # Permittivity of free space
CGAMMA = float(_jl.CGAMMA)
__E0 = float(_jl.__E0)  # Electron rest mass energy (GeV)
__HBAR_C = float(_jl.__HBAR_C)  # ħc in MeV⋅fm
Cq = float(_jl.Cq)  # Quantum constant (m)  

# Coordinate limits
CoordLimit = float(_jl.CoordLimit)
AngleLimit = float(_jl.AngleLimit)

def randn_approx(n: int, m: int) -> np.ndarray:
    """Generate N×M matrix of approximately normal random numbers
    
    Args:
        n: Number of rows
        m: Number of columns
    
    Returns:
        NumPy array of shape (n, m) with approximately normal distribution
    """
    result = _jl.randn_approx(n, m)
    return _to_numpy_array(result)

def symplectic(matrix: np.ndarray) -> float:
    """Check symplecticity of matrix"""
    julia_matrix = _to_julia_matrix(matrix)
    return float(_jl.symplectic(julia_matrix))

def matrix_to_array(matrix):
    """Convert Julia matrix to array"""
    return _to_numpy_array(_jl.matrix_to_array(matrix))

def array_to_matrix(array: np.ndarray):
    """Convert array to Julia matrix"""
    return _jl.array_to_matrix(_to_julia_matrix(array))

# ADVANCED: Direct Julia Access
def get_julia_module():
    """Get direct access to Julia Main module (for advanced use only)"""
    return _jl

def include(filename: str):
    """
    Include and execute a Julia file.
    
    Args:
        filename: Path to Julia file (can use forward or back slashes)
    
    Example:
        >>> jt.include("../src/demo/SPEAR3/spear3.jl")
        >>> ring = jt.Lattice(jt.call_julia("spear3"))
    """
    filename = filename.replace('\\', '/')
    return _jl.include(filename)

def call_julia(func_name: str, *args, **kwargs):
    """
    Call a Julia function by name.
    
    Args:
        func_name: Name of Julia function
        *args: Positional arguments
        **kwargs: Keyword arguments
    
    Returns:
        Result from Julia function
    
    Example:
        >>> jt.include("../src/demo/SPEAR3/spear3.jl")
        >>> elements = jt.call_julia("spear3")
        >>> ring = jt.Lattice(elements)
    """
    julia_func = getattr(_jl, func_name)
    return julia_func(*args, **kwargs)

def seval(code: str):
    """Evaluate Julia code string (for advanced use)"""
    return _jl.seval(code)

# MODULE METADATA
__version__ = "1.0.0"
__author__ = "JuTrack.jl Integration"

__all__ = [
    # Classes
    'Lattice', 'Beam',
    
    # Element types (for findelem by type)
    'element_types', 'ElementTypes',
    
    # Basic elements
    'DRIFT', 'MARKER',
    
    # Magnets
    'QUAD', 'KQUAD', 'KSEXT', 'KOCT',
    'SBEND', 'RBEND', 'LBEND', 'ESBEND', 'ERBEND',
    'SOLENOID', 'thinMULTIPOLE',
    
    # Correctors
    'CORRECTOR', 'HKICKER', 'VKICKER',
    
    # RF and special
    'RFCA', 'CRABCAVITY', 'CRABCAVITY_K2', 'easyCRABCAVITY',
    'LongitudinalRLCWake',
    
    # Space charge
    'SPACECHARGE', 'QUAD_SC', 'DRIFT_SC', 'KQUAD_SC', 'KSEXT_SC',
    'KOCT_SC', 'SBEND_SC', 'RBEND_SC',
    
    # Transformations
    'TRANSLATION', 'YROTATION',
    
    # Tracking
    'linepass', 'ringpass', 'plinepass', 'pringpass', 'check_lost',
    
    # TPSA
    'set_tps_dim', 'DTPSAD', 'NVAR', 'zeros', 'to_dtpsad_matrix', 'pow',
    'Number2TPSAD', 'TPSAD2Number', 'convert_to_DTPSAD', 'Gradient', 'Jacobian',
    
    # Twiss and optics
    'twissring', 'twissline', 'ADtwissring', 
    'periodicEdwardsTengTwiss', 'ADperiodicEdwardsTengTwiss',
    
    # Transfer matrices
    'findm66', 'fastfindm66', 'findm66_refpts',
    
    # Closed orbit
    'find_closed_orbit', 'find_closed_orbit_4d', 'find_closed_orbit_6d',
    
    # Lattice utilities
    'total_length', 'spos', 'findelem', 'buildlatt',
    
    # Serialization
    'load_lattice', 'save_lattice',
    
    # Plotting
    'plot_fma', 'plot_lattice',
    
    # Beam initialization
    'initilize_6DGaussiandist', 'initilize_zslice', 'histogram1DinZ',
    'get_emittance', 'get_centroid', 'get_2nd_moment',
    
    # Dynamic aperture and FMA
    'dynamic_aperture', 'FMA',
    
    # RDT
    'computeRDT', 'ADcomputeRDT',
    
    # Tune and chromaticity
    'gettune', 'getchrom',
    
    # Radiation
    'rad_on', 'rad_off', 'tracking_U0', 'integral_U0', 'ringpara',
    
    # Constants
    'm_e', 'm_p', 'm_goldion', 'charge_e', 'speed_of_light',
    'epsilon_0', 'CGAMMA', '__E0', '__HBAR_C', 'Cq',
    'CoordLimit', 'AngleLimit',
    
    # Utilities
    'randn_approx', 'symplectic', 'matrix_to_array', 'array_to_matrix',
    
    # Advanced
    'get_julia_module', 'seval',
]
