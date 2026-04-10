using Printf
using LinearAlgebra

"""
    translate_Julia_to_Python(elements::Vector, filename::String)

Translate a Julia sequence of lattice elements to a Python script for pyJuTrack.
Skips default values (zeros for T1/T2, identity for R1/R2) to keep output clean.

# Arguments
- elements::Vector: Vector of JuTrack lattice elements
- filename::String: Output Python filename

# Example
```julia
using JuTrack
d1 = DRIFT(name="D1", len=0.5)
q1 = KQUAD(name="Q1", len=0.3, k1=1.2)
elements = [d1, q1]
translate_Julia_to_Python(elements, "lattice.py")
```
"""
function translate_Julia_to_Python(elements::Vector, filename::String)
    open(filename, "w") do io
        # Write header
        println(io, "# Auto-generated Python lattice file for pyJuTrack")
        println(io, "# Generated from JuTrack.jl elements")
        println(io, "")
        println(io, "import pyJuTrack as jt")
        println(io, "import numpy as np")
        println(io, "")
        println(io, "# Define lattice elements")
        println(io, "")
        
        element_names = String[]
        
        for ele in elements
            # Get element type
            eletype = ele.eletype
            name = ele.name
            
            # Build Python constructor call
            py_str = element_to_python(ele, eletype, name)
            
            if !isnothing(py_str)
                println(io, py_str)
                push!(element_names, name)
            end
        end
        
        # Write element list
        println(io, "")
        println(io, "# Create element list")
        println(io, "elements = [", join(element_names, ", "), "]")
    end
    
    println("Python lattice file created: $filename")
end

"""
Convert a single Julia element to Python pyJuTrack constructor string.
Only includes non-default parameters.
"""
function element_to_python(ele, eletype::String, name::String)
    # Handle different element types
    if eletype == "DRIFT"
        return drift_to_python(ele, name)
    elseif eletype == "QUAD"
        return quad_to_python(ele, name)
    elseif eletype == "KQUAD"
        return kquad_to_python(ele, name)
    elseif eletype == "KSEXT"
        return ksext_to_python(ele, name)
    elseif eletype == "KOCT"
        return koct_to_python(ele, name)
    elseif eletype == "SBEND"
        return sbend_to_python(ele, name)
    elseif eletype == "RBEND"
        return rbend_to_python(ele, name)
    elseif eletype == "LBEND"
        return lbend_to_python(ele, name)
    elseif eletype == "ESBEND"
        return esbend_to_python(ele, name)
    elseif eletype == "ERBEND"
        return erbend_to_python(ele, name)
    elseif eletype == "SOLENOID"
        return solenoid_to_python(ele, name)
    elseif eletype == "thinMULTIPOLE"
        return thinmultipole_to_python(ele, name)
    elseif eletype == "HKICKER"
        return hkicker_to_python(ele, name)
    elseif eletype == "VKICKER"
        return vkicker_to_python(ele, name)
    elseif eletype == "WIGGLER"
        return wiggler_to_python(ele, name)
    elseif eletype == "RFCA"
        return rfca_to_python(ele, name)
    elseif eletype == "MARKER"
        return marker_to_python(ele, name)
    elseif eletype == "CORRECTOR"
        return corrector_to_python(ele, name)
    elseif eletype == "DRIFT_SC"
        return drift_sc_to_python(ele, name)
    elseif eletype == "DRIFT_SC2P5D"
        return drift_sc2p5d_to_python(ele, name)
    elseif eletype == "QUAD_SC"
        return quad_sc_to_python(ele, name)
    elseif eletype == "QUAD_SC2P5D"
        return quad_sc2p5d_to_python(ele, name)
    elseif eletype == "KQUAD_SC"
        return kquad_sc_to_python(ele, name)
    elseif eletype == "KQUAD_SC2P5D"
        return kquad_sc2p5d_to_python(ele, name)
    elseif eletype == "KSEXT_SC"
        return ksext_sc_to_python(ele, name)
    elseif eletype == "KSEXT_SC2P5D"
        return ksext_sc2p5d_to_python(ele, name)
    elseif eletype == "KOCT_SC"
        return koct_sc_to_python(ele, name)
    elseif eletype == "KOCT_SC2P5D"
        return koct_sc2p5d_to_python(ele, name)
    elseif eletype == "SBEND_SC"
        return sbend_sc_to_python(ele, name)
    elseif eletype == "SBEND_SC2P5D"
        return sbend_sc2p5d_to_python(ele, name)
    elseif eletype == "RBEND_SC"
        return rbend_sc_to_python(ele, name)
    elseif eletype == "RBEND_SC2P5D"
        return rbend_sc2p5d_to_python(ele, name)
    elseif eletype == "SPACECHARGE2P5D"
        return spacecharge2p5d_to_python(ele, name)
    else
        @warn "Unknown element type: $eletype, skipping element $name"
        return nothing
    end
end

# Helper functions to check if arrays are at default values
is_zero_vector(v::Vector) = all(x -> x == 0, v)
is_identity_matrix(m::Matrix) = m == I(size(m, 1))

function format_vector_python(v::Vector)
    """Format a vector for Python numpy array"""
    return "np.array([" * join([@sprintf("%.10e", x) for x in v], ", ") * "])"
end

function format_matrix_python(m::Matrix)
    """Format a matrix for Python numpy array"""
    rows = []
    for i in 1:size(m, 1)
        row_str = "[" * join([@sprintf("%.10e", m[i, j]) for j in 1:size(m, 2)], ", ") * "]"
        push!(rows, row_str)
    end
    return "np.array([" * join(rows, ", ") * "])"
end

function add_misalignment_params(params::Vector{String}, ele)
    """Add T1, T2, R1, R2 parameters if non-default"""
    if hasfield(typeof(ele), :T1) && !is_zero_vector(ele.T1)
        push!(params, "T1=$(format_vector_python(ele.T1))")
    end
    if hasfield(typeof(ele), :R1) && !is_identity_matrix(ele.R1)
        push!(params, "R1=$(format_matrix_python(ele.R1))")
    end
    if hasfield(typeof(ele), :T2) && !is_zero_vector(ele.T2)
        push!(params, "T2=$(format_vector_python(ele.T2))")
    end
    if hasfield(typeof(ele), :R2) && !is_identity_matrix(ele.R2)
        push!(params, "R2=$(format_matrix_python(ele.R2))")
    end
end

# Element-specific converters
function drift_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len)]
    add_misalignment_params(params, ele)
    return "$name = jt.DRIFT(" * join(params, ", ") * ")"
end

function quad_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("k1=%.10e", ele.k1)]
    add_misalignment_params(params, ele)
    return "$name = jt.QUAD(" * join(params, ", ") * ")"
end

function kquad_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("k1=%.10e", ele.k1)]
    add_misalignment_params(params, ele)
    return "$name = jt.KQUAD(" * join(params, ", ") * ")"
end

function ksext_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("k2=%.10e", ele.k2)]
    add_misalignment_params(params, ele)
    return "$name = jt.KSEXT(" * join(params, ", ") * ")"
end

function koct_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("k3=%.10e", ele.k3)]
    add_misalignment_params(params, ele)
    return "$name = jt.KOCT(" * join(params, ", ") * ")"
end

function sbend_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("%.10e", ele.angle)]
    
    # Add optional bend parameters
    if ele.e1 != 0.0
        push!(params, @sprintf("e1=%.10e", ele.e1))
    end
    if ele.e2 != 0.0
        push!(params, @sprintf("e2=%.10e", ele.e2))
    end
    
    # Check PolynomB for non-zero multipole components
    if hasfield(typeof(ele), :PolynomB) && any(x -> x != 0, ele.PolynomB)
        push!(params, "PolynomB=$(format_vector_python(ele.PolynomB))")
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.SBEND(" * join(params, ", ") * ")"
end

function rbend_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("%.10e", ele.angle)]
    
    # Add optional bend parameters
    if hasfield(typeof(ele), :e1) && ele.e1 != 0.0
        push!(params, @sprintf("e1=%.10e", ele.e1))
    end
    if hasfield(typeof(ele), :e2) && ele.e2 != 0.0
        push!(params, @sprintf("e2=%.10e", ele.e2))
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.RBEND(" * join(params, ", ") * ")"
end

function lbend_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("%.10e", ele.angle)]
    
    # Add optional bend parameters
    if hasfield(typeof(ele), :e1) && ele.e1 != 0.0
        push!(params, @sprintf("e1=%.10e", ele.e1))
    end
    if hasfield(typeof(ele), :e2) && ele.e2 != 0.0
        push!(params, @sprintf("e2=%.10e", ele.e2))
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.LBEND(" * join(params, ", ") * ")"
end

function esbend_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("%.10e", ele.angle)]
    
    # Add optional bend parameters
    if hasfield(typeof(ele), :e1) && ele.e1 != 0.0
        push!(params, @sprintf("e1=%.10e", ele.e1))
    end
    if hasfield(typeof(ele), :e2) && ele.e2 != 0.0
        push!(params, @sprintf("e2=%.10e", ele.e2))
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.ESBEND(" * join(params, ", ") * ")"
end

function erbend_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("%.10e", ele.angle)]
    
    # Add optional bend parameters
    if hasfield(typeof(ele), :e1) && ele.e1 != 0.0
        push!(params, @sprintf("e1=%.10e", ele.e1))
    end
    if hasfield(typeof(ele), :e2) && ele.e2 != 0.0
        push!(params, @sprintf("e2=%.10e", ele.e2))
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.ERBEND(" * join(params, ", ") * ")"
end

function solenoid_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("ks=%.10e", ele.ks)]
    add_misalignment_params(params, ele)
    return "$name = jt.SOLENOID(" * join(params, ", ") * ")"
end

function thinmultipole_to_python(ele, name::String)
    params = ["\"$name\""]
    
    # Add PolynomA and PolynomB if non-zero
    if hasfield(typeof(ele), :PolynomA) && any(x -> x != 0, ele.PolynomA)
        push!(params, "PolynomA=$(format_vector_python(ele.PolynomA))")
    end
    if hasfield(typeof(ele), :PolynomB) && any(x -> x != 0, ele.PolynomB)
        push!(params, "PolynomB=$(format_vector_python(ele.PolynomB))")
    end
    
    return "$name = jt.thinMULTIPOLE(" * join(params, ", ") * ")"
end

function hkicker_to_python(ele, name::String)
    params = ["\"$name\""]
    
    if hasfield(typeof(ele), :len) && ele.len != 0.0
        push!(params, @sprintf("%.10e", ele.len))
    end
    if hasfield(typeof(ele), :xkick) && ele.xkick != 0.0
        push!(params, @sprintf("xkick=%.10e", ele.xkick))
    end
    
    return "$name = jt.HKICKER(" * join(params, ", ") * ")"
end

function vkicker_to_python(ele, name::String)
    params = ["\"$name\""]
    
    if hasfield(typeof(ele), :len) && ele.len != 0.0
        push!(params, @sprintf("%.10e", ele.len))
    end
    if hasfield(typeof(ele), :ykick) && ele.ykick != 0.0
        push!(params, @sprintf("ykick=%.10e", ele.ykick))
    end
    
    return "$name = jt.VKICKER(" * join(params, ", ") * ")"
end

function wiggler_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("len=%.10e", ele.len)]
    
    # Add wiggler-specific parameters
    if hasfield(typeof(ele), :b0) && ele.b0 != 0.0
        push!(params, @sprintf("b0=%.10e", ele.b0))
    end
    if hasfield(typeof(ele), :period_len) && ele.period_len != 0.0
        push!(params, @sprintf("period_len=%.10e", ele.period_len))
    end
    if hasfield(typeof(ele), :nperiod)
        push!(params, "nperiod=$(ele.nperiod)")
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.WIGGLER(" * join(params, ", ") * ")"
end

function rfca_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("%.10e", ele.volt), @sprintf("%.10e", ele.freq)]
    
    # Add optional RF parameters
    if hasfield(typeof(ele), :h) && ele.h != 0
        push!(params, @sprintf("h=%.10e", ele.h))
    end
    if hasfield(typeof(ele), :energy) && ele.energy != 0
        push!(params, @sprintf("energy=%.10e", ele.energy))
    end
    
    return "$name = jt.RFCA(" * join(params, ", ") * ")"
end

function marker_to_python(ele, name::String)
    return "$name = jt.MARKER(\"$name\")"
end

function corrector_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len)]
    
    if hasfield(typeof(ele), :xkick) && ele.xkick != 0.0
        push!(params, @sprintf("xkick=%.10e", ele.xkick))
    end
    if hasfield(typeof(ele), :ykick) && ele.ykick != 0.0
        push!(params, @sprintf("ykick=%.10e", ele.ykick))
    end
    
    return "$name = jt.CORRECTOR(" * join(params, ", ") * ")"
end

# Space charge variants
function drift_sc_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len)]
    
    # Add space charge parameters
    if hasfield(typeof(ele), :a)
        push!(params, @sprintf("a=%.10e", ele.a))
    end
    if hasfield(typeof(ele), :b)
        push!(params, @sprintf("b=%.10e", ele.b))
    end
    if hasfield(typeof(ele), :Nl)
        push!(params, "Nl=$(ele.Nl)")
    end
    if hasfield(typeof(ele), :Nm)
        push!(params, "Nm=$(ele.Nm)")
    end
    if hasfield(typeof(ele), :Nsteps) && ele.Nsteps != 1
        push!(params, "Nsteps=$(ele.Nsteps)")
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.DRIFT_SC(" * join(params, ", ") * ")"
end

function quad_sc_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("k1=%.10e", ele.k1)]
    
    # Add space charge parameters
    if hasfield(typeof(ele), :a)
        push!(params, @sprintf("a=%.10e", ele.a))
    end
    if hasfield(typeof(ele), :b)
        push!(params, @sprintf("b=%.10e", ele.b))
    end
    if hasfield(typeof(ele), :Nl)
        push!(params, "Nl=$(ele.Nl)")
    end
    if hasfield(typeof(ele), :Nm)
        push!(params, "Nm=$(ele.Nm)")
    end
    if hasfield(typeof(ele), :Nsteps) && ele.Nsteps != 1
        push!(params, "Nsteps=$(ele.Nsteps)")
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.QUAD_SC(" * join(params, ", ") * ")"
end

function kquad_sc_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("k1=%.10e", ele.k1)]
    
    # Add space charge parameters
    if hasfield(typeof(ele), :a)
        push!(params, @sprintf("a=%.10e", ele.a))
    end
    if hasfield(typeof(ele), :b)
        push!(params, @sprintf("b=%.10e", ele.b))
    end
    if hasfield(typeof(ele), :Nl)
        push!(params, "Nl=$(ele.Nl)")
    end
    if hasfield(typeof(ele), :Nm)
        push!(params, "Nm=$(ele.Nm)")
    end
    if hasfield(typeof(ele), :Nsteps) && ele.Nsteps != 1
        push!(params, "Nsteps=$(ele.Nsteps)")
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.KQUAD_SC(" * join(params, ", ") * ")"
end

function ksext_sc_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("k2=%.10e", ele.k2)]
    
    # Add space charge parameters
    if hasfield(typeof(ele), :a)
        push!(params, @sprintf("a=%.10e", ele.a))
    end
    if hasfield(typeof(ele), :b)
        push!(params, @sprintf("b=%.10e", ele.b))
    end
    if hasfield(typeof(ele), :Nl)
        push!(params, "Nl=$(ele.Nl)")
    end
    if hasfield(typeof(ele), :Nm)
        push!(params, "Nm=$(ele.Nm)")
    end
    if hasfield(typeof(ele), :Nsteps) && ele.Nsteps != 1
        push!(params, "Nsteps=$(ele.Nsteps)")
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.KSEXT_SC(" * join(params, ", ") * ")"
end

function koct_sc_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("k3=%.10e", ele.k3)]
    
    # Add space charge parameters
    if hasfield(typeof(ele), :a)
        push!(params, @sprintf("a=%.10e", ele.a))
    end
    if hasfield(typeof(ele), :b)
        push!(params, @sprintf("b=%.10e", ele.b))
    end
    if hasfield(typeof(ele), :Nl)
        push!(params, "Nl=$(ele.Nl)")
    end
    if hasfield(typeof(ele), :Nm)
        push!(params, "Nm=$(ele.Nm)")
    end
    if hasfield(typeof(ele), :Nsteps) && ele.Nsteps != 1
        push!(params, "Nsteps=$(ele.Nsteps)")
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.KOCT_SC(" * join(params, ", ") * ")"
end

function sbend_sc_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("%.10e", ele.angle)]
    
    # Add optional bend parameters
    if ele.e1 != 0.0
        push!(params, @sprintf("e1=%.10e", ele.e1))
    end
    if ele.e2 != 0.0
        push!(params, @sprintf("e2=%.10e", ele.e2))
    end
    
    # Add space charge parameters
    if hasfield(typeof(ele), :a)
        push!(params, @sprintf("a=%.10e", ele.a))
    end
    if hasfield(typeof(ele), :b)
        push!(params, @sprintf("b=%.10e", ele.b))
    end
    if hasfield(typeof(ele), :Nl)
        push!(params, "Nl=$(ele.Nl)")
    end
    if hasfield(typeof(ele), :Nm)
        push!(params, "Nm=$(ele.Nm)")
    end
    if hasfield(typeof(ele), :Nsteps) && ele.Nsteps != 1
        push!(params, "Nsteps=$(ele.Nsteps)")
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.SBEND_SC(" * join(params, ", ") * ")"
end

function rbend_sc_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("%.10e", ele.angle)]
    
    # Add optional bend parameters
    if hasfield(typeof(ele), :e1) && ele.e1 != 0.0
        push!(params, @sprintf("e1=%.10e", ele.e1))
    end
    if hasfield(typeof(ele), :e2) && ele.e2 != 0.0
        push!(params, @sprintf("e2=%.10e", ele.e2))
    end
    
    # Add space charge parameters
    if hasfield(typeof(ele), :a)
        push!(params, @sprintf("a=%.10e", ele.a))
    end
    if hasfield(typeof(ele), :b)
        push!(params, @sprintf("b=%.10e", ele.b))
    end
    if hasfield(typeof(ele), :Nl)
        push!(params, "Nl=$(ele.Nl)")
    end
    if hasfield(typeof(ele), :Nm)
        push!(params, "Nm=$(ele.Nm)")
    end
    if hasfield(typeof(ele), :Nsteps) && ele.Nsteps != 1
        push!(params, "Nsteps=$(ele.Nsteps)")
    end
    
    add_misalignment_params(params, ele)
    return "$name = jt.RBEND_SC(" * join(params, ", ") * ")"
end

function add_sc2p5d_params(params::Vector{String}, ele)
    if hasfield(typeof(ele), :xsize)
        push!(params, "xsize=$(ele.xsize)")
    end
    if hasfield(typeof(ele), :ysize)
        push!(params, "ysize=$(ele.ysize)")
    end
    if hasfield(typeof(ele), :zsize)
        push!(params, "zsize=$(ele.zsize)")
    end
    if hasfield(typeof(ele), :pipe_radius)
        push!(params, @sprintf("pipe_radius=%.10e", ele.pipe_radius))
    end
    if hasfield(typeof(ele), :xy_ratio) && ele.xy_ratio != 1.0
        push!(params, @sprintf("xy_ratio=%.10e", ele.xy_ratio))
    end
    if hasfield(typeof(ele), :long_avg_n) && ele.long_avg_n != 3
        push!(params, "long_avg_n=$(ele.long_avg_n)")
    end
    if hasfield(typeof(ele), :Nsteps) && ele.Nsteps != 1
        push!(params, "Nsteps=$(ele.Nsteps)")
    end
end

function drift_sc2p5d_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len)]
    add_sc2p5d_params(params, ele)
    add_misalignment_params(params, ele)
    return "$name = jt.DRIFT_SC2P5D(" * join(params, ", ") * ")"
end

function quad_sc2p5d_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("k1=%.10e", ele.k1)]
    add_sc2p5d_params(params, ele)
    add_misalignment_params(params, ele)
    return "$name = jt.QUAD_SC2P5D(" * join(params, ", ") * ")"
end

function kquad_sc2p5d_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("k1=%.10e", ele.k1)]
    if hasfield(typeof(ele), :NumIntSteps) && ele.NumIntSteps != 10
        push!(params, "NumIntSteps=$(ele.NumIntSteps)")
    end
    add_sc2p5d_params(params, ele)
    add_misalignment_params(params, ele)
    return "$name = jt.KQUAD_SC2P5D(" * join(params, ", ") * ")"
end

function ksext_sc2p5d_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("k2=%.10e", ele.k2)]
    if hasfield(typeof(ele), :NumIntSteps) && ele.NumIntSteps != 10
        push!(params, "NumIntSteps=$(ele.NumIntSteps)")
    end
    add_sc2p5d_params(params, ele)
    add_misalignment_params(params, ele)
    return "$name = jt.KSEXT_SC2P5D(" * join(params, ", ") * ")"
end

function koct_sc2p5d_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("k3=%.10e", ele.k3)]
    if hasfield(typeof(ele), :NumIntSteps) && ele.NumIntSteps != 10
        push!(params, "NumIntSteps=$(ele.NumIntSteps)")
    end
    add_sc2p5d_params(params, ele)
    add_misalignment_params(params, ele)
    return "$name = jt.KOCT_SC2P5D(" * join(params, ", ") * ")"
end

function sbend_sc2p5d_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("%.10e", ele.angle)]
    if hasfield(typeof(ele), :e1) && ele.e1 != 0.0
        push!(params, @sprintf("e1=%.10e", ele.e1))
    end
    if hasfield(typeof(ele), :e2) && ele.e2 != 0.0
        push!(params, @sprintf("e2=%.10e", ele.e2))
    end
    if hasfield(typeof(ele), :NumIntSteps) && ele.NumIntSteps != 10
        push!(params, "NumIntSteps=$(ele.NumIntSteps)")
    end
    add_sc2p5d_params(params, ele)
    add_misalignment_params(params, ele)
    return "$name = jt.SBEND_SC2P5D(" * join(params, ", ") * ")"
end

function rbend_sc2p5d_to_python(ele, name::String)
    params = ["\"$name\"", @sprintf("%.10e", ele.len), @sprintf("%.10e", ele.angle)]
    if hasfield(typeof(ele), :e1) && ele.e1 != 0.0
        push!(params, @sprintf("e1=%.10e", ele.e1))
    end
    if hasfield(typeof(ele), :e2) && ele.e2 != 0.0
        push!(params, @sprintf("e2=%.10e", ele.e2))
    end
    if hasfield(typeof(ele), :NumIntSteps) && ele.NumIntSteps != 10
        push!(params, "NumIntSteps=$(ele.NumIntSteps)")
    end
    add_sc2p5d_params(params, ele)
    add_misalignment_params(params, ele)
    return "$name = jt.RBEND_SC2P5D(" * join(params, ", ") * ")"
end

function spacecharge2p5d_to_python(ele, name::String)
    params = ["\"$name\""]

    if hasfield(typeof(ele), :len) && ele.len != 0.0
        push!(params, @sprintf("length=%.10e", ele.len))
    end
    if hasfield(typeof(ele), :effective_len) && ele.effective_len != 0.0
        push!(params, @sprintf("effective_length=%.10e", ele.effective_len))
    end
    if hasfield(typeof(ele), :xsize)
        push!(params, "xsize=$(ele.xsize)")
    end
    if hasfield(typeof(ele), :ysize)
        push!(params, "ysize=$(ele.ysize)")
    end
    if hasfield(typeof(ele), :zsize)
        push!(params, "zsize=$(ele.zsize)")
    end
    if hasfield(typeof(ele), :pipe_radius)
        push!(params, @sprintf("pipe_radius=%.10e", ele.pipe_radius))
    end
    if hasfield(typeof(ele), :xy_ratio) && ele.xy_ratio != 1.0
        push!(params, @sprintf("xy_ratio=%.10e", ele.xy_ratio))
    end
    if hasfield(typeof(ele), :long_avg_n) && ele.long_avg_n != 3
        push!(params, "long_avg_n=$(ele.long_avg_n)")
    end

    return "$name = jt.SPACECHARGE2P5D(" * join(params, ", ") * ")"
end
