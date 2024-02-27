using JSON

function generate_line_code(line_elements, max_line_length=20)
    line_code = "line = ["
    current_line_length = length(line_code)

    for (index, element) in enumerate(line_elements)
        element_code = "\"$element\""
        if index < length(line_elements)
            element_code *= ", "
        end
        
        if current_line_length + length(element_code) > max_line_length
            line_code *= "\n"  # Start a new line
            current_line_length = 0  # Reset current line length
        end
        
        line_code *= element_code
        current_line_length += length(element_code)
    end

    line_code *= "]"
    return line_code
end

function generate_julia_from_json(json_filename)
    elements = JSON.parsefile(json_filename)

    # Initialize strings to store the generated code and line definition
    generated_code = ""
    line_elements = []

    # Iterate over each element dictionary in the list
    for elem in elements
        name = elem["name"]
        eletype = elem["eletype"]
        len = elem["len"]
        
        if eletype == "KQUAD"
            k1 = elem["k1"]
            element_code = "$(name) = $(eletype)(name=\"$(name)\", len=$(len), k1=$(k1))\n"
        elseif eletype == "DRIFT"
            element_code = "$(name) = $(eletype)(name=\"$(name)\", len=$(len))\n"
        elseif eletype == "MARKER"
            element_code = "$(name) = $(eletype)(name=\"$(name)\")\n"
        elseif eletype == "KSEXT"
            k2 = elem["k2"]
            element_code = "$(name) = $(eletype)(name=\"$(name)\", len=$(len), k2=$(k2))\n"
        elseif eletype == "KOCT"
            k3 = elem["k3"]
            element_code = "$(name) = $(eletype)(name=\"$(name)\", len=$(len), k3=$(k3))\n"
        elseif eletype == "SBEND"
            angle = elem["angle"]
            e1 = elem["e1"]
            e2 = elem["e2"]
            element_code = "$(name) = $(eletype)(name=\"$(name)\", len=$(len), angle=$(angle), e1=$(e1), e2=$(e2))\n"
        elseif eletype == "RBEND"
            angle = elem["angle"]
            element_code = "$(name) = $(eletype)(name=\"$(name)\", len=$(len), angle=$(angle))\n"
        elseif eletype == "RFCA"
            volt = elem["volt"]
            freq = elem["freq"]
            h = elem["h"]
            lag = elem["lag"]
            philag = elem["philag"]
            energy = elem["energy"]
            element_code = "$(name) = $(eletype)(name=\"$(name)\", len=$(len), volt=$(volt), freq=$(freq), h=$(h), lag=$(lag), philag=$(philag), energy=$(energy))\n"
        elseif eletype == "SOLENOID"
            ks = elem["ks"]
            element_code = "$(name) = $(eletype)(name=\"$(name)\", len=$(len), ks=$(ks))\n"
        elseif eletype == "CORRECTOR" || eletype == "HKICKER" || eletype == "VKICKER"
            xkick = elem["xkick"]
            ykick = elem["ykick"]
            element_code = "$(name) = $(eletype)(name=\"$(name)\", len=$(len), xkick=$(xkick), ykick=$(ykick))\n"
        elseif eletype == "CrabCavity"
            volt = elem["volt"]
            freq = elem["freq"]
            k = elem["k"]
            phi = elem["phi"]
            element_code = "$(name) = $(eletype)(name=\"$(name)\", len=$(len), volt=$(volt), freq=$(freq), k=$(k), phi=$(phi))\n"
        elseif eletype == "thinMULTIPOLE"
            PolynomB = elem["PolynomB"]
            element_code = "$(name) = $(eletype)(name=\"$(name)\", len=$(len), PolynomB=$(PolynomB))\n"
        else
            error("Element type not recognized: $(eletype)")
        end
        
        generated_code *= element_code
        push!(line_elements, name)
    end
    
    # Generate the line definition code
    # line_code = generate_line_code(line_elements)

    line_code = "line = [" * join(line_elements, ", ") * "]"
    
    # Combine all generated code
    full_generated_code = generated_code * "\n" * line_code
    
    return full_generated_code
end

# Use the function with your JSON filename
json_filename = "ESR_main.json"
julia_code = generate_julia_from_json(json_filename)

# Display the generated Julia code
println(julia_code)

# Optionally, save the generated Julia code to a file
open("generated_elements.jl", "w") do file
    write(file, julia_code)
end
