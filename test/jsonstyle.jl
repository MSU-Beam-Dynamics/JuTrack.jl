using JSON
using JuTrack

function serialize_element(element::AbstractElement)
    attributes = Dict()
    
    for field in fieldnames(typeof(element))
        val = getfield(element, field)
        attributes[string(field)] = val
    end
    
    element_dict = Dict("Type" => typeof(element).name.name, "Attributes" => attributes)
    return JSON.json(element_dict)
end

function save_elements(elements::Vector{AbstractElement}, filename::String)
    json_elements = "[" * join([serialize_element(el) for el in elements], ",\n") * "]"
    
    open(filename, "w") do file
        write(file, json_elements)
    end
end
elements = [MARKER(), DRIFT(), KQUAD(), KSEXT()]
save_elements(elements, "test/lattice_json_example.json")
