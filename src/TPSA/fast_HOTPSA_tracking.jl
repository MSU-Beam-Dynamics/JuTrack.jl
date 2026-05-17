function _strip(param::T) where {T<:AbstractHOTPSA}
    return jet_primal(param)
end

function _strip(param::Vector{T}) where {T<:AbstractHOTPSA}
    return [_strip(p) for p in param]
end

function _strip(param::Matrix{T}) where {T<:AbstractHOTPSA}
    Mat = zeros(Float64, size(param))
    for i in eachindex(param)
        Mat[i] = _strip(param[i])
    end
    return Mat
end

function to_Number(te::AbstractElement{T}) where {T<:AbstractHOTPSA}
    base_type = typeof(te).name.wrapper
    vals = map(f -> _strip(getfield(te, f)), fieldnames(typeof(te)))
    return base_type(vals...)
end

function TPSAD2Number(line::Vector{<:AbstractElement{T}}) where {T<:AbstractHOTPSA}
    return [to_Number(ele) for ele in line]
end

function _strip_inverse(param::Float64, fieldname::Symbol, ::Type{S}) where {S<:AbstractHOTPSA}
    if fieldname === :green_dx || fieldname === :green_dy
        return param
    end
    return convert(S, param)
end

function _strip_inverse(param::Vector{Float64}, fieldname::Symbol, ::Type{S}) where {S<:AbstractHOTPSA}
    if fieldname == :RApertures || fieldname == :EApertures ||
       fieldname == :z_grid || fieldname == :z_deriv_grid
        return param
    end
    return [convert(S, p) for p in param]
end

function _strip_inverse(param::Matrix{Float64}, fieldname::Symbol, ::Type{S}) where {S<:AbstractHOTPSA}
    if fieldname == :rho_grid || fieldname == :phi_grid
        return param
    end
    Mat = zeros(S, size(param)...)
    for i in eachindex(param)
        Mat[i] = convert(S, param[i])
    end
    return Mat
end

function _strip_inverse(param::Union{String, Int64, Int32, Vector{Int64}, Vector{Int32}, Matrix{Int64}, Matrix{Int32}}, fieldname::Symbol, ::Type{S}) where {S<:AbstractHOTPSA}
    return param
end

function _strip_inverse(param, fieldname::Symbol, ::Type{S}) where {S<:AbstractHOTPSA}
    return param
end

function to_TPSAD(ele::AbstractElement{Float64}, ::Type{S}) where {S<:AbstractHOTPSA}
    base_type = typeof(ele).name.wrapper
    vals = map(f -> _strip_inverse(getfield(ele, f), f, S), fieldnames(typeof(ele)))
    return base_type(vals...)
end

function Number2TPSAD(line::Vector{<:AbstractElement{Float64}}, ::Type{S}) where {S<:AbstractHOTPSA}
    return [to_TPSAD(ele, S) for ele in line]
end
