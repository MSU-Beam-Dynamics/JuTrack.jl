include("../TPSA/TPSA.jl")
using LinearAlgebra

function TransferMap(ele::Drift, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    len = ele.len
    x_new = r[1]  + r[2]*len/(1 + r[5])
    px_new = r[2]
    y_new = r[3] + r[4]*len/(1 + r[5])
    py_new = r[4]
    delta_new = r[5]
    z_new = r[6] + len*(r[2]^2 + r[4]^2)/2

    return [x_new, px_new, y_new, py_new, delta_new, z_new]
end

function TransferMap(ele::Quad, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # the same as QuadLinearPass from AT for defocusing-quad
    # map = Matrix{Float64}(I,6,6)
    len = ele.len
    k = ele.k1
    if k>0
        g = k/(1 + r[5])
        t = sqrt(g)
        phi = t*len
        map = [cos(phi) sin(phi)/t 0 0 0 0;
            -sin(phi)/t*g cos(phi) 0 0 0 0;
            0 0 cosh(phi) sinh(phi)/t 0 0;
            0 0 sinh(phi)/t*g cosh(phi) 0 0;
            0 0 0 0 1.0 0;
            0 0 0 0 0 1.0]
    else
        g = -k/(1 + r[5])
        t = sqrt(g)
        phi = t*len
        map = [cosh(phi) sinh(phi)/t 0 0 0 0;
            sinh(phi)/t*g cosh(phi) 0 0 0 0;
            0 0 cos(phi) sin(phi)/t 0 0;
            0 0 -sin(phi)/t*g cos(phi) 0 0;
            0 0 0 0 1.0 0;
            0 0 0 0 0 1.0]
    end
    x = r[1]
    xpr = r[2]/(1 + r[5])
    y = r[3]
    ypr = r[4]/(1 + r[5])

    x_new = map[1,1]*x + map[1,2]*xpr
    px_new = (map[2,1]*x + map[2,2]*xpr)*(1 + r[5])
    y_new = map[3,3]*y + map[3,4]*ypr
    py_new = (map[4,3]*y + map[4,4]*ypr)*(1 + r[5])
    delta_new = r[5]
    z_new = r[6] + g*(x*x*(len-map[1,1]*map[1,2]) - y*y*(len - map[3,3]*map[3,4]))/4 +
                         (xpr*xpr*(len + map[1,1]*map[1,2]) + ypr*ypr*(len + map[3,3]*map[3,4]))/4 +
                         (x*xpr*map[1,2]*map[2,1] + y*ypr*map[3,4]*map[4,3])/2

    return [x_new, px_new, y_new, py_new, delta_new, z_new]
end

function TransferMap(ele::Bend, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # the same as BendLinearPass from AT. ignore fringe
    len = ele.len
    angle = ele.angle
    e1 = ele.e1
    e2 = ele.e2
    grd = ele.grad

    Kx = angle/len
    G1 = (Kx*Kx + grd)/(1 + r[5]) # G1 must > 0
    G2 = -grd/(1 + r[5]) # G2 must = 0
    arg1 = sqrt(G1) * len
    ByError = 0.0 # ignore By error
    # Build the map matrix without mutating operations
    map = [cos(arg1) sin(arg1)/sqrt(G1) 0 0 0 0;
    -sin(arg1)*sqrt(G1) cos(arg1) 0 0 0 0;
    0 0 1.0 len 0 0;
    0 0 0 1.0 0 0;
    0 0 0 0 1.0 0;
    0 0 0 0 0 1.0]

    x = r[1]
    xpr = r[2]/(1 + r[5])
    y = r[3]
    ypr = r[4]/(1 + r[5])

    x_new = map[1,1]*x + map[1,2]*xpr
    px_new = (map[2,1]*x + map[2,2]*xpr)*(1 + r[5])
    y_new = map[3,3]*y + map[3,4]*ypr
    py_new = (map[4,3]*y + map[4,4]*ypr)*(1 + r[5])
    delta_new = r[5]
    z_new = r[6] + xpr*xpr*(len + map[1,1]*map[1,2])/4
    z_new = z_new + (L-map[1,1]*map[1,2])*(x*x*G1+(r[5]/(1+r[5])-ByError)*(r[5]/(1+r[5])-ByError)*Kx*Kx/G1-2*x*Kx*(r[5]/(1+r[5])-ByError))/4
    z_new = z_new + map[1,2]*map[2,1]*( x*xpr - xpr*(r[5]/(1+r[5])-ByError)*Kx/G1)/2
    z_new = z_new + Kx*x*map[1,2] + xpr*(1-map[1,1])*Kx/G1 + (r[5]/(1+r[5])-ByError)*(len-map[1,2])*Kx*Kx/G1
    z_new = z_new + ((len-map[3,3]*map[3,4])*y*y*G2 + ypr*ypr*(len+map[3,3]*map[3,4]))/4
    z_new = z_new + map[3,4]*map[4,3]*x*xpr/2

    return [x_new, px_new, y_new, py_new, delta_new, z_new]
end


function track(cell::Vector{AbstractElement}, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    r1 = copy(r)
    for ele in cell
        r1 = TransferMap(ele, r1)
    end
    return r1
end

# function linepass(cell, rin)
#     np = length(rin[1,:])
#     r_buff = Zygote.Buffer(rin)
#     for i in 1:np
#         rin_m = [rin[1,i], rin[2,i]* (1 + rin[5]), rin[3,i], rin[4,i]* (1 + rin[5]), rin[5,i], rin[6,i]]
#         rout = track(cell, rin_m)
#         r_buff[:,i] = rout
#     end
#     return copy(r_buff)
# end

# function track_through_refpt(cell, r::Vector{Float64}, refpts::Vector{Int})
#     r1 = copy(r)
#     results = Zygote.Buffer(r, length(r), length(refpts))

#     ref_index = 1
#     current_refpt = refpts[ref_index]

#     # If the first reference point is 1, store the initial state
#     if current_refpt == 1
#         results[:, ref_index] = r1
#         ref_index += 1
#         if ref_index <= length(refpts)
#             current_refpt = refpts[ref_index]
#         else
#             return copy(results)
#         end
#     end

#     for (i, ele) in enumerate(cell)
#         r1, _ = TransferMap(ele, r1)

#         # Check if current element is a reference point
#         if i == current_refpt - 1
#             results[:, ref_index] = r1

#             ref_index += 1
#             if ref_index <= length(refpts)
#                 current_refpt = refpts[ref_index]
#             else
#                 break
#             end
#         end
#     end

#     return copy(results)
# end

# function linepass1(cell, rin, refpts=nothing)
#     np = size(rin, 2)
#     r_buff = Zygote.Buffer(rin, size(rin, 1), np * (refpts === nothing ? 1 : length(refpts)))


#     for i in 1:np
#         rin_m = [rin[1,i], rin[2,i]* (1 + rin[5,i]), rin[3,i], rin[4,i]* (1 + rin[5,i]), rin[5,i], rin[6,i]]

#         # If refpts is specified, store the state at each refpt.
#         if refpts !== nothing
#             rout = track_through_refpt(cell, rin_m, refpts)
#             r_buff[:, (i-1)*size(rout, 2) .+ (1:size(rout, 2))] = rout
#         else
#             rout = track(cell, rin_m)
#             r_buff[:,i] = rout
#         end
#     end

#     return copy(r_buff)
# end



