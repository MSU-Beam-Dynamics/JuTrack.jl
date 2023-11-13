include("../TPSA/TPSA.jl")
using LinearAlgebra

# function TransferMap(ele::Drift, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
#     len = ele.len
#     x_new = r[1]  + r[2]*len/(1 + r[5])
#     px_new = r[2]
#     y_new = r[3] + r[4]*len/(1 + r[5])
#     py_new = r[4]
#     delta_new = r[5]
#     z_new = r[6] + len*(r[2]^2 + r[4]^2)/2

#     return [x_new, px_new, y_new, py_new, delta_new, z_new]
# end

# function TransferMap(ele::Quad, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
#     # the same as QuadLinearPass from AT for defocusing-quad
#     # map = Matrix{Float64}(I,6,6)
#     len = ele.len
#     k = ele.k1
#     if k>0
#         g = k/(1 + r[5])
#         t = sqrt(g)
#         phi = t*len
#         map = [cos(phi) sin(phi)/t 0 0 0 0;
#             -sin(phi)/t*g cos(phi) 0 0 0 0;
#             0 0 cosh(phi) sinh(phi)/t 0 0;
#             0 0 sinh(phi)/t*g cosh(phi) 0 0;
#             0 0 0 0 1.0 0;
#             0 0 0 0 0 1.0]
#     else
#         g = -k/(1 + r[5])
#         t = sqrt(g)
#         phi = t*len
#         map = [cosh(phi) sinh(phi)/t 0 0 0 0;
#             sinh(phi)/t*g cosh(phi) 0 0 0 0;
#             0 0 cos(phi) sin(phi)/t 0 0;
#             0 0 -sin(phi)/t*g cos(phi) 0 0;
#             0 0 0 0 1.0 0;
#             0 0 0 0 0 1.0]
#     end
#     x = r[1]
#     xpr = r[2]/(1 + r[5])
#     y = r[3]
#     ypr = r[4]/(1 + r[5])

#     x_new = map[1,1]*x + map[1,2]*xpr
#     px_new = (map[2,1]*x + map[2,2]*xpr)*(1 + r[5])
#     y_new = map[3,3]*y + map[3,4]*ypr
#     py_new = (map[4,3]*y + map[4,4]*ypr)*(1 + r[5])
#     delta_new = r[5]
#     z_new = r[6] + g*(x*x*(len-map[1,1]*map[1,2]) - y*y*(len - map[3,3]*map[3,4]))/4 +
#                          (xpr*xpr*(len + map[1,1]*map[1,2]) + ypr*ypr*(len + map[3,3]*map[3,4]))/4 +
#                          (x*xpr*map[1,2]*map[2,1] + y*ypr*map[3,4]*map[4,3])/2

#     return [x_new, px_new, y_new, py_new, delta_new, z_new]
# end

# function TransferMap(ele::Bend, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
#     # the same as BendLinearPass from AT. ignore fringe
#     len = ele.len
#     angle = ele.angle
#     e1 = ele.e1
#     e2 = ele.e2
#     grd = ele.grad

#     Kx = angle/len
#     G1 = (Kx*Kx + grd)/(1 + r[5]) # G1 must > 0
#     G2 = -grd/(1 + r[5]) # G2 must = 0
#     arg1 = sqrt(G1) * len
#     ByError = 0.0 # ignore By error
#     # Build the map matrix without mutating operations
#     map = [cos(arg1) sin(arg1)/sqrt(G1) 0 0 0 0;
#     -sin(arg1)*sqrt(G1) cos(arg1) 0 0 0 0;
#     0 0 1.0 len 0 0;
#     0 0 0 1.0 0 0;
#     0 0 0 0 1.0 0;
#     0 0 0 0 0 1.0]

#     x = r[1]
#     xpr = r[2]/(1 + r[5])
#     y = r[3]
#     ypr = r[4]/(1 + r[5])

#     x_new = map[1,1]*x + map[1,2]*xpr
#     px_new = (map[2,1]*x + map[2,2]*xpr)*(1 + r[5])
#     y_new = map[3,3]*y + map[3,4]*ypr
#     py_new = (map[4,3]*y + map[4,4]*ypr)*(1 + r[5])
#     delta_new = r[5]
#     z_new = r[6] + xpr*xpr*(len + map[1,1]*map[1,2])/4
#     z_new = z_new + (L-map[1,1]*map[1,2])*(x*x*G1+(r[5]/(1+r[5])-ByError)*(r[5]/(1+r[5])-ByError)*Kx*Kx/G1-2*x*Kx*(r[5]/(1+r[5])-ByError))/4
#     z_new = z_new + map[1,2]*map[2,1]*( x*xpr - xpr*(r[5]/(1+r[5])-ByError)*Kx/G1)/2
#     z_new = z_new + Kx*x*map[1,2] + xpr*(1-map[1,1])*Kx/G1 + (r[5]/(1+r[5])-ByError)*(len-map[1,2])*Kx*Kx/G1
#     z_new = z_new + ((len-map[3,3]*map[3,4])*y*y*G2 + ypr*ypr*(len+map[3,3]*map[3,4]))/4
#     z_new = z_new + map[3,4]*map[4,3]*x*xpr/2

#     return [x_new, px_new, y_new, py_new, delta_new, z_new]
# end


# function track(cell::Vector{AbstractElement}, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
#     r1 = copy(r)
#     for ele in cell
#         r1 = TransferMap(ele, r1)
#     end
#     return r1
# end

function TransferMap(ele::Drift, x::CTPS{T, TPS_Dim, Max_TPS_Degree},px::CTPS{T, TPS_Dim, Max_TPS_Degree},
    y::CTPS{T, TPS_Dim, Max_TPS_Degree}, py::CTPS{T, TPS_Dim, Max_TPS_Degree}, delta::CTPS{T, TPS_Dim, Max_TPS_Degree},
    z::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    len = ele.len
    x_new = x + px*len/(1 + delta)
    px_new = px
    y_new = y + py*len/(1 + delta)
    py_new = py
    delta_new = delta
    z_new = z + len*(px^2 + py^2)/2

    return x_new, px_new, y_new, py_new, delta_new, z_new
end

function TransferMap(ele::Quad, x::CTPS{T, TPS_Dim, Max_TPS_Degree},px::CTPS{T, TPS_Dim, Max_TPS_Degree},
    y::CTPS{T, TPS_Dim, Max_TPS_Degree}, py::CTPS{T, TPS_Dim, Max_TPS_Degree}, delta::CTPS{T, TPS_Dim, Max_TPS_Degree},
    z::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    # the same as QuadLinearPass from AT for defocusing-quad
    # map = Matrix{Float64}(I,6,6)
    len = ele.len
    k = ele.k1
    if k>0
        g = k/(1 + delta)
        t = sqrt(g)
        phi = t*len
        map = [cos(phi) sin(phi)/t 0 0 0 0;
            -sin(phi)/t*g cos(phi) 0 0 0 0;
            0 0 cosh(phi) sinh(phi)/t 0 0;
            0 0 sinh(phi)/t*g cosh(phi) 0 0;
            0 0 0 0 1.0 0;
            0 0 0 0 0 1.0]
    else
        g = -k/(1 + delta)
        t = sqrt(g)
        phi = t*len
        map = [cosh(phi) sinh(phi)/t 0 0 0 0;
            sinh(phi)/t*g cosh(phi) 0 0 0 0;
            0 0 cos(phi) sin(phi)/t 0 0;
            0 0 -sin(phi)/t*g cos(phi) 0 0;
            0 0 0 0 1.0 0;
            0 0 0 0 0 1.0]
    end
    xpr = px/(1 + delta)
    ypr = py/(1 + delta)

    x_new = map[1,1]*x + map[1,2]*xpr
    px_new = (map[2,1]*x + map[2,2]*xpr)*(1 + delta)
    y_new = map[3,3]*y + map[3,4]*ypr
    py_new = (map[4,3]*y + map[4,4]*ypr)*(1 + delta)
    delta_new = delta
    z_new = z + g*(x*x*(len-map[1,1]*map[1,2]) - y*y*(len - map[3,3]*map[3,4]))/4 +
                         (xpr*xpr*(len + map[1,1]*map[1,2]) + ypr*ypr*(len + map[3,3]*map[3,4]))/4 +
                         (x*xpr*map[1,2]*map[2,1] + y*ypr*map[3,4]*map[4,3])/2

    return x_new, px_new, y_new, py_new, delta_new, z_new
end

function TransferMap(ele::Bend, x::CTPS{T, TPS_Dim, Max_TPS_Degree},px::CTPS{T, TPS_Dim, Max_TPS_Degree},
    y::CTPS{T, TPS_Dim, Max_TPS_Degree}, py::CTPS{T, TPS_Dim, Max_TPS_Degree}, delta::CTPS{T, TPS_Dim, Max_TPS_Degree},
    z::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    # the same as BendLinearPass from AT. ignore fringe
    len = ele.len
    angle = ele.angle
    e1 = ele.e1
    e2 = ele.e2
    grd = ele.grad

    Kx = angle/len
    G1 = (Kx*Kx + grd)/(1 + delta) # G1 must > 0
    G2 = -grd/(1 +delta) # G2 must = 0
    arg1 = sqrt(G1) * len
    ByError = 0.0 # ignore By error
    # Build the map matrix without mutating operations
    map = [cos(arg1) sin(arg1)/sqrt(G1) 0 0 0 0;
    -sin(arg1)*sqrt(G1) cos(arg1) 0 0 0 0;
    0 0 1.0 len 0 0;
    0 0 0 1.0 0 0;
    0 0 0 0 1.0 0;
    0 0 0 0 0 1.0]

    xpr = px/(1 + delta)
    ypr = py/(1 + delta)

    x_new = map[1,1]*x + map[1,2]*xpr
    px_new = (map[2,1]*x + map[2,2]*xpr)*(1 + delta)
    y_new = map[3,3]*y + map[3,4]*ypr
    py_new = (map[4,3]*y + map[4,4]*ypr)*(1 + delta)
    delta_new = delta
    z_new = z + xpr*xpr*(len + map[1,1]*map[1,2])/4
    z_new = z_new + (L-map[1,1]*map[1,2])*(x*x*G1+(delta/(1+delta)-ByError)*(delta/(1+delta)-ByError)*Kx*Kx/G1-2*x*Kx*(rdelta/(1+delta)-ByError))/4
    z_new = z_new + map[1,2]*map[2,1]*( x*xpr - xpr*(delta/(1+delta)-ByError)*Kx/G1)/2
    z_new = z_new + Kx*x*map[1,2] + xpr*(1-map[1,1])*Kx/G1 + (delta/(1+delta)-ByError)*(len-map[1,2])*Kx*Kx/G1
    z_new = z_new + ((len-map[3,3]*map[3,4])*y*y*G2 + ypr*ypr*(len+map[3,3]*map[3,4]))/4
    z_new = z_new + map[3,4]*map[4,3]*x*xpr/2

    return x_new, px_new, y_new, py_new, delta_new, z_new
end


function track(cell::Vector{AbstractElement}, x::CTPS{T, TPS_Dim, Max_TPS_Degree},px::CTPS{T, TPS_Dim, Max_TPS_Degree},
    y::CTPS{T, TPS_Dim, Max_TPS_Degree}, py::CTPS{T, TPS_Dim, Max_TPS_Degree}, delta::CTPS{T, TPS_Dim, Max_TPS_Degree},
    z::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    for ele in cell
        x, px, y, py, delta, z = TransferMap(ele, x, px, y, py, delta, z)
    end
    return x, px, y, py, delta, z
end


