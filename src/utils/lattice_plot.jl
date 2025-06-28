using PyCall
function elem_to_py(e::AbstractElement)
    jd = Dict(
      "type"   => e.eletype,
      "length" => e.len
    )
    if hasproperty(e, :angle)
        jd["angle"] = getfield(e, :angle)
        e1 = getfield(e, :e1)
        e2 = getfield(e, :e2)
        if abs(e1 - jd["angle"]/2.0) < 1e-10 && abs(e2 - jd["angle"]/2.0) < 1e-10
            jd["type"] = "RBEND"
        end
        jd["e1"] = e1
        jd["e2"] = e2
    end
    if hasproperty(e, :k1)
        jd["k1"] = getfield(e, :k1)
    end
    return PyDict(jd)
end

function __init__()
py"""
import numpy as np
try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Wedge, Arc
    from matplotlib.lines import Line2D
    from matplotlib.transforms import Affine2D
except ImportError:
    raise ImportError("Please install matplotlib associated with PyCall to use the lattice plot function.")
# import matplotlib.pyplot as plt
# from matplotlib.patches import Rectangle, Wedge, Arc, Polygon
# from matplotlib.lines import Line2D
# from matplotlib.transforms import Affine2D

def plot_lattice(lattice, width=0.25, axis=True):
    fig, ax = plt.subplots(figsize=(8, 8))
    pos   = np.array([0.0, 0.0])   # beam start
    theta = 0.0                    # heading (rad)

    for elem in lattice:
        L  = elem.get('length', 0.0)
        et = elem.get('type', '').upper()

        trans = Affine2D().rotate(theta).translate(*pos) + ax.transData

        if et in ('KQUAD', 'QUAD'):
            color = 'lightblue'
            if elem.get('k1', 0.0) >= 0:
                quad = Rectangle((0, 0), L, width,
                                 edgecolor='black', facecolor=color,
                                 transform=trans)
            else:
                quad = Rectangle((0, -width), L, width,
                                 edgecolor='black', facecolor=color,
                                 transform=trans)

            ax.add_patch(quad)
            beam = Line2D([0, L], [0, 0], transform=trans,
                          color='black', linewidth=1)
            ax.add_line(beam)

            pos += L * np.array([np.cos(theta), np.sin(theta)])

        elif et in ('KSEXT'):
            color = 'blue' 
            sext = Rectangle((0, -width*0.7/2), L, width*0.7,
                             edgecolor='black', facecolor=color, transform=trans)
            ax.add_patch(sext)
            beam = Line2D([0, L], [0, 0], transform=trans,
                          color='black', linewidth=1)
            ax.add_line(beam)
            pos += L * np.array([np.cos(theta), np.sin(theta)])
        elif et in ('KOCT'):
            color = 'green'
            octu = Rectangle((0, -width*0.4/2), L, width*0.4,
                             edgecolor='black', facecolor=color, transform=trans)
            ax.add_patch(octu)
            beam = Line2D([0, L], [0, 0], transform=trans,
                          color='black', linewidth=1)
            ax.add_line(beam)
            pos += L * np.array([np.cos(theta), np.sin(theta)])

        elif et == 'RBEND':
            φ  = elem.get('angle', 0.0)
            e1 = elem.get('e1',    0.0)
            e2 = elem.get('e2',    0.0)

            θ_face = theta + φ/2
            body_trans = (
                Affine2D()
                .rotate(θ_face)
                .translate(*pos)
                + ax.transData
            )
            body = Rectangle((0, -width/2), L, width,
                            edgecolor='black', facecolor='orange',
                            transform=body_trans)
            ax.add_patch(body)

            if φ != 0:
                R = L / (2 * np.sin(φ/2))

                center = pos + R * np.array([
                    -np.sin(theta),
                    np.cos(theta)
                ])

                start_deg = np.degrees(theta) - 90 + np.degrees(0)
                end_deg   = np.degrees(theta) - 90 + np.degrees(φ - 0)
                if start_deg > end_deg:
                    start_deg, end_deg = end_deg, start_deg
                beam_arc = Arc(center, 2*R, 2*R,
                            theta1=start_deg,
                            theta2=end_deg,
                            edgecolor='black',
                            linewidth=1)
                ax.add_patch(beam_arc)

                theta += φ
                pos = center + R * np.array([
                    np.sin(theta),
                -np.cos(theta)
                ])

            else:
                pos += L * np.array([np.cos(theta),
                                    np.sin(theta)])


        elif et == 'SBEND':
            phi = elem.get('angle', 0.0)
            e1  = elem.get('e1',    0.0)
            e2  = elem.get('e2',    0.0)

            if phi != 0:
                R      = L / phi
                center = pos + R * np.array([-np.sin(theta), np.cos(theta)])
                start_deg = np.degrees(theta) - 90 + np.degrees(e1)
                end_deg   = np.degrees(theta) - 90 + np.degrees(phi - e2)
                if start_deg > end_deg:
                    start_deg, end_deg = end_deg, start_deg
                wedge = Wedge(center, R + width/2, start_deg, end_deg,
                              width=width, edgecolor='black', facecolor='orange')
                ax.add_patch(wedge)

                arc = Arc(center, 2*R, 2*R,
                          theta1=start_deg, theta2=end_deg,
                          edgecolor='black', linewidth=1)
                ax.add_patch(arc)

                theta += phi
                pos = center + R * np.array([np.sin(theta), -np.cos(theta)])
            else:
                end = pos + L * np.array([np.cos(theta), np.sin(theta)])
                ax.plot([pos[0], end[0]], [pos[1], end[1]],
                        color='black', linewidth=1)
                pos = end


        elif et == 'RFCA':
            cav = Rectangle((0, -width/4), L, width/2,
                            edgecolor='black', facecolor='red',
                            hatch='//', transform=trans)
            ax.add_patch(cav)

            beam = Line2D([0, L], [0, 0], transform=trans,
                          color='black', linewidth=1)
            ax.add_line(beam)

            pos += L * np.array([np.cos(theta), np.sin(theta)])


        else:
            end = pos + L * np.array([np.cos(theta), np.sin(theta)])
            ax.plot([pos[0], end[0]], [pos[1], end[1]],
                    color='black', linewidth=1)
            pos = end

    ax.set_aspect('equal')
    if axis:
        ax.set_xlabel('[m]')
        ax.set_ylabel('[m]')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # remove left and bottom spines if axis is false
    if not axis:
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.tight_layout()
    plt.show()

"""
end

"""
plot_lattice(lattice, scale=0.25, axis=true)
Function wrapper to plot a lattice using matplotlib in Python.
lattice: A list of elements in the lattice.
scale: Scale factor for the plotted elements.
axis: If true, show the axis; otherwise, hide it.
"""
function plot_lattice(lattice, scale=0.25, axis=true)
    py_lattice = [elem_to_py(e) for e in lattice]
    py"plot_lattice"(py_lattice, scale, axis)
end

