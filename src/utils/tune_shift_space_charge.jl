# estimate tune shift due to space charge for coasting beam
function tune_shift_SC(emit, beta, gamma, C, I, r0=1.54e-18, e=1.602176634e-19)
    spd_light = 2.99792458e8
    # r0 is classical particle radius, 1.54e-18 m for protons.
    N = I * C / beta / spd_light / e
    dQ = r0 * N / 2 / pi / emit / beta^2 / gamma^3
end

# include("../JuTrack.jl")
# using .JuTrack
# Dr = DRIFT(name="Dr", len=0.5)
# HalfDr = DRIFT(name="HalfDr", len=0.25)
# p2Dr = DRIFT(name="p2Dr", len=0.2)
# SF = KSEXT(name="SF", len=0.1, k2=1.0)
# SD = KSEXT(name="SD", len=0.1, k2=-1.0)
# B1 = SBEND(name="B", len= 1.0, angle=2*pi/40.0)
# Q1 = QUAD(name="Q1", len=0.5, k1=1.2) # optimized k1 starting from -1.0
# Q2 = QUAD(name="Q2", len=0.5, k1=-1.2)

# cell = [HalfDr, B1, p2Dr, SF, p2Dr, Q1, Dr, B1, p2Dr, SD, p2Dr, Q2, HalfDr]
# ring = [cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell...,
#         cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell...]

# s = 100

# beam = Beam(zeros(10001, 6), energy=1.0e9, current=50.0, mass=m_p, charge=1.0, emittance=[1e-6, 1e-6, 0.0])

# beta = beam.beta
# gamma = beam.gamma
# emit = 1e-8/ (beta * gamma)

# dQ = tune_shift_SC(emit, beta, gamma, s[end], beam.current)