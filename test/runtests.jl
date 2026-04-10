# in Julia REPL, type ] to enter package mode
# activate JuTrack
# test JuTrack

using JuTrack
using Test

function f(x)
    k1 = x[1]
    k2 = x[2]
    x1 = CTPS(1.0, 1, 6, 3)
    x2 = CTPS(2.0, 2, 6, 3)
    y = k1*x1^2 + k2*x2^2
    return y.map
end

@testset "JuTrack.jl" begin
ctps = CTPS(Float64, 6, 3)
ctps1 = CTPS(2.0, 6, 3)
ctps2 = CTPS(3.0, 1, 6, 3)
findindex(ctps2, [1, 0, 0, 0, 0, 0, 0])
ind = findindex(ctps2, [1, 0, 0, 0, 0, 0, 0])
ctps3 = ctps2*ctps2
ctps4 = 1.0 / ctps2
ctps5 = exp(ctps2)
ctps6 = log(ctps2)
ctps7 = sqrt(ctps2)
ctps8 = ctps2^ 2
ctps9 = sin(ctps2)
ctps10 = cos(ctps2)
ctps11 = asin(CTPS(0.5, 1, 6, 3))
ctps12 = acos(CTPS(0.5, 1, 6, 3))
ctps13 = sinh(ctps2)
ctps14 = cosh(ctps2)
println(ctps1.map)
println(ctps2.map)
println(ctps3.map)
println(ctps4.map)
println(ctps5.map)
println(ctps6.map)
println(ctps7.map)
println(ctps8.map)
println(ctps9.map)
println(ctps10.map)
println(ctps11.map)
println(ctps12.map)
println(ctps13.map)
println(ctps14.map)

k1 = 3.0
k2 = 1.5
y = f([k1, k2])
println(y)
grad = jacobian(Forward, f, [k1, k2])
println(grad)

# particle tracking
beam = Beam([0.001 0.0 0.0 0.0 0.0 0.0])
D1 = DRIFT(len=1.0)
Q1 = QUAD(k1=1.2, len=0.5)
Q2 = QUAD(k1=-1.2, len=0.5)
line = [D1, Q1, D1, Q2, D1]
linepass!(line, beam)
println(beam.r)

# fast tpsa
set_tps_dim(6)
x = DTPSAD(0.001, 1)
xp = DTPSAD(0.0, 2)
y = DTPSAD(0.0, 3)
yp = DTPSAD(0.0, 4)
z = DTPSAD(0.0, 5)
dp = DTPSAD(0.0, 6)
beam_tpsa = Beam([x xp y yp z dp])
D1 = DRIFT(len=DTPSAD(1.0))
Q1 = QUAD(k1=DTPSAD(1.2), len=DTPSAD(0.5))
Q2 = QUAD(k1=DTPSAD(-1.2), len=DTPSAD(0.5))
line_ = [D1, Q1, D1, Q2, D1]
linepass!(line_, beam_tpsa)
println(beam_tpsa.r)

@testset "SPACECHARGE2P5D DTPSAD" begin
    set_tps_dim(1)
    pts = [
        -1.0e-3  2.0e-4  -0.8e-3  1.5e-4  -1.2e-4  2.0e-6;
        -5.0e-4 -1.0e-4   2.0e-4 -1.0e-4  -6.0e-5 -1.0e-6;
         2.0e-4  1.0e-4  -3.0e-4  2.5e-4  -2.0e-5  1.5e-6;
         6.0e-4 -2.0e-4   7.0e-4 -2.0e-4   3.0e-5 -2.0e-6;
         9.0e-4  3.0e-4  -6.0e-4  1.0e-4   7.0e-5  1.0e-6;
         1.3e-3 -1.5e-4   1.1e-3 -2.5e-4   1.1e-4 -1.5e-6;
    ]

    beam_f = Beam(copy(pts), 1.0e9; np=1_000_000, charge=1.0, mass=m_p)
    sc_f = SPACECHARGE2P5D(effective_len=0.05, xsize=8, ysize=8, zsize=4,
        pipe_radius=0.013, xy_ratio=1.0, long_avg_n=3)
    pass!(sc_f, beam_f.r, beam_f.nmacro, beam_f)

    beam_d = Beam(DTPSAD.(pts); energy=1.0e9, np=1_000_000, charge=1.0, mass=m_p)
    sc_d = SPACECHARGE2P5D(effective_len=DTPSAD(0.05, 1), xsize=8, ysize=8, zsize=4,
        pipe_radius=DTPSAD(0.013), xy_ratio=DTPSAD(1.0), long_avg_n=3)
    pass!(sc_d, beam_d.r, beam_d.nmacro, beam_d)

    maxdiff = 0.0
    for i in eachindex(beam_f.r)
        maxdiff = max(maxdiff, abs(beam_f.r[i] - beam_d.r[i].val))
    end
    @test maxdiff < 1.0e-18
    @test abs(beam_d.r[1, 2].deriv[1]) > 0.0
    @test abs(beam_d.r[1, 6].deriv[1]) > 0.0

    line_sc = [SPACECHARGE2P5D(effective_len=0.05, xsize=8, ysize=8, zsize=4,
        pipe_radius=0.013, xy_ratio=1.0, long_avg_n=3)]
    line_sc_d = Number2TPSAD(line_sc, DTPSAD{1, Float64})
    @test line_sc_d[1].rho_grid isa Matrix{Float64}
    @test line_sc_d[1].green_dx isa Float64
end

@testset "Thick SC2P5D elements" begin
    set_tps_dim(1)
    pts = [
        -8.0e-4  1.5e-4  -6.0e-4  1.0e-4  -8.0e-5  1.0e-6;
        -2.0e-4 -1.2e-4   1.5e-4 -8.0e-5  -2.0e-5 -5.0e-7;
         3.5e-4  8.0e-5  -2.5e-4  1.8e-4   2.0e-5  7.0e-7;
         7.0e-4 -1.8e-4   5.5e-4 -1.4e-4   7.5e-5 -1.2e-6;
    ]

    make_beam_f() = Beam(copy(pts), 1.0e9; np=1_000_000, charge=1.0, mass=m_p)
    make_beam_d() = Beam(DTPSAD.(pts); energy=1.0e9, np=1_000_000, charge=1.0, mass=m_p)

    elements = [
        QUAD_SC2P5D(len=0.08, k1=0.7, xsize=8, ysize=8, zsize=4, pipe_radius=0.013, xy_ratio=1.0, Nsteps=2),
        KSEXT_SC2P5D(len=0.08, k2=4.5, xsize=8, ysize=8, zsize=4, pipe_radius=0.013, xy_ratio=1.0, Nsteps=2),
        KOCT_SC2P5D(len=0.08, k3=12.0, xsize=8, ysize=8, zsize=4, pipe_radius=0.013, xy_ratio=1.0, Nsteps=2),
        SBEND_SC2P5D(len=0.10, angle=0.01, e1=0.005, e2=0.005, xsize=8, ysize=8, zsize=4,
            pipe_radius=0.013, xy_ratio=1.0, Nsteps=2, NumIntSteps=4),
    ]

    rb = RBEND_SC2P5D(len=0.10, angle=0.01, xsize=8, ysize=8, zsize=4, pipe_radius=0.013,
        xy_ratio=1.0, Nsteps=2, NumIntSteps=4)
    @test rb isa SBEND_SC2P5D
    push!(elements, rb)

    for ele in elements
        beam_f = make_beam_f()
        pass!(ele, beam_f.r, beam_f.nmacro, beam_f)
        @test sum(beam_f.lost_flag) == 0

        ele_d = Number2TPSAD([ele], DTPSAD{1, Float64})[1]
        beam_d = make_beam_d()
        pass!(ele_d, beam_d.r, beam_d.nmacro, beam_d)

        maxdiff = 0.0
        for i in eachindex(beam_f.r)
            maxdiff = max(maxdiff, abs(beam_f.r[i] - beam_d.r[i].val))
        end
        @test maxdiff < 1.0e-12
        @test ele_d.rho_grid isa Matrix{Float64}
        @test ele_d.green_fft isa Matrix{ComplexF64}
    end
end

@testset "SC2P5D corrected linearized drift mode" begin
    set_tps_dim(1)
    pts = [
        -7.5e-4  1.8e-4  -5.0e-4  1.1e-4  -7.0e-5  2.0e-6;
        -1.0e-4 -1.1e-4   1.0e-4 -7.0e-5  -1.0e-5 -9.0e-7;
         4.0e-4  9.0e-5  -2.0e-4  1.6e-4   3.0e-5  1.2e-6;
         8.0e-4 -1.6e-4   6.0e-4 -1.2e-4   9.0e-5 -1.8e-6;
    ]

    elements = [
        QUAD_SC2P5D(len=0.08, k1=0.7, xsize=8, ysize=8, zsize=4, pipe_radius=0.013,
            xy_ratio=1.0, Nsteps=2),
        SBEND_SC2P5D(len=0.10, angle=0.01, e1=0.005, e2=0.005, xsize=8, ysize=8, zsize=4,
            pipe_radius=0.013, xy_ratio=1.0, Nsteps=2, NumIntSteps=4),
    ]

    use_exact_drift(2)
    try
        for ele in elements
            beam_f = Beam(copy(pts), 1.0e9; np=1_000_000, charge=1.0, mass=m_p)
            pass!(ele, beam_f.r, beam_f.nmacro, beam_f)
            @test sum(beam_f.lost_flag) == 0

            ele_d = Number2TPSAD([ele], DTPSAD{1, Float64})[1]
            beam_d = Beam(DTPSAD.(pts); energy=1.0e9, np=1_000_000, charge=1.0, mass=m_p)
            pass!(ele_d, beam_d.r, beam_d.nmacro, beam_d)

            maxdiff = 0.0
            for i in eachindex(beam_f.r)
                maxdiff = max(maxdiff, abs(beam_f.r[i] - beam_d.r[i].val))
            end
            @test maxdiff < 1.0e-12
        end
    finally
        use_exact_drift(1)
    end
end

@testset "Basic tracking kernels" begin
    coords = [
         2.0e-4   3.0e-5  -1.0e-4   2.0e-5   7.0e-4   1.0e-4;
        -3.0e-4  -4.0e-5   2.0e-4  -3.0e-5  -5.0e-4  -8.0e-5;
         1.5e-4  -2.0e-5   3.5e-4   5.0e-5   2.0e-4   6.0e-5;
    ]
    beam = Beam(copy(coords), 1.0e9)
    beti = use_exact_beti == 1 ? 1.0 / beam.beta : 1.0
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / (beam.gamma^2) : 0.0

    zeros6 = zeros(Float64, 6)
    zeros66 = zeros(Float64, 6, 6)
    zeros2 = zeros(Float64, 2)
    fringe0 = zeros(Float64, 5)
    polyA = zeros(Float64, 4)
    polyB = [0.0, 1.1, 0.0, 0.0]
    bendA = zeros(Float64, 4)
    bendB = [0.0, 0.0, 0.0, 0.0]

    function compare_kernel!(serial_runner, threaded_runner; tol = 1.0e-12)
        rs = copy(coords)
        ls = zeros(Int, size(coords, 1))
        serial_runner(rs, ls)
        rp = copy(coords)
        lp = zeros(Int, size(coords, 1))
        threaded_runner(rp, lp)
        @test ls == lp
        @test maximum(abs.(rs .- rp)) < tol
    end

    compare_kernel!(
        (r, lost) -> JuTrack.DriftPass!(r, 0.4, beti, gamma2i, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, size(r, 1), lost),
        (r, lost) -> JuTrack.DriftPass_P!(r, 0.4, beti, gamma2i, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, size(r, 1), lost),
    )

    compare_kernel!(
        (r, lost) -> JuTrack.QuadLinearPass!(r, 0.2, 0.9, beti, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, size(r, 1), lost),
        (r, lost) -> JuTrack.QuadLinearPass_P!(r, 0.2, 0.9, beti, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, size(r, 1), lost),
    )

    compare_kernel!(
        (r, lost) -> JuTrack.StrMPoleSymplectic4Pass!(r, 0.2, beti, polyA, copy(polyB), 1, 8, 0, 0, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, zeros2, size(r, 1), lost, gamma2i),
        (r, lost) -> JuTrack.StrMPoleSymplectic4Pass_P!(r, 0.2, beti, polyA, copy(polyB), 1, 8, 0, 0, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, zeros2, size(r, 1), lost, gamma2i),
    )

    compare_kernel!(
        (r, lost) -> JuTrack.BendSymplecticPass!(r, 0.3, beti, 0.05, bendA, copy(bendB), 0, 8, 0.0, 0.0, 0, 0, 0.0, 0.0, 0.0, 0, 0, fringe0, fringe0, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, zeros2, size(r, 1), lost),
        (r, lost) -> JuTrack.BendSymplecticPass_P!(r, 0.3, beti, 0.05, bendA, copy(bendB), 0, 8, 0.0, 0.0, 0, 0, 0.0, 0.0, 0.0, 0, 0, fringe0, fringe0, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, zeros2, size(r, 1), lost),
    )

    compare_kernel!(
        (r, lost) -> JuTrack.CorrectorPass!(r, 0.1, 1.0e-4, -1.5e-4, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, size(r, 1), lost),
        (r, lost) -> JuTrack.CorrectorPass_P!(r, 0.1, 1.0e-4, -1.5e-4, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, size(r, 1), lost),
    )

    compare_kernel!(
        (r, lost) -> JuTrack.ThinMPolePass!(r, 0.0, polyA, copy(polyB), 1, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, zeros2, size(r, 1), lost),
        (r, lost) -> JuTrack.ThinMPolePass_P!(r, 0.0, polyA, copy(polyB), 1, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, zeros2, size(r, 1), lost),
    )

    t0 = 1.0 / 500.0e6
    compare_kernel!(
        (r, lost) -> JuTrack.RFCavityPass!(r, 0.05, 2.0e-4, 500.0e6, 400.0, 0.0, 0.0, 0, t0, beam.beta, beti, size(r, 1), lost, zeros6, zeros6),
        (r, lost) -> JuTrack.RFCavityPass_P!(r, 0.05, 2.0e-4, 500.0e6, 400.0, 0.0, 0.0, 0, t0, beam.beta, beti, size(r, 1), lost, zeros6, zeros6),
    )

    compare_kernel!(
        (r, lost) -> begin
            b = Beam(copy(r), 1.0e9)
            b.lost_flag .= lost
            pass!(SOLENOID(len = 0.08, ks = 0.3, T1 = zeros6, T2 = zeros6, R1 = zeros66, R2 = zeros66), b.r, b.nmacro, b)
            r .= b.r
            lost .= b.lost_flag
        end,
        (r, lost) -> begin
            b = Beam(copy(r), 1.0e9)
            b.lost_flag .= lost
            pass_P!(SOLENOID(len = 0.08, ks = 0.3, T1 = zeros6, T2 = zeros6, R1 = zeros66, R2 = zeros66), b.r, b.nmacro, b)
            r .= b.r
            lost .= b.lost_flag
        end,
    )

    function quad_kernel_obj(k1)
        r = [1.0e-4 2.0e-5 -1.2e-4 1.0e-5 1.0e-3 2.0e-4]
        lost = [0]
        JuTrack.QuadLinearPass!(r, 0.1, k1, 1.0, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, 1, lost)
        return r[1, 1]
    end
    function kquad_kernel_obj(k1)
        r = [1.0e-4 2.0e-5 -1.2e-4 1.0e-5 1.0e-3 2.0e-4]
        lost = [0]
        A = zeros(4)
        B = zeros(4)
        B[2] = k1
        JuTrack.StrMPoleSymplectic4Pass!(r, 0.1, 1.0, A, B, 1, 8, 0, 0, zeros6, zeros6, zeros66, zeros66, zeros6, zeros6, zeros2, 1, lost, 0.0)
        return r[1, 1]
    end
    quad_grad, quad_primal = autodiff(ForwardWithPrimal, quad_kernel_obj, Duplicated(0.8, 1.0))
    kquad_grad, kquad_primal = autodiff(ForwardWithPrimal, kquad_kernel_obj, Duplicated(0.8, 1.0))
    @test isfinite(quad_grad[1])
    @test isfinite(quad_primal)
    @test isfinite(kquad_grad[1])
    @test isfinite(kquad_primal)
end
end
