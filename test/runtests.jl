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

converted4 = Number2TPSAD([CORRECTOR()], DTPSAD{4, Float64})
@test converted4[1].len isa DTPSAD{4, Float64}
@test converted4[1].R1[1, 1] isa DTPSAD{4, Float64}

@testset "Taylor core infrastructure" begin
    table = monomial_table(Val(2), Val(3))
    proto = CTPS(Float64, 2, 3)
    @test length(table) == binomial(2 + 3, 3)
    @test collect(degree_range(table, 0)) == [1]
    @test length(collect(degree_range(table, 1))) == 2
    @test length(collect(degree_range(table, 2))) == 3
    @test length(collect(degree_range(table, 3))) == 4
    for exponent in table.exponents
        @test monomial_index(table, exponent) == findindex(proto, collect(exponent))
        @test monomial_exponents(table, monomial_index(table, exponent)) == exponent
    end
    @test monomial_product_index(table, (1, 1), (1, 0)) == monomial_index(table, (2, 1))
    @test monomial_product_index(table, (2, 1), (0, 1)) == 0

    f_scalar(x1, x2) = x1^2 * x2 + sin(x2)
    F_vector(x1, x2) = [x1 * x2, x1 + cos(x2)]
    x0 = [1.2, -0.3]

    dtpsad_type = jet_type(DTPSADBackend(), 2)
    @test dtpsad_type == DTPSAD{2, Float64}
    dtpsad_seed = jet_seed(dtpsad_type, x0[1], 1)
    @test jet_backend(dtpsad_seed) isa DTPSADBackend
    @test jet_nvars(dtpsad_seed) == 2
    @test jet_order(dtpsad_seed) == 1
    @test jet_primal(dtpsad_seed) == x0[1]
    @test jet_coefficient(dtpsad_seed, (1, 0)) == 1.0
    @test jet_coefficient(dtpsad_seed, (0, 1)) == 0.0
    @test jet_hessian(f_scalar(seeded_jets(dtpsad_type, x0)...)) == zeros(2, 2)

    grad_dtpsad, primal_dtpsad = TaylorGradient(dtpsad_type, f_scalar, x0, primal=true)
    @test isapprox(primal_dtpsad, x0[1]^2 * x0[2] + sin(x0[2]); atol=1.0e-12)
    @test all(isapprox.(grad_dtpsad, [2 * x0[1] * x0[2], x0[1]^2 + cos(x0[2])]; atol=1.0e-12))
    @test TaylorGradient(dtpsad_type, f_scalar, x0) == Gradient(f_scalar, x0)

    J_dtpsad, vals_dtpsad = TaylorJacobian(dtpsad_type, F_vector, x0, primal=true)
    @test isapprox(vals_dtpsad[1], x0[1] * x0[2]; atol=1.0e-12)
    @test isapprox(vals_dtpsad[2], x0[1] + cos(x0[2]); atol=1.0e-12)
    @test all(isapprox.(J_dtpsad, [x0[2] x0[1]; 1.0 -sin(x0[2])]; atol=1.0e-12))
    @test TaylorJacobian(dtpsad_type, F_vector, x0) == Jacobian(F_vector, x0)

    ctps_type = jet_type(CTPSBackend(), 2, order=3)
    @test ctps_type == CTPS{Float64, 2, 3}
    ctps_seed = jet_seed(ctps_type, x0[2], 2)
    @test jet_backend(ctps_seed) isa CTPSBackend
    @test jet_nvars(ctps_seed) == 2
    @test jet_order(ctps_seed) == 3
    @test jet_primal(ctps_seed) == x0[2]
    @test jet_coefficient(ctps_seed, (0, 1)) == 1.0
    @test jet_coefficient(ctps_seed, (1, 0)) == 0.0

    grad_ctps, primal_ctps = TaylorGradient(ctps_type, f_scalar, x0, primal=true)
    @test isapprox(primal_ctps, primal_dtpsad; atol=1.0e-12)
    @test all(isapprox.(grad_ctps, grad_dtpsad; atol=1.0e-12))

    y_ctps = f_scalar(seeded_jets(ctps_type, x0)...)
    H_ctps = jet_hessian(y_ctps)
    H_expected = [2 * x0[2] 2 * x0[1]; 2 * x0[1] -sin(x0[2])]
    @test all(isapprox.(H_ctps, H_expected; atol=1.0e-12))

    J_ctps, vals_ctps = TaylorJacobian(ctps_type, F_vector, x0, primal=true)
    @test all(isapprox.(vals_ctps, vals_dtpsad; atol=1.0e-12))
    @test all(isapprox.(J_ctps, J_dtpsad; atol=1.0e-12))

    hot_type = jet_type(HOTPSABackend(), 2, order=3)
    @test hot_type == hotpsa_type(2, 3)
    hot_seed = jet_seed(hot_type, x0[1], 1)
    @test jet_backend(hot_seed) isa HOTPSABackend
    @test jet_nvars(hot_seed) == 2
    @test jet_order(hot_seed) == 3
    @test jet_primal(hot_seed) == x0[1]
    @test jet_coefficient(hot_seed, (1, 0)) == 1.0
    @test jet_coefficient(hot_seed, (0, 1)) == 0.0

    grad_hot, primal_hot = TaylorGradient(hot_type, f_scalar, x0, primal=true)
    @test isapprox(primal_hot, primal_dtpsad; atol=1.0e-12)
    @test all(isapprox.(grad_hot, grad_dtpsad; atol=1.0e-12))

    y_hot = f_scalar(seeded_jets(hot_type, x0)...)
    H_hot = jet_hessian(y_hot)
    @test all(isapprox.(H_hot, H_expected; atol=1.0e-12))

    J_hot, vals_hot = TaylorJacobian(hot_type, F_vector, x0, primal=true)
    @test all(isapprox.(vals_hot, vals_dtpsad; atol=1.0e-12))
    @test all(isapprox.(J_hot, J_dtpsad; atol=1.0e-12))
end

@testset "Taylor transport core internals" begin
    R = zeros(6, 6)
    for i in 1:6
        R[i, i] = 1.0
    end
    R[1, 2] = 0.01
    R[3, 4] = -0.02
    T1 = [1.0e-5, 0.0, -2.0e-5, 0.0, 0.0, 0.0]
    T2 = [0.0, 0.0, 0.0, 0.0, 3.0e-6, 0.0]
    line_taylor = [
        DRIFT(len=0.35, T1=copy(T1), T2=copy(T2), R1=copy(R), R2=copy(R)),
        CORRECTOR(len=0.08, xkick=1.2e-4, ykick=-9.0e-5),
        QUAD(len=0.22, k1=0.7),
        MARKER(),
    ]

    r0 = [1.0e-3, 2.5e-4, -8.0e-4, 1.5e-4, 4.0e-5, 2.0e-3]
    beam_single = Beam(reshape(copy(r0), 1, 6); energy=3.0e9)
    linepass!(line_taylor, beam_single)

    r_taylor = copy(r0)
    JuTrack.linepass_Taylor!(line_taylor, r_taylor; E0=3.0e9, m0=m_e)
    @test all(isapprox.(r_taylor, vec(beam_single.r[1, :]); atol=1.0e-12))

    ctconst(v) = CTPS(v, 1, 2)
    len_seed = CTPS(0.4, 1, 1, 2)
    ctzeros = [ctconst(0.0) for _ in 1:6]
    ctidentity = [i == j ? ctconst(1.0) : ctconst(0.0) for i in 1:6, j in 1:6]
    drift_seed = DRIFT("DRIFT", len_seed, copy(ctzeros), copy(ctzeros), copy(ctidentity), copy(ctidentity), zeros(6), zeros(6), "DRIFT")
    r_seed = [ctconst(1.0e-3), ctconst(2.5e-4), ctconst(-8.0e-4), ctconst(1.5e-4), ctconst(4.0e-5), ctconst(2.0e-3)]
    JuTrack.pass_Taylor!(drift_seed, r_seed; E0=3.0e9, m0=m_e)
    if use_exact_Hamiltonian == 1
        denom = sqrt(1.0 + 2.0 * 2.0e-3 + (2.0e-3)^2 - (2.5e-4)^2 - (1.5e-4)^2)
        expected_dx_dlen = 2.5e-4 / denom
    else
        expected_dx_dlen = 2.5e-4 / (1.0 + 2.0e-3)
    end
    @test isapprox(jet_coefficient(r_seed[1], (1,)), expected_dx_dlen; atol=1.0e-12)
end

@testset "HOTPSA standard tracking dispatch" begin
    hot_type = jet_type(HOTPSABackend(), 2, order=3)
    converted_hot = Number2TPSAD([CORRECTOR()], hot_type)
    @test converted_hot[1].len isa hot_type
    @test converted_hot[1].R1[1, 1] isa hot_type

    line = [
        DRIFT(len=0.35),
        CORRECTOR(len=0.08, xkick=1.2e-4, ykick=-9.0e-5),
        QUAD(len=0.22, k1=0.7),
        MARKER(),
    ]
    r0 = [1.0e-3, 2.5e-4, -8.0e-4, 1.5e-4, 4.0e-5, 2.0e-3]

    beam_single = Beam(reshape(copy(r0), 1, 6); energy=3.0e9)
    linepass!(line, beam_single)

    hot_state = [
        jet_seed(hot_type, r0[1], 1),
        jet_seed(hot_type, r0[2], 2),
        hot_type(r0[3]),
        hot_type(r0[4]),
        hot_type(r0[5]),
        hot_type(r0[6]),
    ]
    ctps_state = [
        CTPS(r0[1], 1, 2, 3),
        CTPS(r0[2], 2, 2, 3),
        CTPS(r0[3], 2, 3),
        CTPS(r0[4], 2, 3),
        CTPS(r0[5], 2, 3),
        CTPS(r0[6], 2, 3),
    ]

    hot_beam = Beam(reshape(copy(hot_state), 1, 6); energy=3.0e9)
    linepass!(line, hot_beam)

    linepass!(line, hot_state; E0=3.0e9, m0=m_e)
    linepass_TPSA!(line, ctps_state; E0=3.0e9, m0=m_e)

    table = monomial_table(Val(2), Val(3))
    function assert_hotpsa_matches_existing(line, hot_state, ctps_state;
        float_reference=nothing, E0=3.0e9, m0=m_e, primal_atol=1.0e-12, coeff_atol=1.0e-10)
        linepass!(deepcopy(line), hot_state; E0=E0, m0=m0)
        linepass_TPSA!(deepcopy(line), ctps_state; E0=E0, m0=m0)
        for i in 1:6
            if float_reference === nothing
                @test isapprox(jet_primal(hot_state[i]), jet_primal(ctps_state[i]); atol=primal_atol)
            else
                @test isapprox(jet_primal(hot_state[i]), float_reference[i]; atol=primal_atol)
                @test isapprox(jet_primal(ctps_state[i]), float_reference[i]; atol=primal_atol)
            end
            for exponent in table.exponents
                @test isapprox(jet_coefficient(hot_state[i], exponent),
                    jet_coefficient(ctps_state[i], exponent); atol=coeff_atol)
            end
        end
    end
    for i in 1:6
        @test isapprox(jet_primal(hot_state[i]), beam_single.r[1, i]; atol=1.0e-12)
        @test isapprox(jet_primal(hot_beam.r[1, i]), beam_single.r[1, i]; atol=1.0e-12)
        for exponent in table.exponents
            @test isapprox(jet_coefficient(hot_state[i], exponent),
                jet_coefficient(ctps_state[i], exponent); atol=1.0e-10)
            @test isapprox(jet_coefficient(hot_beam.r[1, i], exponent),
                jet_coefficient(ctps_state[i], exponent); atol=1.0e-10)
        end
    end

    geo_line = [
        TRANSLATION(dx=2.0e-4, dy=-1.5e-4, ds=3.0e-3),
        YROTATION(angle=6.0e-3),
    ]
    hot_geo = [
        jet_seed(hot_type, 6.0e-4, 1),
        hot_type(1.5e-4),
        jet_seed(hot_type, -4.0e-4, 2),
        hot_type(-1.2e-4),
        hot_type(3.0e-5),
        hot_type(1.0e-3),
    ]
    ctps_geo = [
        CTPS(6.0e-4, 1, 2, 3),
        CTPS(1.5e-4, 2, 3),
        CTPS(-4.0e-4, 2, 2, 3),
        CTPS(-1.2e-4, 2, 3),
        CTPS(3.0e-5, 2, 3),
        CTPS(1.0e-3, 2, 3),
    ]
    linepass!(geo_line, hot_geo; E0=3.0e9, m0=m_e)
    linepass_TPSA!(geo_line, ctps_geo; E0=3.0e9, m0=m_e)
    for i in 1:6
        @test isapprox(jet_primal(hot_geo[i]), jet_primal(ctps_geo[i]); atol=1.0e-12)
        for exponent in table.exponents
            @test isapprox(jet_coefficient(hot_geo[i], exponent),
                jet_coefficient(ctps_geo[i], exponent); atol=1.0e-10)
        end
    end

    mpole_line = [
        thinMULTIPOLE(PolynomA=[0.0, 0.0, 0.02, -0.01], PolynomB=[0.0, 0.15, -0.03, 0.01], MaxOrder=3),
        KQUAD(len=0.05, k1=0.7, MaxOrder=1, NumIntSteps=4),
        KSEXT(len=0.04, k2=12.0, MaxOrder=2, NumIntSteps=4),
        KOCT(len=0.03, k3=-80.0, MaxOrder=3, NumIntSteps=4),
    ]
    hot_mpole = [
        jet_seed(hot_type, 5.0e-4, 1),
        hot_type(1.8e-4),
        jet_seed(hot_type, -3.0e-4, 2),
        hot_type(-1.1e-4),
        hot_type(2.5e-5),
        hot_type(8.0e-4),
    ]
    ctps_mpole = [
        CTPS(5.0e-4, 1, 2, 3),
        CTPS(1.8e-4, 2, 3),
        CTPS(-3.0e-4, 2, 2, 3),
        CTPS(-1.1e-4, 2, 3),
        CTPS(2.5e-5, 2, 3),
        CTPS(8.0e-4, 2, 3),
    ]
    linepass!(mpole_line, hot_mpole; E0=3.0e9, m0=m_e)
    linepass_TPSA!(mpole_line, ctps_mpole; E0=3.0e9, m0=m_e)
    for i in 1:6
        @test isapprox(jet_primal(hot_mpole[i]), jet_primal(ctps_mpole[i]); atol=1.0e-11)
        for exponent in table.exponents
            @test isapprox(jet_coefficient(hot_mpole[i], exponent),
                jet_coefficient(ctps_mpole[i], exponent); atol=1.0e-8)
        end
    end

    sbend_line = [
        SBEND(len=0.42, angle=0.08, e1=0.03, e2=-0.025,
            PolynomB=[0.0, 0.35, -0.04, 0.01], MaxOrder=3, NumIntSteps=4,
            rad=0, fint1=0.12, fint2=0.09, gap=0.018,
            FringeBendEntrance=1, FringeBendExit=1,
            FringeQuadEntrance=1, FringeQuadExit=1,
            FringeIntM0=[0.6, 0.2, 0.05, -0.03, 0.01],
            FringeIntP0=[0.4, -0.1, 0.03, 0.02, -0.01]),
    ]
    r0_sbend = [4.0e-4, 1.2e-4, -2.5e-4, -1.0e-4, 0.0, 0.0]
    beam_sbend = Beam(reshape(copy(r0_sbend), 1, 6); energy=3.0e9)
    linepass!(deepcopy(sbend_line), beam_sbend)
    hot_sbend = [
        jet_seed(hot_type, r0_sbend[1], 1),
        hot_type(r0_sbend[2]),
        jet_seed(hot_type, r0_sbend[3], 2),
        hot_type(r0_sbend[4]),
        hot_type(r0_sbend[5]),
        hot_type(r0_sbend[6]),
    ]
    ctps_sbend = [
        CTPS(r0_sbend[1], 1, 2, 3),
        CTPS(r0_sbend[2], 2, 3),
        CTPS(r0_sbend[3], 2, 2, 3),
        CTPS(r0_sbend[4], 2, 3),
        CTPS(r0_sbend[5], 2, 3),
        CTPS(r0_sbend[6], 2, 3),
    ]
    assert_hotpsa_matches_existing(sbend_line, hot_sbend, ctps_sbend;
        float_reference=vec(beam_sbend.r[1, :]), primal_atol=1.0e-11, coeff_atol=1.0e-8)

    sbend_elegant_line = [
        SBEND(len=0.42, angle=0.08, e1=0.03, e2=-0.025,
            PolynomB=[0.0, 0.35, -0.04, 0.01], MaxOrder=3, NumIntSteps=4,
            rad=0, fint1=0.12, fint2=0.09, gap=0.018,
            FringeBendEntrance=1, FringeBendExit=1,
            FringeQuadEntrance=2, FringeQuadExit=2,
            FringeIntM0=[0.6, 0.2, 0.05, -0.03, 0.01],
            FringeIntP0=[0.4, -0.1, 0.03, 0.02, -0.01]),
    ]
    hot_sbend_elegant = [
        jet_seed(hot_type, r0_sbend[1], 1),
        hot_type(r0_sbend[2]),
        jet_seed(hot_type, r0_sbend[3], 2),
        hot_type(r0_sbend[4]),
        hot_type(r0_sbend[5]),
        hot_type(r0_sbend[6]),
    ]
    linepass!(deepcopy(sbend_elegant_line), hot_sbend_elegant; E0=3.0e9, m0=m_e)
    beam_sbend_elegant = Beam(reshape(copy(r0_sbend), 1, 6); energy=3.0e9)
    linepass!(deepcopy(sbend_elegant_line), beam_sbend_elegant)
    for i in 1:6
        @test isapprox(jet_primal(hot_sbend_elegant[i]), beam_sbend_elegant.r[1, i]; atol=1.0e-11)
    end
    fd_step = 1.0e-8
    for (coord, exponent) in ((1, (1, 0)), (3, (0, 1)))
        rp = copy(r0_sbend)
        rm = copy(r0_sbend)
        rp[coord] += fd_step
        rm[coord] -= fd_step
        bp = Beam(reshape(rp, 1, 6); energy=3.0e9)
        bm = Beam(reshape(rm, 1, 6); energy=3.0e9)
        linepass!(deepcopy(sbend_elegant_line), bp)
        linepass!(deepcopy(sbend_elegant_line), bm)
        fd = (vec(bp.r[1, :]) .- vec(bm.r[1, :])) ./ (2fd_step)
        for i in 1:6
            @test isapprox(jet_coefficient(hot_sbend_elegant[i], exponent), fd[i];
                rtol=1.0e-5, atol=1.0e-7)
        end
    end

    esbend_line = [
        ESBEND(len=0.38, angle=0.06, e1=0.025, e2=-0.018,
            PolynomB=[0.0, 0.28, -0.02, 0.008], MaxOrder=3, NumIntSteps=4,
            rad=0, gK=0.35, FringeBendEntrance=1, FringeBendExit=1,
            FringeQuadEntrance=1, FringeQuadExit=1),
    ]
    r0_esbend = [3.0e-4, 1.1e-4, -2.0e-4, -8.0e-5, 2.0e-5, 5.0e-4]
    beam_esbend = Beam(reshape(copy(r0_esbend), 1, 6); energy=3.0e9)
    linepass!(deepcopy(esbend_line), beam_esbend)
    hot_esbend = [
        jet_seed(hot_type, r0_esbend[1], 1),
        hot_type(r0_esbend[2]),
        jet_seed(hot_type, r0_esbend[3], 2),
        hot_type(r0_esbend[4]),
        hot_type(r0_esbend[5]),
        hot_type(r0_esbend[6]),
    ]
    ctps_esbend = [
        CTPS(r0_esbend[1], 1, 2, 3),
        CTPS(r0_esbend[2], 2, 3),
        CTPS(r0_esbend[3], 2, 2, 3),
        CTPS(r0_esbend[4], 2, 3),
        CTPS(r0_esbend[5], 2, 3),
        CTPS(r0_esbend[6], 2, 3),
    ]
    assert_hotpsa_matches_existing(esbend_line, hot_esbend, ctps_esbend;
        float_reference=vec(beam_esbend.r[1, :]), primal_atol=1.0e-11, coeff_atol=1.0e-8)

    rfca_line = [
        RFCA(len=0.06, volt=2.5e6, freq=500.0e6, h=864.0, lag=1.5e-3,
            philag=0.23, energy=3.0e9),
    ]
    r0_rfca = [2.5e-4, 1.1e-4, -1.8e-4, -7.0e-5, 3.5e-4, 8.0e-4]
    beam_rfca = Beam(reshape(copy(r0_rfca), 1, 6); energy=3.0e9)
    linepass!(deepcopy(rfca_line), beam_rfca)
    hot_rfca = [
        jet_seed(hot_type, r0_rfca[1], 1),
        hot_type(r0_rfca[2]),
        hot_type(r0_rfca[3]),
        hot_type(r0_rfca[4]),
        jet_seed(hot_type, r0_rfca[5], 2),
        hot_type(r0_rfca[6]),
    ]
    ctps_rfca = [
        CTPS(r0_rfca[1], 1, 2, 3),
        CTPS(r0_rfca[2], 2, 3),
        CTPS(r0_rfca[3], 2, 3),
        CTPS(r0_rfca[4], 2, 3),
        CTPS(r0_rfca[5], 2, 2, 3),
        CTPS(r0_rfca[6], 2, 3),
    ]
    assert_hotpsa_matches_existing(rfca_line, hot_rfca, ctps_rfca;
        float_reference=vec(beam_rfca.r[1, :]), primal_atol=1.0e-11, coeff_atol=1.0e-9)

    crab_line = [
        CRABCAVITY(len=0.04, volt=1.2e6, freq=499.0e6, phi=0.31, energy=3.0e9),
    ]
    r0_crab = [3.0e-4, 1.4e-4, -2.2e-4, -9.0e-5, 2.8e-4, 6.0e-4]
    beam_crab = Beam(reshape(copy(r0_crab), 1, 6); energy=3.0e9)
    linepass!(deepcopy(crab_line), beam_crab)
    hot_crab = [
        jet_seed(hot_type, r0_crab[1], 1),
        hot_type(r0_crab[2]),
        hot_type(r0_crab[3]),
        hot_type(r0_crab[4]),
        jet_seed(hot_type, r0_crab[5], 2),
        hot_type(r0_crab[6]),
    ]
    ctps_crab = [
        CTPS(r0_crab[1], 1, 2, 3),
        CTPS(r0_crab[2], 2, 3),
        CTPS(r0_crab[3], 2, 3),
        CTPS(r0_crab[4], 2, 3),
        CTPS(r0_crab[5], 2, 2, 3),
        CTPS(r0_crab[6], 2, 3),
    ]
    assert_hotpsa_matches_existing(crab_line, hot_crab, ctps_crab;
        float_reference=vec(beam_crab.r[1, :]), primal_atol=1.0e-11, coeff_atol=1.0e-9)

    crab_k2_line = [
        CRABCAVITY_K2(len=0.05, volt=1.0e6, freq=501.0e6, phi=0.27, k2=18.0, energy=3.0e9),
    ]
    r0_crab_k2 = [2.7e-4, 1.2e-4, -2.4e-4, -8.5e-5, 3.1e-4, 5.5e-4]
    beam_crab_k2 = Beam(reshape(copy(r0_crab_k2), 1, 6); energy=3.0e9)
    linepass!(deepcopy(crab_k2_line), beam_crab_k2)
    hot_crab_k2 = [
        jet_seed(hot_type, r0_crab_k2[1], 1),
        hot_type(r0_crab_k2[2]),
        jet_seed(hot_type, r0_crab_k2[3], 2),
        hot_type(r0_crab_k2[4]),
        hot_type(r0_crab_k2[5]),
        hot_type(r0_crab_k2[6]),
    ]
    ctps_crab_k2 = [
        CTPS(r0_crab_k2[1], 1, 2, 3),
        CTPS(r0_crab_k2[2], 2, 3),
        CTPS(r0_crab_k2[3], 2, 2, 3),
        CTPS(r0_crab_k2[4], 2, 3),
        CTPS(r0_crab_k2[5], 2, 3),
        CTPS(r0_crab_k2[6], 2, 3),
    ]
    assert_hotpsa_matches_existing(crab_k2_line, copy(hot_crab_k2), ctps_crab_k2;
        float_reference=vec(beam_crab_k2.r[1, :]), primal_atol=1.0e-11, coeff_atol=1.0e-8)
    linepass!(deepcopy(crab_k2_line), hot_crab_k2; E0=3.0e9, m0=m_e)
    for i in 1:6
        @test isapprox(jet_primal(hot_crab_k2[i]), beam_crab_k2.r[1, i]; atol=1.0e-11)
    end
    fd_step_crab_k2 = 1.0e-8
    for (coord, exponent) in ((1, (1, 0)), (3, (0, 1)))
        rp = copy(r0_crab_k2)
        rm = copy(r0_crab_k2)
        rp[coord] += fd_step_crab_k2
        rm[coord] -= fd_step_crab_k2
        bp = Beam(reshape(rp, 1, 6); energy=3.0e9)
        bm = Beam(reshape(rm, 1, 6); energy=3.0e9)
        linepass!(deepcopy(crab_k2_line), bp)
        linepass!(deepcopy(crab_k2_line), bm)
        fd = (vec(bp.r[1, :]) .- vec(bm.r[1, :])) ./ (2fd_step_crab_k2)
        for i in 1:6
            @test isapprox(jet_coefficient(hot_crab_k2[i], exponent), fd[i];
                rtol=1.0e-5, atol=1.0e-7)
        end
    end

    solenoid_R1 = zeros(6, 6)
    solenoid_R2 = zeros(6, 6)
    for i in 1:6
        solenoid_R1[i, i] = 1.0
        solenoid_R2[i, i] = 1.0
    end
    solenoid_R1[1, 2] = 0.015
    solenoid_R2[3, 4] = -0.012
    solenoid_line = [
        SOLENOID(len=0.18, ks=0.9,
            T1=[2.0e-5, 0.0, -1.5e-5, 0.0, 0.0, 0.0],
            T2=[0.0, 0.0, 0.0, 0.0, 4.0e-6, 0.0],
            R1=solenoid_R1, R2=solenoid_R2),
    ]
    r0_solenoid = [4.0e-4, 1.6e-4, -3.0e-4, -1.2e-4, 3.0e-5, 9.0e-4]
    beam_solenoid = Beam(reshape(copy(r0_solenoid), 1, 6); energy=3.0e9)
    linepass!(deepcopy(solenoid_line), beam_solenoid)
    hot_solenoid = [
        jet_seed(hot_type, r0_solenoid[1], 1),
        hot_type(r0_solenoid[2]),
        jet_seed(hot_type, r0_solenoid[3], 2),
        hot_type(r0_solenoid[4]),
        hot_type(r0_solenoid[5]),
        hot_type(r0_solenoid[6]),
    ]
    ctps_solenoid = [
        CTPS(r0_solenoid[1], 1, 2, 3),
        CTPS(r0_solenoid[2], 2, 3),
        CTPS(r0_solenoid[3], 2, 2, 3),
        CTPS(r0_solenoid[4], 2, 3),
        CTPS(r0_solenoid[5], 2, 3),
        CTPS(r0_solenoid[6], 2, 3),
    ]
    assert_hotpsa_matches_existing(solenoid_line, copy(hot_solenoid), ctps_solenoid;
        float_reference=vec(beam_solenoid.r[1, :]), primal_atol=1.0e-11, coeff_atol=1.0e-8)
    linepass!(deepcopy(solenoid_line), hot_solenoid; E0=3.0e9, m0=m_e)
    for i in 1:6
        @test isapprox(jet_primal(hot_solenoid[i]), beam_solenoid.r[1, i]; atol=1.0e-11)
    end
    fd_step_solenoid = 1.0e-8
    for (coord, exponent) in ((1, (1, 0)), (3, (0, 1)))
        rp = copy(r0_solenoid)
        rm = copy(r0_solenoid)
        rp[coord] += fd_step_solenoid
        rm[coord] -= fd_step_solenoid
        bp = Beam(reshape(rp, 1, 6); energy=3.0e9)
        bm = Beam(reshape(rm, 1, 6); energy=3.0e9)
        linepass!(deepcopy(solenoid_line), bp)
        linepass!(deepcopy(solenoid_line), bm)
        fd = (vec(bp.r[1, :]) .- vec(bm.r[1, :])) ./ (2fd_step_solenoid)
        for i in 1:6
            @test isapprox(jet_coefficient(hot_solenoid[i], exponent), fd[i];
                rtol=1.0e-5, atol=1.0e-7)
        end
    end

    c_light = 2.99792458e8
    hotpsa_extra_lines = [
        ("easyCRABCAVITY",
            [easyCRABCAVITY(halfthetac=1.0e-3, k=2.0 * pi * 500.0e6 / c_light, phi=0.2)]),
        ("AccelCavity",
            [AccelCavity(volt=2.0e6, freq=500.0e6, phis=0.15)]),
        ("LBEND",
            [LBEND(len=0.36, angle=0.045, e1=0.012, e2=-0.01,
                K=0.08, ByError=2.0e-4, fint1=0.08, fint2=0.07, FullGap=0.018)]),
        ("LorentzBoost",
            [LorentzBoost(0.015)]),
        ("InvLorentzBoost",
            [InvLorentzBoost(0.015)]),
    ]
    r0_extra = [4.0e-4, 1.2e-4, -2.5e-4, -1.0e-4, 3.0e-5, 7.0e-4]
    fd_step_extra = 1.0e-8
    for (_, extra_line) in hotpsa_extra_lines
        beam_extra = Beam(reshape(copy(r0_extra), 1, 6); energy=3.0e9)
        linepass!(deepcopy(extra_line), beam_extra)
        hot_extra = [
            jet_seed(hot_type, r0_extra[1], 1),
            hot_type(r0_extra[2]),
            jet_seed(hot_type, r0_extra[3], 2),
            hot_type(r0_extra[4]),
            hot_type(r0_extra[5]),
            hot_type(r0_extra[6]),
        ]
        linepass!(deepcopy(extra_line), hot_extra; E0=3.0e9, m0=m_e)
        for i in 1:6
            @test isapprox(jet_primal(hot_extra[i]), beam_extra.r[1, i]; atol=1.0e-11)
        end
        for (coord, exponent) in ((1, (1, 0)), (3, (0, 1)))
            rp = copy(r0_extra)
            rm = copy(r0_extra)
            rp[coord] += fd_step_extra
            rm[coord] -= fd_step_extra
            bp = Beam(reshape(rp, 1, 6); energy=3.0e9)
            bm = Beam(reshape(rm, 1, 6); energy=3.0e9)
            linepass!(deepcopy(extra_line), bp)
            linepass!(deepcopy(extra_line), bm)
            fd = (vec(bp.r[1, :]) .- vec(bm.r[1, :])) ./ (2fd_step_extra)
            for i in 1:6
                @test isapprox(jet_coefficient(hot_extra[i], exponent), fd[i];
                    rtol=1.0e-5, atol=1.0e-7)
            end
        end
    end
end
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

    base_line = [DRIFT(len=1.0), KQUAD(len=0.5, k1=0.0), DRIFT(len=1.0), KQUAD(len=0.5, k1=0.0)]
    params = [LatticeParameter(base_line, 2, :k1), LatticeParameter(base_line, 4, :k1)]
    function overlay_tracking(x)
        beam = Beam([1.0e-3 0.0 0.0 0.0 0.0 0.0])
        ADlinepass!(base_line, beam, params, x)
        return beam.r
    end
    function direct_tracking(x)
        line = [DRIFT(len=1.0), KQUAD(len=0.5, k1=x[1]), DRIFT(len=1.0), KQUAD(len=0.5, k1=x[2])]
        beam = Beam([1.0e-3 0.0 0.0 0.0 0.0 0.0])
        linepass!(line, beam)
        return beam.r
    end
    overlay_grad, overlay_primal = jacobian(set_runtime_activity(ForwardWithPrimal), Const(overlay_tracking), [0.8, -0.7])
    direct_grad, direct_primal = jacobian(set_runtime_activity(ForwardWithPrimal), Const(direct_tracking), [0.8, -0.7])
    @test all(isapprox.(overlay_primal, direct_primal; atol=1.0e-12))
    @test all(isapprox.(overlay_grad[1], direct_grad[1]; atol=1.0e-12))
    @test base_line[2].k1 == 0.0
    @test base_line[4].k1 == 0.0
end
end
