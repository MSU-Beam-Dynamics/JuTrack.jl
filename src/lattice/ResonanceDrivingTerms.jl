
function avedata(ring, dpp)
    twi = twissring(ring, dpp, 1)
    long = findall(x -> x.len > 0, ring)
    beta, alpha, gamma, mu, dp = array_optics(twi)

    # make the style consistent with AT
    beta_has_begin = [beta[end:end,:]; beta]
    alpha_has_begin = [alpha[end:end,:]; alpha]
    mu_has_begin = [[0.0 0.0]; mu]
    dp_has_begin = [dp[end:end,:]; dp]

    avebeta = beta_has_begin[1:end-1, :]
    avealpha = alpha_has_begin[1:end-1, :]
    avemu = mu_has_begin[1:end-1, :]
    avedisp = dp_has_begin[1:end-1, :]

    initial = long[1:end]
    final = [long[2:end]; long[1]]
    beta0 = avebeta[initial, :]
    alpha0 = avealpha[initial, :]
    mu0 = avemu[initial, :]
    disp0 = avedisp[initial, :]

    beta1 = avebeta[final, :]
    alpha1 = avealpha[final, :]
    mu1 = mu[long, :]
    disp1 = avedisp[final, :]

    L = zeros(length(long))
    for i in eachindex(long)
        L[i] = ring[long[i]].len
    end
    L2 = [L L]
    # println(size(beta0), size(beta1), size(alpha0), size(L2))
    avebeta[long, :] = betadrift(beta0, beta1, alpha0, L2)
    avemu[long, :] = 0.5 * (mu0 + mu1)
    avedisp[long, [1, 3]] = 0.5 * (disp0[:, [1, 3]] + disp1[:, [1, 3]])


    foc = Int[]
    for i in eachindex(long)
        if :PolynomB in fieldnames(typeof(ring[long[i]])) && ring[long[i]].PolynomB[2] != 0
            push!(foc, i)
        end
    end

    if length(foc) > 0
        K = zeros(size(L))
        K[foc] = [ring[long[foc[i]]].PolynomB[2] for i in eachindex(foc)]
        K2 = [K -K]

        avebeta[long[foc], :] = betafoc(beta1[foc, :], alpha0[foc, :], alpha1[foc, :], K2[foc, :], L2[foc, :])
        avedisp[long[foc], [1, 3]] = dispfoc(disp0[foc, [2, 4]], disp1[foc, [2, 4]], K2[foc, :], L2[foc, :])
    end
    avedisp[long, [2, 4]] = (disp1[:, [1, 3]] - disp0[:, [1, 3]]) ./ L2

    return avebeta, avemu, avedisp, beta, alpha, mu, dp
end

function betadrift(beta0, beta1, alpha0, L)
    gamma0 = (1 .+ alpha0 .* alpha0) ./ beta0
    avebeta = 0.5 * (beta0 + beta1) - gamma0 .* L .* L / 6
    return avebeta
end
function betafoc(beta1, alpha0, alpha1, K, L)
    gamma1 = (1 .+ alpha1 .* alpha1) ./ beta1
    avebeta = 0.5 * ((gamma1 + K.*beta1) .* L + alpha1 - alpha0) ./ K ./ L
    return avebeta
end

dispfoc(dispp0, dispp1, K, L) = (dispp0 .- dispp1) ./ K ./ L

# avebeta, avemu, avedisp = avedata(ring, 0.0)

function cal_FromindexDQSO(indDQSO, indexii)
    c = 0
    for i in eachindex(indDQSO)
        if indDQSO[i] < indexii
            c += 1
        else
            break
        end
    end
    return c + 1
end



struct DrivingTerms
    h21000::Vector{Float64}
    h30000::Vector{Float64}
    h10110::Vector{Float64}
    h10020::Vector{Float64}
    h10200::Vector{Float64}
    h11001::Vector{Float64}
    h00111::Vector{Float64}
    h20001::Vector{Float64}
    h00201::Vector{Float64}
    h10002::Vector{Float64}
    h10010::Vector{Float64}
    h10100::Vector{Float64}
    h22000::Vector{Float64}
    h11110::Vector{Float64}
    h00220::Vector{Float64}
    h31000::Vector{Float64}
    h40000::Vector{Float64}
    h20110::Vector{Float64}
    h11200::Vector{Float64}
    h20020::Vector{Float64}
    h20200::Vector{Float64}
    h00310::Vector{Float64}
    h00400::Vector{Float64}
    dnux_dJx::Float64
    dnux_dJy::Float64
    dnuy_dJy::Float64
end

struct ElementData
    betax::Float64
    betay::Float64
    rbetax::Float64
    rbetay::Float64
    betax2::Float64
    betay2::Float64
    phix::Float64
    phiy::Float64
    px::Vector{ComplexF64}
    py::Vector{ComplexF64}
    b2L::Float64
    b3L::Float64
    s::Float64
end

function computeDrivingTerms(s, betax, betay, phix, phiy, etax, Lista2L, Listb2L, Listb3L, Listb4L, tune, flags, nPeriods)
    NumElem = length(s)
    ed = Vector{ElementData}(undef, NumElem)
    ii = im
    # Initialize complex variables
    h11001 = h00111 = h20001 = h00201 = h10002 = h10010 = h10100 = 0.0 + 0.0im
    h21000 = h30000 = h10110 = h10020 = h10200 = 0.0+ 0.0im
    h22000 = h11110 = h00220 = h31000 = h40000 = 0.0 + 0.0im
    h20110 = h11200 = h20020 = h20200 = h00310 = h00400 = 0.0 + 0.0im

    periodicFactor = ones(9, 9)
    # periodicFactor = Matrix{ComplexF64}(undef, 9, 9)
    # if nPeriods != 1
    #     for i in 1:9
    #         for j in 1:9
    #             a1 = PIx2 * (tune[1] * (i-5) + tune[2] * (j-5))  # Adjusted indices for 1-based indexing
    #             a2 = a1 / nPeriods
    #             periodicFactor[i, j] = (exp(ii * a1) - 1) / (exp(ii * a2) - 1)
    #         end
    #     end
    # else
    #     for i in 1:9
    #         for j in 1:9
    #             periodicFactor[i, j] = 1.0 + 0im
    #         end
    #     end
    # end
    dnux_dJx = dnux_dJy = dnuy_dJy = 0.0
    # Loop through each element
    for nE in 1:NumElem
        a2L = Lista2L[nE]
        b2L = Listb2L[nE]
        b3L = Listb3L[nE]
        b4L = Listb4L[nE]

        betax1 = betax[nE]
        betay1 = betay[nE]
        phix1 = phix[nE]
        phiy1 = phiy[nE]
        etax1 = etax[nE]
        # Fill the element data
        ed[nE] = ElementData(betax1, betay1, sqrt(betax1), sqrt(betay1), betax1^2, betay1^2,
                            phix1, phiy1, [exp(ii * phix1 * j) for j in 1:5], [exp(ii * phiy1 * j) for j in 1:5],
                            b2L, b3L, s[nE])

        if abs(a2L) > 1e-6
            if flags.Coupling1
                h10010 += (a2L / 4.0) * ed[nE].rbetax * ed[nE].rbetay * ed[nE].px[1] / exp(ii * phiy1) * periodicFactor[5, 3]
                h10100 += (a2L / 4.0) * ed[nE].rbetax * ed[nE].rbetay * ed[nE].px[1] * exp(ii * phiy1) * periodicFactor[5, 5]
            end
        end

        if abs(b2L) > 1e-6 || abs(b3L) > 1e-6
            if flags.Chromatic1
                h11001::ComplexF64 += (b3L * betax1 * etax1 / 2 - b2L * betax1 / 4) * nPeriods
                h00111::ComplexF64 += (b2L * betay1 / 4 - b3L * betay1 * etax1 / 2) * nPeriods
                h20001::ComplexF64 += (b3L * betax1 * etax1 / 2 - b2L * betax1 / 4) / 2 * (exp(ii * 2 * phix1)) * periodicFactor[6, 4]
                h00201::ComplexF64 += (b2L * betay1 / 4 - b3L * betay1 * etax1 / 2) / 2 * (ed[nE].py[2]) * periodicFactor[4, 6]
                h10002::ComplexF64 += (b3L*ed[nE].rbetax*etax1^2-b2L*ed[nE].rbetax*etax1)/2*ed[nE].px[1]*periodicFactor[5, 4]
            end
        end

        if abs(ed[nE].b3L) > 1e-6
            if flags.Geometric1
                h21000::ComplexF64 += b3L * ed[nE].rbetax * betax1 / 8 * ed[nE].px[1] * periodicFactor[5, 4]
                h30000::ComplexF64 += b3L * ed[nE].rbetax * betax1 / 24 * ed[nE].px[3] * periodicFactor[7, 4]
                h10110::ComplexF64 += -b3L * ed[nE].rbetax * betay1 / 4 * ed[nE].px[1] * periodicFactor[5, 4]
                h10020::ComplexF64 += -b3L * ed[nE].rbetax * betay1 / 8 * ed[nE].px[1] * conj(ed[nE].py[2]) * periodicFactor[5, 2]
                h10200::ComplexF64 += -b3L * ed[nE].rbetax * betay1 / 8 * ed[nE].px[1] * (ed[nE].py[2]) * periodicFactor[5, 6]
            end
        end

        if abs(b4L) > 1e-6
            if flags.TuneShifts
                dnux_dJx += 3 * b4L * (betax1^2) / (8 * π) * nPeriods
                dnux_dJy -= 3 * b4L * betax1 * betay1 / (4 * π) * nPeriods
                dnuy_dJy += 3 * b4L * (betay1^2) / (8 * π) * nPeriods
            end
            if flags.Geometric2
                h22000::ComplexF64 += 3 * b4L * betax1^2 / 32 * nPeriods
                h11110::ComplexF64 -= 3 * b4L * betax1 * betay1 / 8 * nPeriods
                h00220::ComplexF64 += 3 * b4L * betay1^2 / 32 * nPeriods
                h31000::ComplexF64 += b4L * betax1^2 / 16 * exp(ii * phix1 * 2) * periodicFactor[6, 4]
                h40000::ComplexF64 += b4L * betax1^2 / 64 * exp(ii * phix1 * 4) * periodicFactor[8, 4]
                h20110::ComplexF64 -= 3 * b4L * betax1 * betay1 / 16 * exp(ii * phix1 * 2) * periodicFactor[6, 4]
                h11200::ComplexF64 -= 3 * b4L * betax1 * betay1 / 16 * exp(ii * phiy1 * 2) * periodicFactor[4, 6]
                h20020::ComplexF64 -= 3 * b4L * betax1 * betay1 / 32 * exp(ii * phix1 * 2) * conj(exp(ii * phiy1 * 2)) * periodicFactor[6, 2]
                h20200::ComplexF64 -= 3 * b4L * betax1 * betay1 / 32 * exp(ii * phix1 * 2) * exp(ii * phiy1 * 2) * periodicFactor[6, 6]
                h00310::ComplexF64 += b4L * betay1^2 / 16 * exp(ii * phiy1 * 2) * periodicFactor[4, 6]
                h00400::ComplexF64 += b4L * betay1^2 / 64 * exp(ii * phiy1 * 4) * periodicFactor[4, 8]
            end
            
        end
    end

    if flags.Geometric2 || flags.TuneShifts
        if nPeriods != 1
            println("Warning: Tune shifts and geometric terms are only valid for a single period.")
        end
        nux = tune[1]
        nuy = tune[2]
        for iE in 1:NumElem
            for jE in 1:NumElem
                if abs(ed[iE].b3L) > 1e-6 && abs(ed[jE].b3L) > 1e-6
                    if flags.TuneShifts
                        dnux_dJx += ed[iE].b3L * ed[jE].b3L / (-16 * pi) * (ed[iE].betax * ed[jE].betax)^1.5 *
                            (3 * cos(abs(ed[iE].phix - ed[jE].phix) - pi * nux) / sin(pi * nux) + cos(abs(3 * (ed[iE].phix - ed[jE].phix)) - 3 * pi * nux) / sin(3 * pi * nux))
                        dnux_dJy += ed[iE].b3L * ed[jE].b3L / (8 * pi) * sqrt(ed[iE].betax * ed[jE].betax) * ed[iE].betay *
                            (2 * ed[jE].betax * cos(abs(ed[iE].phix - ed[jE].phix) - pi * nux) / sin(pi * nux) -
                             ed[jE].betay * cos(abs(ed[iE].phix - ed[jE].phix) + 2 * abs(ed[iE].phiy - ed[jE].phiy) - pi * (nux + 2 * nuy)) / sin(pi * (nux + 2 * nuy)) +
                             ed[jE].betay * cos(abs(ed[iE].phix - ed[jE].phix) - 2 * abs(ed[iE].phiy - ed[jE].phiy) - pi * (nux - 2 * nuy)) / sin(pi * (nux - 2 * nuy)))
                        dnuy_dJy += ed[iE].b3L * ed[jE].b3L / (-16 * pi) * sqrt(ed[iE].betax * ed[jE].betax) * ed[iE].betay * ed[jE].betay *
                            (4 * cos(abs(ed[iE].phix - ed[jE].phix) - pi * nux) / sin(pi * nux) +
                             cos(abs(ed[iE].phix - ed[jE].phix) + 2 * abs(ed[iE].phiy - ed[jE].phiy) - pi * (nux + 2 * nuy)) / sin(pi * (nux + 2 * nuy)) +
                             cos(abs(ed[iE].phix - ed[jE].phix) - 2 * abs(ed[iE].phiy - ed[jE].phiy) - pi * (nux - 2 * nuy)) / sin(pi * (nux - 2 * nuy)))
                    end
    
                    if flags.Geometric2
                        termSign = sign(ed[iE].s - ed[jE].s)
                        if termSign != 0
                            h22000::ComplexF64 += (1.0 / 64.0) * termSign * ii * ed[iE].b3L * ed[jE].b3L *
                            sqrt(ed[iE].betax) * sqrt(ed[jE].betax) * ed[iE].betax * ed[jE].betax *
                            (ed[iE].px[3] * conj(ed[jE].px[3]) + 3 * ed[iE].px[1] * conj(ed[jE].px[1]))
                            h31000::ComplexF64 += (1.0 / 32.0) * termSign * ii * ed[iE].b3L * ed[jE].b3L *
                            sqrt(ed[iE].betax) * sqrt(ed[jE].betax) * ed[iE].betax * ed[jE].betax *
                            ed[iE].px[3] * conj(ed[jE].px[1])
                        
                            t1 = conj(ed[iE].px[1]) * ed[jE].px[1]
                            t2 = ed[iE].px[1] * conj(ed[jE].px[1])
                            h11110::ComplexF64 += (1.0 / 16.0) * termSign * ii * ed[iE].b3L * ed[jE].b3L *
                            sqrt(ed[iE].betax) * sqrt(ed[jE].betax) * ed[iE].betay *
                            (ed[jE].betax * (t1 - conj(t1)) + ed[jE].betay * ed[iE].py[2] * conj(ed[jE].py[2]) * (conj(t1) + t1))
    
                            t1 = exp(-ii * (ed[iE].phix - ed[jE].phix))
                            t2 = conj(t1)
                            h11200::ComplexF64 += (1.0 / 32.0) * termSign * ii * ed[iE].b3L * ed[jE].b3L *
                            sqrt(ed[iE].betax) * sqrt(ed[jE].betax) * ed[iE].betay * exp(ii * (2 * ed[iE].phiy)) *
                            (ed[jE].betax * (t1 - t2) + 2 * ed[jE].betay * (t2 + t1))
                        
                            h40000::ComplexF64 += (1.0 / 64.0) * termSign * ii * ed[iE].b3L * ed[jE].b3L *
                            sqrt(ed[iE].betax) * sqrt(ed[jE].betax) * ed[iE].betax * ed[jE].betax *
                            ed[iE].px[3] * ed[jE].px[1]
                            h20020::ComplexF64 += (1.0 / 64.0) * termSign * ii * ed[iE].b3L * ed[jE].b3L *
                            sqrt(ed[iE].betax) * sqrt(ed[jE].betax) * ed[iE].betay *
                            (ed[jE].betax * conj(ed[iE].px[1] * ed[iE].py[2]) * ed[jE].px[3] -
                             (ed[jE].betax + 4 * ed[jE].betay) * ed[iE].px[1] * ed[jE].px[1] * conj(ed[iE].py[2]))

                            h20110::ComplexF64 += (1.0 / 32.0) * termSign * ii * ed[iE].b3L * ed[jE].b3L *
                             sqrt(ed[iE].betax) * sqrt(ed[jE].betax) * ed[iE].betay *
                             (ed[jE].betax * (conj(ed[iE].px[1]) * ed[jE].px[3] - ed[iE].px[1] * ed[jE].px[1]) +
                              2 * ed[jE].betay * ed[iE].px[1] * ed[jE].px[1] * ed[iE].py[2] * conj(ed[jE].py[2]))
     
                            # Calculate h20200
                            h20200::ComplexF64 += (1.0 / 64.0) * termSign * ii * ed[iE].b3L * ed[jE].b3L *
                             sqrt(ed[iE].betax) * sqrt(ed[jE].betax) * ed[iE].betay *
                             (ed[jE].betax * conj(ed[iE].px[1]) * ed[jE].px[3] * ed[iE].py[2] -
                              (ed[jE].betax - 4 * ed[jE].betay) * ed[iE].px[1] * ed[jE].px[1] * ed[iE].py[2])
     
                            # Calculate h00220
                            h00220::ComplexF64 += (1.0 / 64.0) * termSign * ii * ed[iE].b3L * ed[jE].b3L *
                             sqrt(ed[iE].betax) * sqrt(ed[jE].betax) * ed[iE].betay * ed[jE].betay *
                             (ed[iE].px[1] * ed[iE].py[2] * conj(ed[jE].px[1] * ed[jE].py[2]) +
                              4 * ed[iE].px[1] * conj(ed[jE].px[1]) -
                              conj(ed[iE].px[1] * ed[jE].py[2]) * ed[jE].px[1] * ed[iE].py[2])
     
                            # Calculate h00310
                            h00310::ComplexF64 += (1.0 / 32.0) * termSign * ii * ed[iE].b3L * ed[jE].b3L *
                             sqrt(ed[iE].betax) * sqrt(ed[jE].betax) * ed[iE].betay * ed[jE].betay * ed[iE].py[2] *
                             (ed[iE].px[1] * conj(ed[jE].px[1]) - conj(ed[iE].px[1]) * ed[jE].px[1])
     
                            # Calculate h00400
                            h00400::ComplexF64 += (1.0 / 64.0) * termSign * ii * ed[iE].b3L * ed[jE].b3L *
                             sqrt(ed[iE].betax) * sqrt(ed[jE].betax) * ed[iE].betay * ed[jE].betay *
                             ed[iE].px[1] * conj(ed[jE].px[1]) * ed[iE].py[2] * ed[jE].py[2]
                        end
                    end
                end
            end
        end
    end
    # Create and return the structure of driving terms
    DrivingTerms([abs(h21000), real(h21000), imag(h21000)], [abs(h30000), real(h30000), imag(h30000)],
                  [abs(h10110), real(h10110), imag(h10110)], [abs(h10020), real(h10020), imag(h10020)],
                  [abs(h10200), real(h10200), imag(h10200)], [abs(h11001), real(h11001), imag(h11001)],
                  [abs(h00111), real(h00111), imag(h00111)], [abs(h20001), real(h20001), imag(h20001)],
                  [abs(h00201), real(h00201), imag(h00201)], [abs(h10002), real(h10002), imag(h10002)],
                  [abs(h10010), real(h10010), imag(h10010)], [abs(h10100), real(h10100), imag(h10100)],
                  [abs(h22000), real(h22000), imag(h22000)], [abs(h11110), real(h11110), imag(h11110)],
                  [abs(h00220), real(h00220), imag(h00220)], [abs(h31000), real(h31000), imag(h31000)],
                  [abs(h40000), real(h40000), imag(h40000)], [abs(h20110), real(h20110), imag(h20110)],
                  [abs(h11200), real(h11200), imag(h11200)], [abs(h20020), real(h20020), imag(h20020)],
                  [abs(h20200), real(h20200), imag(h20200)], [abs(h00310), real(h00310), imag(h00310)],
                  [abs(h00400), real(h00400), imag(h00400)], dnux_dJx, dnux_dJy, dnuy_dJy)
end


struct RDTflags
    Geometric1::Bool
    Geometric2::Bool
    Chromatic1::Bool
    Coupling1::Bool
    TuneShifts::Bool
end
function juliaRDT(s, betax, betay, phix, phiy, etax, Lista2L, Listb2L, Listb3L, Listb4L, Tunex, Tuney, NumElem, Geometric1, Geometric2, Chromatic1, Coupling1, TuneShifts, nPeriods)    
    flags = RDTflags(Geometric1, Geometric2, Chromatic1, Coupling1, TuneShifts)
    tune = [Tunex, Tuney]
    d = computeDrivingTerms(s, betax, betay, phix, phiy, etax, Lista2L, Listb2L, Listb3L, Listb4L, tune, flags, nPeriods)
    return d
end

function computeRDT(ring, index; chromatic=false, coupling=false, geometric1=false, geometric2=false, tuneshifts=false)
    # Compute Hamiltonian resonance driving terms (RDTs)
    # ring: lattice sequence
    # index: index of the element to compute the RDTs
    # This code is based on the elegant function computeDrivingTerms and the AT function computeRDT.
    # The function is modified to work with JuTrack. The formulas are not changed.
    # The results have been validated against elegant and AT.
    if chromatic == false && coupling == false && geometric1 == false && geometric2 == false && tuneshifts == false
        chromatic, coupling, geometric1, geometric2, tuneshifts = true, true, true, true, true
    end

    indB = findelem(ring, SBEND)
    indQ = findelem(ring, KQUAD)
    if length(indQ) == 0
        indQ = findelem(ring, QUAD)
    end
    indS = findelem(ring, KSEXT)
    indO = findelem(ring, KOCT)
    indDQSO = sort(union(indB, indQ, indS, indO))

    AVEBETA, AVEMU, AVEDISP, beta, alpha, mu, dp = avedata(ring, 0.0)  
    sIndex = spos(ring, indDQSO.-1)
    s = spos(ring)
    # make the style consistent with AT
    s = [0.0; s[1:end-1]]
    sEnd = spos(ring, [length(ring)])

    betax, betay = AVEBETA[indDQSO, 1], AVEBETA[indDQSO, 2]
    etax = AVEDISP[indDQSO, 1]
    phix = AVEMU[indDQSO, 1]
    phiy = AVEMU[indDQSO, 2]

    PolyA2 = [ring[indDQSO[i]].PolynomA[2] for i in eachindex(indDQSO)]
    PolyB2 = [ring[indDQSO[i]].PolynomB[2] for i in eachindex(indDQSO)]
    PolyB3 = [ring[indDQSO[i]].PolynomB[3] for i in eachindex(indDQSO)] ./ 2 # AT style
    PolyB4 = [ring[indDQSO[i]].PolynomB[4] for i in eachindex(indDQSO)] ./ 6 # AT style
    len_list = [ring[indDQSO[i]].len for i in eachindex(indDQSO)]
    a2L = PolyA2 .* len_list
    b2L = PolyB2 .* len_list
    b3L = PolyB3 .* len_list
    b4L = PolyB4 .* len_list

    Mux = mu[end-1, 1]
    Muy = mu[end-1, 2]
    Tunex = Mux / (2.0 * pi)
    Tuney = Muy / (2.0 * pi)
    nElem = length(indDQSO)

    dlist = Vector{DrivingTerms}(undef, length(index))
    for ii in eachindex(index)
        FromindexDQSO = cal_FromindexDQSO(indDQSO, index[ii])
        betax_Fromindex = [betax[FromindexDQSO:end]; betax[1:FromindexDQSO-1]]
        betay_Fromindex = [betay[FromindexDQSO:end]; betay[1:FromindexDQSO-1]]
        etax_Fromindex = [etax[FromindexDQSO:end]; etax[1:FromindexDQSO-1]]
        phix_Fromindex = [phix[FromindexDQSO:end] .- AVEMU[index[ii], 1]; phix[1:FromindexDQSO-1] .+ Mux .- AVEMU[index[ii], 1]]
        phiy_Fromindex = [phiy[FromindexDQSO:end] .- AVEMU[index[ii], 2]; phiy[1:FromindexDQSO-1] .+ Muy .- AVEMU[index[ii], 2]]
        s_Fromindex = [sIndex[FromindexDQSO:end] .- s[index[ii]]; sIndex[1:FromindexDQSO-1] .+ sEnd .- s[index[ii]]]
        a2L_Fromindex = [a2L[FromindexDQSO:end]; a2L[1:FromindexDQSO-1]]
        b2L_Fromindex = [b2L[FromindexDQSO:end]; b2L[1:FromindexDQSO-1]]
        b3L_Fromindex = [b3L[FromindexDQSO:end]; b3L[1:FromindexDQSO-1]]
        b4L_Fromindex = [b4L[FromindexDQSO:end]; b4L[1:FromindexDQSO-1]]
        d = juliaRDT(s_Fromindex, betax_Fromindex, betay_Fromindex,
                                        phix_Fromindex, phiy_Fromindex, etax_Fromindex,
                                        a2L_Fromindex, b2L_Fromindex,
                                        b3L_Fromindex, b4L_Fromindex,
                                        Tunex, Tuney, nElem,
                                        geometric1, geometric2, chromatic, coupling, tuneshifts, 1)
        dlist[ii] = d
    end
    s1 = spos(ring, index)
    return dlist, s1
end

function get_polynom(ele, AB, n)
    if AB == 1
        result = get_polynom_value(ele.PolynomA, n)
    elseif AB == 2
        result = get_polynom_value(ele.PolynomB, n)
    else
        result = 0.0
    end
    return result
end
function get_polynom_value(poly, n)
    return poly[n]
end

function ADavedata(ring, dpp, changed_ids, changed_elems)
    refpts = [i for i in 1:length(ring)]
    twi = ADtwissring(ring, dpp, 1, refpts, changed_ids, changed_elems)
    nlong = 0
    for i in eachindex(ring)
        if get_len(ring[i]) > 0
            nlong += 1
        end
    end
    long = zeros(Int, nlong)
    nlong = 0
    for i in eachindex(ring)
        if get_len(ring[i]) > 0
            nlong += 1
            long[nlong] = i
        end
    end
    beta, alpha, gamma, mu, dp = array_optics(twi)

    # make the style consistent with AT
    beta_has_begin = [beta[end:end,:]; beta]
    alpha_has_begin = [alpha[end:end,:]; alpha]
    mu_has_begin = [[0.0 0.0]; mu]
    dp_has_begin = [dp[end:end,:]; dp]

    avebeta = beta_has_begin[1:end-1, :]
    avealpha = alpha_has_begin[1:end-1, :]
    avemu = mu_has_begin[1:end-1, :]
    avedisp = dp_has_begin[1:end-1, :]

    initial = long[1:end]
    final = [long[2:end]; long[1]]
    beta0 = avebeta[initial, :]
    alpha0 = avealpha[initial, :]
    mu0 = avemu[initial, :]
    disp0 = avedisp[initial, :]

    beta1 = avebeta[final, :]
    alpha1 = avealpha[final, :]
    mu1 = mu[long, :]
    disp1 = avedisp[final, :]

    L = zeros(length(long))
    for i in eachindex(long)
        if long[i] in changed_ids
            for j in eachindex(changed_ids)
                if long[i] == changed_ids[j]
                    L[i] = changed_elems[j].len
                end
                break
            end
            # L[i] = changed_elems[findfirst(x -> x == long[i], changed_ids)].len
        else
            L[i] = get_len(ring[long[i]])
        end
    end
    L2 = [L L]
    # println(size(beta0), size(beta1), size(alpha0), size(L2))
    avebeta[long, :] = betadrift(beta0, beta1, alpha0, L2)
    avemu[long, :] = 0.5 * (mu0 + mu1)
    avedisp[long, [1, 3]] = 0.5 * (disp0[:, [1, 3]] + disp1[:, [1, 3]])

    nfoc = 0
    for i in eachindex(long)
        if ring[long[i]] isa KQUAD || ring[long[i]] isa KSEXT || ring[long[i]] isa KOCT || 
            ring[long[i]] isa SBEND || ring[long[i]] isa QUAD || ring[long[i]] isa thinMULTIPOLE
            if get_polynom(ring[long[i]], 2, 2) != 0
                nfoc += 1
            end
        end
    end
    foc = zeros(Int, nfoc)
    nfoc = 0
    for i in eachindex(long)
        if ring[long[i]] isa KQUAD || ring[long[i]] isa KSEXT || ring[long[i]] isa KOCT || 
            ring[long[i]] isa SBEND || ring[long[i]] isa QUAD || ring[long[i]] isa thinMULTIPOLE
            if get_polynom(ring[long[i]], 2, 2) != 0
                nfoc += 1
                foc[nfoc] = i
            end
        end
    end
    
    K = zeros(length(long))
    if length(foc) > 0
        for i in eachindex(foc)
            if long[foc[i]] in changed_ids
                for j in eachindex(changed_ids)
                    if long[foc[i]] == changed_ids[j]
                        K[foc[i]] = changed_elems[j].PolynomB[2]
                    end
                    break
                end
                ###### K[foc[i]] = changed_elems[findfirst(x -> x == long[foc[i]], changed_ids)].PolynomB[2]
            else
                K[foc[i]] = get_polynom(ring[long[foc[i]]], 2, 2) 
            end
        end
        K2 = [K -K]

        avebeta[long[foc], :] = betafoc(beta1[foc, :], alpha0[foc, :], alpha1[foc, :], K2[foc, :], L2[foc, :])
        avedisp[long[foc], [1, 3]] = dispfoc(disp0[foc, [2, 4]], disp1[foc, [2, 4]], K2[foc, :], L2[foc, :])
    end
    avedisp[long, [2, 4]] = (disp1[:, [1, 3]] - disp0[:, [1, 3]]) ./ L2
    return avebeta, avemu, avedisp, beta, alpha, mu, dp
end

function findBQSO(ring::Vector)
    c = 0
    for i in eachindex(ring)
        if typeof(ring[i]) == SBEND || typeof(ring[i]) == KQUAD || typeof(ring[i]) == KSEXT || 
            typeof(ring[i]) == KOCT || typeof(ring[i]) == QUAD || typeof(ring[i]) == thinMULTIPOLE
            c += 1
        end
    end
    indDQSO = zeros(Int, c)
    c = 0
    for i in eachindex(ring)
        if typeof(ring[i]) == SBEND || typeof(ring[i]) == KQUAD || typeof(ring[i]) == KSEXT || 
            typeof(ring[i]) == KOCT || typeof(ring[i]) == QUAD || typeof(ring[i]) == thinMULTIPOLE
            c += 1
            indDQSO[c] = i
        end
    end
    return indDQSO
end

function ADcomputeRDT(ring, index, changed_ids, changed_elems; chromatic=true, coupling=true, geometric1=true, geometric2=true, tuneshifts=true)
    # Compute Hamiltonian resonance driving terms (RDTs)
    # ring: lattice sequence
    # index: index of the element to compute the RDTs
    # This code is based on the elegant function computeDrivingTerms and the AT function computeRDT.
    # The function is modified to work with JuTrack. The formulas are not changed.
    # The results have been validated against elegant and AT.

    # if chromatic == false && coupling == false && geometric1 == false && geometric2 == false && tuneshifts == false
    #     chromatic, coupling, geometric1, geometric2, tuneshifts = true, true, true, true, true
    # end

    indDQSO = findBQSO(ring)
    AVEBETA, AVEMU, AVEDISP, beta, alpha, mu, dp = ADavedata(ring, 0.0, changed_ids, changed_elems)  
    sIndex = spos(ring, indDQSO.-1)
    s = spos(ring)
    # make the style consistent with AT
    sEnd = s[end]
    s = [0.0; s[1:end-1]]
     

    betax, betay = AVEBETA[indDQSO, 1], AVEBETA[indDQSO, 2]
    etax = AVEDISP[indDQSO, 1]
    phix = AVEMU[indDQSO, 1]
    phiy = AVEMU[indDQSO, 2]

    a2L = zeros(length(indDQSO))
    b2L = zeros(length(indDQSO))
    b3L = zeros(length(indDQSO))
    b4L = zeros(length(indDQSO))
    for i in eachindex(indDQSO)
        PolyA2 = 0.0
        PolyB2 = 0.0
        PolyB3 = 0.0
        PolyB4 = 0.0
        len_list = 0.0
        if indDQSO[i] in changed_ids
            for j in eachindex(changed_ids)
                if indDQSO[i] == changed_ids[j]
                    PolyA2 = get_polynom(changed_elems[j], 1, 2)
                    PolyB2 = get_polynom(changed_elems[j], 2, 2)
                    PolyB3 = get_polynom(changed_elems[j], 2, 3) / 2.0
                    PolyB4 = get_polynom(changed_elems[j], 2, 4) / 6.0
                    len_list = changed_elems[j].len
                end
            end
        else
            PolyA2 = get_polynom(ring[indDQSO[i]], 1, 2)
            PolyB2 = get_polynom(ring[indDQSO[i]], 2, 2)
            PolyB3 = get_polynom(ring[indDQSO[i]], 2, 3) / 2.0
            PolyB4 = get_polynom(ring[indDQSO[i]], 2, 4) / 6.0
            len_list = get_len(ring[indDQSO[i]])
        end
        a2L[i] = PolyA2 * len_list
        b2L[i] = PolyB2 * len_list
        b3L[i] = PolyB3 * len_list
        b4L[i] = PolyB4 * len_list
    end
    Mux = mu[end-1, 1]
    Muy = mu[end-1, 2]
    Tunex = Mux / (2.0 * pi)
    Tuney = Muy / (2.0 * pi)
    nElem = length(indDQSO)

    dlist = Vector{DrivingTerms}(undef, length(index))
    for i in eachindex(index)
        FromindexDQSO = cal_FromindexDQSO(indDQSO, index[i])
        betax_Fromindex = [betax[FromindexDQSO:end]; betax[1:FromindexDQSO-1]]
        betay_Fromindex = [betay[FromindexDQSO:end]; betay[1:FromindexDQSO-1]]
        etax_Fromindex = [etax[FromindexDQSO:end]; etax[1:FromindexDQSO-1]]
        phix_Fromindex = [phix[FromindexDQSO:end] .- AVEMU[index[i], 1]; phix[1:FromindexDQSO-1] .+ Mux .- AVEMU[index[i], 1]]
        phiy_Fromindex = [phiy[FromindexDQSO:end] .- AVEMU[index[i], 2]; phiy[1:FromindexDQSO-1] .+ Muy .- AVEMU[index[i], 2]]
        s_Fromindex = [sIndex[FromindexDQSO:end] .- s[index[i]]; sIndex[1:FromindexDQSO-1] .+ sEnd .- s[index[i]]]
        a2L_Fromindex = [a2L[FromindexDQSO:end]; a2L[1:FromindexDQSO-1]]
        b2L_Fromindex = [b2L[FromindexDQSO:end]; b2L[1:FromindexDQSO-1]]
        b3L_Fromindex = [b3L[FromindexDQSO:end]; b3L[1:FromindexDQSO-1]]
        b4L_Fromindex = [b4L[FromindexDQSO:end]; b4L[1:FromindexDQSO-1]]
        d = juliaRDT(s_Fromindex, betax_Fromindex, betay_Fromindex,
                                        phix_Fromindex, phiy_Fromindex, etax_Fromindex,
                                        a2L_Fromindex, b2L_Fromindex,
                                        b3L_Fromindex, b4L_Fromindex,
                                        Tunex, Tuney, nElem,
                                        geometric1, geometric2, chromatic, coupling, tuneshifts, 1)
        dlist[i] = d
    end
    s1 = spos(ring, index)
    return dlist, s1
end

function getRDTvalues(dlist, term::Symbol)
    # e.g., val_21000 = getRDTvalues(dlist, :h21000)
    values = zeros(ComplexF64, length(dlist))
    for i in eachindex(dlist)
        real_part = getfield(dlist[i], term)[2] 
        imag_part = getfield(dlist[i], term)[3] 
        values[i] = real_part + imag_part * im
    end
    return values
end


