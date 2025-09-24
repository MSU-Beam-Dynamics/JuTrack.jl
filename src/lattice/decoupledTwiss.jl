# under construction
struct DecoupledTwiss{T} <: AbstractTwiss
	betax::T
	betay::T
	alphax::T
	alphay::T
	gammax::T
	gammay::T
	dx::T
	dpx::T
	dy::T
	dpy::T
	mux::T
	muy::T
end
function DecoupledTwiss(betax::Float64, betay::Float64,
                        alphax::Float64, alphay::Float64;
                        dx::Float64=0.0, dpx::Float64=0.0,
                        dy::Float64=0.0, dpy::Float64=0.0,
                        mux::Float64=0.0, muy::Float64=0.0)

    gammax=(1.0+alphax*alphax)/betax
    gammay=(1.0+alphay*alphay)/betay
    DecoupledTwiss{Float64}(betax, betay, alphax, alphay, gammax, gammay, dx, dpx, dy, dpy, mux, muy)
end
function DecoupledTwiss(betax::DTPSAD{N, T}, betay::DTPSAD{N, T},
                        alphax::DTPSAD{N, T}, alphay::DTPSAD{N, T};
                        dx::DTPSAD{N, T}=zero(DTPSAD{N, T}),
                        dpx::DTPSAD{N, T}=zero(DTPSAD{N, T}),
                        dy::DTPSAD{N, T}=zero(DTPSAD{N, T}),
                        dpy::DTPSAD{N, T}=zero(DTPSAD{N, T}),
                        mux::DTPSAD{N, T}=zero(DTPSAD{N, T}),
                        muy::DTPSAD{N, T}=zero(DTPSAD{N, T})) where {N, T}

    gammax=(1.0+alphax*alphax)/betax
    gammay=(1.0+alphay*alphay)/betay
    DecoupledTwiss{DTPSAD{N, T}}(betax, betay, alphax, alphay, gammax, gammay, dx, dpx, dy, dpy, mux, muy)
end

function uncoupledtwissring(seq::Vector{<:AbstractElement{Float64}}, dp::Float64=0.0;
                        order::Int=0, E0::Float64=0.0, m0::Float64=0.0)
    refpts = [i for i in 1:length(seq)]
    orb, M0 = find_closed_orbit(seq, dp, mass=m0, energy=E0)
    orb1, M1 = find_closed_orbit(seq, dp+1e-8, mass=m0, energy=E0)
    rout1 = linepass!(seq, Beam([orb[1] orb[2] orb[3] orb[4] orb[5] orb[6]], energy=E0, mass=m0), [i for i in 1:length(seq)])
    rout2 = linepass!(seq, Beam([orb1[1] orb1[2] orb1[3] orb1[4] orb1[5] orb1[6]], energy=E0, mass=m0), [i for i in 1:length(seq)])
    disp = zeros(Float64, length(seq), 4)
    for i in 1:length(seq)
        disp[i, 1] = (rout2[i][1] - rout1[i][1]) / 1e-8
        disp[i, 2] = (rout2[i][2] - rout1[i][2]) / 1e-8
        disp[i, 3] = (rout2[i][3] - rout1[i][3]) / 1e-8
        disp[i, 4] = (rout2[i][4] - rout1[i][4]) / 1e-8
    end
    DX = disp[:, 1]
    DPX = disp[:, 2]
    DY = disp[:, 3]
    DPY = disp[:, 4]
    if order == 0
		MS = fastfindm66_refpts(seq, dp, refpts, E0=E0, m0=m0, orb=orb)
        M44 = fastfindm66(seq, dp, E0=E0, m0=m0, orb=orb)[1:4, 1:4]
	else
		MS = findm66_refpts(seq, dp, order, refpts, E0=E0, m0=m0, orb=orb)
        M44 = findm66(seq, dp, order, refpts, E0=E0, m0=m0, orb=orb)[1:4, 1:4]
	end
    for i in 2:length(refpts)
        MS[:,:,i] = MS[:,:,i] * MS[:,:,i-1]
    end
    cos_mu_x = (M44[1,1] + M44[2,2]) / 2
    cos_mu_y = (M44[3,3] + M44[4,4]) / 2

    sin_mu_x = sign(M44[1,2]) * sqrt(-M44[1,2]*M44[2,1] - (M44[1,1] - M44[2,2])^2 / 4)
    sin_mu_y = sign(M44[3,4]) * sqrt(-M44[3,4]*M44[4,3] - (M44[3,3] - M44[4,4])^2 / 4)

    ax = (M44[1,1] - M44[2,2]) / (2 * sin_mu_x)
    ay = (M44[3,3] - M44[4,4]) / (2 * sin_mu_y)

    bx = M44[1,2] / sin_mu_x
    by = M44[3,4] / sin_mu_y

    MS11 = vec(MS[1,1,:]); MS12 = vec(MS[1,2,:])
    MS21 = vec(MS[2,1,:]); MS22 = vec(MS[2,2,:])
    MS33 = vec(MS[3,3,:]); MS34 = vec(MS[3,4,:])
    MS43 = vec(MS[4,3,:]); MS44 = vec(MS[4,4,:])

    BX = ((MS11 .* bx .- MS12 .* ax).^2 .+ MS12.^2) ./ bx
    BY = ((MS33 .* by .- MS34 .* ay).^2 .+ MS34.^2) ./ by

    AX = -(((MS11 .* bx .- MS12 .* ax) .* (MS21 .* bx .- MS22 .* ax)) .+ MS12 .* MS22) ./ bx
    AY = -(((MS33 .* by .- MS34 .* ay) .* (MS43 .* by .- MS44 .* ay)) .+ MS34 .* MS44) ./ by

    GX = (1 .+ AX.^2) ./ BX
    GY = (1 .+ AY.^2) ./ BY

    MX = atan.(MS12, MS11 .* bx .- MS12 .* ax)
    MY = atan.(MS34, MS33 .* by .- MS34 .* ay)
    tune = mod.(atan.([sin_mu_x, sin_mu_y], [cos_mu_x, cos_mu_y]) ./ (2π), 1.0)

    NR = length(refpts)

    twiss = Vector{DecoupledTwiss{Float64}}(undef, NR)
    for i in 1:NR
        twiss[i] = DecoupledTwiss{Float64}(
            BX[i], BY[i], AX[i], AY[i], GX[i], GY[i],
            DX[i], DPX[i], DY[i], DPY[i], MX[i], MY[i]
        )
    end
end


function uncoupledtwissring(seq::Vector{<:AbstractElement{DTPSAD{N, T}}}, dp::Float64=0.0;
                        order::Int=0, E0::Float64=0.0, m0::Float64=0.0) where {N, T}
    refpts = [i for i in 1:length(seq)]
    orb, M0 = find_closed_orbit(seq, dp, mass=m0, energy=E0)
    orb1, M1 = find_closed_orbit(seq, dp+1e-8, mass=m0, energy=E0)
    rout1 = linepass!(seq, Beam([orb[1] orb[2] orb[3] orb[4] orb[5] orb[6]], energy=E0, mass=m0), [i for i in 1:length(seq)])
    rout2 = linepass!(seq, Beam([orb1[1] orb1[2] orb1[3] orb1[4] orb1[5] orb1[6]], energy=E0, mass=m0), [i for i in 1:length(seq)])
    disp = zeros(DTPSAD{N, T}, length(seq), 4)
    for i in 1:length(seq)
        disp[i, 1] = (rout2[i][1] - rout1[i][1]) / 1e-8
        disp[i, 2] = (rout2[i][2] - rout1[i][2]) / 1e-8
        disp[i, 3] = (rout2[i][3] - rout1[i][3]) / 1e-8
        disp[i, 4] = (rout2[i][4] - rout1[i][4]) / 1e-8
    end
    DX = disp[:, 1]
    DPX = disp[:, 2]
    DY = disp[:, 3]
    DPY = disp[:, 4]
    if order == 0
		MS = fastfindm66_refpts(seq, dp, refpts, E0=E0, m0=m0, orb=orb)
        M44 = fastfindm66(seq, dp, E0=E0, m0=m0, orb=orb)[1:4, 1:4]
	else
		MS = findm66_refpts(seq, dp, order, refpts, E0=E0, m0=m0, orb=orb)
        M44 = findm66(seq, dp, order, refpts, E0=E0, m0=m0, orb=orb)[1:4, 1:4]
	end
    for i in 2:length(refpts)
        MS[:,:,i] = MS[:,:,i] * MS[:,:,i-1]
    end
    cos_mu_x = (M44[1,1] + M44[2,2]) / 2
    cos_mu_y = (M44[3,3] + M44[4,4]) / 2

    sin_mu_x = sign(M44[1,2]) * sqrt(-M44[1,2]*M44[2,1] - (M44[1,1] - M44[2,2])^2 / 4)
    sin_mu_y = sign(M44[3,4]) * sqrt(-M44[3,4]*M44[4,3] - (M44[3,3] - M44[4,4])^2 / 4)

    ax = (M44[1,1] - M44[2,2]) / (2 * sin_mu_x)
    ay = (M44[3,3] - M44[4,4]) / (2 * sin_mu_y)

    bx = M44[1,2] / sin_mu_x
    by = M44[3,4] / sin_mu_y

    MS11 = vec(MS[1,1,:]); MS12 = vec(MS[1,2,:])
    MS21 = vec(MS[2,1,:]); MS22 = vec(MS[2,2,:])
    MS33 = vec(MS[3,3,:]); MS34 = vec(MS[3,4,:])
    MS43 = vec(MS[4,3,:]); MS44 = vec(MS[4,4,:])

    BX = ((MS11 .* bx .- MS12 .* ax).^2 .+ MS12.^2) ./ bx
    BY = ((MS33 .* by .- MS34 .* ay).^2 .+ MS34.^2) ./ by

    AX = -(((MS11 .* bx .- MS12 .* ax) .* (MS21 .* bx .- MS22 .* ax)) .+ MS12 .* MS22) ./ bx
    AY = -(((MS33 .* by .- MS34 .* ay) .* (MS43 .* by .- MS44 .* ay)) .+ MS34 .* MS44) ./ by

    GX = (1 .+ AX.^2) ./ BX
    GY = (1 .+ AY.^2) ./ BY

    MX = atan.(MS12, MS11 .* bx .- MS12 .* ax)
    MY = atan.(MS34, MS33 .* by .- MS34 .* ay)
    tune = mod.(atan.([sin_mu_x, sin_mu_y], [cos_mu_x, cos_mu_y]) ./ (2π), 1.0)

    NR = length(refpts)

    twiss = Vector{DecoupledTwiss{Float64}}(undef, NR)
    for i in 1:NR
        twiss[i] = DecoupledTwiss{Float64}(
            BX[i], BY[i], AX[i], AY[i], GX[i], GY[i],
            DX[i], DPX[i], DY[i], DPY[i], MX[i], MY[i]
        )
    end
end
