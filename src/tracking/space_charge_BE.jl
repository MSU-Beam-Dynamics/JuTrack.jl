#=mutable struct SC_lens <: AbstractElement
    optics::AbstractOptics4D
    ds::Float64
    nSC::Int64
    turns::Int64
    sigma_xi_SC::Matrix{Float64}
    sigma_yi_SC::Matrix{Float64}
    sigma_x_SC::Vector{Float64}
    sigma_y_SC::Vector{Float64}

    function SC_lens(optics, ds, nSC, turns)
        sigma_xi_SC = zeros(Float64, turns, nSC)
        sigma_yi_SC = zeros(Float64, turns, nSC)

        sigma_x_SC = zeros(Float64, nSC)
        sigma_y_SC = zeros(Float64, nSC)
        new(optics, ds, nSC, turns, sigma_xi_SC, sigma_yi_SC, sigma_x_SC, sigma_y_SC)
    end
end
=#

function track_SC!(rin, σx, σy, σz, temp1, bbSC::Beam, factorSC::Float64)

    dpx_dz = zeros(length(bbSC.nmacro))
    dpy_dz = zeros(length(bbSC.nmacro))
   
    fieldvec_thread=[MVector{3}(0.0, 0.0, 0.0)  for j = 1:Threads.nthreads()]
    @inbounds Threads.@threads :static for j in eachindex(bbSC.nmacro)
        
        r6 = @view rin[(j-1)*6+1:j*6]

        Bassetti_Erskine!(fieldvec_thread[Threads.threadid()], r6[1], r6[3], σx, σy)

        # gaussian function for particle distribution, lambda_z
        temp1[j] = 1.0/sqrt(2*π)/σz*exp((-0.5)*r6[5]^2/σz^2)

        dpx_dz[j] = factorSC*temp1[j] *fieldvec_thread[Threads.threadid()][1]* r6[5]/σz^2
        dpy_dz[j] = factorSC*temp1[j] *fieldvec_thread[Threads.threadid()][2]* r6[5]/σz^2


        r6[2] += factorSC* temp1[j]* fieldvec_thread[Threads.threadid()][1]
        r6[4] += factorSC* temp1[j]* fieldvec_thread[Threads.threadid()][2]
        r6[6] += dpx_dz[j]*r6[1] + dpy_dz[j]*r6[3]

    end
end


function SC_kick!(SC::SC_lens, bbSC::Beam, rin)
    factorSC = 2*bbSC.classrad0/bbSC.beta^2/bbSC.gamma^3*SC.ds*bbSC.np #bbSC.np = 5e8, to have same Ji's I
    
    get_emittance!(bbSC)
    
    σz = bbSC.beamsize[5]
    σx = bbSC.beamsize[1] #sqrt(bbSC.emittance[1] * betax)  # beamsize x at SC_point
    σy = bbSC.beamsize[3]

    track_SC!(rin, σx, σy, σz, bbSC.temp1, bbSC, factorSC)
    
end


function pass!(SC::SC_lens, r_in::Array{Float64,1}, num_particles::Int64, bbSC::Beam)

    #phi_advx = LinRange(0, 29.228, SC.nSC+1)
    #phi_advy = LinRange(0, 30.210, SC.nSC+1)
    """
    d = Dirichlet(SC.nSC+1, 1.0)
    phi_advx = sort(rand(d)*29.228)
    phi_advy = sort(rand(d)*30.210)


    TM1 = Array{TransferMap4D, 1}(undef, SC.nSC)
    TM1_inv = Array{Inverse_TransferMap4D, 1}(undef, SC.nSC)

    opIPp = optics4DUC(0.8, 0.0, 0.072, 0.0)

    for j in 1:(SC.nSC)
        
        TM1[j] = TransferMap4D(opIPp, SC.optics, phi_advx[j], phi_advy[j])
        pass!(TM1[j], r_in, num_particles, bbSC) #line 148 Transfermaps
        
        SC_kick!(SC, bbSC, r_in)
        
        TM1_inv[j] = Inverse_TransferMap4D(opIPp, SC.optics, phi_advx[j], phi_advy[j])
        pass!(TM1_inv[j], r_in, num_particles, bbSC)
    end
    """
    SC_kick!(SC, bbSC, r_in)

end


# function for smoothing, replacing beamsize
#=
function SC_kick!(SC::SC_lens, bbSC::BunchedBeam, sx, sy)
    factorSC = 2*bbSC.particle.classrad0/bbSC.beta^2/bbSC.gamma^3*bbSC.num_particle*SC.ds
    
    get_emittance!(bbSC)
    σz = bbSC.beamsize[5]

    σx = sx *weight + bbSC.beamsize[1]*(1.0-weight)
    σy = sy *weight + bbSC.beamsize[3]*(1.0-weight)

    track!(σx, σy, σz, bbSC.temp1, bbSC, factorSC)
    
end=#