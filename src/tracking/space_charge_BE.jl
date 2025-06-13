function track_SC!(rin, σx, σy, σz, temp1, bbSC::Beam, factorSC::Float64, num_macro)

    #dpx_dz = zeros(length(bbSC.nmacro))
    #dpy_dz = zeros(length(bbSC.nmacro))
   
    fieldvec = zeros(3)

    @inbounds for j in 1:num_macro

        r6 = @view rin[(j-1)*6+1:j*6]

        #println("Before: px = ", r6[2], ", py = ", r6[4])

        Bassetti_Erskine!(fieldvec, r6[1], r6[3], σx, σy)
        #println("Field: ", fieldvec[1], " ;", fieldvec[2])

        # gaussian function for particle distribution, lambda_z
        temp1[j] = 1.0/sqrt(2*π)/σz*exp((-0.5)*r6[5]^2/σz^2)

        r6[2] += factorSC* temp1[j]* fieldvec[1]
        r6[4] += factorSC* temp1[j]* fieldvec[2]
        
        #dpx_dz[j] = factorSC*temp1[j] *fieldvec[1]* r6[5]/σz^2
        #dpy_dz[j] = factorSC*temp1[j] *fieldvec][2]* r6[5]/σz^2
        #r6[6] += dpx_dz[j]*r6[1] + dpy_dz[j]*r6[3]

        #println("After: px = ", r6[2], ", py = ", r6[4]) 
        #println("total is:", factorSC* temp1[j]* fieldvec[1])

    end
end


function SC_kick!(SC::SC_lens, bbSC::Beam, r_in, num_macro)
    
    factorSC = 2*bbSC.classrad0/bbSC.beta^2/bbSC.gamma^3*bbSC.np*SC.ds
    
    get_emittance!(bbSC)
    
    σx = bbSC.beamsize[1] #sqrt(bbSC.emittance[1] * betax)  # beamsize x at SC_point
    σy = bbSC.beamsize[3]
    σz = bbSC.beamsize[5]

    track_SC!(r_in, σx, σy, σz, bbSC.temp1, bbSC, factorSC, num_macro)
    
end

function pass!(SC::SC_lens, r_in::Array{Float64,1}, num_macro::Int64, bbSC::Beam)
    SC_kick!(SC, bbSC, r_in, num_macro)
end


"""
# imperfect - add phi advx, TM as struct element in SC lens??
############### pass with TM+inv_TM #############
function SC_kick_ring!(SC::SC_lens, bbSC::Beam, r_in, num_macro)
    
    ### maps ###
    #phi_advx = LinRange(0.0, 29.228, SC.nSC+1)
    #phi_advy = LinRange(0.0, 30.210, SC.nSC+1)

    phi_advx = sort(29.228 .* rand(SC.nSC+1))
    phi_advy = sort(30.210 .* rand(SC.nSC+1))

    opIPp = optics4DUC(0.8, 0.0, 0.072, 0.0) #beam.optics
    
    TM = Array{TransferMap4D, 1}(undef, SC.nSC)
    TM_inv = Array{Inverse_TransferMap4D, 1}(undef, SC.nSC)

    ### sc kick ###
    for j in 1:SC.nSC
        
        TM[j] = TransferMap4D(opIPp, SC.optics, phi_advx[j], phi_advy[j])
        pass!(TM[j], r_in, num_macro, bbSC) #line 148 Transfermaps
        
        SC_kick!(SC, bbSC, r_in, num_macro)
        
        TM_inv[j] = Inverse_TransferMap4D(opIPp, SC.optics, phi_advx[j], phi_advy[j])
        pass!(TM_inv[j], r_in, num_macro, bbSC)
    end
    
end

function pass!(SC::SC_lens, r_in::Array{Float64,1}, num_macro::Int64, bbSC::Beam)
    SC_kick_ring!(SC, bbSC, r_in, num_macro)
end
"""


"""
###### parallelizing ########
function track_SC_P!(rin, σx, σy, σz, temp1, bbSC::Beam, factorSC::Float64, num_macro)

    #dpx_dz = zeros(length(bbSC.nmacro))
    #dpy_dz = zeros(length(bbSC.nmacro))
   
    fieldvec_thread=[MVector{3}(0.0, 0.0, 0.0)  for j = 1:Threads.nthreads()]
    @inbounds Threads.@threads :static for j in eachindex(num_macro)
    
        r6 = @view rin[(j-1)*6+1:j*6]

        Bassetti_Erskine!(fieldvec_thread[Threads.threadid()], r6[1], r6[3], σx, σy)

        # gaussian function for particle distribution, lambda_z
        temp1[j] = 1.0/sqrt(2*π)/σz*exp((-0.5)*r6[5]^2/σz^2)

        r6[2] += factorSC* temp1[j]* fieldvec_thread[Threads.threadid()][1]
        r6[4] += factorSC* temp1[j]* fieldvec_thread[Threads.threadid()][2]
        
        #dpx_dz[j] = factorSC*temp1[j] *fieldvec_thread[Threads.threadid()][1]* r6[5]/σz^2
        #dpy_dz[j] = factorSC*temp1[j] *fieldvec_thread[Threads.threadid()][2]* r6[5]/σz^2
        #r6[6] += dpx_dz[j]*r6[1] + dpy_dz[j]*r6[3]

        #println("total is:", factorSC* temp1[j]* fieldvec_thread[Threads.threadid()][1])


    end
end

function SC_kick_P!(SC::SC_lens, bbSC::Beam, r_in, num_macro)
    
    factorSC = 2*bbSC.classrad0/bbSC.beta^2/bbSC.gamma^3*bbSC.np*SC.ds*
    
    get_emittance!(bbSC)
    
    σx = bbSC.beamsize[1] #sqrt(bbSC.emittance[1] * betax)  # beamsize x at SC_point
    σy = bbSC.beamsize[3]
    σz = bbSC.beamsize[5]

    track_SC_P!(r_in, σx, σy, σz, bbSC.temp1, bbSC, factorSC, num_macro)
    
end

function pass_P!(SC::SC_lens, r_in::Array{Float64,1}, num_macro::Int64, bbSC::Beam)
    SC_kick_P!(SC, bbSC, r_in, num_macro)
end

"""



"""
# function for smoothing, replacing beamsize

function SC_kick!(SC::SC_lens, bbSC::BunchedBeam, sx, sy)
    factorSC = 2*bbSC.particle.classrad0/bbSC.beta^2/bbSC.gamma^3*bbSC.num_particle*SC.ds
    
    get_emittance!(bbSC)
    σz = bbSC.beamsize[5]

    σx = sx *weight + bbSC.beamsize[1]*(1.0-weight)
    σy = sy *weight + bbSC.beamsize[3]*(1.0-weight)

    track!(σx, σy, σz, bbSC.temp1, bbSC, factorSC)
    
end

"""


