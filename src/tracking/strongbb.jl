# strong beam-beam is still under development
# don't use this file for now
using SpecialFunctions
using Statistics
# using StaticArrays

function normal_mat(op2d::optics2D)
    return [1.0/sqrt(op2d.beta) 0.0; op2d.alpha/sqrt(op2d.beta)  sqrt(op2d.beta)]
end

function invnormal_mat(op2d::optics2D)
    return [sqrt(op2d.beta) 0.0; -op2d.alpha*sqrt(op2d.beta)  1.0/sqrt(op2d.beta)]
end


function normal_mat(op4d::optics4DUC)
    nfmx=[1.0/sqrt(op4d.optics_x.beta) 0.0; op4d.optics_x.alpha/sqrt(op4d.optics_x.beta)  sqrt(op4d.optics_x.beta)]
    nfmy=[1.0/sqrt(op4d.optics_y.beta) 0.0; op4d.optics_y.alpha/sqrt(op4d.optics_y.beta)  sqrt(op4d.optics_y.beta)] 
    return [nfmx zeros(2,2); zeros(2,2) nfmy]
end

function invnormal_mat(op4d::optics4DUC)
    invnfmx=[sqrt(op4d.optics_x.beta) 0.0; -op4d.optics_x.alpha*sqrt(op4d.optics_x.beta)  1.0/sqrt(op4d.optics_x.beta)]
    invnfmy=[sqrt(op4d.optics_y.beta) 0.0; -op4d.optics_y.alpha*sqrt(op4d.optics_y.beta)  1.0/sqrt(op4d.optics_y.beta)]
    return[invnfmx zeros(2,2); zeros(2,2) invnfmy]
end

function new_collect(zmin, zmax, n)
    v = zeros(n)
    for i in 1:n
        v[i] = zmin + (i-1)*(zmax-zmin)/(n-1)
    end
    return v
end
function initilize_zslice!(beam::StrongGaussianBeam, profile::Symbol, slice_type::Symbol, zrange::Float64=5.0)
    zmin=-zrange*beam.beamsize[3]
    zmax=zrange*beam.beamsize[3]
    if profile == :gaussian
        if slice_type == :evenzsep
            zedge = new_collect(zmin, zmax, beam.nzslice+1)
            beam.zslice_center .= 0.5.*(zedge[1:end-1]+zedge[2:end])
            beam.zslice_npar .= exp.(-0.5.*(beam.zslice_center.^2)./beam.beamsize[3]^2)
            beam.zslice_npar .= beam.zslice_npar./sum(beam.zslice_npar).*beam.num_particle
        end
        if slice_type == :evennpar   # Here the zslicecenter is where split the npar in the slice, not the center of zposition
            npartedge = new_collect(0.0, 1.0, beam.nzslice+1)
            npartcenter=0.5.*(npartedge[1:end-1]+npartedge[2:end])  
            beam.zslice_center .= (sqrt(2.0)*beam.beamsize[3]).*erfinv.(2.0.*npartcenter.-1.0)
            beam.zslice_npar .= zeros(beam.nzslice).+1.0/beam.nzslice*beam.num_particle
        end
    end

    if profile == :uniform
        zedge = new_collect(zmin, zmax, beam.nzslice+1)
        beam.zslice_center .= 0.5.*(zedge[1:end-1]+zedge[2:end])
        beam.zslice_npar .= zeros(beam.nzslice).+1.0/beam.nzslice*beam.num_particle
    end
    return nothing
end

function initilize_zslice!(beam::StrongGaussianBeam, zlist::Vector{Float64}, slice_type::Symbol) # Convert from distribution
    zmin=minimum(zlist)
    zmax=maximum(zlist)
    if slice_type == :evenzsep
        zedge=collect(range(zmin, stop=zmax, length=beam.nzslice+1))
        beam.zslice_center=0.5.*(zedge[1:end-1]+zedge[2:end])
        beam.zslice_npar=zeros(beam.nzslice)
        for i in 1:beam.nzslice
            beam.zslice_npar[i]=sum((zlist.>=zedge[i]).&(zlist.<zedge[i+1]))
        end
        beam.zslice_npar=beam.zslice_npar./sum(beam.zslice_npar).*beam.num_particle
    end
    if slice_type == :evennpar
        sort_zlist=sort(zlist)
        nzlist=length(sort_zlist)
        npartedge=Int.(floor.(collect(range(0, nzlist, length=beam.nzslice+1))))
        beam.zslice_npar .= (npartedge[2:end]-npartedge[1:end-1])/nzlist*beam.num_particle
        beam.zslice_center .= zeros(beam.nzslice)
        for i in 1:beam.nzslice
            beam.zslice_center[i]=mean(sort_zlist[npartedge[i]+1:npartedge[i+1]])
        end
    end
    return nothing
end


function Bassetti_Erskine_xgty!(res::AbstractVector, x::Float64, y::Float64, sigmax::Float64, sigmay::Float64) # x size greater than y
    # Only positive y is valid for this function
    # for y<0, Ex = Ex, Ey = -Ey
    if y < 0.0
        Bassetti_Erskine_xgty!(res, x, -y, sigmax, sigmay)
        res[2] = -res[2]
        return nothing
    end
    termexp=exp(-x*x/2/sigmax/sigmax-y*y/2/sigmay/sigmay)
	sqrtδsigma2=sqrt(Complex(2*(sigmax*sigmax-sigmay*sigmay)))
	term1=erfcx(-1.0im*(x+1.0im*y)/sqrtδsigma2)
	term2=erfcx(-1im*(x*sigmay/sigmax+1im*y*sigmax/sigmay)/sqrtδsigma2)
	
	complex_e=-1im*2*sqrt(pi)/sqrtδsigma2*(term1-termexp*term2)
	res[1]=real(complex_e)
    res[2]=-imag(complex_e)
    res[3]=termexp/2.0/π/sigmax/sigmay
    return nothing
end

function Bassetti_Erskine_ygtx!(res::AbstractVector, x::Float64, y::Float64, sigmax::Float64, sigmay::Float64) # x size greater than y
    # Only negative x is valid for this function
    # for x>0, Ex = -Ex, Ey = Ey
    if x > 0.0
        Bassetti_Erskine_ygtx!(res, -x, y, sigmax, sigmay)
        res[1] = -res[1]
        return nothing
    end
    termexp=exp(-x*x/2/sigmax/sigmax-y*y/2/sigmay/sigmay)
	sqrtδsigma2=sqrt(Complex(2*(sigmax*sigmax-sigmay*sigmay)))
	term1=erfcx(-1.0im*(x+1.0im*y)/sqrtδsigma2)
	term2=erfcx(-1im*(x*sigmay/sigmax+1im*y*sigmax/sigmay)/sqrtδsigma2)
	
	complex_e=-1im*2*sqrt(pi)/sqrtδsigma2*(term1-termexp*term2)
    res[1]=real(complex_e)
    res[2]=-imag(complex_e)
    res[3]=termexp/2.0/π/sigmax/sigmay
	return nothing
end

function Bassetti_Erskine!(res::AbstractVector, x::Float64, y::Float64, sigmax::Float64, sigmay::Float64)
    if sigmax > sigmay
        Bassetti_Erskine_xgty!(res, x, y, sigmax, sigmay)
        return nothing
    else
        Bassetti_Erskine_ygtx!(res, x, y, sigmax, sigmay) 
        return nothing
    end
end


function track_sbb!(rin, num_macro, temp1, temp2, temp3, temp4, temp5, sgb::StrongGaussianBeam, factor::Float64) 
    #factor=wb.particle.classrad0/wb.gamma*wb.particle.charge*sgb.particle.charge
    
    lumi=0.0
    fieldvec = zeros(3)

    @inbounds for i in 1:sgb.nzslice
        slicelumi=0.0
        @inbounds for j in 1:num_macro
            r6 = @view rin[(j-1)*6+1:j*6]
            # temp1: collision zlocation, temp2: beamsize x, temp3: beamsize y, temp4: beta x, temp5: beta y
            temp1[j] = (r6[5] .+ sgb.zslice_center[i])./2.0
            temp4[j] = sgb.optics.optics_x.beta .+ sgb.optics.optics_x.gamma .* temp1[j] .* temp1[j] .- 2.0 .* sgb.optics.optics_x.alpha .* temp1[j]
            temp2[j] = sgb.beamsize[1] .* sqrt.(temp4[j] ./ sgb.optics.optics_x.beta)
            temp5[j] = sgb.optics.optics_y.beta .+ sgb.optics.optics_y.gamma .* temp1[j] .* temp1[j] .- 2.0 .* sgb.optics.optics_y.alpha .* temp1[j]
            temp3[j] = sgb.beamsize[2] .* sqrt.(temp5[j] ./ sgb.optics.optics_y.beta)
        
            # temp4 and temp5 are free to change now.
            r6[1] += (r6[2] .* temp1[j])
            r6[3] += (r6[4] .* temp1[j])
            Bassetti_Erskine!(fieldvec, r6[1], r6[3], temp2[j], temp3[j])
            r6[2] += (sgb.zslice_npar[i]*factor) * fieldvec[1]
            r6[4] += (sgb.zslice_npar[i]*factor) * fieldvec[2]
            slicelumi += fieldvec[3]

            r6[1] -= (r6[2] .* temp1[j])
            r6[3] -= (r6[4] .* temp1[j])
        end
       
        lumi += slicelumi * sgb.zslice_npar[i] #  Will do it outside* wb.num_particle / wb.num_macro
    end

    return lumi

end


function pass!(sgb::StrongGaussianBeam, r_in::Array{Float64,1}, num_macro::Int, wb::Beam)
    factor=wb.classrad0/wb.gamma*wb.charge*sgb.charge
    lumi=track_sbb!(r_in, num_macro, wb.temp1, wb.temp2, wb.temp3, wb.temp4, wb.temp5, sgb, factor)
    lumi *= wb.np / wb.nmacro
    return nothing
end

