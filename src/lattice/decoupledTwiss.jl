# under construction
# struct DecoupledTwiss
# 	betx::RealType
# 	alfx::RealType
# 	gamx::RealType
# 	dx::RealType
# 	dpx::RealType
# 	smux::RealType
# 	cmux::RealType
# 	bety::RealType
# 	alfy::RealType
# 	gamy::RealType
# 	dy::RealType
# 	dpy::RealType
# 	smuy::RealType
# 	cmuy::RealType
# end

# DecoupledTwiss(;betx::RealType,bety::RealType,
# 			   alfx::RealType=RealType(0),alfy::RealType=RealType(0),
# 			   dx::RealType=RealType(0),dy::RealType=RealType(0),
# 			   dpx::RealType=RealType(0),dpy::RealType=RealType(0),
# 			   mux::RealType=RealType(0),muy::RealType=RealType(0))=DecoupledTwiss(betx,alfx,(RealType(1)+alfx*alfx)/betx,dx,dpx,sincos(mux)...,
# 																				   bety,alfy,(RealType(1)+alfy*alfy)/bety,dy,dpy,sincos(muy)...)


# twissTransform(M)=begin
#     m11,m21,m12,m22=M
#     [m11*m11 -2m11*m12 m12*m12
#     -m11*m21 1.0+2m12*m21 -m12*m22
#     m21*m21 -2m21*m22 m22*m22]
# end

# function decoupledTwissPropagate(tin::DecoupledTwiss,M::Matrix{RealType})
# 	Mx=@view M[1:2,1:2]	
# 	Nx=twissTransform(Mx)
# 	My=@view M[3:4,3:4]	
# 	Ny=twissTransform(My)
# 	betx,alfx,gamx=Nx*[tin.betx;tin.alfx;tin.gamx]
# 	bety,alfy,gamy=Ny*[tin.bety;tin.alfy;tin.gamy]
# 	dx,dpx=Mx*[tin.dx;tin.dpx]+(@view M[1:2,6])
# 	dy,dpy=My*[tin.dy;tin.dpy]+(@view M[3:4,6])
# 	sin_dmux=M[1,2]/sqrt(betx*tin.betx)
# 	cos_dmux=M[1,1]*sqrt(tin.betx/betx)-tin.alfx*sin_dmux
# 	sin_dmuy=M[3,4]/sqrt(bety*tin.bety)
# 	cos_dmuy=M[3,3]*sqrt(tin.bety/bety)-tin.alfy*sin_dmuy
# 	smux=sin_dmux*tin.cmux+cos_dmux*tin.smux
# 	cmux=cos_dmux*tin.cmux-sin_dmux*tin.smux
# 	smuy=sin_dmuy*tin.cmuy+cos_dmuy*tin.smuy
# 	cmuy=cos_dmuy*tin.cmuy-sin_dmuy*tin.smuy
# 	#mux=tin.mux+atan(tin.gamx*M[1,2]/M[1,1]-tin.alfx)+atan(tin.alfx)
# 	#muy=tin.muy+atan(tin.gamy*M[3,4]/M[3,3]-tin.alfy)+atan(tin.alfy)
# 	return DecoupledTwiss(betx,alfx,gamx,dx,dpx,smux,cmux,bety,alfy,gamy,dy,dpy,smuy,cmuy)
# end

# function decoupledTwissPropagate(tin::DecoupledTwiss,seq::Sequence,beam::AbstractBeam=_beam[])
# 	ret=Vector{DecoupledTwiss}(undef,1+length(seq))
# 	ss=zeros(RealType,length(ret))
# 	names=Vector{String}(undef,length(ret))
# 	ret[1]=tin
# 	names[1]="Start"
# 	for (index,mag) in enumerate(seq.Line)
# 		M=transferMatrix(mag,beam)
# 		ret[index+1]=decoupledTwissPropagate(ret[index],M)
# 		ss[index+1]=refS(mag)+ss[index]
# 		names[index+1]=mag.Name
# 	end
# 	return ss,names,StructArray(ret)
# end