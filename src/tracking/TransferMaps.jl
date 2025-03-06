#TransferMapd4D
using StaticArrays
struct TransferMap4D <: AbstractTransverseMap
    dim::Int64
    tune::SVector{2,Float64}
    linearmap::SMatrix{4,4,Float64}
    TransferMap4D(tune::AbstractVector, mat::AbstractMatrix)=new(4, tune, mat)
end

function TransferMap4D(o1::AbstractOptics4D, o2::AbstractOptics4D, tunex::Float64, tuney::Float64)
    s1,c1=sincos(2π*tunex)
    r1=@SMatrix [c1 s1; -s1 c1]
    s2,c2=sincos(2π*tuney)
    r2=@SMatrix [c2 s2; -s2 c2]
    zm=@SMatrix [0.0 0.0; 0.0 0.0]
    rotation=SMatrix{4,4}([r1 zm; zm r2])
    tune=SVector{2,Float64}(tunex, tuney)
    return TransferMap4D(tune, invnormal_mat(o2)*rotation*normal_mat(o1))
end
function TransferMap4D(o1::AbstractOptics4D, tunex::Float64, tuney::Float64)
    return TransferMap4D(o1, o1, tunex, tuney)
end


struct Inverse_TransferMap4D <: AbstractTransverseMap
    dim::Int64
    tune::SVector{2,Float64}
    linearmap::SMatrix{4,4,Float64}
    Inverse_TransferMap4D(tune::AbstractVector, mat::AbstractMatrix)=new(4, tune, mat)
end
function Inverse_TransferMap4D(o1::AbstractOptics4D, o2::AbstractOptics4D, tunex::Float64, tuney::Float64)
    s1,c1=sincos(2π*(tunex))
    r1=@SMatrix [c1 s1; -s1 c1]
    s2,c2=sincos(2π*(tuney))
    r2=@SMatrix [c2 s2; -s2 c2]
    zm=@SMatrix [0.0 0.0; 0.0 0.0]
    rotation=(SMatrix{4,4}([r1 zm; zm r2]))
    tune=SVector{2,Float64}(tunex, tuney)
    S_inv = inv((invnormal_mat(o2)*rotation*normal_mat(o1)))
    return Inverse_TransferMap4D(tune,S_inv)
end

function Inverse_TransferMap4D(o1::AbstractOptics4D, tunex::Float64, tuney::Float64)
    return Inverse_TransferMap4D(o1, o1, tunex, tuney)
end


#TransferMap4DChrom
struct TransferMap4DChrom <: AbstractTransverseMap
    dim::Int64
    tune::SVector{2,Float64}
    chrom::SVector{2,Float64}
    umat::SMatrix{4,4,Float64}
    invumat::SMatrix{4,4,Float64}
    TransferMap4DChrom(tune::AbstractVector, chrom::AbstractVector, umat::AbstractMatrix, invumat::AbstractMatrix)=new(4, tune, chrom, umat, invumat)
end
function TransferMap4DChrom(o1::AbstractOptics4D, o2::AbstractOptics4D, tunex::Float64, tuney::Float64, chromx::Float64, chromy::Float64)
    tune=SVector{2,Float64}(tunex, tuney)
    chrom=SVector{2,Float64}(chromx, chromy)
    return TransferMap4DChrom(tune, chrom, normal_mat(o1), invnormal_mat(o2))
end

function TransferMap4DChrom(o1::AbstractOptics4D, tunex::Float64, tuney::Float64, chromx::Float64, chromy::Float64)
    return TransferMap4DChrom(o1, o1, tunex, tuney, chromx, chromy)
end



#start track! function to apply TransferMap matrix
function track!(rin, b::Beam, tm::TransferMap4D, tempx, temppx, tempy, num_macro)

    for j in 1:num_macro
        r6 = @view rin[(j-1)*6+1:j*6]
        tempx[j]  = tm.linearmap[1,1] * r6[1] + tm.linearmap[1,2] * r6[2]  + tm.linearmap[1,3] * r6[3] + tm.linearmap[1,4] * r6[4] 
        temppx[j] = tm.linearmap[2,1] * r6[1] + tm.linearmap[2,2] * r6[2] + tm.linearmap[2,3] * r6[3] + tm.linearmap[2,4] * r6[4]
        tempy[j]  = tm.linearmap[3,1] * r6[1] + tm.linearmap[3,2] * r6[2] + tm.linearmap[3,3] * r6[3] + tm.linearmap[3,4] * r6[4]
        r6[4] = tm.linearmap[4,1] * r6[1] + tm.linearmap[4,2] * r6[2] + tm.linearmap[4,3] * r6[3] + tm.linearmap[4,4] * r6[4]
        r6[1] = tempx[j]
        r6[2] = temppx[j]
        r6[3] = tempy[j]
    end
    return nothing
end


function pass!(tm::TransferMap4D, r_in::Array{Float64,1}, num_macro::Int, beam::Beam)
#function track!(beam::Beam, tm::TransferMap4D)
    track!(r_in, beam, tm, beam.temp1, beam.temp2, beam.temp3, num_macro)
end


#se commento tempy e r6[3] è tutto okay
function track!(rin, b::Beam, tm::Inverse_TransferMap4D, tempx, temppx, tempy, num_macro)
    
    for j in num_macro
        r6 = @view rin[(j-1)*6+1:j*6]
        tempx[j]  = tm.linearmap[1,1] * r6[1] + tm.linearmap[1,2] * r6[2] + tm.linearmap[1,3] * r6[3] + tm.linearmap[1,4] * r6[4] 
        temppx[j] = tm.linearmap[2,1] * r6[1] + tm.linearmap[2,2] * r6[2] + tm.linearmap[2,3] * r6[3] + tm.linearmap[2,4] * r6[4]
        tempy[j] = tm.linearmap[3,1] * r6[1] + tm.linearmap[3,2] * r6[2] + tm.linearmap[3,3] * r6[3] + tm.linearmap[3,4] * r6[4]
        r6[4] = tm.linearmap[4,1] * r6[1] + tm.linearmap[4,2] * r6[2] + tm.linearmap[4,3] * r6[3] + tm.linearmap[4,4] * r6[4]
        r6[1] = tempx[j]
        r6[2] = temppx[j]
        r6[3] = tempy[j]
    end
    return nothing
end


function pass!(tm::Inverse_TransferMap4D, r_in::Array{Float64,1}, num_macro::Int, beam::Beam)
#function track!(beam::Beam, tm::Inverse_TransferMap4D)    
    track!(r_in, beam, tm, beam.temp1, beam.temp2, beam.temp3, num_macro)
end


function track!(rin, b::Beam, tm::TransferMap4DChrom, temp1, temp2, temp3, sinphi, cosphi, num_macro) #where T
    
    for j in num_macro
        r6 = @view rin[(j-1)*6+1:j*6]
        temp1[j]  = tm.umat[1,1] * r6[1] + tm.umat[1,2] * r6[2] + tm.umat[1,3] * r6[3] + tm.umat[1,4] * r6[4] 
        temp2[j] = tm.umat[2,1] * r6[1] + tm.umat[2,2] * r6[2] + tm.umat[2,3] * r6[3] + tm.umat[2,4] * r6[4]
        temp3[j]  = tm.umat[3,1] * r6[1] + tm.umat[3,2] * r6[2] + tm.umat[3,3] * r6[3] + tm.umat[3,4] * r6[4]
        r6[4] = tm.umat[4,1] * r6[1] + tm.umat[4,2] * r6[2] + tm.umat[4,3] * r6[3] + tm.umat[4,4] * r6[4]
        r6[1] = temp1[j]
        r6[2] = temp2[j]
        r6[3] = temp3[j]

        temp1[j] = 2π*tm.tune[1] + (2π*tm.chrom[1]) * r6[6]
        sinphi = sin.(temp1[j])
        cosphi = cos.(temp1[j])
        temp3[j]  = cosphi * r6[1] + sinphi * r6[2]
        r6[2] =  cosphi * r6[2] - sinphi * r6[1]
        r6[1] = temp3[j]

        temp2[j] = 2π*tm.tune[2] + (2π*tm.chrom[2]) * r6[5]
        sinphi = sin.(temp2[j])
        cosphi = cos.(temp2[j])
        temp3[j]  = cosphi * r6[3] + sinphi * r6[4]
        r6[4] = cosphi * r6[4] - sinphi * r6[3]
        r6[3] = temp3[j]


        temp1[j]  = tm.invumat[1,1] * r6[1] + tm.invumat[1,2] * r6[2] + tm.invumat[1,3] * r6[3] + tm.invumat[1,4] * r6[4]
        temp2[j] = tm.invumat[2,1] * r6[1] + tm.invumat[2,2] * r6[2] + tm.invumat[2,3] * r6[3] + tm.invumat[2,4] * r6[4]
        temp3[j]  = tm.invumat[3,1] * r6[1] + tm.invumat[3,2] * r6[2] + tm.invumat[3,3] * r6[3] + tm.invumat[3,4] * r6[4]
        r6[4] = tm.invumat[4,1] * r6[1] + tm.invumat[4,2] * r6[2] + tm.invumat[4,3] * r6[3] + tm.invumat[4,4] * r6[4]
        r6[1] = temp1[j]
        r6[2] = temp2[j]
        r6[3] = temp3[j]
    end
    return nothing

end

    
function pass!(tm::TransferMap4DChrom, r_in::Array{Float64,1}, num_macro::Int, beam::Beam)
#function track!(beam::Beam, tm::TransferMap4DChrom)
    track!(r_in, beam, tm,  beam.temp1, beam.temp2, beam.temp3, beam.temp4, beam.temp5, num_macro)
end

