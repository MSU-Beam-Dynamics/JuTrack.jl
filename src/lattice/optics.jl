abstract type AbstractOptics end
abstract type AbstractOptics2D <:AbstractOptics end
abstract type AbstractOptics4D <:AbstractOptics end
mutable struct optics2D <: AbstractOptics2D
    beta::Float64
    alpha::Float64
    gamma::Float64
    phase::Float64
    eta::Float64
    etap::Float64
end
optics2D(b::Float64, a::Float64, ph::Float64, eta::Float64, etap::Float64)=optics2D(b,a,(1.0+a^2)/b,ph,eta,etap)
optics2D(b::Float64, a::Float64)=optics2D(b,a,(1.0+a^2)/b,0.0,0.0,0.0)

mutable struct optics4DUC <: AbstractOptics4D # 4D linear transformation uncoupled
    optics_x::optics2D
    optics_y::optics2D
end
optics4DUC(bx::Float64, ax::Float64, by::Float64, ay::Float64)=optics4DUC(optics2D(bx,ax),optics2D(by,ay))