export domain, rangedomain

abstract AbstractMarkovMap{D<:Domain,R<:Domain,T,FF} <: Function
abstract AbstractDerivativeMarkovMap{D<:Domain,R<:Domain,T,FF} <: AbstractMarkovMap{D,R,T,FF}

Base.summary(m::AbstractMarkovMap) =  string(typeof(m).name.name)*":"*string(domain(m))*"↦"*string(rangedomain(m)) #branches??
Base.eltype(m::AbstractMarkovMap) = eltype(rangedomain(m))
Base.show(io::IO,m::AbstractMarkovMap) = print(io,typeof(m)) #temporary


immutable MarkovMap{D<:Domain,R<:Domain,T,ff} <: AbstractMarkovMap{D,R,T,ff}
    domain::D
    range::R
    f::ff
    bl::Vector{T} #lower edges of branches
    bu::Vector{T} #upper edges of branches
    sgns::Vector{T} #sign of derivative of f on branches
#     λhat::T #maximum derivative
#     C1::T #first distortion bound
    # cache of inverse branches of f
    function FunctionMarkovMap{T}(
        dom::Domain,ran::Domain,f,bl::Vector{T},bu::Vector{T},sgns::Vector{T})
        @assert length(bl) == length(bu)
        @assert length(bl) == length(sgns)
        @assert all(sgns.^2 == 1)
        @assert all(bl.<bu)
        @assert sum(bu-bl) <= arclength(dom)
            # assert branches non-overlapping, in D,...
        new{typeof(dom),typeof(ran),T,typeof(f)}(dom,ran,f,bl,bu,sgns)
    end
end

length(m::MarkovMap) = length(m.bl)


immutable InverseDerivativeMarkovMap{D<:Domain,R<:Domain,T,FF} <: AbstractDerivativeMarkovMap{D,R,T,FF}
    domain::D
    rangedomain::R
    v::Vector{FF} #branches of inverse
    dvdx::Vector{FF} #derivatives of branches of inverse
    sgns::Vector{T} #sign of derivative of f on branches
#     λhat::T #maximum derivative
#     C1::T #first distortion bound
#    InverseDerivativeMarkovMap(v,dvdx,sgns,Ds,Rs,T,FF) = new(v,)

    # assert stuff??
 #   function InverseDerivativeMarkovMap{FF,T}(dom::Domain,ran::Domain,v::Vector{FF},dvdx::Vector{FF},sgns::Vector{T})
 #       @assert length(v) == length(dvdx)
 #       @assert length(v) == length(sgns)
 #       @assert all(sgns.^2 == 1)
 #       new(dom,ran,v,dvdx,sgns)
 #   end
end

function InverseDerivativeMarkovMap{FF}(dom::Domain,ran::Domain,v::Vector{FF},dvdx::Vector{FF})
    sgnpt = rand(ran); #hack
    return InverseDerivativeMarkovMap(dom,ran,v,dvdx,eltype(dom)[sign(dvdxi(sgnpt)) for dvdxi in dvdx])
end

length(m::InverseDerivativeMarkovMap) = length(m.v)
function getindex(m::InverseDerivativeMarkovMap,d::Symbol)
    d==:v && return m.v
    d==:dvdx && return m.dvdx

    error("Symbol not recognised")
end
#getindex(m::InverseDerivativeMarkovMap,:dvdx) = m.dvdx

for TYP in (:MarkovMap,:InverseDerivativeMarkovMap)
    @eval ApproxFun.domain(m::$TYP) = m.domain
    @eval rangedomain(m::$TYP) = m.rangedomain
end

abstract MarkovRegularityWrapper{D<:Domain,R<:Domain,T,FF} <: AbstractMarkovMap{D,R,T,FF}

type AnalyticMarkovRegularityWrapper{D<:Domain,R<:Domain,T,FF} <: MarkovRegularityWrapper{D,R,T,FF}
  mm::AbstractMarkovMap{D,R,T,FF}
  delta::T
  p::T
  K::T
end

function getanalyticp{D<:Domain,R<:PeriodicInterval,T,FF}(mm::AbstractDerivativeMarkovMap{D,R,T,FF},delta::T)
  try
    maxabs(Fun(x->imag(mm[:dvdx](x+im*delta)),space(R)))
  catch
    mmf = Fun(mm[:dvdx],space(R))
    maxabs(Fun(x->imag(mmf(x+im*delta)),space(R)))
  end
end


type DifferentiableMarkovRegularityWrapper{D<:Domain,R<:Domain,T,FF} <: MarkovRegularityWrapper{D,R,T,FF}
  mm::AbstractMarkovMap{D,R,T,FF}
  r::Int
  W::Array{T}
end

