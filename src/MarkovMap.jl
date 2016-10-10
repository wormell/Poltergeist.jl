abstract AbstractMarkovMap{D<:Space,R<:Space,ff,T} <: Function

immutable FunctionMarkovMap{D<:Space,R<:Space,ff<:Function,T} <: AbstractMarkovMap{D,R,ff,T}
    domainspace::D
    rangespace::R
    f::ff
    bp::Vector{T} #break points in f
    sgns::Vector{T} #sign of derivative of f on branches
#     λhat::T #maximum derivative
#     C1::T #first distortion bound
    # cache of inverse branches of f
end

immutable InverseDerivativeMarkovMap{D<:Space,R<:Space,T,FF} <: AbstractMarkovMap{D,R,FF,T}
    domainspace::D
    rangespace::R
    v::Vector{FF} #branches of inverse
    dvdx::Vector{FF} #derivatives of branches of inverse
    sgns::Vector{T} #sign of derivative of f on branches
#     λhat::T #maximum derivative
#     C1::T #first distortion bound
#    InverseDerivativeMarkovMap(v,dvdx,sgns,Ds,Rs,T,FF) = new(v,)
end

function InverseDerivativeMarkovMap{D<:Space,R<:Space,FF}(domainspace::D,rangespace::R,v::Vector{FF},dvdx::Vector{FF})
    sgnpt = rand(rangespace); #hack
    return InverseDerivativeMarkovMap(domainspace,rangespace,v,dvdx,eltype(D)[sign(dvdxi(sgnpt)) for dvdxi in dvdx])
end

for OP in (:domainspace,:rangespace), TYP in (:FunctionMarkovMap,:InverseDerivativeMarkovMap)
    @eval ApproxFun.$OP{MM<:$TYP}(M::MM) = M.$OP
end

Base.summary(M::AbstractMarkovMap) =  string(typeof(M).name.name)*":"*string(domainspace(M))*"↦"*string(rangespace(M)) #branches??
Base.show(io::IO,M::AbstractMarkovMap) = print(io,typeof(M)) #temporary

