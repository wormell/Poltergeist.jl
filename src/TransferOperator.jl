abstract AbstractTransferOperator{T} <: Operator{T}

immutable TransferOperator{D<:Space,R<:Space,T,M<:AbstractMarkovMap} <: AbstractTransferOperator{T}
    m::M
    domainspace::D #Domain space
    rangespace::R
#    f::Array{Vector{T},1}
end

TransferOperator{T,MM<:AbstractMarkovMap}(::Type{T},m::MM)=
TransferOperator{typeof(domainspace(m)),typeof(rangespace(m)),T,MM}(m,domainspace(m),rangespace(m))


for OP in (:domainspace,:rangespace)
    @eval ApproxFun.$OP{Tr<:AbstractTransferOperator}(L::Tr) = L.$OP
end

Base.show(io::IO,L::TransferOperator) = print(io,typeof(L)) #temporary

