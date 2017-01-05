export SchurInv

# Schur wrappers
type SchurInvWrapper{QR<:ApproxFun.QROperator,F<:Fun,T} <: Operator{T}
  op::QR
  u::F
end
SchurInvWrapper(op::ApproxFun.QROperator,f::Fun) = SchurInvWrapper{typeof(op),typeof(f),eltype(op)}(op,f)
ApproxFun.qrfact!(op::SchurInvWrapper) = op
ApproxFun.qrfact(op::SchurInvWrapper) = op
#SchurInvWrapper(op::Operator) = SchurInvWrapper(op,uniform(domainspace(op)))
ApproxFun.@wrapper SchurInvWrapper

#ApproxFun.linsolve{S,T,DD,dim}(A::SchurInvWrapper,b::Fun{MatrixSpace{S,T,DD,dim}};kwds...) = linsolve(L.op,b;kwds...) # avoid method ambiguity
ApproxFun.linsolve(L::SchurInvWrapper,b;kwds...) = linsolve(L.op,b;kwds...)

uniform(S::Space) = Fun(1.,S)/sum(Fun(1.,S))
uniform(D::Domain) = Fun(1.,D)/sum(Fun(1.,D))
function SchurInv(L::Operator,u::Fun=uniform(domainspace(L)))
  @assert domain(L) == rangedomain(L)
  if isa(domainspace(L),ApproxFun.ProductDomain) #TODO: put this into ApproxFun
    di = DefiniteIntegral(domainspace(L).spaces[1])⊗DefiniteIntegral(domainspace(L).spaces[2])
    for i = 3:length(domainspace(L).spaces)
      di = di⊗DefiniteIntegral(domainspace(L).spaces[2])
    end
  else
    di = DefiniteIntegral(domainspace(L))
  end


  SchurInvWrapper(qrfact(I-L + cache(u/sum(u) * di)),u)
end
SchurInv(M::MarkovMap,u::Fun=uniform(Space(domain(M)))) = SchurInv(Transfer(M,space(u)),u)


# function DefiniteIntegral(sp::ProductSpace) # likely a hack
#   di = DefiniteIntegral(sp.spaces[1])
#   for i = 2:length(sp.spaces)
#     di = di⊗DefiniteIntegral(sp.spaces[i])
#   end
#   di
# end
