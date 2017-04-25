export SolutionInv

# Solution wrappers
type SolutionInvWrapper{QR<:ApproxFun.QROperator,F<:Fun,T} <: Operator{T}
  op::QR
  u::F
end
SolutionInvWrapper(op::ApproxFun.QROperator,f::Fun) = SolutionInvWrapper{typeof(op),typeof(f),eltype(op)}(op,f)
ApproxFun.qrfact!(op::SolutionInvWrapper) = op
ApproxFun.qrfact(op::SolutionInvWrapper) = op
#SolutionInvWrapper(op::Operator) = SolutionInvWrapper(op,uniform(domainspace(op)))
ApproxFun.@wrapper SolutionInvWrapper

#ApproxFun.\{S,T,DD,dim}(A::SolutionInvWrapper,b::Fun{MatrixSpace{S,T,DD,dim}};kwds...) = \(L.op,b;kwds...) # avoid method ambiguity
(\)(L::SolutionInvWrapper,b::ApproxFun.Fun;kwds...) = \(L.op,b;kwds...)

function uniform(S::Space)
  u = Fun(one,S)
  scale!(u.coefficients,1/sum(u))
  u
end
uniform(D::Domain) = uniform(Space(D))

function SolutionInv(L::Operator,u::Fun=uniform(domainspace(L)))
  @assert domain(L) == rangedomain(L)
  if isa(domainspace(L),ApproxFun.TensorSpace) #TODO: put this into ApproxFun
    di = DefiniteIntegral(domainspace(L).spaces[1])⊗DefiniteIntegral(domainspace(L).spaces[2])
    for i = 3:length(domainspace(L).spaces)
      di = di⊗DefiniteIntegral(domainspace(L).spaces[2])
    end
  else
    di = DefiniteIntegral(domainspace(L))
  end


  SolutionInvWrapper(qrfact(I-L + cache(u/sum(u) * di)),u)
end
SolutionInv(M::MarkovMap,u::Fun=uniform(Space(domain(M)))) = SolutionInv(Transfer(M,space(u)),u)


# function DefiniteIntegral(sp::ProductSpace) # likely a hack
#   di = DefiniteIntegral(sp.spaces[1])
#   for i = 2:length(sp.spaces)
#     di = di⊗DefiniteIntegral(sp.spaces[i])
#   end
#   di
# end
