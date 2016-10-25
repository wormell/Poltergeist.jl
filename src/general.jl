rangedomain(op::Operator) = domain(rangespace(op))

containstransfer(M::Operator) = ApproxFun.iswrapper(M::Operator) && containstransfer(M.op)
for optype in (:(ApproxFun.InterlaceOperator),:(ApproxFun.PlusOperator),:(ApproxFun.MultiplicationOperator))
  eval(:(containstransfer(M::$optype) = any(containstransfer,M.ops)))
end



#Base.isless(a::ApproxFun.Infinity{Bool},b::ApproxFun.Infinity{Bool}) = a.angle && ~b.angle

#function ApproxFun.colstop{T,RM<:ApproxFun.RaggedMatrix,M<:ApproxFun.Operator,DS<:Any,RS<:Any,BI<:Any}(op::CachedOperator{T,RM,M,DS,RS,BI},j::Integer)

#function ApproxFun.resizedata!(ApproxFun.c)
