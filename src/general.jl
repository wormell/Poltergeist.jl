rangedomain(op::Operator) = domain(rangespace(op))

containstransfer(M::Operator) = ApproxFun.iswrapper(M::Operator) && containstransfer(M.op)
for optype in (:(ApproxFun.InterlaceOperator),:(ApproxFun.PlusOperator),:(ApproxFun.TimesOperator))
  eval(:(containstransfer(M::$optype) = any(containstransfer,M.ops)))
end

ApproxFun.colstop(L::ApproxFun.LowRankPertOperator,k::Integer) = max(colstop(L.op,k::Integer),colstop(L.pert,k::Integer))

#Base.isless(a::ApproxFun.Infinity{Bool},b::ApproxFun.Infinity{Bool}) = a.angle && ~b.angle

#function ApproxFun.colstop{T,RM<:ApproxFun.RaggedMatrix,M<:ApproxFun.Operator,DS<:Any,RS<:Any,BI<:Any}(op::CachedOperator{T,RM,M,DS,RS,BI},j::Integer)

#function ApproxFun.resizedata!(ApproxFun.c)
