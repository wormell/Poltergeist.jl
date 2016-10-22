rangedomain(op::Operator) = domain(rangespace(op))
#Base.isless(a::ApproxFun.Infinity{Bool},b::ApproxFun.Infinity{Bool}) = a.angle && ~b.angle

#function ApproxFun.colstop{T,RM<:ApproxFun.RaggedMatrix,M<:ApproxFun.Operator,DS<:Any,RS<:Any,BI<:Any}(op::CachedOperator{T,RM,M,DS,RS,BI},j::Integer)
