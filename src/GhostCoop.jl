module GhostCoop
  using Base, ApproxFun#, Compat, FastTransforms

import Base: values,getindex,setindex!,*,.*,+,.+,-,.-,==,<,<=,>,|,
                >=,./,/,.^,^,\,∪,transpose, size
import ApproxFun: domainspace, rangespace


include("MarkovMap.jl")
include("TransferOperator.jl")

end # module
