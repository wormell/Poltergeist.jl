export Transfer, transfer#, colstop, resizedata!

@compat abstract type AbstractTransfer{T} <: Operator{T} end

ApproxFun.israggedbelow(L::AbstractTransfer) = true
#Base.issymmetric(L::AbstractTransfer) = false

#for OP in (:(ApproxFun.domain),:rangedomain)
#    @eval $OP(L::AbstractTransfer) = $OP(markovmap(L))
#end

@compat struct ConcreteTransfer{T,D<:Space,R<:Space,M<:AbstractMarkovMap} <: AbstractTransfer{T}
  m::M
  domainspace::D #Domain space
  rangespace::R
  colstops::Array{Int,1}
  #    f::Array{Vector{T},1}
  # function ConcreteTransfer{T,D<:Space,R<:Space,M<:AbstractMarkovMap}(m::AbstractMarkovMap,domainspace::Space,rangespace::Space,colstops::Array{Int,1}=Int[])
  #   @assert domain(m)==domain(domainspace)
  #   @assert rangedomain(m)==domain(rangespace)
  #   nneutral(m) != 0 && error("Neutral fixed points not supported for transfer operator")
  #   new{T,typeof(domainspace),typeof(rangespace),typeof(m)}(m,domainspace,rangespace,colstops)
  # end
end

"""
    Transfer(m::AbstractMarkovMap)

Create a `CachedOperator` of a `ConcreteTransfer` operator encoding the transfer operator of `m`.

Caching is used for speed, as entries of the transfer operator are most efficiently calculated whole columns at a time.
"""
Transfer(stuff...;padding=false) = cache(ConcreteTransfer(stuff...),padding=padding)

# ConcreteTransfer(::Type{T},m::AbstractMarkovMap,dom::Space=Space(domain(m)),
#                     ran::Space=(domain(m)==rangedomain(m) ? dom : Space(rangedomain(m))),
#                     colstops::Array{Int,1}=Int[]) where T =
#   ConcreteTransfer{T,typeof(dom),typeof(ran),typeof(m)}(m,dom,ran,colstops)#,domainspace(m),rangespace(m))
ConcreteTransfer(m,dom::Space=Space(domain(m)),ran=(domain(m)==rangedomain(m) ? dom : Space(rangedomain(m))),colstops::Array{Int,1}=Int[]) =
  ConcreteTransfer{prectype(rangedomain(m)),typeof(dom),typeof(ran),typeof(m)}(m,dom,ran,colstops)

for OP in (:domainspace,:rangespace)
  @eval ApproxFun.$OP(L::ConcreteTransfer) = L.$OP
end

function Base.convert(::Type{Operator{T}},D::ConcreteTransfer) where T
    if T==eltype(D)
        D
    else
        ConcreteTransfer{T,typeof(D.domainspace),typeof(D.rangespace),typeof(D.m)}(D.m,D.domainspace,D.rangespace,D.colstops)
    end
end


markovmap(L::ConcreteTransfer) = L.m
markovmap(L::CachedOperator) = markovmap(L.op)

function resizecolstops!(L::ConcreteTransfer,n)
  colstopslength = length(L.colstops)
  if n > colstopslength
    append!(L.colstops,-ones(Int,n-colstopslength)) #-1 is clunky but w/e
  end
end

function setcolstops!(L::ConcreteTransfer,inds,vals) #should be generic range/collection I guess?
  resizecolstops!(L,maximum(inds))
  L.colstops[inds] = vals
end



# Transfer a fun - TODO
function transfer(m::AbstractIntervalMap,fn)
  @inline tf(x) = transferfunction(x,m,fn)
  Fun(tf,rangespace(m)) # TODO: MarkovMaps don't have spaces
end
transfer(m::AbstractIntervalMap,fn,x) = transferfunction(x,m,fn)

# Indexing

transferfunction_nodes(L::ConcreteTransfer{TT,D,R,M},n::Integer,kk) where {TT,D,R,M<:AbstractIntervalMap} =
  prectype(domainspace(L))[transferfunction(p,markovmap(L),BasisFun(domainspace(L),kk)) for p in points(rangespace(L),n)]
# transferfunction_nodes{TT,D,R,M<:MarkovInverseCache}(L::ConcreteTransfer{TT,D,R,M},n::Integer,kk) =
#   T[transferfunction(InterpolationNode(rangespace(markovmap(L)),k,n),markovmap(L),BasisFun(domainspace(L),kk)) for k = 1:n]


function transfer_getindex(L::ConcreteTransfer{T},
    jdat::Tuple{Integer,Integer,Union{Integer,Infinity{Bool}}},
    k::AbstractRange,padding::Bool=false) where T
  Tr = real(T)
  @compat dat = Array{T}(undef,0)
  @compat cols = Array{eltype(k)}(undef,Base.length(k)+1)
  cols[1] = 1
  K = isfinite(jdat[3]) ? length(jdat[1]:jdat[2]:jdat[3]) : 0 # largest colstop
  padding && start(k) > 1 && (K = max(maximum([colstop(L,i) for i =1:start(k)-1]),K))
  rs = rangespace(L)

  resizecolstops!(L,maximum(k))
  mc = maximum(L.colstops[1:first(k)])

  for (kind,kk) in enumerate(k)
    kind > 1 && (mc = max(mc,maximum(L.colstops[k[kind-1]+1:kk])))

    tol = (Tr==Any) ? 200eps() : 200eps(Tr)

    if L.colstops[kk] >= 1
      @compat coeffs = ApproxFun.transform(rs,transferfunction_nodes(L,max(16,nextpow(2,L.colstops[kk])),kk))[1:L.colstops[kk]]
      maxabsc = max(maximum(abs.(coeffs)),one(Tr))
      chop!(coeffs,tol*maxabsc*log2(length(coeffs)))
    elseif L.colstops[kk]  == 0
      coeffs = zeros(T,1)
    else

      #       if mc ≤ 2^4
      #         coeffs = Fun(f,rs).coefficients
      #       else

      # code from ApproxFun v0.5 (src/Fun/constructors.jl) under BSD license (same as this package)
      # ( the difference is we start at n \approx mc )

      r=ApproxFun.checkpoints(rs)
      fr=[transferfunction(rr,markovmap(L),BasisFun(domainspace(L),kk)) for rr in r]
      maxabsfr=norm(fr,Inf)

      logn = min(round(Int,log2(max(mc,16)),RoundUp),20)
      while logn < 21
        coeffs = ApproxFun.transform(rs,transferfunction_nodes(L,2^logn,kk))

        maxabsc = max(one(Tr),maximum(abs.(coeffs)))
        if maxabsc == 0 && maxabsfr == 0
          coeffs = zeros(T,1)
          break
        else
          maxabsfr = max(maxabsfr,one(Tr))
          b = ApproxFun.block(rs,length(coeffs))
          bs = ApproxFun.blockstart(rs,max(div(2Integer(b),3),1))
          if length(coeffs) > 8 && maximum(abs.(coeffs[bs:end])) < 20tol*maxabsc*logn &&
              all(kkk->norm(Fun(rs,coeffs)(r[kkk])-fr[kkk],1)<tol*length(coeffs)*maxabsfr*500,1:length(r))
            chop!(coeffs,tol*maxabsc*logn)
            break
          end
        end
        logn += 1
      end

      if logn == 21
        warn("Maximum number of coefficients "*string(2^20+1)*" reached in constructing $(kk)th column.")
        coeffs = ApproxFun.transform(rs,transferfunction_nodes(L,2^logn,kk))
      end

    end

    #    tol = 4ceil(log2(length(tk.coefficients)))*eps(T) # 4 is dummy should be like 1 in Fourier, 2C+1 in Cheby
    #    chop!(tk,tol)
    lcfc = length(coeffs)
    maxabsc = lcfc > 0 ? maximum(abs.(coeffs)) : zero(eltype(coeffs))

    for i = 1:lcfc # set small enough coefficients (in middle of vector) to zero
      abs(coeffs[i])<tol*maxabsc*log2(lcfc)/10 && (coeffs[i] = 0)
    end

    cutcfc = coeffs[(jdat[1]:jdat[2]:min(lcfc,jdat[3]))::AbstractRange]

    K = max(K,length(cutcfc))
    padding && ApproxFun.pad!(cutcfc,K)
    append!(dat,cutcfc)
    cols[kind+1] = cols[kind]+length(cutcfc)

    setcolstops!(L,kk,lcfc)
    #   kk == 1 && display(stacktrace())
  end
  RaggedMatrix(dat,cols,K)
end

Base.getindex(L::ConcreteTransfer,j::AbstractRange,k::AbstractRange) = transfer_getindex(L,(start(j),step(j),last(j)),k)
Base.getindex(L::ConcreteTransfer,j::Colon,k::AbstractRange) = transfer_getindex(L,(1,1,ApproxFun.∞),k)
Base.getindex(L::ConcreteTransfer,j::ApproxFun.AbstractCount,k::AbstractRange) = transfer_getindex(L,(start(j),step(j),ApproxFun.∞),k)
Base.getindex(L::ConcreteTransfer,j::Integer,k::Integer) = Base.getindex(L,j:j,k:k)[1,1]
Base.getindex(L::ConcreteTransfer,j::Integer,k::AbstractRange) = Base.getindex(L,j:j,k).data
Base.getindex(L::ConcreteTransfer,j::AbstractRange,k::Integer) = Base.getindex(L,j,k:k).data

Base.convert(::Type{RaggedMatrix},S::ApproxFun.SubOperator{T,LL,Tuple{R1,R2}}) where
        {T,LL<:AbstractTransfer,R1<:Union{AbstractRange,ApproxFun.AbstractCount},R2<:AbstractRange} =
  Base.getindex(parent(S),parentindexes(S)[1],parentindexes(S)[2])



function ApproxFun.colstop(L::ConcreteTransfer,k::Integer)
  k <= length(L.colstops) && L.colstops[k] != -1 && return L.colstops[k]
  transfer_getindex(L,(1,1,ApproxFun.∞),k:k)
  L.colstops[k]
end

function ApproxFun.resizedata!(co::ApproxFun.CachedOperator{T,RaggedMatrix{T},CT},
    ::Colon,n::Integer) where {T<:Number,CT<:ConcreteTransfer}
  if n > co.datasize[2]
    RO = transfer_getindex(co.op,(1,1,ApproxFun.∞),(co.datasize[2]+1):n,co.padding)
    if co.datasize[2] == 0
      co.data = RO
      co.datasize = (co.data.m,n)
    else
      append!(co.data.data,RO.data)
      append!(co.data.cols,RO.cols[2:end].+(co.data.cols[end]-1))
      co.data.m = max(co.data.m,RO.m)
      co.datasize = (co.data.m,n)
    end
  end

  co
end

function ApproxFun.colstop(co::ApproxFun.CachedOperator{T,DM,CT},n::Integer) where
    {T<:Number,DM<:AbstractMatrix,CT<:ConcreteTransfer}
  ApproxFun.resizedata!(co,:,n)
  ApproxFun.colstop(co.data,n)
end

# # fast colstop for cached operators
# function ApproxFun.colstop{T,M<:ConcreteTransfer,DS<:Space,RS<:Space}(B::ApproxFun.CachedOperator{T,RaggedMatrix{T},M,DS,RS,Tuple{Infinity{Bool}}},k::Integer)
#   ApproxFun.resizedata!(B,:,k)
#   ApproxFun.colstop(co.op,k)
# end
