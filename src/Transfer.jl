abstract AbstractTransfer{T} <: Operator{T}

ApproxFun.israggedbelow(L::AbstractTransfer) = true
Base.issymmetric(L::AbstractTransfer) = false

#for OP in (:(ApproxFun.domain),:rangedomain)
#    @eval $OP(L::AbstractTransfer) = $OP(getmap(L))
#end

export Transfer
Transfer(stuff...) = cache(ConcreteTransfer(stuff...))

# fast colstop for transfer operators
function ApproxFun.colstop{T<:Number,M<:AbstractTransfer,DS,RS}(B::ApproxFun.CachedOperator{T,RaggedMatrix{T},M,DS,RS,Tuple{Infinity{Bool}}},k::Integer)
  ApproxFun.resizedata!(B,:,k)
  B.data.cols[k+1] - B.data.cols[k]
end

immutable ConcreteTransfer{T,D<:Space,R<:Space,M<:AbstractMarkovMap} <: AbstractTransfer{T}
  m::M
  domainspace::D #Domain space
  rangespace::R
  colstops::Array{Int,1}
  #    f::Array{Vector{T},1}
  function ConcreteTransfer(m::AbstractMarkovMap,domainspace::Space,rangespace::Space)
    @assert domain(m)==domain(domainspace)
    @assert rangedomain(m)==domain(rangespace)
    new{eltype(m),typeof(domainspace),typeof(rangespace),typeof(m)}(m,domainspace,rangespace,eltype(M)[])
  end
end

ConcreteTransfer{T}(::Type{T},m::AbstractMarkovMap,dom::Space=Space(domain(m)),ran::Space=Space(rangedomain(m)))=
  ConcreteTransfer{T,typeof(dom),typeof(ran),typeof(m)}(m,dom,ran)#,domainspace(m),rangespace(m))

for OP in (:domainspace,:rangespace)
  @eval ApproxFun.$OP(L::ConcreteTransfer) = L.$OP
end

getmap(L::ConcreteTransfer) = L.m

Base.show(io::IO,L::AbstractTransfer) = print(io,typeof(L)) #temporary

Transfer{T}(::Type{T},m::AbstractMarkovMap) = cache(ConcreteTransfer(T,m),padding=true)

function extendcolstops!(L::ConcreteTransfer,inds,vals) #should be generic range/collection I guess?
  colstopslength = length(L.colstops)
  max_ind = maximum(inds)
  if max_ind > colstopslength
    append!(L.colstops,-ones(Int,max_ind-colstopslength)) #-1 is clunky but w/e
  end

  L.colstops[inds] = vals
end

function transferfunction{TT,D}(x,L::ConcreteTransfer{TT,Chebyshev{D}},kk::Integer,T)
  y = zero(x);
  fn(xx) = cos(kk*acos(tocanonical(domainspace(L),xx))) #this probably introduces lots of roundoff error??? oh well

  for i = 1:Base.length(L.m)
    y += getmap(L)[:dvdx][i](x).*fn(getmap(L)[:v][i](x))
  end;
  abs(y) < 200eps(one(y)) ? zero(y) : y
end

function transferfunction{TT,D}(x,L::ConcreteTransfer{TT,Fourier{D}},kk::Integer,T)
  y = zero(x);
  if rem(kk,2) == 1
    fnh = cos
  else
    fnh = sin
  end

  #     if rem(kk,2) == 1
  #         fn(xx) = cos(div(kk-1,2)*tocanonical(domainspace(L),xx))
  #     else
  #         fn(xx) = sin(div(kk,2)*tocanonical(domainspace(L),xx))
  #     end

  for i = 1:Base.length(L.m)
    y += getmap(L)[:dvdx][i](x).*fnh(fld(kk,2)*tocanonical(domainspace(L),getmap(L)[:v][i](x)))
  end;
  abs(y) < 20kk*eps(one(y)) ? zero(y) : y

end

function default_transferfunction{TT,DD}(x,L::ConcreteTransfer{TT,DD},kk::Integer,T)
  y = zero(x);
  fn = Fun([zeros(T,kk-1);one(T)],domainspace(L));

  for i = 1:Base.length(L.m)
    y += getmap(L)[:dvdx][i](x).*fn(getmap(L)[:v][i](x))
  end;
  abs(y) < 200eps(one(y)) ? zero(y) : y
end
transferfunction(x,L,kk,T) = default_transferfunction(x,L,kk,T)

# Indexing

function transfer_getindex{T,D,R,M<:InverseDerivativeMarkovMap}(L::ConcreteTransfer{T,D,R,M},jdat::Tuple{Integer,Integer,Union{Integer,Infinity{Bool}}},k::Range)
  #T = eltype(getmap(L))
  dat = Array(T,0)
  cols = Array(eltype(k),Base.length(k)+1)
  cols[1] = 1
  K = isfinite(jdat[3]) ? length(jdat[1]:jdat[2]:jdat[3]) : 0 # largest colstop
  colstops = Array(Int,length(k))

  for (kind,kk) in enumerate(k)
    tk = Fun(x->transferfunction(x,L,kk,T),rangespace(L))
    tol = 4ceil(log2(length(tk.coefficients)))*eps(T) # 4 is dummy should be like 1 in Fourier, 2C+1 in Cheby
    chop!(tk,tol)

    ltkc = length(tk.coefficients)
    for i = 1:ltkc
      abs(tk.coefficients[i])<tol && (tk.coefficients[i] = 0)
    end

    tkc = tk.coefficients[(jdat[1]:jdat[2]:min(ltkc,jdat[3]))::Range]
    append!(dat,tkc)
    cols[kind+1] = cols[kind]+length(tkc)
    K = max(K,length(tkc))
    colstops[kind] = ltkc
  end
  extendcolstops!(L,k,colstops)
  RaggedMatrix(dat,cols,K)
end

Base.getindex(L::ConcreteTransfer,j::Range,k::Range) = transfer_getindex(L,(start(j),step(j),last(j)),k)
Base.getindex(L::ConcreteTransfer,j::Colon,k::Range) = transfer_getindex(L,(1,1,ApproxFun.∞),k)
Base.getindex(L::ConcreteTransfer,j::ApproxFun.AbstractCount,k::Range) = transfer_getindex(L,(start(j),step(j),ApproxFun.∞),k)
Base.getindex(L::ConcreteTransfer,j::Integer,k::Integer) = Base.getindex(L,j:j,k:k)[1,1]
Base.getindex(L::ConcreteTransfer,j::Integer,k::Range) = Base.getindex(L,j:j,k).data
Base.getindex(L::ConcreteTransfer,j::Range,k::Integer) = Base.getindex(L,j,k:k).data

function ApproxFun.colstop(L::ConcreteTransfer,k::Integer)
  k <= length(L.colstops) && L.colstops[k] != -1 && return L.colstops[k]
  length(transfer_getindex(L,(1,1,ApproxFun.∞),k:k).data)
end
