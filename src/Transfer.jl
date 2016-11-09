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

chebyTk(x,d,k::Integer) = cos((k-1)*acos(tocanonical(d,x))) #roundoff error grows linearly(??) with k may not be bad wrt x too
function transferfunction{TT,D}(x,L::ConcreteTransfer{TT,Chebyshev{D}},k::Integer,T)
  y = zero(x);

  for i = 1:length(L.m)
    y += getmap(L)[:dvdx][i](x).*chebyTk(getmap(L)[:v][i](x),domain(L),k)
  end;
  abs(y) < 20k*eps(one(abs(y))) ? zero(y) : y
end

fourierCSk(x,d,k::Integer) = rem(k,2) == 1 ? cos(fld(k,2)*tocanonical(d,x)) : sin(fld(k,2)*tocanonical(d,x))
function transferfunction{TT,D}(x,L::ConcreteTransfer{TT,Fourier{D}},k::Integer,T)
  y = zero(x);

  for i = 1:Base.length(L.m)
    y += getmap(L)[:dvdx][i](x).*fourierCSk(getmap(L)[:v][i](x),domain(L),k)
  end;
  abs(y) < 20k*eps(one(abs(y))) ? zero(y) : y

end


function default_transferfunction{TT,DD}(x,L::ConcreteTransfer{TT,DD},kk::Integer,T)
  y = zero(x);
  fn = Fun([zeros(T,kk-1);one(T)],domainspace(L));

  for i = 1:Base.length(L.m)
    y += getmap(L)[:dvdx][i](x).*fn(getmap(L)[:v][i](x))
  end;
  abs(y) < 200eps(one(abs(y))) ? zero(y) : y
end
transferfunction(x,L,kk,T) = default_transferfunction(x,L,kk,T)

# Indexing

function transfer_getindex{T,D,R,M<:InverseDerivativeMarkovMap}(L::ConcreteTransfer{T,D,R,M},jdat::Tuple{Integer,Integer,Union{Integer,Infinity{Bool}}},k::Range)
  #T = eltype(getmap(L))
  dat = Array(T,0)
  cols = Array(eltype(k),Base.length(k)+1)
  cols[1] = 1
  K = isfinite(jdat[3]) ? length(jdat[1]:jdat[2]:jdat[3]) : 0 # largest colstop

  rs = rangespace(L)

  resizecolstops!(L,maximum(k))
  mc = maximum(L.colstops[1:first(k)])

  for (kind,kk) in enumerate(k)
    kind > 1 && (mc = max(mc,maximum(L.colstops[k[kind-1]+1:kk])))

    f(x) = transferfunction(x,L,kk,T)
    tol =T==Any?200eps():200eps(T)

    if L.colstops[kk] >= 1
      cf = Fun(f,rs,convert(Int,2^max(4,ceil(log2(L.colstops[kk])))))
    elseif L.colstops[kk]  == 0
      cf = zeros(rs)

    else

      if mc ≤ 2^4
        cf = Fun(f,rs)
      else
        # code from ApproxFun (src/Fun/constructors.jl) - the difference is we start at n \approx mc

        r=ApproxFun.checkpoints(rs)
        fr=map(f,r)
        maxabsfr=norm(fr,Inf)

        logn = min(convert(Int,ceil(log2(mc))),20)
        while logn < 21
          cf = Fun(f,rs,2^logn)

          maxabsc = maxabs(cf.coefficients)
          if maxabsc == 0 && maxabsfr == 0
            cf = (zeros(rs))
            break
          else
            b = ApproxFun.block(rs,length(cf.coefficients))
            bs = ApproxFun.blockstart(rs,max(b-2,1))
            if ncoefficients(cf) > 8 && maxabs(cf.coefficients[bs:end]) < tol*maxabsc &&
                all(kkk->norm(cf(r[kkk])-fr[kkk],1)<tol*length(cf.coefficients)*maxabsfr*1000,1:length(r))
              cf =  chop!(cf,tol*maxabsc/10)
              break
            end
          end
          logn += 1
        end

        if logn == 21
          warn("Maximum number of coefficients "*string(2^20+1)*" reached in constructing Fun.")
          Fun(f,rs,2^21)
        end

      end
    end

    #    tol = 4ceil(log2(length(tk.coefficients)))*eps(T) # 4 is dummy should be like 1 in Fourier, 2C+1 in Cheby
    #    chop!(tk,tol)
    maxabsc = maxabs(cf.coefficients)
    lcfc = length(cf.coefficients)

    for i = 1:lcfc
      abs(cf.coefficients[i])<tol*maxabsc/10 && (cf.coefficients[i] = 0)
    end

    cutcfc = cf.coefficients[(jdat[1]:jdat[2]:min(lcfc,jdat[3]))::Range]

    append!(dat,cutcfc)
    cols[kind+1] = cols[kind]+length(cutcfc)
    K = max(K,length(cutcfc))

    setcolstops!(L,kk,lcfc)
 #   kk == 1 && display(stacktrace())
  end
  RaggedMatrix(dat,cols,K)
end

Base.getindex(L::ConcreteTransfer,j::Range,k::Range) = transfer_getindex(L,(start(j),step(j),last(j)),k)
Base.getindex(L::ConcreteTransfer,j::Colon,k::Range) = transfer_getindex(L,(1,1,ApproxFun.∞),k)
Base.getindex(L::ConcreteTransfer,j::ApproxFun.AbstractCount,k::Range) = transfer_getindex(L,(start(j),step(j),ApproxFun.∞),k)
Base.getindex(L::ConcreteTransfer,j::Integer,k::Integer) = Base.getindex(L,j:j,k:k)[1,1]
Base.getindex(L::ConcreteTransfer,j::Integer,k::Range) = Base.getindex(L,j:j,k).data
Base.getindex(L::ConcreteTransfer,j::Range,k::Integer) = Base.getindex(L,j,k:k).data

function ApproxFun.default_raggedmatrix{T,LL<:AbstractTransfer,R1<:Union{Range,ApproxFun.AbstractCount},R2<:Range}(
  S::ApproxFun.SubOperator{T,LL,Tuple{R1,R2}})
 Base.getindex(parent(S),parentindexes(S)[1],parentindexes(S)[2])
end



function ApproxFun.colstop(L::ConcreteTransfer,k::Integer)
  k <= length(L.colstops) && L.colstops[k] != -1 && return L.colstops[k]
  transfer_getindex(L,(1,1,ApproxFun.∞),k:k)
  L.colstops[k]
end
