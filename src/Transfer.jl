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
#    f::Array{Vector{T},1}
    function ConcreteTransfer(m::AbstractMarkovMap,domainspace::Space,rangespace::Space)
        @assert domain(m)==domain(domainspace)
        @assert rangedomain(m)==domain(rangespace)
        new{eltype(m),typeof(domainspace),typeof(rangespace),typeof(m)}(m,domainspace,rangespace)
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

function transferfunction(x,L::ConcreteTransfer,kk::Number,T)
    y = zero(x);
    fn = Fun([zeros(T,kk-1);one(T)],domainspace(L));
    for i = 1:Base.length(L.m)
        y += getmap(L)[:dvdx][i](x).*fn(getmap(L)[:v][i](x))
    end;
    abs(y) < 200eps(one(y)) ? zero(y) : y
 #   abs(y) < 200eps(one(y)) && zero(y) : y
end
# Indexing
#function gind{T,D<:Space,R<:Space}(L::ConcreteTransfer{T,D,R,InverseDerivativeMarkovMap},j::Colon,k::Range)
function transfer_getindex{T,D,R,M<:InverseDerivativeMarkovMap}(L::ConcreteTransfer{T,D,R,M},jdat::Tuple{Integer,Integer,Union{Integer,Infinity{Bool}}},k::Range)
    #T = eltype(getmap(L))
    dat = Array(T,0)
    cols = Array(eltype(k),Base.length(k)+1)
    cols[1] = 1
    K = isfinite(jdat[3]) ? length(jdat[1]:jdat[2]:jdat[3]) : 0 # largest colstop
    for (kind,kk) in enumerate(k)
        #            tk = Fun(zero(T),rangespace(L))
        #            for i = 1:Base.length(L.m)
        #                tk += Fun(x->L.m.dvdx[i](x)*Fun([zeros(T,kk-1);one(T)],domainspace(L))(L.m.v[i](x)),rangespace(L))
        #            end
        tk = Fun(x->transferfunction(x,L,kk,T),rangespace(L))
        ##            @profile (tk = Fun(x->transferfunction(x,L,kk,T),rangespace(L)))
        tkc = tk.coefficients[(jdat[1]:jdat[2]:min(length(tk.coefficients),jdat[3]))::Range]
        append!(dat,tkc)
        cols[kind+1] = cols[kind]+length(tkc)
        K = max(K,length(tkc))
    end
    RaggedMatrix(dat,cols,K)
end

Base.getindex(L::ConcreteTransfer,j::Range,k::Range) = transfer_getindex(L,(start(j),step(j),last(j)),k)
Base.getindex(L::ConcreteTransfer,j::Colon,k::Range) = transfer_getindex(L,(1,1,ApproxFun.∞),k)
Base.getindex(L::ConcreteTransfer,j::ApproxFun.AbstractCount,k::Range) = transfer_getindex(L,(start(j),step(j),ApproxFun.∞),k)
Base.getindex(L::ConcreteTransfer,j::Integer,k::Integer) = Base.getindex(L,j:j,k:k)[1,1]
Base.getindex(L::ConcreteTransfer,j::Integer,k::Range) = Base.getindex(L,j:j,k).data
Base.getindex(L::ConcreteTransfer,j::Range,k::Integer) = Base.getindex(L,j,k:k).data

ApproxFun.colstop(L::ConcreteTransfer,k::Integer) = length(transfer_getindex(L,(1,1,ApproxFun.∞),k:k).data)
