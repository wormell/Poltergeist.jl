struct InducedMap{M<:AbstractIntervalMap, D<:Domain, R<:Domain, HD, I} <: AbstractMarkovMap{D,R}
    he::HofbauerExtension{I,HD,M}
    domain::D
    rangedomain::R
    d_ind::Int
    r_ind::Int

    function InducedMap(he::HofbauerExtension{I,HD,M}, d::D, rd::R, d_ind::Int, r_ind::Int) where
            {M<:AbstractIntervalMap,D<:Domain,R<:Domain,HD,I}
    @assert d == Domain(hdomains(he)[d_ind])
    @assert rd == Domain(hdomains(he)[r_ind])
    new{M,D,R,HD,I}(he,d,rd,d_ind,r_ind)
    end
end
function InducedMap(he::HofbauerExtension, domain::Domain, rangedomain::Domain)
    @assert domain ∈ domains.(he)
    @assert rangedomain ∈ domains.(he)
    InducedMap(he, domain, rangedomain, findfirst(domains.(he),domain),
        findfirst(domains.(he),rangedomain))
end

InducedMap(he::HofbauerExtension,d_ind::Int,r_ind::Int) =
    InducedMap(he,Domain(hdomains(he)[d_ind]),Domain(hdomains(he)[r_ind]),d_ind,r_ind)
InducedMap(he::HofbauerExtension,d,r) = InducedMap(he,Domain(d),Domain(r))
InducedMap(he::HofbauerExtension,d) = InducedMap(he,d,d)
InducedMap(he::HofbauerExtension) = InducedMap(he,1,1)

InducedMap(m::AbstractIntervalMap;maxdepth=100,forcereturn=true) = InducedMap(hofbauerextension(m;maxdepth=maxdepth,forcereturn=forcereturn))
InducedMap(m::AbstractIntervalMap,d;maxdepth=100,forcereturn=trues(length(d))) = InducedMap(hofbauerextension(m,d;maxdepth=maxdepth,forcereturn=forcereturn),d)
#length(d::Interval) not defined so:
InducedMap(m::AbstractIntervalMap,d::ClosedInterval;maxdepth=100,forcereturn=true) = InducedMap(hofbauerextension(m,Domain(d);maxdepth=maxdepth,forcereturn=forcereturn),d)

induce(m::AbstractIntervalMap,args...) = InducedMap(m)

Base.show(io::IO, i::InducedMap) = print(io, "Induced map on $(i.domain)→$(i.rangedomain) of $(i.he.m)")

map_inducing(b::AbstractBranch,x) = b(x)
map_inducing(b::NeutralBranch,x) = _NI("map_inducing for NeutralBranch")

function (im::InducedMap)(x;maxiter=10000)
    @assert x ∈ domain(im)
    m = getmap(im.he)
    graphind = im.d_ind
    for i = 1:maxiter
        bind = getbranchind(m,x)
        x = map_inducing(branches(m)[bind],x)
        graphinds = searchsorted(im.he.fedgelist[graphind],HofbauerEdge(0,0,bind,false),by=branchind)
        isempty(graphinds) && error("reached end of generated Hofbauer extension at depth $(hdomains(im.he)[graphind].depth)")
        for g in graphinds
            dst_ind = dst(im.he.fedgelist[graphind][g])
            if x ∈ Domain(hdomains(im.he)[dst_ind])
                graphind = dst_ind
                break
            end
        end
        graphind == im.r_ind && return x
    end
    error("reached max number of iterations: $maxiter")
end
