"""
    HofbauerDomain

A `Domain` object annotated with a depth index for `HofbauerExtensions`.
It is not a `Domain` subtype itself as it should only be used in the context of Hofbauer extensions.
"""
struct HofbauerDomain{D<:Domain}
    domain::D
    depth::Int
end
show(io::IO, d::HofbauerDomain) = print(io, "$(d.domain) at depth $(d.depth)")
convert(Domain,d::HofbauerDomain) = d.domain

maxrangedomain(m::AbstractIntervalMap) =
    rangedomain(branches(m)[findmax(arclength.(rangedomain.(branches(m))))[2]])

"""
    HofbauerPoint(x, graphind)

Constructs a `HofbauerPoint`, representing point `x` on the `graphind`th domain
of a Hofbauer extension.
"""
struct HofbauerPoint{T}
    x::T
    graphind::Int
end
show(io::IO,hp::HofbauerPoint) = print(io,"$(hp.x) on domain #$(hp.graphind)")

# mapping Hofbauer extension

"""
    maphb(b::AbstractBranch, x)

Maps point `x` according to `b`. For `NeutralBranches` calculates the return map out of `domain(b)`.
"""
maphb(b::AbstractBranch,x) = b(x)
maphb(b::NeutralBranch,x) = _NI("maphb for NeutralBranch") # use domain(b) to fix.
for (M,MHB) = ((:mapD,:maphbD),(:mapP,:maphbP),(:mapinv,:maphbinv),
        (:mapinvD,:maphbinvD),(:mapinvP,:maphbinvP))
    @eval $MHB(b::AbstractBranch,x) = $M(b,x)
    @eval $MHB(b::NeutralBranch,x) = _NI("$MHB for NeutralBranch")
end


function _hofbauerdestination(he::HofbauerExtension, bind, x, graphind)
    graphinds = searchsorted(he.fedgelist[graphind], HofbauerEdge(0,0,bind,false), by=branchind)
    isempty(graphinds) && error("reached end of generated Hofbauer extension at depth $(hdomains(im.he)[graphind].depth)")
    for g in graphinds
        dst_ind = dst(im.he.fedgelist[graphind][g])
        if x ∈ convert(Domain,hdomains(he)[dst_ind])
            graphind = dst_ind
            break
        end
    end
    graphind
end
@inline _hofbauerdestination(he::HofbauerExtension, bind, p) =
    _hofbauerdestination(he, bind, p.x, p.graphind)

function (he::HofbauerExtension)(p::HofbauerPoint)
    bind = getbranchind(he.m, p.x)
    fx = maphb(branches(he.m)[bind], p.x)
    HofbauerPoint(x, _hofbauerdestination(he, bind, fx))
end
function mapP(he::HofbauerExtension, p::HofbauerPoint)
    bind = getbranchind(he.m, p.x)
    fx, dfx = maphbP(branches(he.m)[bind], p.x)
    HofbauerPoint(fx,graphind), dfx
end
mapD(he::HofbauerExtension, p::HofbauerPoint) = maphbD(he.m, p.x)

_hofbauersource(he::HofbauerExtension, bind, graphind) =
    searchsortedfirst(he.bedgelist[graphind], HofbauerEdge(0,0,bind,false), by=branchind)

# TODO: assert in rangedomain
mapinv(he::HofbauerExtension, i, p::HofbauerPoint) =
    HofbauerPoint(maphbinv(branches(he.m)[i],p.x), _hofbauersource(he,i,p.graphind))
mapinvD(he::HofbauerExtension, i, p::HofbauerPoint) =
        maphbinvD(branches(he.m)[i], p.x)
function mapinvP(he::HofbauerExtension, i, p::HofbauerPoint)
    (vx, dvx) = maphbinvP(branches(he.m)[i], p.x)
    HofbauerPoint(vx, _hofbauersource(he,i,p.graphind)), dvx
end

function mapinv(he::HofbauerExtension,e::HofbauerEdge,p::HofbauerPoint)
    @assert dst(e) == p.graphind
    HofbauerPoint(maphbinv(branches(he.m)[branchind(e)],p.x),src(e))
end
function mapinvD(he::HofbauerExtension,e::HofbauerEdge,p::HofbauerPoint)
    @assert dst(e) == p.graphind
    maphbinvD(branches(he.m)[branchind(e)],p.x)
end
function mapinvP(he::HofbauerExtension,e::HofbauerEdge,p::HofbauerPoint)
    @assert dst(e) == p.graphind
    vx,dvx = maphbinvP(branches(he.m)[branchind(e)],p.x)
    HofbauerPoint(vx,src(e)), dvx
end


"""
    hofbauerextension(m,basedomains;maxdepth=100,forcereturn=trues(length(basedomains)))

Generates a Hofbauer extension of `m` starting from `basedomains`, which may be
a single domain.

The keyword `maxdepth` says how deep the Hofbauer extension should go, and
`forcereturn` sets whether a return to given members of `basedomains` should
 forced if possible.

 For example, in the case of `f(x) = mod(2x,1)`, if `basedomains` is set to `Segment(0.,0.3)` then
`Segment(0.,0.5)` might map to `Segment(0.,1.)`` for a given map if `forcereturn=false` but would map to
Segment(0,0.3)` and `Segment(0.3,1)` if `forcereturn=true`.
"""
function hofbauerextension(m::AbstractIntervalMap,basedoms=(maxrangedomain(m));maxdepth=100,forcereturn=trues(length(basedoms)))
    @assert domain(m) == rangedomain(m)

    G = HofbauerExtension{Int,HofbauerDomain}(m)
    for (n,bd) in enumerate(basedoms)
        basehdom = Hofbauerconvert(Domain,convert(Domain,bd),1)
        add_vertex!(G,basehdom,forcereturn[n])
    end

    i = 1
    while nv(G) >= i && hdomains(G)[i].depth < maxdepth
        towerwalk!(G,i,hdomains(G)[i])
        i += 1
    end

    D = promote_type(collect(typeof(h.domain) for h in hdomains(G))...)
    # HofbauerExtension(ne(G),fedgelist(G),bedgelist(G),getmap(G),
    #     convert(Vector{HofbauerDomain{D}},hdomains(G)))
    G
end


# walks through domains in HofbauerExtension
function towerwalk!(G::HofbauerExtension, graphind, hdom)
    @assert hdomains(G)[graphind] == hdom # putting in hdom for type inference
    m = getmap(G)
    for branchind in eachbranchindex(m)
        b = branches(m)[branchind]
        bdom = domain(b)
        capdom = bdom ∩ hdom.domain
        if ~isempty(capdom) # hdom is partially contained in the domain of b
            if isa(b,NeutralBranch) #&& (neutralfixedpoint(b) ∈ capdom) # inducing at a neutral fixed point
                isnfp = true
                newdom = rangedomain(b) \ ((bdom ≈ capdom) ? domain(b) : b(capdom))
            else # normal bounce inducing
                isnfp = false
                newdom = (bdom ≈ capdom) ? rangedomain(b) : b(capdom)
            end
            for i in G.returnto
                forcedom = convert(Domain,hdomains(G)[i])
                if issubset(forcedom,newdom)
                    towerextend!(G, graphind, hdom, forcedom, branchind, isnfp)
                    newdom = newdom \ forcedom
                end
            end
            towerextend!(G, graphind, hdom, newdom, branchind, isnfp)
        end
    end
end

function towerextend!(G::HofbauerExtension, graphind, hdom, newdom, branchind, isnfp)
    arclength(newdom) < 300eps(arclength(convert(Domain,hdom))) && return true
    newgraphind = findfirst(collect(h.domain≈newdom for h in hdomains(G)))
    addgraphind = (newgraphind == 0) # does the new domain already exist

    # add new domain if necessary
    if addgraphind
        newhdom = Hofbauerconvert(Domain,newdom,hdom.depth+1)
        add_vertex!(G,newhdom)
        newgraphind = nv(G)
    end

    # adding new edge
    newedge = HofbauerEdge(graphind, newgraphind, branchind, isnfp)
    add_edge!(G,newedge)
end
function towerextend!(G::HofbauerExtension, graphind, hdom, newdom::UnionDomain, branchind, isnfp)
    for d in newdom.domains
        towerextend!(G, graphind, hdom, d, branchind, isnfp)
    end
end
towerextend!(G::HofbauerExtension, graphind, hdom, newdom::ApproxFun.EmptyDomain, branchind, isnfp) = true
