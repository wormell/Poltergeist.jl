"""
    HofbauerDomain
An ApproxFun Domain object annotated with a depth index for HofbauerExtensions.
It is not a Domain subtype itself as it should only be used in the context of Hofbauer extensions.
"""
struct HofbauerDomain{D<:Domain}
    domain::D
    depth::Int
end
show(io::IO, d::HofbauerDomain) = print(io, "$(d.domain) at depth $(d.depth)")
Domain(d::HofbauerDomain) = d.domain

maxrangedomain(m::AbstractIntervalMap) =
    rangedomain(branches(m)[findmax(arclength.(rangedomain.(branches(m))))[2]])

function hofbauerextension(m::AbstractIntervalMap,basedoms=(maxrangedomain(m));maxdepth=100,forcereturn=trues(length(basedoms)))
    @assert domain(m) == rangedomain(m)

    G = HofbauerExtension{Int,HofbauerDomain}(m)
    for (n,bd) in enumerate(basedoms)
        basehdom = HofbauerDomain(Domain(bd),1)
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
                forcedom = Domain(hdomains(G)[i])
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
    arclength(newdom) < 300eps(arclength(Domain(hdom))) && return true
    newgraphind = findfirst(collect(h.domain≈newdom for h in hdomains(G)))
    addgraphind = (newgraphind == 0) # does the new domain already exist

    # add new domain if necessary
    if addgraphind
        newhdom = HofbauerDomain(newdom,hdom.depth+1)
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
