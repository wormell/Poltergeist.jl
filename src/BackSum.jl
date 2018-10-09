struct BackSum{I<:InducedMap}
    im::I
end

struct BackSumEval{I<:InducedMap,Q}
    S::BackSum{I}
    q::Q
    domainonly::Bool
end
(S::BackSum)(Q, domainonly::Bool=false) = BackSumEval(S, Q, domainonly)

(SQ::BackSumEval)(hx::HofbauerPoint{T}; maxdepth=1000, dvmin=eps(T)) where T =
    backsumrecursion(SQ, hx, one(T), 0, maxdepth, dvmin)
(SQ::BackSumEval)(x; kwargs...) =
    SQ(HofbauerPoint(x,SQ.S.im.r_ind); kwargs...)

function backsumrecursion(SQ::BackSumEval, x::HofbauerPoint{T}, dv, n, depthleft, dvmin) where T
    im = SQ.S.im
    backsum = zero(T)
    # TODO: accomodate NeutralBranch
    for e in bedgelist(im.he, x.graphind)
        (nx, dvprod) = mapinvP(im.he, e, x)
        ndv = abs(dv)*dvprod # check this is OK for neutral branches
        is_d_ind = nx.graphind == im.d_ind
        (~SQ.domainonly || is_d_ind) && (backsum += SQ.q(nx.x,ndv,n))
        abs(ndv) ≥ dvmin && nx.graphind != im.r_ind && ~is_d_ind && depthleft > 0 &&
            (backsum += backsumrecursion(SQ, nx, ndv, n+1, depthleft-1, dvmin))
    end
    backsum
end
