convert(::Type{MarkovMap},t::IntervalMap) = MarkovMap(t.branches,t.domain,t.rangedomain)
convert(::Type{IntervalMap},t::MarkovMap) = IntervalMap(t.branches,t.domain,t.rangedomain)

convert(::Type{MarkovMap},t::FwdCircleMap) = forwardmodulomap(t.f,t.domain,t.rangedomain,t.dfdx)
convert(::Type{MarkovMap},t::RevCircleMap) = reversemodulomap(t.v,t.domain,t.rangedomain,t.dvdx)
convert(::Type{IntervalMap},t::AbstractCircleMap) = IntervalMap(MarkovMap(t))
