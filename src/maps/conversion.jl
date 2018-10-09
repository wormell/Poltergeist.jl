MarkovMap(t::IntervalMap) = MarkovMap(t.branches,t.domain,t.rangedomain)
IntervalMap(t::MarkovMap) = IntervalMap(t.branches,t.domain,t.rangedomain)

MarkovMap(t::FwdCircleMap) = forwardmodulomap(t.f,t.domain,t.rangedomain,t.dfdx)
MarkovMap(t::RevCircleMap) = reversemodulomap(t.v,t.domain,t.rangedomain,t.dvdx)
IntervalMap(t::AbstractCircleMap) = IntervalMap(MarkovMap(t))
