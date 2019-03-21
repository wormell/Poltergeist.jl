MarkovMap(t::IntervalMap) = MarkovMap(t.branches,t.domain,t.rangedomain)
IntervalMap(t::MarkovMap) = IntervalMap(t.branches,t.domain,t.rangedomain)

MarkovMap(t::FwdCircleMap) = forwardmodulomap(t.f,Interval(t.domain),Interval(t.rangedomain),t.dfdx)
MarkovMap(t::RevCircleMap) = reversemodulomap(t.v,Interval(t.domain),Interval(t.rangedomain),t.dvdx)
IntervalMap(t::AbstractCircleMap) = IntervalMap(MarkovMap(t))
