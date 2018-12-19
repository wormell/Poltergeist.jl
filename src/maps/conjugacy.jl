# MarkovMaps
function inv(m::MarkovMap)
  ~isa(domain(m),AbstractInterval) && error("Non-interval domains not supported")
  nbranches(m) > 1 && error("Cannot invert map with multiple branches")
  domain(m.branches[1]) != domain(m) && error("Cannot invert map that isn't defined everywhere")
  unsafe_inv(m)
end

function unsafe_inv(m::MarkovMap{D,R,B}) where {D,R} where B<:FwdExpandingBranch
  b = m.branches[1]
  b2 = RevExpandingBranch(b.f,b.dfdx,b.rangedomain,b.domain)
  MarkovMap([b2],m.rangedomain,m.domain)
end

function unsafe_inv(m::MarkovMap{D,R,B}) where {D,R} where B<:RevExpandingBranch
  b = m.branches[1]
  b2 = FwdExpandingBranch(b.v,b.dvdx,b.rangedomain,b.domain)
  MarkovMap([b2],m.rangedomain,m.domain)
end

# ComposedMarkovMaps
inv(m::ComposedMarkovMap) = ComposedMarkovMap([inv(mm) for mm in reverse(m.maps)],m.rangedomain,m.domain)

# CircleMaps
function inv(m::AbstractCircleMap)
  ncover(m) > 1 && error("Cannot invert map with covering index > 1")
  unsafe_inv(m)
end
unsafe_inv(m::FwdCircleMap) = RevCircleMap(m.f,m.dfdx,m.rangedomain,m.domain)
unsafe_inv(m::RevCircleMap) = FwdCircleMap(m.v,m.dvdx,m.rangedomain,m.domain)
