export timeseries, timehist
function map_n(m::AbstractMarkovMap,x::Number,n::Integer)
  y = copy(x)
  for i = 1:n
    y = m(x)*(1+randn()*eps(prectype(m)))
  end
  y
end

function timeseries_initialvalue(m::AbstractMarkovMap,n::Integer,xinit::T) where T
  @assert domain(m) == rangedomain(m)
  @compat x_ts = Array{T}(undef,n)#randn(n)*eps(T)
  x_ts[1] = xinit
  x = copy(xinit)
  for i = 2:n
    x = m(x)#*(1+x_ts[i])
    x_ts[i] = x
  end
  x_ts
end
timeseries(m::AbstractMarkovMap,n::Integer;x0=map_n(m,rand(domain(m)),10^4)) = timeseries_initialvalue(m,n,x0)
samplable(ρ::Fun) = isa(domain(ρ),AbstractInterval) ? ρ : Fun(ρ,Space(Interval(domain(ρ))))
        # probably Interval(..) should be a convert routine but that's how ApproxFun has it...
timeseries(m::AbstractMarkovMap,n::Integer,ρ::Fun) = timeseries_initialvalue(m,n,sample(samplable(ρ)))


function timehist_initialvalue(m::AbstractMarkovMap,n::Integer,nbins::Integer,xinit::T) where T
  @assert domain(m) == rangedomain(m)
  bin_start = leftendpoint(domain(m))
  bin_step = (rightendpoint(domain(m))-bin_start)/nbins
  x_hist = zeros(Int,nbins)
  x = copy(xinit)
  for i = 1:n
    x = m(x)#*(1+randn()*eps(T))
    x_hist[min(nbins,max(1,convert(Int,ceil((x-bin_start)/bin_step))))] += 1
  end

  range(bin_start,stop=rightendpoint(domain(m)),length=nbins+1),x_hist
end
timehist(m::AbstractMarkovMap,n::Integer,nbins::Integer;x0=map_n(m,rand(domain(m)),10^4)) = timehist_initialvalue(m,n,nbins,x0)
timehist(m::AbstractMarkovMap,n::Integer,nbins::Integer,ρ::Fun) = timehist_initialvalue(m,n,nbins,sample(samplable(ρ)))
