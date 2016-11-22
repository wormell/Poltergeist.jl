function map_n(m::AbstractMarkovMap,x::Number,n::Integer)
  y = copy(x)
  for i = 1:n
    y = m(x)
  end
  y
end

function timeseries_init{T}(m::AbstractMarkovMap,n::Integer,xinit::T)
  @assert domain(m) == range(m)
  x_ts = Array(T,n)
  x_ts[1] = xinit
  x = copy(xinit)
  for i = 2:n
    x = m(x)
    x_ts[i] = x
  end
  x_ts
end
timeseries(m::AbstractMarkovMap,n::Integer;x0=map(m,rand(domain(m)),10^4)) = timeseries_init(m,n,x0)
timeseries{D}(m::AbstractMarkovMap,n::Integer,rho::Fun) = timeseries_init(m,n,sample(rho))
#convert(Int,floor(sqrt(n)/2))


function timehist_init{T}(m::AbstractMarkovMap,n::Integer,nbins,xinit::T)
  @assert domain(m) == range(m)
  bin_start = domain(m).a
  bin_step = (domain(m).b-domain(m).a)/nbins
  x_hist = Array(T,n)
  x = copy(xinit)
  for i = 2:n
    x = m(x)
    x_ts[min(nbins,max(1,convert(Int,ceil((x-bin_start)/bin_step))))] += 1
  end

  bin_start:bin_step:domain(m).b,x_ts
end
timehist(m::AbstractMarkovMap,n::Integer,n;x0=map(m,rand(domain(m)),10^4)) = timehist_init(m,n,x0)
timehist{D}(m::AbstractMarkovMap,n::Integer,rho::Fun) = timehist_init(m,n,sample(rho))
