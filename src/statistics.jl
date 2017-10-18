# Overloads
export acim, linearresponse, correlationsum, birkhoffcov, birkhoffvar
export lyapunov, covariancefunction, MA_process, AR_process


acim(K::SolutionInvWrapper) = K.op\K.u
acim(L::Operator,u::Fun=uniform(domainspace(L))) = SolutionInv(L,u)\u
acim(M::AbstractMarkovMap) = acim(Transfer(M))

function zero_to(f::Fun,r::Fun=uniform(space(f)))
  rf = r*f
  rf -= r*sum(rf)/sum(r)
end

function linearresponse(K::SolutionInvWrapper,X::Fun)
  @assert all(isapprox.(X.(∂(domain(X))),0;atol=sqrt(eps(eltype(X))*arclength(domain(X)))))
  K\(-acim(K)*X)'
end
correlationsum(K::SolutionInvWrapper,A::Fun,r=acim(K)) = K \ zero_to(A,r)

function birkhoffcov(K::SolutionInvWrapper,A::Fun,B::Fun)
  r = acim(K); Az = zero_to(A,r); Bz = zero_to(B,r)
  sum(B*(K\Az)) + sum(A*(K\Bz)) - sum(B*Az)
end

function birkhoffvar(K::SolutionInvWrapper,A::Fun)
  Az = zero_to(A,acim(K))
  2sum(A*(K\Az)) - sum(A*Az)
end

for OP in (:linearresponse,:correlationsum,:birkhoffvar)
  @eval $OP(L::Operator,X::Fun) = $OP(SolutionInv(L),X)
  @eval $OP(M::AbstractMarkovMap,X::Fun) = $OP(Transfer(M),X)
  @eval $OP(obj::Operator,X) = $OP(obj,Fun(X,domainspace(obj)))
  @eval $OP(obj::AbstractMarkovMap,X) = $OP(obj,Fun(X,domain(obj)))
end
birkhoffcov(L::Operator,X::Fun,Y::Fun) = birkhoffcov(SolutionInv(L),X,Y)
birkhoffcov(M::AbstractMarkovMap,X::Fun,Y::Fun) = birkhoffcov(Transfer(M),X,Y)

# lag-correlation function
for TYP=(:AbstractMarkovMap,:SolutionInvWrapper)
  @eval covariancefunction(K::$TYP,A::Fun,n::Int) = covariancefunction(Transfer(K),A,n)
  @eval covariancefunction(K::$TYP,A::Fun,B::Fun,n::Int) = covariancefunction(Transfer(K),A,B,n)
  @eval covariancefunction(K::$TYP,A::Fun;tol=eps(norm(A)^2)) = covariancefunction(Transfer(K),A;tol=tol)
  @eval covariancefunction(K::$TYP,A::Fun,B::Fun;tol=eps(norm(A)*norm(B))) = covariancefunction(Transfer(K),A,B;tol=tol)
end
  """
      cfA = covariancefunction(L, A, n)
  Compute the lag-covariance function against transfer operator L for observable A for n time steps
  in each direction. Specifically, the (k+1)th entry is the expectation of A A∘f^(k).
  """
function covariancefunction(L::Operator,A::Fun,n::Int;r=acim(L))
  cf = Array{eltype(A)}(n+1) # TODO: update to coefficienttype
  Az = zero_to(A,r)
  cf[1] = sum(A*Az)
  for i = 1:n
    Az = L*Az
    Az -= r*sum(Az)
    chop!(Az)
    cf[i+1] = sum(A*Az)
  end
  cf
end

"""
    cf = covariancefunction(L, A; tol = eps(A))
Choose n so that the covariance declines to a given tolerance
"""
function covariancefunction(L::Operator,A::Fun;r=acim(L),tol=eps(norm(A.coefficients)^2))
  n = 16
  cf = Array{eltype(A)}(n+1)
  Az = zero_to(A,r)
  cf[1] = sum(A*Az)
  for i = 1:n
    Az = L*Az
    Az -= r*sum(Az)
    chop!(Az)
    cf[i+1] = sum(A*Az)
  end
  while maximum(abs.(cf[div(2n,3):end]))>tol && n <= 2^20
    n_old = n
    n = 2n
    pad!(cf,n+1)
    for i = n_old+1:n
      Az = L*Az
      Az -= r*sum(Az)
      chop!(Az)
      cf[i+1] = sum(A*Az)
    end
    # chop!(Az)
  end

  if n==2^20
    warn("Maximum number of coefficients reached (n=$n)")
    cf
  else
    chop!(cf)
  end
end

"""
    cfA, cfB = covariancefunction(L, A, B, n)
Compute the lag-covariance function against transfer operator L between observables A and B for n time steps
in each direction. Specifically, the (k+1)th entry of cfA is the expectation of B A∘f^(k) and vice versa for cfB.
"""
function covariancefunction(L::Operator,A::Fun,B::Fun,n::Int;r=acim(L))
  cfA = Array{eltype(A)}(n+1)  # TODO: update to coefficienttype
  cfB = Array{eltype(A)}(n+1)
  Az = zero_to(A,r); Bz = zero_to(B,r)
  cfA[1] = sum(A*Bz)
  cfB[1] = cfA[1]
  for i = 1:n
    Az = L*Az;Bz = L*Bz
    Az -= r*sum(Az); Bz -= r*sum(Bz)
    chop!(Az); chop!(Bz)
    cfA[i+1] = sum(A*Bz); cfB[i+1] = sum(B*Az)
  end
  cfA,cfB
end

"""
    cfA, cfB = covariancefunction(L, A, B; tol=eps(A*B))
Chooses n so that the covariance declines to a given tolerance.
"""
function covariancefunction(L::Operator,A::Fun,B::Fun;r=acim(L),tol=eps(norm(coefficients(A))*norm(coefficients(B))))
  n = 16
  cfA,cfB = covariancefunction(L,A,B,n;r=r)
  while max(maximum(abs.(cfA[div(2n,3):end])),maximum(abs.(cfB[div(2n,3):end])))>tol && n <= 2^20
    n = 2n
    cfA,cfB = covariancefunction(L,A,B,n;r=r)   # chop!(Az)
  end

  if n==2^20
    warn("Maximum number of coefficients reached (n=$n)")
    cfA,cfB
  else
    chop!(cfA),chop!(cfB)
  end
end

# Lyapunov exponent
function lyapunov(f::AbstractMarkovMap,r=acim(f),sp=Space(rangedomain(f)))
  sum(Fun(x->transferfunction(x,f,x->log(abs(f'(x)))*r(x),eltype(f)),sp))
end
function lyapunov(K::SolutionInvWrapper)
  L = Transfer(K)
  lyapunov(markovmap(L),acim(K),rangespace(L))
end
lyapunov(m::Operator) = lyapunov(markovmap(m),acim(m),rangespace(m))

struct LyapContainer{M<:AbstractMarkovMap,rr}
  m::M
  tr::rr
end
(lc::LyapContainer)(x) = log(abs(lc.m'(x)))*lc.tr(x)

function lyapunov(f::ComposedMarkovMap,r=acim(f))#,sp=Space(rangedomain(f)))
  T = eltype(f)
  tr = copy(r)
  lyap = sum(Fun(x->transferfunction(x,f.maps[end],LyapContainer(f.maps[end],tr),T),rangedomain(f.maps[end])))
  for k = complength(f)-1:-1:1
    tr = Fun(x->transferfunction(x,f.maps[k+1],tr,T),rangedomain(f.maps[k+1]))
    lyap += sum(Fun(x->transferfunction(x,f.maps[k],LyapContainer(f.maps[k],tr),T),rangedomain(f.maps[k])))
  end
  lyap
end

# Gaussian process parametrisation
function logMA_process(L,A::Fun)
  @assert isreal(A)
  corf = covariancefunction(L,A)
  # d = Circle(); sp = Taylor(d)
  Corfn = Fun(CosSpace(),[corf[1];2corf[2:end]])
  # lCf = log(Corfn)/2
  # SlogMA = Fun(SinSpace(),coefficients(lCf)[2:end])
  # MA = cos(SlogMA)*exp(lCf)
  #  MA = 1+expm1(Fun(sp,chop!(coefficients(log(Corfn)/2))))
  #  # TODO: because the special functions in ApproxFun are buggy
  # chop!(real.(coefficients(MA)[1:2:end]))

  ## OR
  # Corfn = Fun(Laurent(d),[corf[1];repeat(corf[2:end],inner=2)])
  # logMA = Fun(sp,
  # coefficients(Fun(x->log(abs(Corfn(e^(im*x)))),CosSpace())))/2

  logMA = Fun(Taylor(),coefficients(log(Corfn)))/2
end

function MA_process(L,A::Fun)
  logMA = logMA_process(L,A)
  chop!(real.(coefficients(Fun(x->exp(logMA(x)),space(logMA)))))
end

function AR_process(L,A::Fun)
  logMA = logMA_process(L,A)
  chop!(real.(coefficients(Fun(x->exp(-logMA(x)),space(logMA)))))
end
