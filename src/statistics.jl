# Overloads
export acim, linearresponse, correlationsum, birkhoffcov, birkhoffvar
export lyapunov, covariancefunction, MA_process, AR_process

"""
    acim(L)

Output a Fun object giving the acim of the associated map"""
acim(K::SolutionInvWrapper) = K.op\K.u
acim(L::Operator,u::Fun=uniform(domainspace(L))) = SolutionInv(L,u)\u
acim(M::AbstractMarkovMap) = acim(Transfer(M))

"""
    linearresponse(L, X::Fun)

Output a Fun object giving the first-order change in the acim of the map under perturbation X"""
function linearresponse(K::SolutionInvWrapper,X::Fun)
  @assert all(isapprox.(X.(∂(domain(X))),0;atol=sqrt(eps(cfstype(X))*arclength(domain(X)))))
  K\(-acim(K)*X)'
end

"""
    zero_to(A::Fun, ρ::Fun=uniform(space(A)))


Output ρ*(A - ∫A dρ)"""
function zero_to(A::Fun,ρ::Fun=uniform(space(A)))
  ρA = ρ*A
  ρA -= ρ*sum(ρA)/sum(ρ)
end

"""
    correlationsum(L, A::Fun, ρ=acim(L))

Output resolvent of transfer operator applied to ρ*(A - ∫A dρ)"""
correlationsum(K::SolutionInvWrapper,A::Fun,ρ=acim(K)) = K \ zero_to(A,ρ)

"""
    birkhoffcov(L, A::Fun, B::Fun)

Output covariance of CLT-normalised Birkhoff sums of A and B under the map"""
function birkhoffcov(K::SolutionInvWrapper,A::Fun,B::Fun)
  ρ = acim(K); Az = zero_to(A,ρ); Bz = zero_to(B,ρ)
  sum(B*(K\Az)) + sum(A*(K\Bz)) - sum(B*Az)
end

"""
    birkhoffvar(L, A::Fun)

Output diffusion coefficient of A under the map using Green-Kubo formula"""
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
  @compat cf = Array{cfstype(A)}(undef,n+1) # TODO: update to coefficienttype
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
  @compat cf = Array{cfstype(A)}(undef,n+1)
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
  @compat cfA = Array{cfstype(A)}(undef,n+1)  # TODO: update to coefficienttype
  @compat cfB = Array{cfstype(A)}(undef,n+1)
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
"""
    lyapunov(f, r=acim(f), sp=Space(rangedomain(f)))

Calculate Lyapunov exponent associated with f"""
function lyapunov(f::AbstractMarkovMap,r=acim(f),sp=Space(rangedomain(f)))
<<<<<<< HEAD
  sum(Fun(x->transferfunction(x,f,x->log(abs(f'(x)))*r(x),cfstype(f)),sp))
=======
  sum(Fun(x->transferfunction(x,f,x->log(abs(f'(x)))*r(x)),sp))
>>>>>>> 4c33d1554aa02b6a20399f0090f695de8ad3e517
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
<<<<<<< HEAD
  T = cfstype(f)
=======
>>>>>>> 4c33d1554aa02b6a20399f0090f695de8ad3e517
  tr = copy(r)
  lyap = sum(Fun(x->transferfunction(x,f.maps[end],LyapContainer(f.maps[end],tr)),rangedomain(f.maps[end])))
  for k = complength(f)-1:-1:1
    tr = Fun(x->transferfunction(x,f.maps[k+1],tr),rangedomain(f.maps[k+1]))
    lyap += sum(Fun(x->transferfunction(x,f.maps[k],LyapContainer(f.maps[k],tr)),rangedomain(f.maps[k])))
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

  logMA = Fun(Taylor(),coefficients(log(Corfn))/2)
end

"""
    MA_process(f,A::Fun)

Calculate coefficients of Gaussian MA process with the same autocorrelation function as A under the map f.

See Appendix A.2 of Wormell, C.L. & Gottwald, G.A. 'On the Validity of Linear Response Theory in High-Dimensional Deterministic Dynamical Systems' (2018).
"""
function MA_process(L,A::Fun)
  logMA = logMA_process(L,A)
  chop!(real.(coefficients(Fun(x->exp(logMA(x)),space(logMA)))))
end

"""
   AR_process(f,A::Fun)

Calculate coefficients of Gaussian AR process with the same autocorrelation function as A under the map f.

See Appendix A.2 of Wormell, C.L. & Gottwald, G.A. 'On the Validity of Linear Response Theory in High-Dimensional Deterministic Dynamical Systems' (2018).
"""
function AR_process(L,A::Fun)
  logMA = logMA_process(L,A)
  chop!(real.(coefficients(Fun(x->exp(-logMA(x)),space(logMA)))))
end
