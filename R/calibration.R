ac_projection <- function(F, intercept=FALSE) {
  F <- as.matrix(F); storage.mode(F) <- "double"
  if (any(!is.finite(F))) stop("Nonfinite factor design")
  T <- nrow(F)
  X <- if (intercept) cbind(1,F) else F
  if (!ncol(X)) return(diag(T))
  qx <- qr(X, tol=1e-10)
  if (qx$rank != ncol(X)) stop("Rank-deficient factor design")
  Q <- qr.Q(qx)[,seq_len(ncol(X)),drop=FALSE]
  M <- diag(T)-tcrossprod(Q)
  (M+t(M))/2
}

ac_nnls <- function(A,b,max_iter=NULL) {
  # Lawson-Hanson active-set NNLS. Tolerances concern numerical KKT
  # precision only; no positivity floor or statistical regularizer is used.
  A <- as.matrix(A); b <- as.numeric(b)
  if (nrow(A)!=length(b) || any(!is.finite(A)) || any(!is.finite(b)))
    stop("Invalid NNLS system")
  n <- ncol(A)
  if (is.null(max_iter)) max_iter <- 30L*max(1L,n)
  x <- numeric(n); passive <- rep(FALSE,n)
  w <- as.numeric(crossprod(A,b))
  tol <- 64*.Machine$double.eps*max(dim(A))*max(1,max(abs(w)))
  count <- 0L
  while (any(!passive & w>tol)) {
    candidates <- which(!passive)
    chosen <- candidates[which.max(w[candidates])]
    passive[chosen] <- TRUE
    repeat {
      count <- count+1L
      if (count>max_iter) stop("NNLS active set did not converge")
      z <- numeric(n)
      if (any(passive))
        z[passive] <- as.numeric(qr.solve(A[,passive,drop=FALSE],b,tol=1e-12))
      bad <- which(passive & z<=0)
      if (!length(bad)) break
      den <- x[bad]-z[bad]
      usable <- den>0
      if (!any(usable)) stop("Degenerate NNLS active-set step")
      step <- min(x[bad][usable]/den[usable])
      x <- x+step*(z-x)
      zero_tol <- 64*.Machine$double.eps*max(dim(A))*max(1,max(abs(x)))
      removed <- passive & x<=zero_tol
      x[removed] <- 0
      passive[removed] <- FALSE
    }
    x <- z
    w <- as.numeric(crossprod(A,b-as.numeric(A%*%x)))
  }
  x[x<0 & x>=-tol] <- 0
  if (any(x<0)) stop("NNLS returned a negative component")
  residual <- as.numeric(A%*%x)-b
  list(x=x,residual_norm=sqrt(sum(residual^2)),iterations=count,
       kkt_max_inactive=max(c(0,w[!passive])),
       kkt_max_active=max(c(0,abs(w[passive]))),kkt_tolerance=tol)
}

ac_time_variance <- function(E,F,d) {
  # E must be the unrestricted OLS residual matrix M_[1,F] Y.
  # With fixed coordinate scales, the common-radial working model gives
  # E{rowMeans((E/d)^2)|F}=(M_X elementwise-squared) lambda, up to a common
  # cross-sectional scale. Fitted d introduces a plug-in approximation.
  E <- as.matrix(E); d <- as.numeric(d)
  T <- nrow(E)
  if (length(d)!=ncol(E) || any(!is.finite(E)) ||
      any(!is.finite(d)) || any(d<=0)) stop("Invalid residuals/scales")
  MX <- ac_projection(F,intercept=TRUE)
  if (nrow(MX)!=T) stop("Residual/factor dimension mismatch")
  s <- rowMeans(sweep(E,2,d,"/")^2)
  A <- MX^2
  unconstrained <- tryCatch(as.numeric(solve(A,s)),error=function(e) NULL)
  tolerance <- if (is.null(unconstrained)) NA_real_ else
    64*.Machine$double.eps*T*max(1,max(abs(unconstrained)))
  if (!is.null(unconstrained) && min(unconstrained)>=-tolerance) {
    lambda <- pmax(unconstrained,0)
    residual_norm <- sqrt(sum((as.numeric(A%*%lambda)-s)^2))
    solver <- "nonnegative_direct_solution"
    iterations <- 0L
    kkt <- as.numeric(crossprod(A,as.numeric(A%*%lambda)-s))
    kkt_active <- max(c(0,abs(kkt[lambda>0])))
    kkt_inactive <- max(c(0,-kkt[lambda==0]))
  } else {
    ans <- ac_nnls(A,s)
    lambda <- ans$x; residual_norm <- ans$residual_norm
    solver <- "NNLS"; iterations <- ans$iterations
    kkt_active <- ans$kkt_max_active
    kkt_inactive <- ans$kkt_max_inactive
  }
  common_scale <- mean(lambda)
  if (!is.finite(common_scale) || common_scale<=0)
    stop("Degenerate estimated time covariance")
  list(lambda=lambda,residual_norm=residual_norm,
       zeros=sum(lambda==0),zero_count=sum(lambda==0),
       lambda_removed_scale=common_scale,solver=solver,iterations=iterations,
       unconstrained_lambda_min=if (is.null(unconstrained)) NA_real_ else min(unconstrained),
       maximum_design_leverage=max(1-diag(MX)),
       kkt_max_active=kkt_active,kkt_max_inactive=kkt_inactive,
       row_square_mean=s)
}

ac_angular_gamma <- function(E,F,d,lambda) {
  # Gaussian-whitened residual trace-ratio estimator. The Wishart moment
  # subtraction removes its leading covariance-estimation noise under the
  # working Gaussian model. Estimated scales/row variances are not fixed
  # nuisances in the actual data, so exact Wishart unbiasedness is not claimed.
  E <- as.matrix(E); d <- as.numeric(d); lambda <- as.numeric(lambda)
  T <- nrow(E); N <- ncol(E)
  if (length(lambda)!=T || any(!is.finite(lambda)) || any(lambda<0) ||
      mean(lambda)<=0 || length(d)!=N || any(d<=0))
    stop("Invalid angular trace inputs")
  MX <- ac_projection(F,intercept=TRUE)
  Ct <- sweep(MX,2,lambda/mean(lambda),"*")%*%MX
  Ct <- (Ct+t(Ct))/2
  ee <- eigen(Ct,symmetric=TRUE)
  tolerance <- 64*.Machine$double.eps*T*max(1,max(abs(ee$values)))
  keep <- ee$values>tolerance
  rank <- sum(keep)
  if (rank<=2L) stop("Insufficient positive time-covariance rank")
  E0 <- sweep(E,2,d,"/")
  W <- sweep(crossprod(ee$vectors[,keep,drop=FALSE],E0),
             1,sqrt(ee$values[keep]),"/")
  Gram <- tcrossprod(W)
  g1 <- sum(diag(Gram)); g2 <- sum(Gram^2)
  if (!is.finite(g1) || g1<=0) stop("Degenerate whitened residual covariance")
  unbounded <- (g2-g1^2/rank)/((1+1/rank-2/rank^2)*g1^2)
  # trace(R^2)/N^2 >= 1/N for a correlation matrix. This is a parameter-space
  # constraint, not a fitted size correction.
  list(gamma=max(1/N,unbounded),gamma_unbounded=unbounded,
       rank=rank,whitening_rank=rank,g1=g1,g2=g2,
       eigen_tolerance=tolerance)
}

ac_ss_fit <- function(Z,tol=1e-7,max_iter=500L) {
  # Current converged spatial median + diagonal-sign fit; no legacy formula.
  Z <- as.matrix(Z); T <- nrow(Z); N <- ncol(Z)
  if (T<3L || N<2L || any(!is.finite(Z))) stop("Invalid spatial sample")
  theta <- colMeans(Z); d <- apply(Z,2,sd)
  if (any(d<=0)) stop("Nonpositive initial spatial scale")
  d <- d/exp(mean(log(d))); converged <- FALSE
  for (iteration in seq_len(max_iter)) {
    V <- sweep(sweep(Z,2,theta,"-"),2,d,"/")
    r <- sqrt(rowSums(V^2))
    if (min(r)<=100*.Machine$double.eps*median(r))
      stop("Collapsed spatial median radius")
    U <- V/r; ir <- 1/r
    theta_new <- colSums(Z*(ir/sum(ir)))
    d_new <- d*sqrt(N*colMeans(U^2))
    d_new <- d_new/exp(mean(log(d_new)))
    change <- max(max(abs((theta_new-theta)/d))/(median(r)/sqrt(N)),
                  max(abs(log(d_new/d))))
    theta <- theta_new; d <- d_new
    if (change<tol) {
      V <- sweep(sweep(Z,2,theta,"-"),2,d,"/")
      r <- sqrt(rowSums(V^2)); U <- V/r
      eq1 <- sqrt(N)*max(abs(colMeans(U)))
      eq2 <- max(abs(N*colMeans(U^2)-1))
      if (max(eq1,eq2)<tol) {converged <- TRUE;break}
    }
  }
  if (!converged) stop("Spatial median did not converge")
  list(theta=theta,d=d,r=r,iterations=iteration,
       equation=max(eq1,eq2),location_equation=eq1,scale_equation=eq2)
}

ac_ss_components <- function(F,r) {
  r <- as.numeric(r); F <- as.matrix(F); T <- length(r)
  if (nrow(F)!=T || any(!is.finite(r)) || any(r<=0))
    stop("Invalid frozen SS radii")
  M <- ac_projection(F)
  h <- as.numeric(M%*%rep(1,T)); H <- sum(h^2); eta <- 1-H/T
  if (H<=0) stop("Unidentified intercept")
  a <- 1/r; a <- a/mean(a)
  v <- as.numeric(M%*%a)
  C <- M-tcrossprod(rep(1,T),v)/sum(a)
  K <- crossprod(C,a^2*C); K <- (K+t(K))/2
  m_1 <- mean(r); m_minus1 <- mean(1/r); m_2 <- mean(r^2)
  B <- 1-2*eta*m_minus1*m_1+eta*m_2*m_minus1^2
  if (!is.finite(B) || B<=0) stop("Nonpositive formal SS radial denominator")
  list(M=M,h=h,H=H,eta=eta,nu=T-ncol(F)-1L,a=a,v=v,C=C,K=K,
       B=B,m_minus1=m_minus1,m_1=m_1,m_2=m_2,
       radial_J_coefficient=(1-eta*m_minus1*m_1)/B)
}

ac_ss_quadratic_tail <- function(x,components,lambda,
                                 rel_tol=1e-10,abs_tol=1e-15) {
  # Exact tail in the frozen Gaussian quadratic-ratio model. At most one
  # eigenvalue is positive. The integral is nonoscillatory and deterministic.
  x <- as.numeric(x); lambda <- as.numeric(lambda)
  T <- length(components$v)
  if (length(x)!=1L || is.na(x) || x<0 ||
      length(lambda)!=T || any(!is.finite(lambda)) ||
      any(lambda<0) || mean(lambda)<=0) stop("Invalid quadratic-tail inputs")
  if (x==0) return(list(p=1,log_p=0,positive_eigenvalue=NA_real_,
                        negative_count=NA_integer_,integral_error=0))
  if (is.infinite(x)) return(list(p=0,log_p=-Inf,positive_eigenvalue=0,
                                 negative_count=NA_integer_,integral_error=0))
  sq <- sqrt(lambda/mean(lambda))
  A <- (tcrossprod(components$v)-x*components$B*components$K)*tcrossprod(sq)
  ev <- eigen((A+t(A))/2,symmetric=TRUE,only.values=TRUE)$values
  tolerance <- 64*.Machine$double.eps*T*max(1,max(abs(ev)))
  positive <- ev[ev>tolerance]; negative <- -ev[ev< -tolerance]
  if (length(positive)>1L) stop("More than one positive quadratic eigenvalue")
  if (!length(positive)) return(list(p=0,log_p=-Inf,positive_eigenvalue=0,
                                    negative_count=length(negative),
                                    integral_error=0,eigen_tolerance=tolerance))
  if (!length(negative)) return(list(p=1,log_p=0,positive_eigenvalue=positive,
                                    negative_count=0L,integral_error=0,
                                    eigen_tolerance=tolerance))
  ratio <- negative/positive
  integrand <- function(theta) vapply(theta,function(z) {
    s2 <- sin(z)^2
    if (s2==0) return(0)
    exp(-.5*sum(log1p(ratio/s2)))
  },numeric(1))
  val <- integrate(integrand,0,pi/2,subdivisions=1000L,
                   rel.tol=rel_tol,abs.tol=abs_tol,stop.on.error=FALSE)
  p <- (2/pi)*val$value; error <- (2/pi)*val$abs.error
  if (!is.finite(p) || p<0 || p>1+1e-12 ||
      (val$message!="OK" && error>max(abs_tol,1e-8*p)))
    stop(paste("Quadratic tail integration failed:",val$message))
  p <- min(1,p)
  list(p=p,log_p=log(p),positive_eigenvalue=as.numeric(positive),
       negative_count=length(negative),integral_error=error,
       integral_message=val$message,eigen_tolerance=tolerance)
}

ac_max_from_marginal <- function(p,N,calibration="sidak") {
  calibration <- tolower(calibration)
  if (length(p)!=1L || !is.finite(p) || p<0 || p>1 || N<2)
    stop("Invalid marginal tail or dimension")
  if (p==0) return(0)
  if (p==1) return(1)
  if (calibration=="sidak") return(-expm1(N*log1p(-p)))
  if (calibration!="gumbel") stop("Unknown maximum calibration")
  z <- qnorm(p/2,lower.tail=FALSE)
  centered <- z^2-2*log(N)+log(log(N))
  log_rate <- -.5*(centered+log(pi))
  if (log_rate>log(.Machine$double.xmax)) return(1)
  -expm1(-exp(log_rate))
}

ac_ss_analytic <- function(Y,F,calibration="sidak",
                           fit_tol=1e-7,max_iter=500L) {
  Y <- as.matrix(Y); F <- as.matrix(F); T <- nrow(Y); N <- ncol(Y)
  if (nrow(F)!=T || any(!is.finite(Y))) stop("Invalid panel")
  M <- ac_projection(F); MX <- ac_projection(F,intercept=TRUE)
  Z <- M%*%Y; E <- MX%*%Y
  fit <- ac_ss_fit(Z,fit_tol,max_iter)
  co <- ac_ss_components(F,fit$r)
  zeta <- N*co$m_minus1^2/co$B
  coordinate_raw <- T*zeta*(fit$theta/fit$d)^2
  raw <- max(coordinate_raw)
  centered <- raw-2*log(N)+log(log(N))
  log_rate <- -.5*(centered+log(pi))
  p_formal <- if (log_rate>log(.Machine$double.xmax)) 1 else
    -expm1(-exp(log_rate))
  radial <- ac_time_variance(E,F,fit$d)
  angular <- ac_angular_gamma(E,F,fit$d,radial$lambda)
  cJ <- 1-angular$gamma*co$radial_J_coefficient
  if (!is.finite(cJ) || cJ<=0)
    stop("Nonpositive first-order SS angular-Jacobian factor")
  # The original raw statistic is preserved and returned unchanged.
  # cJ modifies only the working reference-tail mapping.
  tail <- ac_ss_quadratic_tail(cJ^2*raw,co,radial$lambda)
  p_sidak <- ac_max_from_marginal(tail$p,N,"sidak")
  p_gumbel <- ac_max_from_marginal(tail$p,N,"gumbel")
  selected <- ac_max_from_marginal(tail$p,N,calibration)
  list(p=selected,original_p=p_formal,p_formal=p_formal,p_calibrated=selected,
       p_calibrated_sidak=p_sidak,p_calibrated_gumbel=p_gumbel,
       marginal_p=tail$p,marginal_tail=tail$p,raw=raw,raw_formal=raw,
       centered_formal=centered,gamma=angular$gamma,lambda=radial$lambda,
       coordinate_raw=coordinate_raw,zeta_formal=zeta,cJ=cJ,
       reference_raw=cJ^2*raw,fit=fit,radial=radial,angular=angular,
       components=co,tail_diagnostics=tail,
       calibration=tolower(calibration),
       model_scope=paste("Frozen Gaussian/common-time-variance working model;",
                         "fitted nuisance parameters and IC-t3 are not exact conditioning."))
}

ac_pqr <- function(E) {
  E <- as.matrix(E); T <- nrow(E); N <- ncol(E)
  pairs <- which(upper.tri(matrix(FALSE,T,T)),arr.ind=TRUE)
  k <- ceiling(nrow(pairs)/2)
  q <- vapply(seq_len(N),function(j) {
    a <- abs(E[pairs[,1],j]-E[pairs[,2],j])
    sort(a,partial=k)[k]
  },numeric(1))
  if (any(!is.finite(q)) || any(q<=0)) stop('Nonpositive PQR scale')
  q/sqrt(mean(q^2))
}

ac_in_fit <- function(Z,F,d,tol=1e-8,max_iter=1000L) {
  Z <- as.matrix(Z); F <- as.matrix(F); T <- nrow(Z); N <- ncol(Z)
  M <- ac_projection(F); h <- as.numeric(M%*%rep(1,T)); H <- sum(h^2)
  nu <- T-ncol(F)-1L
  if (H<=1 || nu<=4) stop('Insufficient intercept information or residual df')
  alpha <- as.numeric(crossprod(h,Z))/H; converged <- FALSE
  for (iteration in seq_len(max_iter)) {
    E <- sweep(Z-tcrossprod(h,alpha),2,d,'/')
    r <- sqrt(rowSums(E^2))
    if (min(r)<=100*.Machine$double.eps*median(r)) stop('Collapsed IN radius')
    ir2 <- 1/r^2
    updated <- colSums(Z*(h*ir2))/sum(h^2*ir2)
    change <- max(abs((updated-alpha)/d))/(median(r)/sqrt(N)); alpha <- updated
    if (change<tol) {
      E <- sweep(Z-tcrossprod(h,alpha),2,d,'/')
      r <- sqrt(rowSums(E^2)); ir2 <- 1/r^2
      eq <- sqrt(N)*median(r)*max(abs(colMeans(E*(h*ir2))))
      if (eq<tol) {converged <- TRUE; break}
    }
  }
  if (!converged) stop('IN root did not converge')
  A <- mean(h^2*ir2); a <- h*ir2
  projection <- as.numeric((diag(T)-M)%*%a)
  Bc <- A+mean(r^2*projection^2)
  z <- sqrt(nu*N)*A*alpha/(d*sqrt(Bc))
  list(alpha=alpha,d=d,r=r,A=A,Bc=Bc,z=z,raw=max(z^2),
       h=h,H=H,nu=nu,iterations=iteration,equation=eq)
}

ac_pqr_heterogeneous_fast <- function(lambda,tol=2e-9) {
  # Exact aggregation of the same independent-row Hajek variances. One
  # centered integral in log physical coordinates resolves all positive
  # variance scales without altering a zero or small positive lambda.
  lambda <- as.numeric(lambda); T <- length(lambda)
  if (T<3 || any(!is.finite(lambda)) || any(lambda<0) || !any(lambda>0))
    stop('Invalid nonnegative time-variance profile')
  global_scale <- mean(lambda); lambda <- lambda/global_scale
  pairs <- which(upper.tri(matrix(FALSE,T,T)),arr.ind=TRUE)
  pair_sd <- sqrt(lambda[pairs[,1]]+lambda[pairs[,2]])
  mass <- mean(pair_sd==0); positive_sd <- pair_sd[pair_sd>0]
  if (mass>=.5) stop('Degenerate zero population PQR median')
  cdf <- function(q) mass+sum(2*pnorm(q/positive_sd)-1)/length(pair_sd)
  q <- uniroot(function(q)cdf(q)-.5,c(0,max(positive_sd)*qnorm(.75)*1.01),tol=1e-13)$root
  u <- sum(2*dnorm(q/positive_sd)/positive_sd)/length(pair_sd)
  pair_matrix <- sqrt(outer(lambda,lambda,'+'))
  probability <- matrix(1,T,T)
  positive <- pair_matrix>0
  probability[positive] <- 2*pnorm(q/pair_matrix[positive])-1
  diag(probability) <- 0
  pj <- rowSums(probability)/(T-1)
  keep <- lambda>0; sd <- sqrt(lambda[keep]); nzero <- sum(!keep)
  coefficient <- 16/(T*(q*u)^2)
  integrand <- function(log_x,below_q) vapply(log_x,function(yy) {
    xx <- exp(yy)
    a <- pnorm((xx+q)/sd)-pnorm((xx-q)/sd)
    total <- sum(a)+if (below_q) nzero else 0
    centered <- (total-a)/(T-1)-pj[keep]
    z <- xx/sd
    density_jacobian <- numeric(length(z)); active <- z<40
    density_jacobian[active] <- dnorm(z[active])*z[active]
    2*sum(density_jacobian*centered^2)
  },numeric(1))
  lower <- log(min(sd))-48; upper <- log(12*max(sd)); log_q <- log(q)
  absolute_tolerance <- min(1e-10,tol/coefficient/10)
  left <- integrate(function(y)integrand(y,TRUE),lower,log_q,
                    abs.tol=absolute_tolerance,rel.tol=2e-11,
                    subdivisions=300L,stop.on.error=TRUE)
  right <- integrate(function(y)integrand(y,FALSE),log_q,upper,
                     abs.tol=absolute_tolerance,rel.tol=2e-11,
                     subdivisions=300L,stop.on.error=TRUE)
  variance_sum <- left$value+right$value
  kappa <- coefficient*variance_sum
  tail_bound <- 2*dnorm(0)*exp(-48)*sum(min(sd)/sd)+
                2*sum(keep)*pnorm(12,lower.tail=FALSE)
  list(kappa=kappa,q=q*sqrt(global_scale),u=u/sqrt(global_scale),
       variance_sum=variance_sum,mean_pair_probability=mean(pj),
       integration_error=coefficient*(left$abs.error+right$abs.error+tail_bound),
       positive_profile_rank=sum(keep),zero_count=nzero,
       integration_method='aggregated centered log-coordinate quadrature')
}

ac_pqr_t_match <- function(nu,kappa) {
  if (!is.finite(kappa) || kappa<2-1e-6 || nu<=4)
    stop('PQR moment approximation is outside its stated range')
  rho <- 1+(kappa-2)/nu
  m <- 4+2/(rho*(nu-2)/(nu-4)-1)
  s2 <- (nu/(nu-2))/(m/(m-2))
  list(df=m,scale_squared=s2,rho=rho,kappa=kappa)
}

ac_in_analytic <- function(Y,F,calibration='sidak',pqr_profile='heterogeneous') {
  Y <- as.matrix(Y); F <- as.matrix(F); T <- nrow(Y); N <- ncol(Y)
  M <- ac_projection(F); MX <- ac_projection(F,intercept=TRUE)
  Z <- M%*%Y; E <- MX%*%Y; d <- ac_pqr(E)
  fit <- ac_in_fit(Z,F,d)
  radial <- ac_time_variance(E,F,d)
  lambda <- radial$lambda
  angular <- ac_angular_gamma(E,F,d,lambda)
  v <- as.numeric(M%*%(fit$h/fit$r^2))
  kappa_variance <- (fit$nu*N/T)*sum(lambda*v^2)/sum(fit$r^2*v^2)
  if (!is.finite(kappa_variance) || kappa_variance<=0) stop('Invalid working variance')
  pqr <- switch(pqr_profile,
    heterogeneous=ac_pqr_heterogeneous_fast(lambda),
    gaussian=list(kappa=2.316186819891847),
    stop('pqr_profile must be heterogeneous or gaussian'))
  match <- ac_pqr_t_match(fit$nu,pqr$kappa)
  cJ <- 1-2*angular$gamma
  if (!is.finite(cJ) || cJ<=0) stop('Invalid leading angular correction')
  argument <- cJ*sqrt(fit$raw/(kappa_variance*match$scale_squared))
  marginal_p <- 2*pt(argument,df=match$df,lower.tail=FALSE)
  p <- ac_max_from_marginal(marginal_p,N,calibration)
  original_p <- ac_max_from_marginal(2*pnorm(sqrt(fit$raw),lower.tail=FALSE),N,'gumbel')
  c(fit,list(p=p,original_p=original_p,marginal_p=marginal_p,
             kappa_variance=kappa_variance,gamma=angular$gamma,cJ=cJ,
             pqr_kappa=pqr$kappa,effective_df=match$df,t_scale_squared=match$scale_squared,
             lambda=lambda,radial=radial,angular=angular,pqr=pqr,
             calibration=calibration,pqr_profile=pqr_profile))
}
