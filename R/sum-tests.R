INMAX_METHODS <- c("LS-SUM","LS-MAX","LS-CC","SS-SUM","SS-MAX","SS-CC",
                   "IN-SUM","IN-MAX","IN-CC")

inmax_projection <- function(F) {
  F <- as.matrix(F); storage.mode(F) <- "double"
  T <- nrow(F); p <- ncol(F)
  qrf <- qr(F)
  if (qrf$rank != p) stop("Rank-deficient factor design")
  Q <- qr.Q(qrf)
  h <- as.numeric(1 - Q %*% crossprod(Q, rep(1,T)))
  h2 <- sum(h^2); df <- T-p-1L
  if (h2 <= 1 || df <= 4) stop("Insufficient intercept information or residual df")
  list(F=F,Q=Q,h=h,h2=h2,eta=1-h2/T,df=df)
}

inmax_project <- function(Y, pr) Y-pr$Q %*% crossprod(pr$Q,Y)

inmax_gumbel <- function(maximum, N) {
  M <- maximum-2*log(N)+log(log(N))
  log_rate <- -(M+log(pi))/2
  p <- if (log_rate > 6.6) 1 else -expm1(-exp(log_rate))
  c(p=p, statistic=M)
}

inmax_cauchy <- function(p, truncated=FALSE) {
  if (any(!is.finite(p)) || any(p<0 | p>1)) return(c(p=NA_real_,statistic=NA_real_))
  if (any(p==0)) return(c(p=0,statistic=Inf))
  if (!truncated && any(p==1)) return(c(p=1,statistic=-Inf))
  active <- if (truncated) p<.5 else rep(TRUE,length(p))
  terms <- numeric(length(p)); tiny <- active & p<1e-8
  terms[tiny] <- 1/(pi*p[tiny])
  regular <- active & !tiny
  terms[regular] <- tan(pi*(.5-p[regular]))
  C <- mean(terms)
  tail <- if (C>1) atan(1/C)/pi else if (C < -1) 1+atan(1/C)/pi else .5-atan(C)/pi
  c(p=tail,statistic=C)
}

inmax_ols <- function(Z,pr) {
  T <- nrow(Z); N <- ncol(Z)
  alpha <- as.numeric(crossprod(pr$h,Z)/pr$h2)
  E <- Z-tcrossprod(pr$h,alpha)
  sse <- colSums(E^2)
  if (any(sse<=0)) stop("Nonpositive OLS variance")
  sd <- sqrt(sse/pr$df)
  t <- sqrt(pr$h2)*alpha/sd
  C <- crossprod(sweep(E,2,sqrt(sse),"/"))
  rho <- C[upper.tri(C)]
  retained <- abs(sqrt(pr$df)*rho)>qnorm(1-.1/(2*N))
  a2 <- 1+2*sum(rho[retained]^2)/N
  center <- pr$df/(pr$df-2)
  scale <- center*sqrt(2*(pr$df-1)/(pr$df-4)*a2)
  S <- sum(t^2-center)/sqrt(N)/scale
  pm <- inmax_gumbel(max(t^2),N)
  ps <- pnorm(S,lower.tail=FALSE)
  pc <- inmax_cauchy(c(ps,pm["p"]))
  list(alpha=alpha,residual=E,sd=sd,t=t,S=S,M=unname(pm["statistic"]),
       C=unname(pc["statistic"]),p=c(ps,unname(pm["p"]),unname(pc["p"])))
}

inmax_diag_sign <- function(Z,tol=1e-7,max_iter=500L) {
  T <- nrow(Z); N <- ncol(Z); sq <- Z^2
  D <- colMeans(sq)
  if (any(D<=0)) stop("Degenerate initial sign scale")
  D <- D/mean(D); converged <- FALSE
  for (iteration in seq_len(max_iter)) {
    radius2 <- rowSums(sweep(sq,2,D,"/"))
    if (any(!is.finite(radius2)) || any(radius2<=0)) stop("Zero restricted sign radius")
    updated <- N*colMeans(sq/radius2)
    moments <- updated/D
    updated <- updated/mean(updated)
    eq <- max(abs(moments-1)); change <- max(abs(log(updated/D)))
    D <- updated
    if (max(eq,change)<=tol) {converged <- TRUE;break}
  }
  if (!converged) stop("Restricted sign scale did not converge")
  E <- sweep(Z,2,sqrt(D),"/"); r <- sqrt(rowSums(E^2)); U <- E/r
  eq <- max(abs(N*colMeans(U^2)-1))
  if (eq>tol) stop("Restricted sign equation residual exceeds tolerance")
  list(D=D,U=U,r=r,iterations=iteration,equation=eq)
}

inmax_sign_q <- function(U,h) {
  weighted <- as.numeric(crossprod(h,U))
  ncol(U)/sum(h^2)*(sum(weighted^2)-sum(h^2*rowSums(U^2)))
}

inmax_split_trace_literal <- function(Y,pr,D) {
  T <- nrow(Y); N <- ncol(Y); total <- 0
  for (s in seq_len(T-1L)) for (t in seq.int(s+1L,T)) {
    rest <- setdiff(seq_len(T),c(s,t)); k <- floor(length(rest)/2)
    first <- rest[seq_len(k)]; second <- rest[-seq_len(k)]
    F1 <- pr$F[first,,drop=FALSE]; F2 <- pr$F[second,,drop=FALSE]
    B1 <- solve(crossprod(F1),crossprod(F1,Y[first,,drop=FALSE]))
    B2 <- solve(crossprod(F2),crossprod(F2,Y[second,,drop=FALSE]))
    u <- function(i,B) {e <- as.numeric(Y[i,]-pr$F[i,]%*%B)/sqrt(D);e/sqrt(sum(e^2))}
    us1 <- u(s,B1);ut2 <- u(t,B2);ut1 <- u(t,B1);us2 <- u(s,B2)
    total <- total+pr$h[s]^2*pr$h[t]^2*(sum(us1*ut2)^2+sum(ut1*us2)^2)
  }
  N^2/(pr$h2*(pr$h2-1))*total
}

inmax_split_trace <- function(Y,pr,D) {
  T <- nrow(Y); N <- ncol(Y)
  if (T%%2L) return(inmax_split_trace_literal(Y,pr,D))
  pairs <- which(upper.tri(matrix(FALSE,T,T)),arr.ind=TRUE)
  s_all <- pairs[,1];t_all <- pairs[,2];H <- T%/%2L
  specs <- list(
    list(keep=t_all<=H,rows1=seq_len(H+1L),rows2=seq.int(H+2L,T),del1=c(1L,2L),del2=integer()),
    list(keep=s_all<=H & t_all>H,rows1=seq_len(H),rows2=seq.int(H+1L,T),del1=1L,del2=2L),
    list(keep=s_all>H,rows1=seq_len(H-1L),rows2=seq.int(H,T),del1=integer(),del2=c(1L,2L)))
  total <- 0
  for (spec in specs) {
    s <- s_all[spec$keep];t <- t_all[spec$keep]
    if (!length(s)) next
    ff <- list(pr$F[s,,drop=FALSE],pr$F[t,,drop=FALSE])
    yy <- list(Y[s,,drop=FALSE],Y[t,,drop=FALSE])
    half_residuals <- function(rows,deleted) {
      Fr <- pr$F[rows,,drop=FALSE]
      invA <- solve(crossprod(Fr))
      B <- crossprod(Fr,Y[rows,,drop=FALSE])
      v <- lapply(ff,function(x) x%*%invA)
      answer <- vector("list",2L)
      for (target in 1:2) {
        w <- v[[target]]
        if (length(deleted)==1L) {
          j <- deleted[1]
          denominator <- 1-rowSums(ff[[j]]*v[[j]])
          if (any(denominator<=0)) stop("Singular pair-deleted half fit")
          w <- w+v[[j]]*(rowSums(w*ff[[j]])/denominator)
        } else if (length(deleted)==2L) {
          aa <- 1-rowSums(ff[[1]]*v[[1]])
          dd <- 1-rowSums(ff[[2]]*v[[2]])
          bb <- -rowSums(ff[[1]]*v[[2]])
          determinant <- aa*dd-bb^2
          if (any(determinant<=0)) stop("Singular pair-deleted half fit")
          b1 <- rowSums(w*ff[[1]]);b2 <- rowSums(w*ff[[2]])
          w <- w+v[[1]]*((b1*dd-b2*bb)/determinant)+
            v[[2]]*((b2*aa-b1*bb)/determinant)
        }
        E <- yy[[target]]-w%*%B
        for (j in deleted) E <- E+yy[[j]]*rowSums(w*ff[[j]])
        E <- sweep(E,2,sqrt(D),"/")
        rr <- sqrt(rowSums(E^2))
        if (any(rr<=0)) stop("Zero leave-two-out radius")
        answer[[target]] <- E/rr
      }
      answer
    }
    a <- half_residuals(spec$rows1,spec$del1)
    b <- half_residuals(spec$rows2,spec$del2)
    forward <- rowSums(a[[1]]*b[[2]])
    reverse <- rowSums(a[[2]]*b[[1]])
    total <- total+sum(pr$h[s]^2*pr$h[t]^2*(forward^2+reverse^2))
  }
  ans <- N^2/(pr$h2*(pr$h2-1))*total
  if (!is.finite(ans) || ans<=0) stop("Nonpositive split trace")
  ans
}

inmax_bootstrap_bias <- function(E,pr,B=100L,multipliers=NULL,tol=1e-7,
                                 max_iter=500L,batch_size=20L) {
  T <- nrow(E); N <- ncol(E); values <- numeric(B)
  if (!is.null(multipliers) && !identical(dim(multipliers),as.integer(c(B,T,N))))
    stop("Bootstrap multipliers must be a B by T by N array")
  for (start in seq.int(1L,B,by=batch_size)) {
    finish <- min(B,start+batch_size-1L); bb <- finish-start+1L
    group <- rep(seq_len(bb),each=T)
    if (is.null(multipliers)) {
      xi <- matrix(sample(c(-1,1),bb*T*N,replace=TRUE),bb*T,N)
    } else {
      xi <- matrix(aperm(multipliers[start:finish,,,drop=FALSE],c(2,1,3)),bb*T,N)
    }
    Z <- E[rep(seq_len(T),bb),,drop=FALSE]*xi
    Zwide <- matrix(Z,nrow=T)
    Z <- matrix(inmax_project(Zwide,pr),nrow=bb*T,ncol=N)
    sq <- Z^2
    D <- rowsum(sq,group,reorder=FALSE)/T
    D <- D/rowMeans(D); converged <- FALSE
    for (iteration in seq_len(max_iter)) {
      r2 <- rowSums(sq/D[group,,drop=FALSE])
      if (any(r2<=0)) stop("Zero bootstrap radius")
      updated <- N*rowsum(sq/r2,group,reorder=FALSE)/T
      moments <- updated/D
      updated <- updated/rowMeans(updated)
      eq <- max(abs(moments-1));change <- max(abs(log(updated/D)))
      D <- updated
      if (max(eq,change)<=tol) {converged <- TRUE;break}
    }
    if (!converged) stop("Bootstrap sign scale did not converge")
    U <- Z/sqrt(D[group,,drop=FALSE]);U <- U/sqrt(rowSums(U^2))
    sums <- rowsum(U*rep(pr$h,bb),group,reorder=FALSE)
    diagonal <- as.numeric(rowsum(rowSums(U^2)*rep(pr$h^2,bb),group,reorder=FALSE))
    values[start:finish] <- N/pr$h2*(rowSums(sums^2)-diagonal)
  }
  list(bias=mean(values),values=values)
}

inmax_stat_pair_scales <- function(E,tol=1e-7,max_iter=500L,block_size=2048L) {
  T <- nrow(E);N <- ncol(E);sq <- E^2
  pairs <- which(upper.tri(matrix(FALSE,T,T)),arr.ind=TRUE)
  left <- pairs[,1];right <- pairs[,2];number <- nrow(pairs)
  initial <- inmax_diag_sign(E,tol,max_iter)$D
  diagonals <- matrix(NA_real_,number,N);iterations <- equations <- numeric()
  for (start in seq.int(1L,number,by=block_size)) {
    finish <- min(number,start+block_size-1L);idx <- seq.int(start,finish);m <- length(idx)
    Di <- matrix(initial,m,N,byrow=TRUE);converged <- FALSE
    exclude1 <- cbind(seq_len(m),left[idx]);exclude2 <- cbind(seq_len(m),right[idx])
    for (iteration in seq_len(max_iter)) {
      radius2 <- tcrossprod(1/Di,sq)
      if (any(!is.finite(radius2)) || any(radius2<=0)) stop("Zero leave-two-out scale radius")
      inverse <- 1/radius2;inverse[exclude1] <- 0;inverse[exclude2] <- 0
      updated <- (N/(T-2))*(inverse%*%sq)
      if (any(!is.finite(updated)) || any(updated<=0)) stop("Nonpositive pair-diagonal update")
      eq <- max(abs(updated/Di-1));updated <- updated/rowMeans(updated)
      change <- max(abs(log(updated/Di)));Di <- updated
      if (max(eq,change)<=tol) {
        radius2 <- tcrossprod(1/Di,sq);inverse <- 1/radius2
        inverse[exclude1] <- 0;inverse[exclude2] <- 0
        final_eq <- max(abs((N/(T-2))*(inverse%*%sq)/Di-1))
        if (final_eq<=tol) {converged <- TRUE;break}
      }
    }
    if (!converged) stop("Leave-two-out pair scales did not converge")
    diagonals[idx,] <- Di;iterations <- c(iterations,iteration);equations <- c(equations,final_eq)
  }
  list(D=diagonals,left=left,right=right,iterations=max(iterations),equation=max(equations),
       minimum_diagonal=min(diagonals))
}

inmax_stat_sum <- function(Z,pr,E=NULL,cache=NULL) {
  T <- nrow(Z);N <- ncol(Z)
  if (is.null(cache)) {
    if (is.null(E)) stop("Unrestricted residuals required for pair-scale estimation")
    cache <- inmax_stat_pair_scales(E)
  }
  left <- cache$left;right <- cache$right;inverseD <- 1/cache$D
  Zi <- Z[left,,drop=FALSE];Zj <- Z[right,,drop=FALSE]
  ri2 <- rowSums(Zi^2*inverseD);rj2 <- rowSums(Zj^2*inverseD)
  if (any(ri2<=0) || any(rj2<=0)) stop("Zero restricted pair radius")
  inner <- rowSums(Zi*Zj*inverseD);kernel <- inner/(ri2*rj2)
  weighted <- pr$h[left]*pr$h[right]*kernel
  total <- sum(weighted);squares <- sum(weighted^2)
  if (!is.finite(total+squares) || squares<=0) stop("Degenerate IN pair-kernel variance")
  H <- pr$h2;Q <- 2*N/H*total;variance <- 4*N^2/H^2*squares
  S <- Q/sqrt(variance);published <- S*sqrt(T*(T-1))/H
  list(S=S,p=pnorm(S,lower.tail=FALSE),Q=Q,variance=variance,
       S_published=published,p_published=pnorm(published,lower.tail=FALSE),
       H=H,design_denominator=H^2,cache=cache)
}
