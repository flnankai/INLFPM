ac_simulation_grid <- function() {
  out <- expand.grid(T=c(60L,120L), N=c(100L,200L,400L),
    distribution=c("normal","t3","mixture","IC-t3"),
    KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  code <- match(out$distribution,c("normal","t3","mixture","IC-t3"))
  out <- out[order(out$T,out$N,code),,drop=FALSE]
  rownames(out) <- NULL
  out
}

ac_driver_positive_integer <- function(x,name) {
  if (length(x)!=1L || !is.finite(x) || x<1 || x!=floor(x) ||
      x>.Machine$integer.max) stop(paste(name,"must be a positive integer"))
  as.integer(x)
}

ac_symmetric_error_root <- local({
  cache <- new.env(parent=emptyenv())
  function(N) {
    N <- ac_driver_positive_integer(N,"N")
    key <- as.character(N)
    if (exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache))
    Sigma <- .5^abs(outer(seq_len(N),seq_len(N),"-"))
    eig <- eigen(Sigma,symmetric=TRUE)
    if (any(eig$values<=0)) stop("The AR correlation matrix is not positive definite")
    root <- eig$vectors%*%(sqrt(eig$values)*t(eig$vectors))
    root <- (root+t(root))/2
    assign(key,root,envir=cache)
    root
  }
})

ac_draw_factors <- function(T) {
  T <- ac_driver_positive_integer(T,"T")
  intercept <- c(.53,.19,.19)
  autoregression <- c(.06,.19,.05)
  variance_intercept <- c(.89,.62,.80)
  variance_persistence <- c(.85,.74,.76)
  innovation_coefficient <- c(.11,.19,.15)
  f_previous <- rep(0,3L)
  variance_previous <- rep(1,3L)
  z_previous <- rep(0,3L)
  factors <- matrix(NA_real_,T+50L,3L)
  for (index in seq_len(T+50L)) {
    # The recursion uses the previous standard-normal innovation squared,
    # exactly as max_core.py does; it is not multiplied by previous variance.
    variance <- variance_intercept+variance_persistence*variance_previous+
      innovation_coefficient*z_previous^2
    z <- rnorm(3L)
    current <- intercept+autoregression*f_previous+sqrt(variance)*z
    factors[index,] <- current
    f_previous <- current
    variance_previous <- variance
    z_previous <- z
  }
  factors[seq.int(51L,T+50L),,drop=FALSE]
}

ac_draw_data <- function(T,N,distribution) {
  T <- ac_driver_positive_integer(T,"T")
  N <- ac_driver_positive_integer(N,"N")
  distribution <- match.arg(distribution,c("normal","t3","mixture","IC-t3"))
  F <- ac_draw_factors(T)
  if (distribution=="IC-t3") {
    independent_components <- matrix(rt(T*N,df=3),T,N,byrow=TRUE)
    # The symmetric positive-definite square root is consequential for this
    # nonelliptical model.  A Cholesky root defines a different joint law.
    epsilon <- independent_components%*%t(ac_symmetric_error_root(N))
  } else {
    z <- matrix(rnorm(T*N),T,N,byrow=TRUE)
    epsilon <- matrix(0,T,N)
    epsilon[,1L] <- z[,1L]
    if (N>1L) for (j in seq.int(2L,N)) {
      epsilon[,j] <- .5*epsilon[,j-1L]+sqrt(.75)*z[,j]
    }
    if (distribution=="t3") {
      # One chi-square multiplier per time row: the coordinate marginal is
      # Student t_3/sqrt(3), with unit variance before cross-sectional mixing.
      epsilon <- sweep(epsilon,1L,sqrt(rchisq(T,df=3)),"/")
    } else if (distribution=="mixture") {
      multiplier <- ifelse(runif(T)<.8,1,3)/sqrt(2.6)
      epsilon <- sweep(epsilon,1L,multiplier,"*")
    }
  }
  beta <- matrix(1,N,3L)
  Y <- epsilon+tcrossprod(F,beta)
  list(Y=Y,F=F)
}

ac_simulation_seed <- function(T,N,distribution,replicate,
                               base_seed=20261109L) {
  t_index <- match(T,c(60L,120L))
  n_index <- match(N,c(100L,200L,400L))
  d_index <- match(distribution,c("normal","t3","mixture","IC-t3"))
  if (anyNA(c(t_index,n_index,d_index))) stop("Cell is outside the stated 24-cell grid")
  replicate <- ac_driver_positive_integer(replicate,"replicate")
  base_seed <- ac_driver_positive_integer(base_seed,"base_seed")
  cell_index <- (t_index-1L)*12L+(n_index-1L)*4L+d_index
  # Double arithmetic avoids R integer overflow.  The mapping is independent
  # of cell execution order and is unique over the stated 24 x 2,000 design.
  as.integer(1+((as.double(base_seed)-1)+(cell_index-1)*1000003+
                 (as.double(replicate)-1)*104729)%%2147483646)
}

inmax_settings <- function(T,N) {
  ans <- data.frame(design="null",s=0L,delta=0)
  if (T==60L && N %in% c(100L,200L))
    ans <- rbind(ans,data.frame(design="sparsity",s=seq(1,46,by=5),delta=.5))
  if (T==60L && N==200L) for (s in c(2L,10L,25L))
    ans <- rbind(ans,data.frame(design="strength",s=s,delta=seq(.1,1,by=.1)))
  ans
}
