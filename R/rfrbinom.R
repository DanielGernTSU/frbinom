rfrbinom <- function(n, max1, p0, h0, c0) {
  if (1<p0 || p0<0) {
    stop("Invalid value for p0 parameter.")
  }

  if (1<h0 || h0<0) {
    stop("Invalid value for h0 parameter.")
  }

  if (min(.5*(-2*p0+2^(2*h0-2)+(4*p0-p0*2^(2*h0)+2^(4*h0-4))^(1/2)),1-p0) < c0
      || c0 < 0) {
    stop("Inavlid value for c0 parameter.")
  }

  r<-c()

  r<-runif(n,0,1)

  r1<-c()

  cum.dist <- pfrbinom(max1, p0, h0, c0)

  for(i in 1:n) {
    r1[i]<-min(which(r[i] <= cum.dist))-1
  }

  r1
}
