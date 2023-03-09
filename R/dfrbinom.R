dfrbinom<-function(n, p0, h0, c0) {
  if (1<p0 || p0<0) {
    stop("Invalid value for p0 parameter.")
  }

  if (1<h0 || h0<0) {
    stop("Invalid value for h0 parameter.")
  }

  if (min(.5*(-2*p0+2^(2*h0-2)+(4*p0-p0*2^(2*h0)+2^(4*h0-4))^(1/2)),1-p0) > c0
      || c0 < 0) {
    stop("Inavlid value for c0 parameter.")
  }

  max2<-n-1

  p.fun<-function(X){
    p<-X[1]
    H<-X[2]
    c<-X[3]
    P<-rep(0,n)
    P[1]<-p+c
    d<-0
    for(i in 2:n){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(p+c*(i-j)^(2*H-2))*P[j]}
    P[i]<-p+c*i^(2*H-2)-d }; P;   }

  p.fun_0<-function(X){
    p<-X[1]
    H<-X[2]
    c<-X[3]
    P<-rep(0,n)
    P[1]<-p
    d<-0
    for(i in 2:n){ d<-0; i1<-i-1; for(j in 1:i1){d<-d+(p+c*(i-j)^(2*H-2))*P[j]}
    P[i]<-p-d }; P;   }

  theo.p<-theo.p_0<-NULL
  PPm<-matrix(0, ncol=n, nrow=n); PP1<-c()
  theo.p_00<-p.fun_0(c(p0, h0, c0))
  theo.p<-p.fun(c(p0, h0, c0))
  theo.p_0<-1-cumsum(theo.p)
  PPm[1,]<-theo.p_00
  for(k in 1:max2){  k1<-k+1; PPm[2:k1,k1]<-PPm[1:k,1:k]%*%rev(theo.p[1:k])  }
  PP1<-PPm[,n]+c(PPm[1:max2, 1:max2]%*%rev(theo.p_0[1:max2]),0)

  c(1-sum(theo.p_00),PP1)   ## it gives you pmf ##
}
