mds.stress.raw.eq <- function(X,Delta) {
  #
  # raw stress criterion with equal weights
  #
  print(sum((as.vector(dist(X))-as.vector(as.dist(Delta)))^2))
}


mds.guttman.eq <- function(X,Delta) {
  #
  # Guttman transform with equal weights
  #
  d <- as.vector(as.dist(Delta))/as.vector(dist(X))
  n <- nrow(X)
  Delta <- matrix(1:n,nrow=n,ncol=n)
  i <- as.vector(as.dist(Delta))
  j <- as.vector(as.dist(t(Delta)))
  Y <- matrix(0,nrow=n,ncol=ncol(X))
  for (k in 1:5) {
    v <- X[i[k],]-X[j[k],]
    v <- d[k]*v
    Y[i[k],] <- Y[i[k],]+v
    Y[j[k],] <- Y[j[k],]-v
    mds.stress.raw.eq(Y,Delta)
  }
  return(Y/n)
}

d<-as.matrix(read.table("http://mypage.iu.edu/~mtrosset/Courses/675/colors.dat"))
delta<-d+t(d)
X.rand <- matrix(rnorm(14*2,sd=40),nrow=14,ncol=2)
d<-mds.stress.raw.eq(X.rand,delta)
apply(c,2,sum)
X.cmds<-cmdscale(delta)
newstress<-mds.guttman.eq(X.cmds,delta)







