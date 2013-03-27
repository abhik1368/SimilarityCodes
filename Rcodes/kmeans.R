kmeans.start <- function(N,k) {
  j <- sample(x=1:N,size=k,replace=F)
  return(j)
}


#  The following functions implement three algorithms for
#  kmeans clustering.  As written, they do not return the
#  clusters or their representatives, only the corresponding
#  value of the error criterion.  However, they can easily
#  be modified to return clusters and/or representatives.

# The code is taken from Michael Trosset course on statistical learning STATS-675
#http://mypage.iu.edu/~mtrosset/Courses/675/kmeans.r



kmeans.macqueen <- function(X,j,niter=2000) {
  N <- nrow(X)
  p <- ncol(X)
  k <- length(j)
  means <- X[j,]
  X <- cbind(rep(0,N),X)
  for (j in 1:N) {
    Y <- means-matrix(rep(X[j,-1],k),byrow=T,ncol=p)
    dist2 <- apply(Y^2,1,sum)
    X[j,1] <- order(dist2)[1] 
  }
  X <- X[order(X[,1]),]
  n <- as.vector(summary(as.factor(X[,1])))
  
  r <- 0
  for (i in 1:k) {
    Y <- X[(r+1):(r+n[i]),-1]
    means[i,] <- apply(matrix(Y,ncol=p),2,mean)
    r <- r+n[i]
  }
  
  for (iter in 1:niter) {
    j <- sample(x=1:N,size=1)
    x <- X[j,-1]
    currC <- X[j,1]
    dist2 <- sum((x-means[currC,])^2)
    bestC <- currC
    for (i in (1:k)[-currC]) {
      if (n[i] > 0) {
        newdist2 <- sum((x-means[i,])^2)
      }
      if (newdist2 < dist2) {
        bestC <- i
        dist2 <- newdist2
      }
    }
    if (bestC != currC) {
      X[j,1] <- bestC
      n1 <- n[currC]
      n[currC] <- n1-1
      if (n[currC] > 0) {
        means[currC,] <- (n1*means[currC,]-x)/n[currC]
      }
      n2 <- n[bestC]
      n[bestC] <- n2+1
      means[bestC,] <- (n2*means[bestC,]+x)/n[bestC]
    }
    
  }
  
  W <- sum((X[order(X[,1]),-1]-matrix(rep(means,rep(n,p)),nrow=N))^2)
  return(list(W=W,k=k))
}


kmeans.macqueen2 <- function(X,j,niter=2000) {
  
  N <- nrow(X)
  p <- ncol(X)
  k <- length(j)
  means <- X[j,]
  X <- cbind(rep(0,N),X)
  for (j in 1:N) {
    Y <- means-matrix(rep(X[j,-1],k),byrow=T,ncol=p)
    dist2 <- apply(Y^2,1,sum)
    X[j,1] <- order(dist2)[1] 
  }
  X <- X[order(X[,1]),]
  n <- as.vector(summary(as.factor(X[,1])))
  
  r <- 0
  for (i in 1:k) {
    Y <- X[(r+1):(r+n[i]),-1]
    means[i,] <- apply(matrix(Y,ncol=p),2,mean)
    r <- r+n[i]
  }
  
  jj <- sample(x=1:N,size=min(c(N,niter)),replace=F)
  for (iter in 1:niter) {
    j <- iter %% N
    if (j<1) j <- N
    j <- jj[j]
    x <- X[j,-1]
    currC <- X[j,1]
    dist2 <- sum((x-means[currC,])^2)
    bestC <- currC
    for (i in (1:k)[-currC]) {
      if (n[i] > 0) {
        newdist2 <- sum((x-means[i,])^2)
      }
      if (newdist2 < dist2) {
        bestC <- i
        dist2 <- newdist2
      }
    }
    if (bestC != currC) {
      X[j,1] <- bestC
      n1 <- n[currC]
      n[currC] <- n1-1
      if (n[currC] > 0) {
        means[currC,] <- (n1*means[currC,]-x)/n[currC]
      }
      n2 <- n[bestC]
      n[bestC] <- n2+1
      means[bestC,] <- (n2*means[bestC,]+x)/n[bestC]
    }
    
  }
  
 W <- sum((X[order(X[,1]),-1]-matrix(rep(means,rep(n,p)),nrow=N))^2)
  return(list(W=W,k=k))
}


kmeans.som <- function(X,j) {
  N <- nrow(X)
  p <- ncol(X)
  k <- length(j)
  means <- X[j,]
  alpha <- 1
  step <- c(.96,.04)/1000
  for (block in 1:2) {
    for (iter in 1:1000) {
      j <- sample(x=1:N,size=1)
      D <- means-matrix(rep(X[j,],k),byrow=T,ncol=p)
      dist2 <- apply(D^2,1,sum)
      i <- order(dist2)[1]
      means[i,] <- (1-alpha)*means[i,]+alpha*X[j,]
      alpha <- alpha-step[block]
    }
  }
  
  
  W <- 0
  for (j in 1:N) {
    D <- means-matrix(rep(X[j,],k),byrow=T,ncol=p)
    dist2 <- apply(D^2,1,sum)
    W <- W+min(dist2)
  }
  return(W)
}
