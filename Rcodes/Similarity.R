##########################################################################
# -----------------------------------------------------------------------#
# -------------------- Similarity Metrics (CF) --------------------------#
# ---------------------Author: Shaohua Zhang  ---------------------------#
##########################################################################


# sample data 1
Mov1 <- c(4,4,3,4,2)
Mov2 <- c(NA,2,NA,4,1)
Mov3 <- c(5,1,2,NA,3)
Mov4 <- c(5,NA,4,NA,5)

username = paste("U",1:5,sep="")
t <- paste("Item",1:7, sep="")

item <- data.frame(Mov1,Mov2,Mov3,Mov4, row.names=username)
user <- t(item)

# sample data 2 (Mahout in Action)
U1 <- c(5.0, 3.0, 2.5, NA, NA, NA, NA)
U2 <- c(2.0, 2.5, 5.0, 2.0, NA, NA, NA)
U3 <- c(2.5, NA, NA, 4.0, 4.5, NA, 5.0)
U4 <- c(5.0, NA, 3.0, 4.5, NA, 4.0, NA)
U5 <- c(4.0, 3.0, 2.0, 4.0, 3.5, 4.0, NA)

data2 <- data.frame(U1,U2,U3,U4,U5, row.names=paste("Item",1:7,sep=""))



###############################
# 1 - Pearson Correlation
###############################

my.user.pearson <- cor(user,use="pairwise.complete.obs",method="pearson")

# meandiff <- function(data,i){
#   if(is.matrix(data)){
#     data[,i] - mean(data[,i])
#   }
#   else{
#     data[i] - mean(data[i])
#   }
# }

meandiff <- function(data,i){
  data[,i] - mean(data[,i])
}

pearsonsim <- function(x) {
  m <- matrix(NA, nrow=ncol(x),ncol=ncol(x),dimnames=list(colnames(x),colnames(x)))
  pc <- as.data.frame(m)
  for(i in 1:ncol(x)) {
    for(j in i:ncol(x)) {
      co_rate_1 <- x[which(x[,i] & x[,j]),i]
      co_rate_2 <- x[which(x[,i] & x[,j]),j]
      co_rate <- data.frame(co_rate_1,co_rate_2)
      pc[i,j]= sum(meandiff(co_rate,1) * meandiff(co_rate,2)) / 
        (sqrt(sum(meandiff(co_rate,1)^2)) * 
        sqrt(sum(meandiff(co_rate,2)^2)))
      pc[j,i]=pc[i,j]        
    }
  }
  return(pc)
}

my.pc.user <- pearsonsim(user)
my.pc.data2 <- pearsonsim(data2)


###############################
# 2 - Spearman Correlation
###############################
my.user.spearman <- cor(user,use="pairwise.complete.obs",method="spearman")

meandiff <- function(data,i){
  data[,i] - mean(data[,i])
}

spearmansim <- function(x) {
  # initialize similarity matrix
  m <- matrix(NA, nrow=ncol(x),ncol=ncol(x),dimnames=list(colnames(x),colnames(x)))
  sp <- as.data.frame(m)
  
  for(i in 1:ncol(x)) {
    for(j in i:ncol(x)) {
      co_rate_1 <- rank(x[which(x[,i] & x[,j]),i],ties.method="average")
      co_rate_2 <- rank(x[which(x[,i] & x[,j]),j],ties.method="average")
      co_rate <- data.frame(co_rate_1,co_rate_2)      
      sp[i,j]= sum(meandiff(co_rate,1) * meandiff(co_rate,2)) / 
        (sqrt(sum(meandiff(co_rate,1)^2)) * 
        sqrt(sum(meandiff(co_rate,2)^2)))
      sp[j,i]=sp[i,j]        
    }
  }
  return(sp)
}

my.sp.user <- spearmansim(user)
my.sp.data2 <- spearmansim(data2)

###############################
# 3 - Cosine Similarity
###############################

mycosine <- function(x,y){
  c <- sum(x*y) / (sqrt(sum(x*x)) * sqrt(sum(y*y)))
  return(c)
}

cosinesim <- function(x) {
  # initialize similarity matrix
  m <- matrix(NA, nrow=ncol(x),ncol=ncol(x),dimnames=list(colnames(x),colnames(x)))
  cos <- as.data.frame(m)
  
  for(i in 1:ncol(x)) {
    for(j in i:ncol(x)) {
      co_rate_1 <- x[which(x[,i] & x[,j]),i]
      co_rate_2 <- x[which(x[,i] & x[,j]),j]  
      cos[i,j]= mycosine(co_rate_1,co_rate_2)
      cos[j,i]=cos[i,j]        
    }
  }
  return(cos)
}

my.cosine.data2 <- cosinesim(data2)
my.cosine.user <- cosinesim(user)

# mycosine(c(1,5,10),c(190,1,1))

###############################
# 3 - Euclidean Distance
###############################
myeuclidean <- function(x,y){
  tmp <- sqrt(sum((x-y)^2))
  return(1/(1+tmp))
}

# myeuclidean(c(1,3,5), c(2,4,6))

euclideansim <- function(x) {
  # initialize similarity matrix
  m <- matrix(NA, nrow=ncol(x),ncol=ncol(x),dimnames=list(colnames(x),colnames(x)))
  eucl <- as.data.frame(m)
  
  for(i in 1:ncol(x)) {
    for(j in i:ncol(x)) {
      co_rate_1 <- x[which(x[,i] & x[,j]),i]
      co_rate_2 <- x[which(x[,i] & x[,j]),j]  
      eucl[i,j]= myeuclidean(co_rate_1,co_rate_2)
      eucl[j,i]=eucl[i,j]        
    }
  }
  return(eucl)
}

my.eucl.data2 <- euclideansim(data2)
my.eucl.user <- euclideansim(user)


###############################
# 4 - Tanimoto Similarity
###############################

tanimotosim <- function(x) {
  # initialize similarity matrix
  m <- matrix(NA, nrow=ncol(x),ncol=ncol(x),dimnames=list(colnames(x),colnames(x)))
  tanimoto <- as.data.frame(m)
  
  for(i in 1:ncol(x)) {
    for(j in i:ncol(x)) {
      tanimoto[i,j]= length(which(x[,i] & x[,j])) / length(which(x[,i] | x[,j]))
      tanimoto[j,i]=tanimoto[i,j]        
    }
  }
  return(tanimoto)
}

my.tanimoto.data2 <- tanimotosim(data2)
my.tanimoto.user <- tanimotosim(user)

###############################
# 5 - LogLikelihood Similarity (LLR)
###############################

# R automatically handles the 0*log(0) situation, no need to create safeLog() function
my.entropy <- function(x){
  cellsize <- x/sum(x,na.rm=TRUE)
  shannon.entropy <- -sum(cellsize*log(cellsize),na.rm=TRUE)
  return(shannon.entropy)
}

my.gtest <- function(n,entropy_xnoNA,entropy_ynoNA,entropy_N){
  return(2*n*(entropy_xnoNA + entropy_ynoNA - entropy_N))
}

llrsim <- function(x) {
  # initialize similarity matrix
  m <- matrix(NA, nrow=ncol(x),ncol=ncol(x),dimnames=list(colnames(x),colnames(x)))
  llr <- as.data.frame(m)
  
  # total item count
  N <- nrow(x)
  
  # LLR main
  for(i in 1:ncol(x)) {
    for(j in i:ncol(x)) {
      # calculate cell count
      xANDy <- length(which(x[,i] & x[,j]))
      xNOy <- length(na.omit(x[,i])) - xANDy
      yNOx <- length(na.omit(x[,j])) - xANDy
      NOxy <- N - length(which(x[,i] | x[,j]))
      x_noNA <- xANDy + xNOy
      y_noNA <- xANDy + yNOx
      
      # calculate shannon's entropy
      entropy.x_noNA <- my.entropy(c(x_noNA,N-x_noNA))
      entropy.y_noNA <- my.entropy(c(y_noNA,N-y_noNA))
      entropy.N <- my.entropy(c(xANDy,xNOy,yNOx,NOxy))
      
      # calculate LLR    
      llr[i,j]= my.gtest(N, entropy.x_noNA, entropy.y_noNA, entropy.N)
      llr[j,i]=llr[i,j]        
    }
  }
  return(1 - 1/(1+llr)) # transform llr so that higher value associate with stronger similarity
}

my.llr.data2 <- llrsim(data2)
my.llr.user <- llrsim(user)