dt<-read.csv("E:/New/rocssim1.csv",header=T,row.names=1)
mat<-as.matrix(dt)

d<-data.frame(col = rep(colnames(m), each = nrow(m)), 
           row = rep(rownames(m), ncol(m)), 
           value = as.vector(m))

write.csv(d,"E:/matrix_edge.csv")
m<-as.data.frame(as.table(mat))
combinations<-combn(colnames(mat), 2 , FUN = function( x ) { paste( x ,collapse = "_" ) } )
m<-m[ m$Var1 != m$Var2 , ]
m = m[ paste( m$Var1 , m$Var2 , sep = "_" ) %in% combinations , ]
write.csv(m,"E:/rocs_edge.csv")