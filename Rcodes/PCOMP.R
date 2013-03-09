#Funcion to calculate principal components

PCOMP<-function(data){
        l<-nrow(data)
        e<-matrix(rep(1,times=l),ncol=1)
        c<-(diag(l)-(e%*%t(e)/l))%*%(data)
        centered.svd<-svd(c,nu=0)
        #print (centered.svd)
        #centered.svd
        var<-nrow(data)-1  
        vectors<-centered.svd$v  #Specify eigenvectors for svd
        #cat("\n Standard Deviation       :  ",values)
        cat ("\n components \n")
        print (vectors)
        sdev<-centered.svd$d/sqrt(var) #Calculate standard deviation of new variables
        eigenval<-sdev^2
        cat ("\n eigen Values            :   ",eigenval)
        cat ("\n Standard Deviation      :   ",sdev)
        scores<-c%*%vectors #Calculate the new scores
        print (scores)
        total.var<-sum(diag(cov(c)))  #Calculate total variance                                                                  
        prop.var<-rep(NA,ncol(data))
        cum.var<-rep(NA,ncol(data))
        for(i in 1:ncol(data)){prop.var[i]<-var(scores[,i])/total.var}
        cat ("\n Proportion of variance  :   " ,prop.var)
        for(i in 1:ncol(data)){cum.var[i]<-sum(prop.var[1:i])}
        cat ("\n cumulative variance     :   ",cum.var)
                
}