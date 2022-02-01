

#Perform anisotropy transformation on location s
anisotransform=function(s,theta,L){
  return(matrix(c(1,0,0,1/L),nrow=2,ncol=2,byrow=T)%*%matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),nrow=2,ncol=2,byrow=T)%*%as.matrix(s))
  
}



spatmodnll_AI=function(par,Exceed.all,Exceed.Inds,c.vec,coords){
  
  
  p=dim(coords)[1]
  
  KA1=par[1]
  KA2=par[2]
  KB1=par[3]
  KB2=par[4]
  KB3=1
  KM1=par[5]
  KM2=par[6]
  KM3=par[7]
  KS1=par[8]
  KS2=par[9]
  KR1=par[10]
  KR2=par[11]
  KD1=par[12]
  KD2=par[13]
  KD3=par[14]
  KD4=1
  theta=par[15]
  L=par[16]
  if(KA1 <=0 || KA2 < 0  || KB1 <0  || KB2 < 0 || KB3 <0  || KB3 > 1  || KM2 < 0 || KM3 <= 0||KS1<=0||KS2<=0|| KR1 <= 0 || KR2 <= 0 ){return(10e10)}
  if(KD1<=0||KD2<0 ||KD3<=0 || KD4 > 1 ){return(10e10)}
  if(theta>0 || theta < -pi/2 || L <= 0){return(10e10)}
  
  coords=apply(coords,1,function(x){
    anisotransform(x,theta=theta,L=L)
  })
  coords=t(coords)
  
  
  nll=0
  n.pairs=dim(Exceed.Inds)[1]
  for( i in 1:n.pairs){
    
    cond.ind=Exceed.Inds[i,1]
    pair.inds=Exceed.Inds[i,2:3]
    h.pairs=rdist.earth(coords[c(cond.ind,pair.inds),],miles=F)
    diag(h.pairs)=0
    alphaSub=exp(-(h.pairs[1,-1]/KA1)^KA2)
    betaSub=KB3*exp(-(h.pairs[1,-1]/KB1)^KB2)
    sigsub=sqrt(2)*(1-exp(-(h.pairs[1,-1]/KS1)^KS2))
    MuSub=KM1*h.pairs[1,-1]^{KM2}*exp(-h.pairs[1,-1]/KM3)
    deltaSub=1+(KD1*h.pairs[1,-1]^{KD2}-KD4)*exp(-h.pairs[1,-1]/KD3)
    
    print(i)
    
    CorSub=Matern(h.pairs,range=KR1,smoothness=KR2)
    
    nllNC<-nllOC<-nllC<-0
    
    CorSub=(CorSub[-1,-1]-CorSub[1,-1]%*%t(CorSub[1,-1]))/((1-CorSub[1,-1]^2)^(1/2)%*%(1-t(CorSub[-1,1])^2)^(1/2))
    
    SigSub=diag(sigsub)%*%CorSub%*%diag(sigsub)
    
    if( !is.finite(sum(sigsub)) ||  !is.finite(sum(alphaSub)) || !is.finite(sum(betaSub)) || !is.finite(sum(MuSub)) || !is.finite(sum(deltaSub))){return(1e11)}
    
    #Case 1 - Neither censored
    
    tempNC=Exceed.all[[i]][[1]]
    if(is.null(dim(tempNC))){
      tempNC=t(as.matrix(tempNC))
    }
 if(length(tempNC!=0)){
    nllNC<-apply(tempNC,1,function(x){
      mu=alphaSub*x[1]+(x[1]^betaSub)*MuSub
      Sigma=diag(c(x[1]^betaSub))%*%SigSub%*%diag(c(x[1]^betaSub))
      dmvdlaplace(x[2:3]-mu,mu= c(0,0),sigmad=sqrt(diag(Sigma)),SigmaChol=chol(Sigma),Sigma=Sigma,
                  delta=deltaSub,log=T)
    })
}
    #Case 2 - Both Censored
    
    tempC=Exceed.all[[i]][[2]]
    censor=c.vec[pair.inds]
    
    if(is.null(dim(tempC))){
      tempC=t(as.matrix(tempC))
    }
    if(length(tempC!=0)){
      nllC<-apply(tempC,1,function(x){
        mu=alphaSub*x[1]+(x[1]^betaSub)*MuSub
        Sigma=diag(c(x[1]^betaSub))%*%SigSub%*%diag(c((x[1]^betaSub)))
        censorN=censor
        censorN[1]=qnorm(pdlaplace(censor[1],mu=mu[1],sigma=sqrt(Sigma[1,1]),delta=deltaSub[1],lower.T=F),mean=mu[1],
                         sd=sqrt(Sigma[1,1]),lower.tail=F)
        censorN[2]=qnorm(pdlaplace(censor[2],mu[2],sigma=sqrt(Sigma[2,2]),delta=deltaSub[2],lower.T=F),mean=mu[2],
                         sd=sqrt(Sigma[2,2]),lower.tail=F)
        
        log(pmnorm(x=censorN,mean=mu,  varcov =Sigma)[1])
      })
    }

    
    #Case 3 - One censored
    tempOC=Exceed.all[[i]][[3]]
    if(is.null(dim(tempOC))){
      tempOC=t(as.matrix(tempOC))
    }
    if(length(tempOC!=0)){
      nllOC<-apply(tempOC,1,function(x){
        Cind=which(is.na(x[2:3]))
        
        muNC= alphaSub*x[1]+(x[1]^betaSub)*MuSub
        sigmaNC=diag(c((x[1]^betaSub)))%*%SigSub%*%diag(c((x[1]^betaSub)))
        
        inv=solve( sigmaNC[-Cind,-Cind])
        xn<-qnorm(pdlaplace(x[2:3][-Cind],mu=muNC[-Cind],sigma=sqrt(sigmaNC[-Cind,-Cind]),delta=deltaSub[-Cind],lower.T=F),
                  mean=muNC[-Cind],sd=sqrt(sigmaNC[-Cind,-Cind]),lower.tail=F)
        
        muC=muNC[Cind]+ sigmaNC[Cind,-Cind]%*%inv%*%c(xn-muNC[-Cind])
        sigC= sigmaNC[Cind,Cind]- sigmaNC[Cind,-Cind]%*%inv%*% sigmaNC[-Cind,Cind]
        
        
        pnorm(qnorm(pdlaplace(censor[Cind],mu=muNC[Cind],sigma=sqrt(sigmaNC[Cind,Cind]),delta=deltaSub[Cind]),
                    mean=muNC[Cind],sd=sqrt(sigmaNC[Cind,Cind])),
              mean=c(muC),sd=sqrt(sigC),log=T)+ddlaplace(x[2:3][-Cind],mu=muNC[-Cind],sigma=sqrt(sigmaNC[-Cind,-Cind]),delta=deltaSub[-Cind],log=T)
        
      }
      )  
    }

    nll<- nll-sum(nllC,nllNC,nllOC)
    if(!is.finite(nll)){break}
    print(nll)
  }
  
  if(is.finite(nll)){
    return(nll)
  }else{return(1e10)}
  
  
}

spatmodnll_AD=function(par,Exceed.all,Exceed.Inds,c.vec,coords){
  
  
  p=dim(coords)[1]
  

  KM1=par[1]
  KM2=par[2]
  KM3=par[3]
  KS1=par[4]
  KS2=par[5]
  KS3=par[6]
  KR1=par[7]
  KR2=par[8]
  KD1=par[9]
  KD2=par[10]
  KD3=par[11]
  KD4=par[12]
  theta=par[13]
  L=par[14]
  if( KM2 < 0 || KM3 <= 0||KS1<=0||KS2<=0 ||KS3<=0 || KR1 <= 0 || KR2 <= 0 ){return(10e10)}
  if(KD1<=0||KD2<0 ||KD3<=0 || KD4 > 1 ){return(10e10)}
  if(theta>0 || theta < -pi/2 || L <= 0){return(10e10)}
  
  coords=apply(coords,1,function(x){
    anisotransform(x,theta=theta,L=L)
  })
  coords=t(coords)
  
  
  nll=0
  n.pairs=dim(Exceed.Inds)[1]
  for( i in 1:n.pairs){
    
    cond.ind=Exceed.Inds[i,1]
    pair.inds=Exceed.Inds[i,2:3]
    h.pairs=rdist.earth(coords[c(cond.ind,pair.inds),],miles=F)
    diag(h.pairs)=0
    alphaSub=c(1,1)
    betaSub=c(0,0)
    sigsub=KS3*(1-exp(-(h.pairs[1,-1]/KS1)^KS2))
    MuSub=KM1*h.pairs[1,-1]^{KM2}*exp(-h.pairs[1,-1]/KM3)
    deltaSub=1+(KD1*h.pairs[1,-1]^{KD2}-KD4)*exp(-h.pairs[1,-1]/KD3)
    
    print(i)
    
    CorSub=Matern(h.pairs,range=KR1,smoothness=KR2)
    
    nllNC<-nllOC<-nllC<-0
    
    CorSub=(CorSub[-1,-1]-CorSub[1,-1]%*%t(CorSub[1,-1]))/((1-CorSub[1,-1]^2)^(1/2)%*%(1-t(CorSub[-1,1])^2)^(1/2))
    
    SigSub=diag(sigsub)%*%CorSub%*%diag(sigsub)
    
    if( !is.finite(sum(sigsub)) ||  !is.finite(sum(alphaSub)) || !is.finite(sum(betaSub)) || !is.finite(sum(MuSub)) || !is.finite(sum(deltaSub))){return(1e11)}
    
    #Case 1 - Neither censored
    
    tempNC=Exceed.all[[i]][[1]]
    if(is.null(dim(tempNC))){
      tempNC=t(as.matrix(tempNC))
    }
 if(length(tempNC!=0)){
    nllNC<-apply(tempNC,1,function(x){
      mu=alphaSub*x[1]+(x[1]^betaSub)*MuSub
      Sigma=diag(c(x[1]^betaSub))%*%SigSub%*%diag(c(x[1]^betaSub))
      dmvdlaplace(x[2:3]-mu,mu= c(0,0),sigmad=sqrt(diag(Sigma)),SigmaChol=chol(Sigma),Sigma=Sigma,
                  delta=deltaSub,log=T)
    })
    }
    #Case 2 - Both Censored
    
    tempC=Exceed.all[[i]][[2]]
    censor=c.vec[pair.inds]
    
    if(is.null(dim(tempC))){
      tempC=t(as.matrix(tempC))
    }
    if(length(tempC!=0)){
      nllC<-apply(tempC,1,function(x){
        mu=alphaSub*x[1]+(x[1]^betaSub)*MuSub
        Sigma=diag(c(x[1]^betaSub))%*%SigSub%*%diag(c((x[1]^betaSub)))
        censorN=censor
        censorN[1]=qnorm(pdlaplace(censor[1],mu=mu[1],sigma=sqrt(Sigma[1,1]),delta=deltaSub[1],lower.T=F),mean=mu[1],
                         sd=sqrt(Sigma[1,1]),lower.tail=F)
        censorN[2]=qnorm(pdlaplace(censor[2],mu[2],sigma=sqrt(Sigma[2,2]),delta=deltaSub[2],lower.T=F),mean=mu[2],
                         sd=sqrt(Sigma[2,2]),lower.tail=F)
        
        log(pmnorm(x=censorN,mean=mu,  varcov =Sigma)[1])
      })
    }
    
    
    #Case 3 - One censored
    tempOC=Exceed.all[[i]][[3]]
    if(is.null(dim(tempOC))){
      tempOC=t(as.matrix(tempOC))
    }
    if(length(tempOC!=0)){
      nllOC<-apply(tempOC,1,function(x){
        Cind=which(is.na(x[2:3]))
        
        muNC= alphaSub*x[1]+(x[1]^betaSub)*MuSub
        sigmaNC=diag(c((x[1]^betaSub)))%*%SigSub%*%diag(c((x[1]^betaSub)))
        
        inv=solve( sigmaNC[-Cind,-Cind])
        xn<-qnorm(pdlaplace(x[2:3][-Cind],mu=muNC[-Cind],sigma=sqrt(sigmaNC[-Cind,-Cind]),delta=deltaSub[-Cind],lower.T=F),
                  mean=muNC[-Cind],sd=sqrt(sigmaNC[-Cind,-Cind]),lower.tail=F)
        
        muC=muNC[Cind]+ sigmaNC[Cind,-Cind]%*%inv%*%c(xn-muNC[-Cind])
        sigC= sigmaNC[Cind,Cind]- sigmaNC[Cind,-Cind]%*%inv%*% sigmaNC[-Cind,Cind]
        
        
        pnorm(qnorm(pdlaplace(censor[Cind],mu=muNC[Cind],sigma=sqrt(sigmaNC[Cind,Cind]),delta=deltaSub[Cind]),
                    mean=muNC[Cind],sd=sqrt(sigmaNC[Cind,Cind])),
              mean=c(muC),sd=sqrt(sigC),log=T)+ddlaplace(x[2:3][-Cind],mu=muNC[-Cind],sigma=sqrt(sigmaNC[-Cind,-Cind]),delta=deltaSub[-Cind],log=T)
        
      }
      )  
    }
    
    nll<- nll-sum(nllC,nllNC,nllOC)
    if(!is.finite(nll)){break}
    print(nll)
  }
  
  if(is.finite(nll)){
    return(nll)
  }else{return(1e10)}
  
  
}
