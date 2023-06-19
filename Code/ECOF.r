# This ECOF package is coded by Feng, Yang (yang.feng@ec.gc.ca) with technical support and scientific review by Qiuzi Han Wen(HanQiuzi.Wen@ec.gc.ca).

#Main functions
#This package contains R functions to conduct detection and attribution analysis using three different algorithms under the optimal fingerprint framework. These include:
#1.ols(), the Ordinary Least Squares method (Allen and Tett, 1999). 
#2.tls.A03 (), the Total Least Squares method (Allen and Scott, 2003). Note that the confidence intervals for the scaling factors are obtained using the method provided in the ROF package by Dr. A. Ribes (Ribes, A., 2012)  
#3.tls.ROF (), the Regularized Optimal Fingerprint method (Ribes et al, 2013a). This function is translated from routines in the ROF package V0.8 coded by Dr. A. Ribes, using SCILAB. (Ribes, A., 2012)
 

setClass(Class='ECOF',
  representation=representation(
    X='matrix',
    Y='matrix',
    noise1='matrix',
    noise2='matrix'
  )
)
# all the required data will be included in class tls, and corresponding 
# validate function: checkOF() will apply checking rules when claim a 
# new class. Later on, user might use different function: tls.A03 or tls.ROF
# using tls class as input data, and also parameters for each function.

checkOF<-function(object){
  if(length(object@X)==0) return('input X is empty')
  if(length(object@Y)==0) return('input Y is empty')
  if(length(object@Y)==0) return('input Y is empty')
  if(length(object@Y)==0) return('input Y is empty')
  if(ncol(object@Y)!=1) return('Y should contains 1 column only')
  np<-nrow(object@X)
  nx<-ncol(object@X)
  if(nrow(object@Y)!=np) return('input Y and X sample size differ')
  if(ncol(object@noise1)!=np) return('input noise1 column number not equal observation number')
  if(ncol(object@noise2)!=np) return('input noise2 column number not equal observation number')
  return(TRUE)
}
setValidity('ECOF',checkOF)

readin<-function(file.obssig,file.noise1,file.noise2){
# read from data files and return class of tls
# sample command lines:
# x<-readin('obs_sig.txt','noise1.txt','noise2.txt')
# x will be input for ols, tls.A03 or tls.ROF
  itmp<-try(obssig<-as.matrix(read.table(file.obssig)),T)
  if(inherits(itmp, "try-error")) stop(paste('read file:',file.obssig,'error'))
  Y<-t(t(obssig[1,]))
  X<-if(nrow(obssig)>2) t(obssig[-1,]) else t(t(obssig[-1,]))
  itmp<-try(noise1<-as.matrix(read.table(file.noise1)),T)
  if(inherits(itmp, "try-error")) stop(paste('read file:',file.noise1,'error'))
  itmp<-try(noise2<-as.matrix(read.table(file.noise2)),T)
  if(inherits(itmp, "try-error")) stop(paste('read file:',file.noise2,'error'))
  np<-length(Y)
  if(ncol(noise1)!=np) stop('dim of noise1 not match obssig')
  if(ncol(noise2)!=np) stop('dim of noise2 not match obssig')
  return(new(Class='ECOF',X=X,Y=Y,noise1=noise1,noise2=noise2))
}

Creg<-function(Cn){
# regularization for noise structure
# input Cn is an n x p matrix, n is number of pieces of noise; 
# nx is point number for each piece
  nx<-ncol(Cn)
  n<-nrow(Cn)
  CE<-t(Cn)%*%Cn/n
  Ip<-diag(1,nx,nx)
  m<-mean(diag(CE))
  XP<-CE-diag(m,nx,nx)
  d2<-mean(diag(XP%*%t(XP)))
  bt<-rep(NA,n)
  for(i in 1:n){
    Mi<-t(t(Cn[i,]))%*%t(Cn[i,]) # Mi is nx x nx matrix
    bt[i]<-mean(diag((Mi-CE)%*%t(Mi-CE)))
  }
  bb2<-sum(bt)/n^2
  b2<-min(bb2,d2)
  a2<-d2-b2
  Creg<-b2*m/d2*Ip+a2/d2*CE
  return(Creg)
}

# given space-temporal structure and p1,p2 as base period start and end point
# location, return EOF(,nt-1) for reduce dimension operator
redop<-function(nt,p1,p2){
  if(p1<1|p1>nt) stop('input p1 error in reduce')
  if(p2<1|p2>nt) stop('input p2 error in reduce')
  M<-diag(1,nt)
  M[,p1:p2]<-M[,p1:p2]-1/(p2-p1+1)
  eM<-eigen(M)
  if(abs(eM$values[nt])>1E-5) stop('ev error in reduce')
  if(any(abs(eM$values[-nt])<1E-5)) stop('ev error2 in reduce')
  u<-eM$vector[,-nt]
  return(u)
}

# given reduce dimension operator u and data vector x and space-temporal structure nt,ns
# return reduced vector (nt-1)*ns
redvec<-function(u,x,nt,ns,timefirst=T){
  if(length(x)!=nt*ns) stop('input dim error in redvec')
  if(all.equal(dim(u),c(nt,nt-1))!=T) stop('input u dim error in redvec')
  x1<-matrix(x,nrow=nt,byrow=timefirst)
  ox<-t(u)%*%x1
  return(as.vector(t(ox)))
}
# usage:
# for example, vector of obs length as 50, listed as 10yrs x 5points (data 1~5 for year 1,
#     6~10 for year 2, so on so for, this kind of sequence should apply to all X and noise
#     pieces), then nt=10, ns=5
#     if the data was anomanies wrt year 6~9, namely baseperiod as 6~9, then p1=6, p2=9
# u<-redop(nt=10,ns=5,p1=6,p2=9)
# newY<-apply(x@Y,2,function(x,u=u,nt=10,ns=5){return(redvec(u,x,nt,ns))},u=u)
# newX<-apply(x@X,2,function(x,u=u,nt=10,ns=5){return(redvec(u,x,nt,ns))},u=u)
# newnoise1<-apply(x@noise1,1,function(x,u=u,nt=10,ns=5){return(redvec(u,x,nt,ns))},u=u)
# newnoise2<-apply(x@noise2,1,function(x,u=u,nt=10,ns=5){return(redvec(u,x,nt,ns))},u=u)
# 

# redECOF is a function to reduce everything (Y,X,noise1,noise2) from a ECOF object and return
# a new ECOF object contains reduced elements (Y,X,noise1,noise2)
redECOF<-function(x,u,nt,ns,timefirst=T){
  newY<-apply(x@Y,2,function(x,u,nt,ns,timefirst){return(redvec(u=u,x,nt=nt,ns=ns,timefirst=timefirst))},u=u,nt=nt,ns=ns,timefirst=timefirst)
  newX<-apply(x@X,2,function(x,u,nt,ns,timefirst){return(redvec(u=u,x,nt=nt,ns=ns,timefirst=timefirst))},u=u,nt=nt,ns=ns,timefirst=timefirst)
  newnoise1<-t(as.matrix(apply(x@noise1,1,function(x,u,nt,ns,timefirst){return(redvec(u=u,x,nt=nt,ns=ns,timefirst=timefirst))},u=u,nt=nt,ns=ns,timefirst=timefirst)))
  newnoise2<-t(as.matrix(apply(x@noise2,1,function(x,u,nt,ns,timefirst){return(redvec(u=u,x,nt=nt,ns=ns,timefirst=timefirst))},u=u,nt=nt,ns=ns,timefirst=timefirst)))
  return(new(Class='ECOF',X=newX,Y=newY,noise1=newnoise1,noise2=newnoise2))
}

tls.A03<-function(Y,X,noise1,noise2,nsig,nsim.CI=1000,df2=NULL,REG=FALSE,plev=.9){
# Mandatory input: Y -- n x 1 matrix
#        X -- n x nx matrix, nx is number of signals
#        noise1 -- nnoise1 x n matrix, each row of this matrix is one piece of noise
#        noise2 -- nnoise2 x n matrix, each row of this matrix is one piece of noise
#        nsig -- vector, length of nx, ensemble size for each signal 
# Optional input:
#        nsim.CI -- simulation size for confidence interval, default as 1000
#        df2  --  degree of freedom for noise2, if no input will treat noise2 pieces as 
#                 independent, namely, df2=nnoise2
#        REG  -- regularization flag, apply regularization on noise1 or not, default as FALSE
#        plev -- confidence level used in residual check, defaul as 0.9
# output as matrix (n-nx) rows (corresponding to EOF number), each row contains EOF#, (beta_low,
#        beta_hat, beta_up) for every signal
# sample command lines:
# x<-readin('obs_sig.txt','noise1.txt','noise2.txt')
# o1<-tls.A03(x@Y,x@X,x@noise1,x@noise2,nsig=5,df2=220,REG=TRUE)

  if(missing(Y)) stop('input Y missing in tls.A03')
  if(missing(X)) stop('input X missing in tls.A03')
  if(missing(noise1)) stop('input noise1 missing in tls.A03')
  if(missing(noise2)) stop('input noise2 missing in tls.A03')
  if(missing(nsig)) stop('input nsig missing in tls.A03')

  n<-length(Y); nx<-ncol(X); nn1<-nrow(noise1); nn2<-nrow(noise2)
  if(is.null(df2)) df2=nn2
  if(length(nsig)!=nx) stop("tls.A03: input nsig length not match input X column number")
  if(REG) C1<-Creg(noise1) else C1<-t(noise1)%*%noise1/nn1
  x<-X
  for(i in 1:nx) x[,i]<-X[,i]*sqrt(nsig[i])
# x is adjusted X with nsig
  eigen1<-eigen(C1)
  P<-eigen1$vectors
  for(i in 1:min(n,nn1)) P[,i]<-P[,i]/sqrt(eigen1$values[i])
  P<-t(P)
# P is pre-whitening operator
  Z<-cbind(P%*%x,P%*%Y)
  betaout<-NULL; betaCI<-NULL
  for(ith in (nx+1):n){
    zz<-Z[1:ith,]
    u<-svd(zz)
    betahat<--u$v[1:nx,(nx+1)]/u$v[(nx+1),(nx+1)]*sqrt(nsig)
    d<-u$d^2
    nd<-length(u$d)
    ui<-t(u$u[,nd])
    n2t<-noise2%*%t(P[1:ith,])
    r1.stat<-d[nd]/((ui%*%t(n2t)%*%n2t%*%t(ui))/nn2)
    dhat<-rep(NA,nx+1)
    for(i in 1:(nx+1)){
      vi<-t(u$u[,i])
      dhat[i]<-d[i]/(vi%*%t(n2t)%*%n2t%*%t(vi)/nn2)
    }
    delta.dhat<-dhat-min(dhat)
# sampling on critical ellipse to construct correspoding beta_simus,
# then get max/min value as upper/lower bound for scaling factor beta_hat
    Crit<-qf(plev,1,df2)
    if(nx>1){
      unit.tmp<-matrix(rnorm(nsim.CI*nx,0,1),nrow=nsim.CI,byrow=T)
      unit<-t(apply(unit.tmp,1,function(x){x/sqrt(sum(x^2))}))
    }
    else unit<-t(t(c(1,-1)))
    ai<-unit*sqrt(Crit)
    bi<-cbind(ai,NA)
    for(i in 1:nx) bi[,i]<-ai[,i]/sqrt(delta.dhat[i])
    bi[,nx+1]<-apply(t(t(bi[,1:nx])),1,function(x) {sqrt(1-sum(x^2))})
    nsim<-dim(bi)[1]
    betaup<-betalow<-betaup.new<-betalow.new<-rep(NA,nx)
    vc.pts<-bi%*%t(u$v)
    for(i in 1:nx) vc.pts[,i]<-vc.pts[,i]*sqrt(nsig[i])
    for(j in 1:nx){
      beta.tmp<--vc.pts[,j]/vc.pts[,(nx+1)]
      betaup[j]<-max(beta.tmp,na.rm=T); betalow[j]<-min(beta.tmp,na.rm=T)
    }
    obeta=NULL
    for(i in 1:nx) obeta<-c(obeta,betalow[i],betahat[i],betaup[i])
    betaout<-rbind(betaout,c(ith,obeta,r1.stat,(ith-nx)*qf(.05,ith-nx,df2),(ith-nx)*qf(.95,ith-nx,df2)))
  }
  cnames<-'#EOF'
  for(i in 1:nx) cnames<-c(cnames,paste('beta',i,'_',c('low','hat','up'),sep=''))
  cnames<-c(cnames,'RCstat','RClow','RCup')
  colnames(betaout)<-cnames
  return(betaout)
}

tls.ROF<-function(Y,X,noise1,noise2,nsig,nsim.CI=1000,nsim.rcc=1000,REG=TRUE,df2=NULL,plev=.9,rcc.flg=TRUE){
# Mandatory input: Y -- n x 1 matrix
#        X -- n x nx matrix, nx is number of signals
#        noise1 -- nnoise1 x n matrix, each row of this matrix is one piece of noise
#        noise2 -- nnoise2 x n matrix, each row of this matrix is one piece of noise
#        nsig -- vector, length of nx, ensemble size for calculating each signal 
# Optinal input: 
#        nsim.CI -- simulation size for confidence interval, default as 1000
#        nsim.rcc -- simulation size for null distirbution of residual check stat
#        df2  --  degree of freedom for noise2, if no input will treat noise2 pieces as 
#                 independent, namely, df2=nnoise2
#        REG  -- regularization flag, apply regularization on noise1 or not
#        rcc.flg -- flag to simulate empirical null distribution for residual test stat, default as TRUE
# sample command lines:
# x<-readin('obs_sig.txt','noise1.txt','noise2.txt')
# o1<-tls.ROF(x@Y,x@X,x@noise1,x@noise2,df2=180,REG=TRUE,nsig=5,rcc.flg=T)

  n<-length(Y); nx<-ncol(X); nn1<-nrow(noise1); nn2<-nrow(noise2)
  if(is.null(df2)) df2=nn2
  if(length(nsig)!=nx) stop("tls.ROF: input nsig length not match input X column number")
  sqrt.matrix<-function(a){ # calculate square root of square matrix
    if(nrow(a)!=ncol(a)) stop('sqrt.matrix: input matrix is not square')
    a.eig <- eigen(a)
    a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
    return(a.sqrt)
  }

  in.sqrt.matrix<-function(a){ # calculate inverse square root of square matrix
    if(nrow(a)!=ncol(a)) stop('in.sqrt.matrix: input matrix is not square')
    a.eig <- eigen(a)
    a.sqrt <- a.eig$vectors %*% diag(1/sqrt(a.eig$values)) %*% solve(a.eig$vectors)
    return(a.sqrt)
  }

  C1=if(REG) Creg(noise1) else t(noise1)%*%noise1/nn1
  C12<-t(in.sqrt.matrix(C1))
  gettlsbeta<-function(X,Y,C12,nsig){
    nx<-ncol(X)
    Z<-cbind(C12%*%(X%*%diag(sqrt(nsig),ncol=nx,nrow=nx)),C12%*%Y)
    u<-svd(Z)
    ns<-ncol(X)
    nd<-length(u$d)
    v<-u$v[,nd]
    beta0<-rep(NA,ns)
    for(i in 1:ns) beta0[i]<- -v[i]*sqrt(nsig[i])/v[ns+1]
    oout<-list()
    oout$u<-u
    oout$beta<-as.matrix(beta0,nrow=nx,ncol=1)
    return(oout)
  }
  o0<-gettlsbeta(X,Y,C12,nsig)
  beta0<-o0$beta # finished estimate beta_hat
  u<-o0$u
  nd<-length(u$d)
  d<-u$d^2
  v<-t(u$v[,nd])
  dhat<-d
  n2t<-t(C12%*%t(noise2))
  ui<-t(u$u[,nd])
  r1.stat<-d[nd]/((ui%*%t(n2t)%*%n2t%*%t(ui))/nn2) # residual test from A03
  r2.stat<-d[nd]/(sum(ui^2%*%t(n2t^2))/nn2)        # residual test from ROF

  # now start estimate CI of beta, multi-signal cases need simulations for nsig sphere
  for(i in 1:(nx+1)){
    vi<-t(u$u[,i])
    dhat[i]<-d[i]/(vi%*%t(n2t)%*%n2t%*%t(vi)/nn2)
  }
  delta.dhat<-dhat-min(dhat)
  Crit<-qf(plev,1,df2)
  if(nx>1){
    unit.tmp<-matrix(rnorm(nsim.CI*nx,0,1),nrow=nsim.CI,byrow=T)
    unit<-t(apply(unit.tmp,1,function(x){x/sqrt(sum(x^2))}))
  }
  else unit<-t(t(c(1,-1)))
  ai<-unit*sqrt(Crit)
  bi<-cbind(ai,NA)
  for(i in 1:nx) bi[,i]<-ai[,i]/sqrt(delta.dhat[i])
  bi[,nx+1]<-apply(t(t(bi[,1:nx])),1,function(x) {sqrt(1-sum(x^2))})
  nsim<-dim(bi)[1]
  betaup<-betalow<-rep(NA,nx)
  vc.pts<-bi%*%t(u$v)
  for(i in 1:nx) vc.pts[,i]<-vc.pts[,i]*sqrt(nsig[i])
  for(j in 1:nx){
    vc.d.pts<-complex(real=vc.pts[,j],imaginary=vc.pts[,(nx+1)])
    vc.d.ref<-complex(real=v[j],imaginary=v[(nx+1)])
    vprod.d<-vc.d.pts/vc.d.ref
    arg<-sort(Im(log(vprod.d)))
    delta.max1<-max(arg[-1]-arg[-length(arg)])
    delta.max<-max(delta.max1,arg[1]-arg[length(arg)]+2*pi)
    arg.ref<-Im(log(vc.d.ref))
    arg.min<-min(arg)+arg.ref; arg.max<-max(arg)+arg.ref
      betalow[j]<--1/tan(arg.min)
      betaup[j]<--1/tan(arg.max)
  }

  Tresi.simu.TLS<-function(Y,X,C,nsim=1000,nn1,nn2,nsig){
    n<-length(Y)
    nx<-ncol(X)
    C12<-t(sqrt.matrix(C))
    C12.in<-t(in.sqrt.matrix(C))
    Z<-cbind(C12.in%*%(X%*%diag(sqrt(nsig),ncol=nx,nrow=nx)),C12.in%*%Y)
    u<-svd(Z)
    ns<-ncol(X)
    nd<-length(u$d)
    v<-u$v[,dim(u$v)[2]]
    beta0<- as.matrix(-v[1:ns]*sqrt(nsig[1:ns])/v[ns+1],ncol=1,nrow=nx)
    r1.stat<-r2.stat<-rep(NA,nsim)
    for(ith in c(1:nsim)){
      Ys<-X%*%beta0+C12%*%rnorm(n,0,1)
      Xs<-X+C12%*%(matrix(rnorm(n*nx),ncol=nx,byrow=T)%*%diag(1/sqrt(nsig),ncol=nx,nrow=nx))
      noise1<-t(C12%*%matrix(rnorm(n*nn1,0,1),nrow=n,byrow=T))
      noise2<-t(C12%*%matrix(rnorm(n*nn2,0,1),nrow=n,byrow=T))
      C1hat<-Creg(noise1)
      Cs12<-sqrt.matrix(C1hat)
      Zs<-cbind(C12.in%*%(Xs%*%diag(sqrt(nsig),ncol=nx,nrow=nx)),C12.in%*%Ys)
      u1<-svd(Zs)
      nd<-length(u1$d)
      ids<-order(u1$d,decreasing=T)
      ui1<-t(u1$u[,ids[nd]])
      v1<-t(u1$v[,ids[nd]])
      n2t<-t(C12.in%*%t(noise2))
      d2<-u1$d[ids[nd]]^2
      r1.stat[ith]<-d2/((ui1%*%t(n2t)%*%n2t%*%t(ui1))/nn2)
      r2.stat[ith]<-d2/(sum(ui1^2%*%t(n2t^2))/nn2)
    }
    oout<-list()
    oout$r1.stat<-r1.stat
    oout$r2.stat<-r2.stat
    return(oout)
  }

  if(rcc.flg){
    noisea<-rbind(noise1,noise2)
    Ca<-t(noisea)%*%noisea/(nn1+nn2)
    rcc<-Tresi.simu.TLS(Y,X,Ca,nsim.rcc,nn1,nn2,nsig)
  }

  oout<-list()
  oout$beta<-beta0
  oout$r1.stat<-r1.stat
  oout$r2.stat<-r2.stat
  oout$n<-n
  oout$beta.CI<-cbind(betalow,beta0,betaup)
  if(rcc.flg) oout$rcc<-rcc
  return(oout)
}

ols<-function(Y,X,noise1,noise2,nsig,df2=NULL,plev=.95){
# sample command lines:
# x<-readin('obs_sig.txt','noise1.txt','noise2.txt')
# o1<-ols(x@Y,x@X,x@noise1,x@noise2,df2=NULL)

  n<-length(Y); nx<-ncol(X); nn1<-nrow(noise1); nn2<-nrow(noise2)
  if(length(nsig)!=nx) stop("ols: input nsig length not match input X column number")
  nmin=min(c(n,nn1,nn2))
  if(is.null(df2)) df2=nn2
  if(df2-nmin+1<0) stop('df2 is to small to carry out residual check')
  C1<-t(noise1)%*%noise1/nn1
  n1eigen<-eigen(C1)
  ev1<-n1eigen$values
  eof1<-n1eigen$vector
  p<-eof1
  for(i in 1:ncol(p)) p[,i]<-eof1[,i]/max(sqrt(ev1[i]),.00005)
  # calculate second covariance matrix (cn2) from n2
  C2<-t(noise2)%*%noise2/nn2
  PX<-t(eof1)%*%X
  PY<-t(eof1)%*%Y
  noise2t<-noise2%*%eof1
  output<-NULL
  for(ith in (nx+1):nmin){
    f<-solve(t(X)%*%p[,1:ith]%*%t(p[,1:ith])%*%X)%*%t(X)%*%p[,1:ith]%*%t(p[,1:ith])
    tbeta<-f%*%Y # value of beta_hat
    vbeta<-f%*%C2%*%t(f) # variance of beta_hat
    betalow<-tbeta-qt(plev,df2)*sqrt(diag(vbeta))*sqrt((nsig+1)/nsig)
    betaup<-tbeta+qt(plev,df2)*sqrt(diag(vbeta))*sqrt((nsig+1)/nsig)
    cn2t<-t(noise2t[,1:ith])%*%noise2t[,1:ith]/nn2
    n2teigen<-eigen(cn2t)
    ev2<-n2teigen$values
    ev2[ev2<0.00005]=0.00005
    eof2<-n2teigen$vector
    u<-PY-PX%*%tbeta # residual base on ith truncation
    tmp<-rep(NA,ith)
    for(j in 1:ith) tmp[j]=sum(u[1:ith]*eof2[1:ith,j])
    stat=sum(tmp^2/ev2[1:ith])
    dfn<-ith-nx
# dfd, and modeifed F-stat proposed in Ribes13
    dfd<-nn2-ith+1
    f1<-dfn*df2*qf(.05,dfn,dfd)/(df2-ith+1)
    f2<-dfn*df2*qf(.95,dfn,dfd)/(df2-ith+1)
    output<-rbind(output,c(ith,as.vector(t(cbind(betalow,tbeta,betaup))),stat,f1,f2))
  }
  cnames<-'#EOF'
  for(i in 1:nx) cnames<-c(cnames,paste('beta',i,'_',c('low','hat','up'),sep=''))
  cnames<-c(cnames,'RCstat','RClow','RCup')
  colnames(output)<-cnames
  return(output)
}

plotbetas<-function(betaout,...){
  cnams<-colnames(betaout)
  sid<-grepl('_hat',cnams)
  sigs<-cnams[sid]
  nsig<-length(sigs)
  sigs.low<-gsub('hat','low',sigs)
  sigs.up<-gsub('hat','up',sigs)
  for(ith in 1:nsig){
    yrange=range(betaout[,c(sigs.low[ith],sigs.up[ith])],na.rm=T)
    yrange[1]=min(0,yrange[1]); yrange[2]<-max(yrange[2],1)
    plot(betaout[,'#EOF'],betaout[,sigs[ith]],type="p",main=paste("best estimates of scaling factors for ",
         strsplit(sigs[ith],'_')[[1]][1],sep=''),ylim=yrange,
         xlab="Number of EOF patterns retained in the truncation",ylab="scaling factors",col="red",...)
    abline(h=0)
    abline(h=1,lty=3)
    for(j in 1:nrow(betaout)){
      lines(rep(betaout[j,'#EOF'],2),c(betaout[j,sigs.low[ith]],betaout[j,sigs.up[ith]]))
      lines(betaout[j,'#EOF']+c(-.1,.1),rep(betaout[j,sigs.low[ith]],2))
      lines(betaout[j,'#EOF']+c(-.1,.1),rep(betaout[j,sigs.up[ith]],2))
    }
  }
}

plotrstat<-function(betaout,...){
  yrange=range(1/betaout[,c('RCstat','RClow','RCup')],na.rm=T)
  plot(betaout[,'#EOF'],1/betaout[,'RCstat'],type="p",col="red",ylim=yrange,
    main="residual consistency test",
    xlab="Number of EOF patterns retained in the truncation",
    ylab="Cumulative ratio model/observation variance",log="y")
  lines(betaout[,'#EOF'],1/betaout[,'RClow'],type="l",lty=3,col="blue")
  lines(betaout[,'#EOF'],1/betaout[,'RCup'],type="l",lty=3,col="blue")
}

