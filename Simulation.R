source("Gen-target.R")
source("Gen-secondary.R")
library(Matrix)
library(refund)
library(nloptr)
library(MASS)
library(Corbi)
library(nlme)

########################################### System Input ##############################################
n1=45; n2=5
n=n1+n2; len=40; pct=0.250; len.ex=40; past=0;   #%=m-1/len-1  # n = Number of output len=length of each output m= observations from p
sigma = 0.3; sigma.ex = 0.05
cut=10 

ploty=30

iterations=100
type=c(rep(1,n1), rep(2,n2))
ex.type=c(rep(1,n1), rep(2,n2))

n.test=1  
test.type=c(2)
textex.type=c(2)
########################################################################################################
######################################### Generate Data #################################################
iman = generation(n=n, len=len, 
                  cutt=cut, sigma=sigma, type=type)
exdata = external_generation(n=n, len=len.ex, cutt=cut, 
                             past=past, sigma=sigma.ex, type=iman$type, 
                             death=iman$death, param1=iman$w1)

test.iman = generation(n=n.test, len=len, 
                       cutt=cut, sigma=sigma, type=textex.type)

test.exdata = external_generation(n=n.test, len=len.ex, cutt=cut, 
                                  past=past, sigma=sigma.ex, type=test.iman$type, 
                                  death=iman$death, param1=test.iman$w1)

data = iman$bdata; 
ex.data = exdata$exdata

test.data = test.iman$bdata 
test.ex.data = test.exdata$exdata

test.w1 = test.iman$w1
test.w2 = test.iman$w2

########################################## Train Test Mean ###############################################
##### unit data
train = data
trains = lapply(1:n,function(i){train$obstime[train$id==i]})
trainy = lapply(1:n,function(i){train$wnoise[train$id==i]})
trainytrue = lapply(1:n,function(i){train$r[train$id==i]})

test = test.data
tests = lapply(1:n.test,function(i){test$obstime[test$id==i]})
testy = lapply(1:n.test,function(i){test$wnoise[test$id==i]})
testytrue = lapply(1:n.test,function(i){test$r[test$id==i]})

##### external factor data
ex.train = ex.data
ex.trains = lapply(1:(n),function(i){ex.train$obstime[ex.train$id==i]})
ex.trainy = lapply(1:(n),function(i){ex.train$wnoise[ex.train$id==i]})
ex.trainytrue = lapply(1:n,function(i){ex.train$ex[ex.train$id==i]})

ex.test = test.ex.data
ex.tests = lapply(1:(n.test),function(i){ex.test$obstime[ex.test$id==i]})
ex.testy = lapply(1:(n.test),function(i){ex.test$wnoise[ex.test$id==i]})
ex.testytrue = lapply(1:n.test,function(i){ex.test$ex[ex.test$id==i]})

########################################## Plot ###############################################
#### unit data
par(mfcol=c(1,2), mai=c(0.45,0.45,0.1,0.1))
for(i in 1:(n)){
  plot(trains[[i]],trainytrue[[i]],xlim=c(0,10),ylim=c(0,40),type="l",pch=2,lwd=1,xlab=NA,ylab=NA,cex.axis=1, col="#C9C8C8", axes=F)
  par(new=T)
}
for(i in 1:(n.test)){
  plot(tests[[i]],testytrue[[i]],xlim=c(0,10),ylim=c(0,40),type="l",lty=2,pch=2,col=2,lwd=2,xlab=NA,ylab=NA,cex.axis=1)
}

for(i in 1:(n)){
  plot(trains[[i]],ex.trainytrue[[i]],xlim=c(0,10),ylim=c(-2,2),type="l",pch=2,lwd=1,xlab=NA,ylab=NA,cex.axis=1, col="#C9C8C8", axes=F)
  par(new=T)
}
for(i in 1:(n.test)){
  plot(tests[[i]],ex.testytrue[[i]],xlim=c(0,10),ylim=c(-2,2),type="l",lty=2,pch=2,col='#2FD320',lwd=2,xlab=NA,ylab=NA,cex.axis=1)
}


par(mfcol=c(2,4), mai=c(0.3,0.3,0.1,0.3))

npclist = c(2,4,4); jj=0
for (pct in c(0.25, 0.5, 0.75)){ jj = jj+1;
obs.n = len*pct
obss = lapply(1:n.test,function(i){test$obstime[test$id==i][c(1:obs.n)]})
obsy = lapply(1:n.test,function(i){test$wnoise[test$id==i][c(1:obs.n)]})
ex.data = exdata$exdata[exdata$exdata$obstime <= pct*cut,]
test.ex.data = test.exdata$exdata[test.exdata$exdata$obstime <= pct*cut, ]
########################################## FPCA for main signal ###############################################
fpca.m = t(matrix(unlist(trainy), ncol=n))
fpca = fpca.sc(Y=fpca.m, argvals=trains[[1]], var=TRUE, npc=npclist[jj])

fpca.mu = as.vector(fpca$mu)
fpca.score = fpca$scores
fpca.eigenf = fpca$efunctions
neigenf = fpca$npc

######################################## FPCA for secondary signal ############################################
fpca2.m = t(matrix(c(unlist(ex.trainy),unlist(ex.testy)), ncol=n+n.test))
fpca2 = fpca.sc(Y=fpca2.m, argvals=ex.tests[[1]], nbasis = 7) #, argvals=ex.trains[[1]])

fpca2.mu = fpca2$mu
fpca2.score = fpca2$scores
neigenf2 = fpca2$npc

fpca.norm = list()
for (j in 1:(n+n.test))
{ 
  fpca.norm[[j]] = unlist(lapply(1:(n+n.test), function(i){sum((fpca2.score[i,]-fpca2.score[j,])^2)})) 
}
fpca.norm = matrix(unlist(fpca.norm), ncol=n+n.test)

XX = submatrix(fpca.norm, rows=c(1:n), cols=c(1:n))

######################################### Covariance Matrix ###################################
covf = function(dmat, H)
{
  H[1]^2*exp(-0.5*dmat/H[2]^2)
}
C=function(dmat,H){
  covf(dmat, H) + H[3]^2*diag(1,nrow=dim(dmat)[1], ncol=dim(dmat)[1])
}

########################################### Optimization ###################################
logL=function(H,fn)
{
  B=C(XX,H)
  deter=det(B)
  if(deter>0) {a=0.5*(log(deter)+t(y)%*%solve(B,y)+log(2*pi)*leny)
  } else {
    ch=chol(B)
    logdeter=2*(sum(log(diag(ch))))
    a=0.5*(logdeter+t(y)%*%solve(B,y)+log(2*pi)*leny)
  }
  return(as.numeric(a))
}
logL_grad=function(H,fn)
{
  return(nl.grad(H,fn))
}

H0=list()
for (j in 1:neigenf)
{
  y = fpca.score[,j]; leny = length(y)
  x0 = c(1,1,1)
  #nloptr::mma #NLOPT_LN_NEWUOA #local_opts <- list( "algorithm" = "NLOPT_LN_NEWUOA","xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 2000)#,  print_level=3) #print_level=3
  t=proc.time()
  one = nloptr(x0=x0, eval_f=logL, eval_grad_f=logL_grad, opts=opts, fn=logL )
  proc.time()-t
  H0[[j]]=one$solution
}

################################### prediction by historical data ####################################
### TRUE
xstar = fpca$argvals; ytrue=list(); par(new=F)
for (i in (1:n.test))
{
  if (textex.type[i] == 1) {ytrue[[i]] = 0.3*xstar^2- 2*sin(test.w1[i]*pi*xstar^0.85) +3*(atan(xstar-5) + pi/2) + test.w2[i]} 
  if (textex.type[i] == 2) {ytrue[[i]] = 0.3*xstar^2- 2*sin(test.w1[i]*pi*xstar) + test.w2[i]}
}

### Prediction Variance and Mean
beta.hat=list();sigma.hat=list()
for (j in 1:neigenf)
{
  MAT = C(fpca.norm, H0[[j]])
  K = submatrix(MAT, rows=c(1:n), cols=c(1:n))
  # K.star = submatrix(MAT, rows=c((n+1):(n+n.test)), cols=c(1:n))
  K.star = matrix(submatrix(MAT, rows=c((n+1):(n+n.test)), cols=c(1:n)), ncol=n)
  K.star.star = submatrix(MAT, rows=c((n+1):(n+n.test)), cols=c((n+1):(n+n.test)))
  
  seok = solve(K, fpca.score[,j])
  beta.hat[[j]] = K.star %*% seok
  sigma.hat[[j]] = diag(K.star.star-K.star%*%solve(K,t(K.star)))
}

################### incorporating new observed data via bayesian linear model  ####################

i=1;
sigma.prior = fpca$sigma2
prior.mu = unlist(beta.hat)
prior.lambda = diag(unlist(sigma.hat))

P = fpca.eigenf[1:obs.n,]
S = obsy[[i]]
mu = fpca.mu[1:obs.n]

C = solve( 1/sigma.prior * t(P) %*% P + solve(prior.lambda) )
d = solve(prior.lambda) %*% prior.mu + 1/sigma.prior * t(P) %*% (S - mu) 

pst.beta.hat = C %*% d
pst.sigma.hat = C


### Prediction
ypred=list(); yvar=list(); ypred.prior=list()
for (i in 1:n.test)
{
  ypred[[i]] = fpca.mu + fpca.eigenf %*% pst.beta.hat
  ypred.prior[[i]] = fpca.mu + fpca.eigenf %*% prior.mu
  yvar[[i]] = sigma.prior + diag(fpca.eigenf %*% pst.sigma.hat %*% t(fpca.eigenf))
}

par(new=F)
if(pct==0.25){
  for(i in 1:(n)){
    plot(trains[[i]],trainytrue[[i]],xlim=c(0,10),ylim=c(0,ploty),type="l",pch=2,lwd=1,xlab=NA,ylab=NA,cex.axis=1, col="#C9C8C8", axes=F)
    par(new=T)
  }
  for(i in 1:(n.test)){
    plot(tests[[i]],testytrue[[i]],xlim=c(0,10),ylim=c(0,ploty),type="l",lty=2,pch=2,col=2,lwd=2,xlab=NA,ylab=NA,cex.axis=1)
  }
  par(new=T)
  plot(xstar,ypred.prior[[i]],xlim=c(0,10),ylim=c(0,ploty),"l",lty=1,lwd=2,cex.axis=1,ann=FALSE,col='blue')
  par(new=T)
}

par(new=F)
if(pct==0.25){
  for(i in 1:(n)){
    plot(trains[[i]],trainytrue[[i]],xlim=c(0,10),ylim=c(0,ploty),type="l",pch=2,lwd=1,xlab=NA,ylab=NA,cex.axis=1, col="#C9C8C8", axes=F)
    par(new=T)
  }
  for(i in 1:(n.test)){
    plot(tests[[i]],testytrue[[i]],xlim=c(0,10),ylim=c(0,ploty),type="l",lty=2,pch=2,col=2,lwd=2,xlab=NA,ylab=NA,cex.axis=1)
  }
  par(new=T)
  plot(xstar,fpca.mu,xlim=c(0,10),ylim=c(0,ploty),"l",lty=1,lwd=2,cex.axis=1,ann=FALSE,col='blue')
  par(new=T)
}

### Plots
par(new=F)
plot(0,0, xlim=c(0,10), ylim=c(0,ploty), type='p', pch=20, col='white' );par(new=T)
polygon(c(xstar, rev(xstar)), c(ypred[[i]]-3*sqrt(yvar[[i]]), rev(ypred[[i]]+3*sqrt(yvar[[i]]))),
        col="#B9E3FF",border = NA);par(new=T)

plot(xstar,ypred[[i]],xlim=c(0,10),ylim=c(0,ploty),"l",lty=1,lwd=2,cex.axis=1,ann=FALSE,col=1)
par(new=T)
plot(xstar,ytrue[[i]],xlim=c(0,10),ylim=c(0,ploty),"l",lty=2,lwd=2.5,cex.axis=1,ann=FALSE,col="red", axes=F)
par(new=T)
plot(obss[[i]],obsy[[i]],xlim=c(0,10),ylim=c(0,ploty),'p',cex=1,pch=16,col=1,lwd=1,xlab=NA,ylab=NA,cex.axis=1, axes=F)

sigma.prior = fpca$sigma2
prior.mu = rep(0,neigenf)
prior.lambda = diag(unlist(fpca$evalues))

C = solve( 1/sigma.prior * t(P) %*% P + solve(prior.lambda) )
d = 1/sigma.prior * t(P) %*% (S - mu) 

pst.beta.hat = C %*% d
pst.sigma.hat = C


### Prediction
ypred=list(); yvar=list(); ypred.prior=list()
for (i in 1:n.test)
{
  ypred[[i]] = fpca.mu + fpca.eigenf %*% pst.beta.hat
  yvar[[i]] = sigma.prior + diag(fpca.eigenf %*% pst.sigma.hat %*% t(fpca.eigenf))
}

### Plots
par(new=F)
plot(0,0, xlim=c(0,10), ylim=c(0,ploty), type='p', pch=20, col='white' );par(new=T)
polygon(c(xstar, rev(xstar)), c(ypred[[i]]-3*sqrt(yvar[[i]]), rev(ypred[[i]]+3*sqrt(yvar[[i]]))),
        col="#B9E3FF",border = NA);par(new=T)
plot(xstar,ypred[[i]],xlim=c(0,10),ylim=c(0,ploty),"l",lty=1,lwd=2,cex.axis=1,ann=FALSE,col=1)
par(new=T)
plot(xstar,ytrue[[i]],xlim=c(0,10),ylim=c(0,ploty),"l",lty=2,lwd=2.5,cex.axis=1,ann=FALSE,col="red", axes=F)
par(new=T)
plot(obss[[i]],obsy[[i]],xlim=c(0,10),ylim=c(0,ploty),'p',cex=1,pch=16,col=1,lwd=1,xlab=NA,ylab=NA,cex.axis=1, axes=F)
}
