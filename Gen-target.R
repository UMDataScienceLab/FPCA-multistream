# rm(list=ls())
# # library(mnormt)
# n=7
# sigma = 0.03
# threshold=14;len=30;cutt=10;type=c(1,1,1,1,2,2,2)
# # 
generation=function(n, len, cutt, sigma, type)
{

k=0
B1=c(); B2=c(); B3=c()
while(1)
 {
  k=k+1
  if (type[k] == 1){
    mean.w1=0.4; sd.w1=0.03
    u.min=0; u.max=5

    # pow=rnorm(1, mean=mean.pow, sd=sd.pow)
    w1=rnorm(1, mean=mean.w1, sd=sd.w1)
    w2=runif(1, min=u.min, max=u.max)
    
    if(w1<mean.w1-1.96*sd.w1|| w1>mean.w1+1.96*sd.w1){k = k-1;next}
    # if(pow<mean.pow-1.96*sd.pow|| pow>mean.pow+1.96*sd.pow){k = k-1;next}
  } else if (type[k] == 2){
    mean.w1=0.4; sd.w1=0.03
    u.min=1.5; u.max=6.5
    
    # pow = NA
    w1=rnorm(1, mean=mean.w1, sd=sd.w1)
    w2=runif(1, min=u.min, max=u.max)
    
    if(w1<mean.w1-1.96*sd.w1|| w1>mean.w1+1.96*sd.w1){k = k-1;next}
  } 
  # w1=rnorm(1, mean=mean, sd=sd)
  # w2=runif(1, min=u.min, max=u.max)
  # pow=rnorm(1, mean=mean.pow, sd=sd.pow)
  
  # if(w1<mean.w1-1.96*sd.w1|| w1>mean.w1+1.96*sd.w1){k = k-1;next}
  # if(pow<mean.pow-1.96*sd.pow|| w1>mean.pow+1.96*sd.pow){k = k-1;next}
  # death[k]=sqrt(B*threshold) #####
  B1[k]=w1
  B2[k]=w2
  # B3[k]=pow
  if (k==n) {cat(" \n","N = ",k,"is reached\n");break}
}

Gt = function(t, ws1, ws2, tt) {
  # if (tt == 1) return(pow*(t-1)^(2) + ws2) + 
  if (tt == 1) return(0.3*t^2- 2*sin(ws1*pi*t^0.85) +3*(atan(t-5) + pi/2) + ws2) 
  if (tt == 2) return(0.3*t^2- 2*sin(ws1*pi*t) + ws2)#
  }

death=c(rep(cutt,10))
# for (i in 1:n){
# Time = function(t) {0.3*t^(2+B3[i]) - 2*sin(B1[i]*pi*t) + B2[i] - threshold}
# root=try(uniroot(Time,c(0,cutt)))
# if(inherits(root, "try-error") || root$root>cutt){k = k-1;next}
# death[i]=root$root
# }


s=seq(0,cutt,length.out=len)
tr=n*len;rNA = rep(NA,tr)
bdata = data.frame(id=rNA, r=rNA, obstime=rNA, deathtime=rNA, wnoise=rNA, type=rNA)
ts=0
for (i in 1:n)
 {
  bdata[(ts+1):(ts+len),1] = rep(i,len)    
  bdata[(ts+1):(ts+len),2] = Gt(s,B1[i], B2[i], type[i]) 
  bdata[(ts+1):(ts+len),3] = s  
  bdata[(ts+1):(ts+len),4] = rep(death[i],len)
  bdata[(ts+1):(ts+len),5] = bdata[(ts+1):(ts+len),2]+rnorm(len,0,sd=sigma)
  bdata[(ts+1):(ts+len),6] = type[i]
  
  ts=ts+len
}

# a=do.call("rbind",lapply(1:n,function(i){  bdata$obstime[bdata$r>threshold & bdata$id==i][1]  }))
# bdata=do.call("rbind",lapply(1:n,function(i){bdata[bdata$obstime <=a[i] & bdata$id==i,]} ))
# 
bdata
return(list(bdata=bdata, w1=B1, w2=B2, death=death, type=type))

}


# dev.off()
# for(i in 1:n)
# {
#   plot(bdata$obstime[bdata$id==i],bdata$r[bdata$id==i],xlim=c(0,10),ylim=c(0,40),xlab=NA,ylab=NA,cex.axis=1.2)
#   par(new=T)
# }
# par(new=T)
# x=seq(0,10,length.out=50)
# # Gt = function(t,Bs) t^2/Bs
# for(i in 1:n)
# {
#   plot(x,Gt(x,B1[i], B2[i], B3[i]),xlim=c(0,10),ylim=c(0,40),type="l",xlab=NA,ylab=NA,cex.axis=1.2)
#   abline(v=death[i])
#   par(new=T)
# }
# abline(h=threshold)






