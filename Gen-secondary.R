# rm(list=ls())
# library(mnormt)
# n=10
# sigma = 0.03
# len=40; cutt=10
# type = c(1,1,1,1,1,2,2,2,2,2)
# past = -5
# death = c(7,8,7,7,7,7,7,7,7,7)

external_generation=function(n, len, cutt, past, sigma, type, death, param1){

#### generate parameter 
k=0; B.list=c(); fr.list=c()
while(1){
  k = k+1
  #### different parameter depending on type 
  B = param1[k]
  if (type[k] == 1){fr=1}
  if (type[k] == 2){fr=0.3}
  B.list[k] = B
  fr.list[k] = fr
  if (k==n) {cat(" \n","N = ",k,"is reached\n");break}
}

#### external factor function
External_signal = function(t, Bs, frs) {sin(t/frs)*(2*Bs)}
s = seq(past, cutt, length.out=len)

#### generate data
tr=n*length(s) ;rNA = rep(NA,tr)
exdata = data.frame(id = rNA, ex = rNA, obstime = rNA, wnoise=rNA, type=rNA)
ts=0
for (i in 1:n){
  exdata[(ts+1):(ts+len),1] = rep(i, length(s))    
  exdata[(ts+1):(ts+len),2] = External_signal(s, B.list[i], fr.list[i])
  exdata[(ts+1):(ts+len),3] = s
  exdata[(ts+1):(ts+len),4] = exdata[(ts+1):(ts+len), 2] + rnorm(length(s), 0, sd=sigma)
  exdata[(ts+1):(ts+len),5] = type[i]
  ts=ts+length(s)
}

#### cut at death
# exdata = do.call("rbind",lapply(1:n,function(i){ exdata[exdata$obstime <=death[i] & exdata$id==i, ] }))

return(list(exdata=exdata, exmat=B.list, exfr=fr.list, type=type))
}

# for(i in 1:n)
# {
#   plot(exdata$obstime[exdata$id==i],exdata$ex[exdata$id==i],xlim=c(-5,10),ylim=c(-5,5),xlab=NA,ylab=NA,cex.axis=1.2)
#   par(new=T)
# }
# par(new=T)
# x=seq(-5,10,length.out=50)
# # Gt = function(t,Bs) t^2/Bs
# for(i in 1:n)
# {
#   plot(x,External_signal(x,B.list[i], fr.list[i]),xlim=c(-5,10),ylim=c(-5,5),type="l",xlab=NA,ylab=NA,cex.axis=1.2)
#   abline(v=death[i])
#   par(new=T)
# }
# abline(h=threshold)




