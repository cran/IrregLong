## ----results='hide',message=FALSE, warning=FALSE------------------------------
library(IrregLong)
library(MEMSS)
library(survival)
library(geepack)
library(frailtypack)

## -----------------------------------------------------------------------------
data(Phenobarb)
Phenobarb$event <- 1-as.numeric(is.na(Phenobarb$conc))
data <- lagfn(Phenobarb, lagvars="dose", id="Subject", time="time", lagfirst = 0)
data <- lagfn(data, lagvars="dose.lag", id="Subject", time="time", lagfirst = 0)
data <- lagfn(data, lagvars="dose.lag.lag", id="Subject", time="time", lagfirst = 0)
data$dose.lag[is.na(data$dose.lag)] <- data$dose.lag.lag[is.na(data$dose.lag)]
data$dose.lag[is.na(data$dose.lag)] <- data$dose.lag.lag.lag[is.na(data$dose.lag)]
data <- data[data$event==1,]
data$id <- as.numeric(data$Subject)
data <- data[data$time<16*24,]
data <- data[order(data$id,data$time),]
head(data)

## -----------------------------------------------------------------------------
summary(tapply(data$event,data$Subject,sum))
abacus.plot(n=59,time="time",id="Subject",data=data,tmin=0,tmax=16*24,
 xlab.abacus="Time in hours",pch=16,col.abacus=gray(0.8))

## -----------------------------------------------------------------------------
i <- iiw.weights(Surv(time.lag,time,event)~Wt + Apgar + I(conc.lag>0) + conc.lag + dose.lag +
      cluster(Subject),id="Subject",time="time",event="event",data=data,
      invariant="Subject",lagvars=c("time","conc"),maxfu=16*24,lagfirst=c(0,0),first=FALSE)
i$m

## -----------------------------------------------------------------------------
i <- iiw.weights(Surv(time.lag,time,event)~Wt + ApgarInd + I(conc.lag>0) + conc.lag + dose.lag +
      cluster(Subject),id="Subject",time="time",event="event",data=data,
      invariant=c("Subject","Wt"),lagvars=c("time","conc"),maxfu=16*24,lagfirst=c(0,0),first=FALSE)
i$m
i <- iiw.weights(Surv(time.lag,time,event)~Wt + I(conc.lag>0) + conc.lag + dose.lag +
      cluster(Subject),id="Subject",time="time",event="event",data=data,
      invariant=c("Subject","Wt"),lagvars=c("time","conc"),maxfu=16*24,lagfirst=c(0,0),first=FALSE)
i$m
i <- iiw.weights(Surv(time.lag,time,event)~Wt + I(conc.lag>0) + conc.lag + dose.lag +          I(dose.lag*Wt^(-1)) + cluster(Subject),
                 id="Subject",time="time",event="event",data=data,
      invariant=c("Subject","Wt"),lagvars=c("time","conc"),maxfu=16*24,lagfirst=c(0,0),first=FALSE)
i$m
i <- iiw.weights(Surv(time.lag,time,event)~Wt + I(conc.lag>0) + conc.lag + sqrt(dose.lag/Wt) +
      cluster(Subject),id="Subject",time="time",event="event",data=data,
      invariant=c("Subject","Wt"),lagvars=c("time","conc"),maxfu=16*24,lagfirst=c(0,0),first=FALSE)
i$m
i <- iiw.weights(Surv(time.lag,time,event)~Wt + I(conc.lag>0) + conc.lag + log(1+dose.lag/Wt) +
      cluster(Subject),id="Subject",time="time",event="event",data=data,
      invariant=c("Subject","Wt"),lagvars=c("time","conc"),maxfu=16*24,lagfirst=c(0,0),first=FALSE)
i$m
i <- iiw.weights(Surv(time.lag,time,event)~Wt *( I(conc.lag>0) + conc.lag) + log(1+dose.lag/Wt) +
      cluster(Subject),id="Subject",time="time",event="event",data=data,
      invariant=c("Subject","Wt"),lagvars=c("time","conc"),maxfu=16*24,lagfirst=c(0,0),first=FALSE)
i$m
i <- iiw.weights(Surv(time.lag,time,event)~Wt + log(1+dose.lag/Wt) + log(1+dose.lag/Wt)*(I(conc.lag>0) + conc.lag) +
      cluster(Subject),id="Subject",time="time",event="event",data=data,
      invariant=c("Subject","wt"),lagvars=c("time","conc"),maxfu=16*24,lagfirst=c(0,0),first=FALSE)
i$m
i <- iiw.weights(Surv(time.lag,time,event)~I(conc.lag>0) + log(1+dose.lag/Wt) + log(1+dose.lag/Wt)*conc.lag +
      cluster(Subject),id="Subject",time="time",event="event",data=data,
      invariant=c("Subject","Wt"),lagvars=c("time","conc"),maxfu=16*24,lagfirst=c(0,0),first=FALSE)
i$m

## ----fig3, fig.height=6, fig.width=6, fig.align="center"----------------------
plot(data$time,data$conc,xlim=c(0,200),pch=16, xlab="Time in hours", ylab="Concentration")

## -----------------------------------------------------------------------------
rsq1 <- array(dim=8)
rsq1[1] <- summary(lm(conc~time,data=data))$adj.r.squared
rsq1[2] <- summary(lm(conc~I((time)^0.5),data=data))$adj.r.squared
rsq1[3] <- summary(lm(conc~I((time)^2),data=data))$adj.r.squared
rsq1[4] <- summary(lm(conc~I((time)^3),data=data))$adj.r.squared
rsq1[5] <- summary(lm(conc~log(time),data=data))$adj.r.squared
rsq1[6] <- summary(lm(conc~I((time)^(-0.5)),data=data))$adj.r.squared
rsq1[7] <- summary(lm(conc~I((time)^(-1)),data=data))$adj.r.squared
rsq1[8] <- summary(lm(conc~I((time)^(-2)),data=data))$adj.r.squared
which.max(rsq1)
rsq1[which.max(rsq1)]

## -----------------------------------------------------------------------------
rsq2 <- array(dim=8)
rsq2[1] <- summary(lm(conc~log(time)+ time,data=data))$adj.r.squared
rsq2[2] <- summary(lm(conc~log(time)+ I((time)^0.5),data=data))$adj.r.squared
rsq2[3] <- summary(lm(conc~log(time) + I((time)^2),data=data))$adj.r.squared
rsq2[4] <- summary(lm(conc~log(time)+ I((time)^3),data=data))$adj.r.squared
rsq2[5] <- summary(lm(conc~log(time)+ time:log(time),data=data))$adj.r.squared
rsq2[6] <- summary(lm(conc~log(time) + I((time)^(-0.5))*log(1+time),data=data))$adj.r.squared
rsq2[7] <- summary(lm(conc~log(time) + I((time)^(-1)),data=data))$adj.r.squared
rsq2[8] <- summary(lm(conc~log(time)+ I((time)^(-2)),data=data))$adj.r.squared
which.max(rsq2)
rsq2[which.max(rsq2)]

## -----------------------------------------------------------------------------
rsq1 <- array(dim=8)
rsq1[1] <- summary(lm(conc~time,data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[2] <- summary(lm(conc~I((time)^0.5),data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[3] <- summary(lm(conc~I((time)^2),data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[4] <- summary(lm(conc~I((time)^3),data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[5] <- summary(lm(conc~log(time),data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[6] <- summary(lm(conc~I((time)^(-0.5)),data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[7] <- summary(lm(conc~I((time)^(-1)),data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[8] <- summary(lm(conc~I((time)^(-2)),data=data,weight=i$iiw.weight))$adj.r.squared
which.max(rsq1)
rsq1[which.max(rsq1)]

## -----------------------------------------------------------------------------
rsq1 <- array(dim=8)
rsq1[1] <- summary(lm(conc~I((time)^3) +time,data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[2] <- summary(lm(conc~I((time)^3) +I((time)^0.5),data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[3] <- summary(lm(conc~I((time)^3) +I((time)^2),data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[4] <- summary(lm(conc~I((time)^3) + log(time):I(time^3),data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[5] <- summary(lm(conc~I((time)^3) +log(time),data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[6] <- summary(lm(conc~I((time)^3) +I((time)^(-0.5)),data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[7] <- summary(lm(conc~I((time)^3) +I((time)^(-1)),data=data,weight=i$iiw.weight))$adj.r.squared
rsq1[8] <- summary(lm(conc~I((time)^3) +I((time)^(-2)),data=data,weight=i$iiw.weight))$adj.r.squared
which.max(rsq1)
rsq1[which.max(rsq1)]

## -----------------------------------------------------------------------------

iiwgee <- iiwgee(conc ~ time + I(time^3) + log(time),Surv(time.lag,time,event)~I(conc.lag>0) + log(1+dose.lag/Wt) + log(1+dose.lag/Wt)*conc.lag + +cluster(id),
        formulanull=NULL,id="id",time="time",event="event",data=data,
        invariant=c("id","Wt"),lagvars=c("time","conc"),maxfu=16*24,lagfirst=c(0,0),first=FALSE)

## -----------------------------------------------------------------------------
summary(iiwgee$geefit)

## ----fig2, fig.height=6, fig.width=6, fig.align="center"----------------------
m <- geeglm(conc ~ time + I(time^3) + log(time) , id=Subject, data=data)
time <- (2:200)
unweighted <- cbind(rep(1,199),time,time^3,log(time))%*%m$coefficients
weighted <- cbind(rep(1,199),time,time^3,log(time))%*%iiwgee$geefit$coefficients
plot(data$time,data$conc,xlim=c(0,199),ylim=c(min(unweighted,weighted,data$conc),max(unweighted,weighted,data$conc)),pch=16,xlab="Time",ylab="Serum phenobarbital concentration")
lines(time,unweighted,type="l")
lines(time,weighted,col=2)
legend (0,60,legend=c("Unweighted","Inverse-intensity weighted"),col=1:2,bty="n",lty=1)

## -----------------------------------------------------------------------------
summary(iiwgee$phfit)

## ----results='hide',message=FALSE, eval=FALSE---------------------------------
#  library(geesmv)
#  
#  
#  reg <- function(data){
#    est <- summary(geeglm(conc~time + I(time^3), id=id,data=data))$coefficients[,1:2]
#    if(max(table(data$id))>1) est[,2] <- GEE.var.md(conc~time + I(time^3) , id=id,data=data)$cov.beta
#    est <- data.matrix(est)
#    return(est)
#  }
#  set.seed(301031)
#  wt <- i$iiw.weight
#  wt[wt>quantile(i$iiw.weight,0.95)] <- quantile(i$iiw.weight,0.95)
#  m.mogee <- mo(20,reg,data,wt, singleobs=FALSE,id="id",time="time",keep.first=FALSE,var=TRUE)
#  m.mogee

## ----results='hide',message=FALSE---------------------------------------------

Liangmo <- function(data,Yname,Xnames,Wnames,maxfu,baseline){
 x <- Liang(data=data,Yname=Yname,Xnames=Xnames,Wnames=Wnames,id="id",time="time",
            maxfu=maxfu,baseline=baseline); print(x); return(x)
}
data$Intercept <- 1
data$time3 <- (data$time)^3
data$logtime <- log(data$time)
data$ApgarInd.time <- as.numeric(data$ApgarInd)*data$time/24
data$ApgarInd.time3 <- as.numeric(data$ApgarInd)*((data$time/24)^3)

set.seed(301031)
ifrailty <- iiw.weights(Surv(time.lag,time,event)~Wt + I(conc.lag>0) + conc.lag +
                      I(log(1+dose.lag/Wt)) + I(log(1+dose.lag/Wt)):conc.lag    
                      +cluster(id),id="id",time="time",event="event",data=data,     
                      invariant=c("id","Wt"),lagvars=c("time","conc"),maxfu=16*24,
                      lagfirst=c(0,0), first=FALSE,frailty=TRUE)

wt <- ifrailty$iiw.weight
wt[wt>quantile(ifrailty$iiw.weight,0.95)] <- quantile(ifrailty$iiw.weight,0.95)
m.moLiang <- mo(20,Liangmo,data,wt, 
         singleobs=FALSE,id="id",time="time",keep.first=FALSE,var=FALSE,Yname="conc",
         Xnames=c("ApgarInd","ApgarInd.time","ApgarInd.time3"), 
         Wnames=c("Intercept"),maxfu=16*24,baseline=0)

m.moLiang$est

## -----------------------------------------------------------------------------
m.moLiang$est

## ----fig4, fig.height=6, fig.width=6, fig.align="center"----------------------
miiwgee <- geeglm(conc ~ ApgarInd *(time + time3),id=id,weight=wt,data=data)

summary(miiwgee)
time <- (1:200)
gp0 <- cbind(rep(1,200),time,time^3)%*%miiwgee$coefficients[c(1,3:4)]
gp1 <- cbind(rep(1,200),time,time^3)%*%(miiwgee$coefficients[c(2,5:6)] + miiwgee$coefficients[c(1,3:4)])
diffLiang <- cbind(rep(1,200),time/24,(time/24)^3)%*%(m.moLiang$est)
plot(time,gp1-gp0,xlab="Time",ylab="Difference in means",type="l",ylim=c(min(gp1-gp0,diffLiang),max(gp1-gp0,diffLiang)))
lines(time,diffLiang,col=2)
legend(10,-3,legend=c("IIW estimate","MO + Liang estimate"),lty=1,col=1:2,bty="n")


