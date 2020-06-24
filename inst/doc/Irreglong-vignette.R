## ----results='hide',message=FALSE, warning=FALSE------------------------------
library(IrregLong)
library(MEMSS)
library(survival)
library(geepack)
library(frailtypack)
library(data.table)

## -----------------------------------------------------------------------------
data(Phenobarb)
Phenobarb$event <- 1-as.numeric(is.na(Phenobarb$conc))
data <- Phenobarb
data <- data[data$event==1,]
data$id <- as.numeric(data$Subject)
data <- data[data$time<16*24,]
data <- data[order(data$id,data$time),]
head(data)

## -----------------------------------------------------------------------------
summary(tapply(data$event,data$Subject,sum))
abacus.plot(n=59,time="time",id="Subject",data=data,tmin=0,tmax=16*24,
 xlab.abacus="Time in hours",pch=16,col.abacus=gray(0.8))

## ----fig.width=8, fig.height=12-----------------------------------------------
counts <- extent.of.irregularity(data,time="time",id="id",
  scheduledtimes=NULL, cutpoints=NULL,ncutpts=50, 
  maxfu=16*24, plot=TRUE,legendx=30,legendy=0.8,
 formula=Surv(time.lag,time,event)~1,tau=16*24)
 counts$auc

## -----------------------------------------------------------------------------
data$Apgar <- as.numeric(data$Apgar)
i <- iiw.weights(Surv(time.lag,time,event)~Wt + Apgar + 
                   I(conc.lag>0 & conc.lag<=20) + 
                I(conc.lag>20 & conc.lag<=30) + I(conc.lag>30)+
      cluster(Subject),id="Subject",time="time",event="event",data=data,
      invariant=c("Subject","Wt","Apgar"),lagvars=c("time","conc"),maxfu=16*24,
      lagfirst=c(0,0),first=FALSE)
i$m
 

## -----------------------------------------------------------------------------
i <- iiw.weights(Surv(time.lag,time,event)~Wt + ApgarInd + I(conc.lag>0 & conc.lag<=20) + 
                I(conc.lag>20 & conc.lag<=30) + I(conc.lag>30)+ 
      cluster(Subject),id="Subject",time="time",event="event",data=data,
      invariant=c("Subject","Wt","ApgarInd"),lagvars=c("time","conc"), maxfu=16*24,lagfirst=c(0,0),first=FALSE)
i$m
i <- iiw.weights(Surv(time.lag,time,event)~Wt + I(conc.lag>0 & conc.lag<=20) + 
                I(conc.lag>20 & conc.lag<=30) + I(conc.lag>30)+ 
      cluster(Subject),id="Subject",time="time",event="event",data=data,
      invariant=c("Subject","Wt"),lagvars=c("time","conc"),maxfu=16*24,lagfirst=c(0,0),first=FALSE)
i$m
i <- iiw.weights(Surv(time.lag,time,event)~I(conc.lag>0 & conc.lag<=20) + 
                I(conc.lag>20 & conc.lag<=30) + I(conc.lag>30) + 
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

iiwgee <- iiwgee(conc ~ time + I(time^3) + log(time),Surv(time.lag,time,event)~I(conc.lag>0 & conc.lag<=20) + 
                I(conc.lag>20 & conc.lag<=30) + I(conc.lag>30) +cluster(id),
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
            maxfu=maxfu,baseline=baseline,Xfn=Xfn,Wfn=Wfn); print(x); return(x)
}
Xfn <- function(id,time){
  # Group is time invariant so just use the first value for each subject
  return(as.numeric(data$ApgarInd[data$id==id][1]))
}
Wfn <- function(id,time){
  return(c(1,time))
}
data$Intercept <- 1
data$time3 <- (data$time)^3
data$logtime <- log(data$time)
data$ApgarInd.time <- as.numeric(data$ApgarInd)*data$time/24
data$ApgarInd.time3 <- as.numeric(data$ApgarInd)*((data$time/24)^3)

set.seed(301031)
ifrailty <- iiw.weights(Surv(time.lag,time,event)~I(conc.lag>0 & conc.lag<=20) + 
                I(conc.lag>20 & conc.lag<=30) + I(conc.lag>30)   
                      +cluster(id),id="id",time="time",event="event",data=data,     
                      invariant=c("id"),lagvars=c("time","conc"),maxfu=16*24,
                      lagfirst=c(0,0), first=FALSE,frailty=TRUE)

wt <- ifrailty$iiw.weight

m.moLiang <- mo(20,Liangmo,data,wt, 
         singleobs=FALSE,id="id",time="time",keep.first=FALSE,var=FALSE,Yname="conc",
         Xnames=c("ApgarInd","ApgarInd.time","ApgarInd.time3"), 
         Wnames=c("Intercept"),maxfu=16*24,baseline=0)

m.moLiang$est

## -----------------------------------------------------------------------------
m.moLiang$est

