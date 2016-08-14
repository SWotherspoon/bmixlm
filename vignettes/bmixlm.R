## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo=TRUE,fig.width=6,fig.height=5)

## ------------------------------------------------------------------------
## Generate 7 uniformly distributed covariates
set.seed(31)
d <- as.data.frame(matrix(runif(1000*7),1000,7))
colnames(d) <- letters[seq_along(d)]

## ------------------------------------------------------------------------
## Generate responses from the two models
sigma <- c(0.4,0.6)
beta1 <- c(runif(1,-0.4,0.4),rnorm(3))
y1 <- model.matrix(~a+b+c,d)%*%beta1+rnorm(nrow(d),0,sigma[1])
beta2 <- c(runif(1,-0.4,0.4),rnorm(3))
y2 <- model.matrix(~c+d+e,d)%*%beta2+rnorm(nrow(d),0,sigma[2])

## ------------------------------------------------------------------------
## Draw the observed response from one model or the other
betap <- c(0,1,-1)
p <- pnorm(model.matrix(~f+g,d)%*%betap)
b <- rbinom(nrow(d),1,p)
d$y <- ifelse(b==0,y1,y2)

## ------------------------------------------------------------------------
## Show the two components and the mixture
library(ggplot2)
ggplot(data.frame(comp=factor(c(b+1,rep("mixture",nrow(d)))),
                  y=c(d$y,d$y)),
       aes(x=y,group=comp,colour=comp))+
  geom_density()

## ------------------------------------------------------------------------
## Draw a burnin sample
library(bmixlm)
fit <- bmixlm(y~a+b+c,y~c+d+e,~f+g,data=d,nsamp=100)

## ------------------------------------------------------------------------
## Traceplots
plot(fit,which="comp1")
plot(fit,which="comp2")
plot(fit,which="probit")
plot(fit,which="error")

## ------------------------------------------------------------------------
## Draw a larger sample
fit <- update(fit,nsamp=2000)
summary(fit)

## ------------------------------------------------------------------------
beta1
beta2
betap
sigma

## ------------------------------------------------------------------------
## Traceplots
plot(fit,which="comp1")
plot(fit,which="comp2")
plot(fit,which="probit")
plot(fit,which="error")

## ------------------------------------------------------------------------
## Fitted values, residuals, and probability of component membership
d.pr <- predictAll(fit)
head(d.pr)

## ------------------------------------------------------------------------
## Residuals vs fitted values
d.rf <- data.frame(comp=rep(1:2,each=nrow(d.pr)),
                   fitted=c(d.pr$y1,d.pr$y2),
                   residual=c(d.pr$r1,d.pr$r2),
                   prob=c(1-d.pr$q,d.pr$q))
library(ggplot2)
ggplot(d.rf[order(d.rf$prob),],
       aes(x=fitted,y=residual,colour=prob)) +
  geom_point(size=1)+
  facet_wrap(~comp,ncol=1)

## ------------------------------------------------------------------------
## Decompose response into components
cl <- classify(fit)
library(ggplot2)
ggplot(data.frame(comp=factor(c(ifelse(cl$q < 0.5,"1","2"),
                                rep("mixture",nrow(d)))),
                  y=rep(d$y,2)),
       aes(x=y,group=comp,colour=comp))+
  geom_density()

## ------------------------------------------------------------------------
## Simulate from posterior predictive distribution
ys <- simulate(fit,nsim=500)

## ------------------------------------------------------------------------
## Show observations against simulations
cl <- classify(fit)
d.pr <- as.data.frame(t(apply(ys,1,quantile,prob=c(0.025,0.5,0.975))))
d.pr <-cbind(d.pr,y=d$y,x=order(order(d.pr$`50%`)),cert=pmax(cl$q,1-cl$q))
library(ggplot2)
ggplot(d.pr,aes(x=x,y=y,ymin=`2.5%`,ymax=`97.5%`,colour=cert))+
  geom_ribbon(col="grey80",fill="grey80")+
  geom_point(size=1)

## ------------------------------------------------------------------------
## Fit a simple two component mixture
fit0 <- bmixlm(y~1,y~1,~1,data=d,nsamp=100)
fit0 <- update(fit0,nsamp=2000)
summary(fit0)

## ------------------------------------------------------------------------
fit1 <- update(fit,formula1=y~a+b,nsamp=100)
fit1 <- update(fit1,nsamp=2000)
summary(fit1)

## ------------------------------------------------------------------------
fit2 <- update(fit,formula1=y~a+b+c+d,nsamp=100)
fit2 <- update(fit2,nsamp=2000)
summary(fit2)

## ------------------------------------------------------------------------
fit3 <- update(fit,formula2=y~d+e,nsamp=100)
fit3 <- update(fit3,nsamp=2000)
summary(fit3)

## ------------------------------------------------------------------------
## Draw a burnin sample
fit0 <- bmixlm(y~a+b+c+e+f+g,y~b+c+d+e+f+g,~a+b+c+d+e+f+g,data=d,nsamp=100)
fit0 <- update(fit0,nsamp=2000)
summary(fit0)

## ------------------------------------------------------------------------
fit0 <- update(fit0,formula1=y~a+b+c,formula2=y~d+e,formulap=~f+g)
summary(fit0)

## ------------------------------------------------------------------------
q <- classify(fit0)$q
summary(step(lm(qnorm(q)~a+b+c+d+e+f+g,data=d),trace=FALSE))
summary(step(lm(y~a+b+c+d+e+f+g,data=d[q<0.4,]),trace=FALSE))
summary(step(lm(y~a+b+c+d+e+f+g,data=d[q>0.6,]),trace=FALSE))

