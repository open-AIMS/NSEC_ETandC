### Load required packages 
require(gslnls)  # Non-linear least-squares package. Used to obtain initial parameter estimates
require(ggplot2)
require(gnlm)
library(R2jags)
library(tidyverse)
library(investr)
library(ggpubr)

### Figure 1 R Code --------------
f1<-function(x,t,b,n){
  y<-t*exp(-b*(x-n)*(x>=n))
  return(y)
}
f2<-function(x,b0,b1,b2){
  y<-b0*exp(-b1*x^b2)
  return (y)
}

t1<-data.frame(x=seq(0,20,0.01),y=f1(seq(0,20,0.01),0.95,0.25,3.5),Equation=factor("1"))
t2<-data.frame(x=seq(0,20,0.01),y=f2(seq(0,20,0.01),0.95,0.005,2.5),Equation=factor("2"))

pdf("Figure1.pdf", width = 7, height = 5)
rbind(t1,t2) |> 
  dplyr::mutate(Curve=ifelse(Equation==1, "Threshold", "No threshold")) |> 
ggplot(aes(x=x,y=y,group=Curve))+ geom_line(aes(color=Curve)) +
  xlab("Concentration (x)") + ylab("Response (y)") +
  theme_classic() +
  theme(legend.title=element_blank(), legend.position = c(0.87, 0.7))
dev.off()

### Figure 2 --------- 
#Note drawn in Powerpoint

### Figure 3 -----------
dat <- read.csv("Table_1.csv")

x<-dat$concentration
n<-dat$n
y<-dat$growth
nsec.dat<-data.frame(x=x,y=y)
k<-nrow(nsec.dat) 

fit <- nls(y ~ b0*exp(-b1*x^b2), data= nsec.dat,start = list(b0=6,b1=0.005,b2=1) )

summary(fit)
params <- summary(fit)$parameters
params

# lower 95th quantile of b0
lwrb0 <- params[1, "Estimate"]-qt(0.95,18)*.31
ols_nsec <- (-log(lwrb0/params[1, "Estimate"])/params[2, "Estimate"])^(1/params[3, "Estimate"])

new.data <- data.frame(x=seq(0, 60, by = 0.01))
p.interval <- as_tibble(predFit(fit, newdata = new.data, interval = "prediction", level= 0.9)) %>% 
  mutate(x = new.data$x)
c.interval <- as_tibble(predFit(fit, newdata = new.data, interval = "confidence", level= 0.9)) %>% 
  mutate(x = new.data$x)

p1 <- ggplot(nsec.dat) +  
  geom_point(aes(x=x, y=y),size=2, 
  colour="black") + xlab("Concentration (x)") + ylab("Response (y)") 

pdf("Figure3.pdf", width = 7, height = 5)
p1 +
  geom_line(data=p.interval, aes(x = x, y = fit ),color="red")+ 
  scale_x_continuous(expand=c(0.01,0),breaks=c(seq(0,60,by=5)))+
  scale_y_continuous(expand=c(0,0))+
  geom_line(data=c.interval,aes(x = x, y = upr ),color="#00BFFF") +
  geom_line(data=c.interval,aes(x = x, y = lwr ),color="#00BFFF") +
  geom_ribbon(data=p.interval, aes(x=x, ymin=lwr, ymax=upr), alpha=0.5, inherit.aes=F, fill="#B0C4DE")+
  geom_segment(aes(x=ols_nsec,y=-2,xend=ols_nsec,yend=lwrb0),inherit.aes = TRUE,linetype=2) + 
  geom_segment(aes(x=0,y=lwrb0,xend=ols_nsec,yend=lwrb0),inherit.aes = TRUE,linetype=2) +
  annotate("text", x = ols_nsec*1.1, y = -1.7, label = "5") +
  theme_classic()
dev.off()

### Figure 4a R Code -------------
dat <- read.csv("Table_2.csv")

x<-dat$concentration
n<-dat$n
y<-dat$response
nsec.dat<-data.frame(x=x,n=n,y=y)
k<-nrow(nsec.dat) 

# Use gsl_nls to obtain vector of initial parameter estimates:
B.cur<-coefficients(gsl_nls(
  fn=y~b0*exp(-b1*x^b2),
  data=data.frame(x=nsec.dat$x,y=nsec.dat$y/nsec.dat$n),
  start=c(b0=0.9,b1=0.005,b2=3)
))
y<-cbind(nsec.dat$y,nsec.dat$n-nsec.dat$y)

#  Equation (2)
mu <- function(B) {
  pi <- B[1] * exp(-B[2] * (nsec.dat$x)^B[3])
  return(pi)
}

pred<-function(x,B) {
  pi <- B[1] * exp(-B[2] * x^B[3])
  return(pi)
}
  
nsec.glim<-gnlr(y, dist = "binomial", mu = mu, pmu = c(0.9, 0.005, 3))
B.new<-nsec.glim$coefficients
H<-nsec.glim$cov


cat("\n","Parameter estimates & standard errors:",
    "\n",paste("b0=",signif(B.cur[1],4)," (",signif(sqrt(H[1,1]),4),")",sep=""),
    "\n",paste("b1=",signif(B.cur[2],4)," (",signif(sqrt(H[2,2]),4),")",sep=""),
    "\n",paste("b2=",signif(B.cur[3],4)," (",signif(sqrt(H[3,3]),4),")",sep=""))
p_1<-B.new[1]+qt(0.01,(k-3))*sqrt(H[1,1])
p_5<-B.new[1]+qt(0.05,(k-3))*sqrt(H[1,1])
p_10<-B.new[1]+qt(0.10,(k-3))*sqrt(H[1,1])
p_20<-B.new[1]+qt(0.20,(k-3))*sqrt(H[1,1])
SEC_99<-((1/B.new[2])*log(B.new[1]/p_1))^(1/B.new[3])
SEC_95<-((1/B.new[2])*log(B.new[1]/p_5))^(1/B.new[3])
SEC_90<-((1/B.new[2])*log(B.new[1]/p_10))^(1/B.new[3])
SEC_80<-((1/B.new[2])*log(B.new[1]/p_20))^(1/B.new[3])
cat("\n",paste("SEC_99>",signif(SEC_99,3)))
cat(paste("SEC_95>",signif(SEC_95,3)))
cat(paste("SEC_90>",signif(SEC_90,3)))
cat(paste("SEC_80>",signif(SEC_80,3)))
t2<-data.frame(x=seq(0,35,0.01),y=pred(seq(0,35,0.01),B.new))
   
pmle<-ggplot(data=nsec.dat,aes(x=x,y=y/n))+ geom_point()+  
  geom_line(data=data.frame(x=seq(0,35,0.01),
                            y=pred(x=seq(0,35,0.01),B.new)),
            aes(x=x,y=y),color="blue") +
  # NSEC 99
  geom_segment(aes(x=SEC_99,y=0,xend=SEC_99,yend=p_1),inherit.aes = TRUE,linetype=2) + 
  geom_segment(aes(x=0,y=p_1,xend=SEC_99,yend=p_1),inherit.aes = TRUE,linetype=2) +
  
  #NSEC 95
  geom_segment(aes(x=SEC_95,y=0,xend=SEC_95,yend=p_5),inherit.aes = TRUE,linetype=2) + 
  geom_segment(aes(x=0,y=p_5,xend=SEC_95,yend=p_5),inherit.aes = TRUE,linetype=2) +

  #NSEC 90
  geom_segment(aes(x=SEC_90,y=0,xend=SEC_90,yend=p_10),inherit.aes = TRUE,linetype=2) + 
  geom_segment(aes(x=0,y=p_10,xend=SEC_90,yend=p_10),inherit.aes = TRUE,linetype=2) +

  #NSEC 80
  geom_segment(aes(x=SEC_80,y=0,xend=SEC_80,yend=p_20),inherit.aes = TRUE,linetype=2) + 
  geom_segment(aes(x=0,y=p_20,xend=SEC_80,yend=p_20),inherit.aes = TRUE,linetype=2) +
  
  annotate("text", x = c(SEC_99, SEC_95, SEC_90, SEC_80), 
           y = -0.02, label = c("1", "5", "10", "20")) +
  
  # ols result
  geom_line(data=data.frame(x=seq(0,35,0.001),y=pred(x=seq(0,35,0.001),B.cur)),
            aes(x=x,y=y),color="red") +
  theme_classic() +
  xlab("Concentration (x)") + ylab("Response (y)") 
pmle


### Figure 4b R Code ----------------------
dat <- read.csv("Table_2.csv")
model_file <- "jags_model.txt"

# create jags model data list
mod.dat <<- list(
  x = dat$concentration, # concentration
  y = dat$response, # response (successes)
  N = nrow(dat), # Sample size
  trials = dat$n # binomial trials
)
response <- dat$response/dat$n

params <- c("b0", "b1", "b2")
set.seed(30)
fit <- R2jags::jags(data = mod.dat, parameters = params, model = model_file)

# get the predicted parameter values
parameter_posteriors <- do.call("cbind", fit$BUGSoutput$sims.list[params])
colnames(parameter_posteriors) <- params
head(parameter_posteriors)

apply(parameter_posteriors, MARGIN = 2, FUN = quantile, probs = 0.5)


x.seq <- seq(0, 35, by=0.01)
y.pred <- apply(parameter_posteriors, MARGIN = 1, FUN = function(r){
  r[1] * exp(-r[2] * (x.seq)^r[3])
})

# predicted values for plotting
m.vals <- apply(y.pred, MARGIN = 1, FUN = quantile, probs = 0.5)
up.vals <- apply(y.pred, MARGIN = 1, FUN = quantile, probs = 0.975)
lw.vals <- apply(y.pred, MARGIN = 1, FUN = quantile, probs = 0.025)

# calculate NSEC
control_quantiles <- quantile(parameter_posteriors[, "b0"], 
                              probs = c(0.01, 0.05, 0.10, 0.20))

# find NSEC using vector of x to predict y
NSEC_vals2 <- apply(y.pred, MARGIN = 2, FUN = function(r){
  sapply(control_quantiles, FUN = function(p){
    x.seq[which.min(abs(r-p))]
  })
})

# find NSEC mathematically
NSEC_vals <- apply(parameter_posteriors, MARGIN = 1, FUN = function(r){
  sapply(control_quantiles, FUN = function(p){
    x <- (-log(p/r[1])/r[2])^(1/r[3])
    #x <- ((1/r[2])*log(r[1]/p))^(1/r[3])
  })
})
NSEC_vals[which(is.na(NSEC_vals))] <- 0

# summarise the posterior estimates
NSECsummary <- apply(NSEC_vals, MARGIN = 1, FUN = quantile, probs = c(0.5, 0.025, 0.975))
NSECsummary

# plot the results
pjags<-ggplot(data=nsec.dat,aes(x=x,y=y/n))+ geom_point()+  
  geom_line(data=data.frame(x=x.seq, y=m.vals),
            aes(x=x,y=y),color="blue") +
  # NSEC 99
  geom_segment(aes(x=NSECsummary[1,1],y=0,xend=NSECsummary[1,1],
                   yend=control_quantiles[1]),inherit.aes = TRUE,linetype=2) + 
  geom_segment(aes(x=0,y=control_quantiles[1],xend=NSECsummary[1,1],
                   yend=control_quantiles[1]),inherit.aes = TRUE,linetype=2) +
  
  # #NSEC 95
  geom_segment(aes(x=NSECsummary[1,2],y=0,xend=NSECsummary[1,2],
                   yend=control_quantiles[2]),inherit.aes = TRUE,linetype=2) + 
  geom_segment(aes(x=0,y=control_quantiles[2],xend=NSECsummary[1,2],
                   yend=control_quantiles[2]),inherit.aes = TRUE,linetype=2) +

  #NSEC 90
  geom_segment(aes(x=NSECsummary[1,3],y=0,xend=NSECsummary[1,3],
                   yend=control_quantiles[3]),inherit.aes = TRUE,linetype=2) +
  geom_segment(aes(x=0,y=control_quantiles[3],xend=NSECsummary[1,3],
                   yend=control_quantiles[3]),inherit.aes = TRUE,linetype=2) +

  #NSEC 80
  geom_segment(aes(x=NSECsummary[1,4],y=0,xend=NSECsummary[1,4],
                   yend=control_quantiles[4]),inherit.aes = TRUE,linetype=2) +
  geom_segment(aes(x=0,y=control_quantiles[4],xend=NSECsummary[1,4],
                   yend=control_quantiles[4]),inherit.aes = TRUE,linetype=2) +
  
  annotate("text", x = NSECsummary[1,], y = -0.02, label = c("1", "5", "10", "20")) +
  
  geom_line(data=data.frame(x=x.seq,y=up.vals),
            aes(x=x,y=y),color="#00BFFF") +
  geom_line(data=data.frame(x=x.seq,y=lw.vals),
            aes(x=x,y=y),color="#00BFFF") +
  theme_classic() +
  xlab("Concentration (x)") + ylab("Response (y)") 


pjags

### Figure 4 single panel plot ----
pdf("Figure4.pdf", width = 7, height = 10)
ggarrange(pmle, pjags, 
          labels = c("(A)", "(B)"),
          ncol = 1, nrow = 2)
dev.off()































