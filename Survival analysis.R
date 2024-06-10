remove(list=ls());cat("\14");
setwd('C:/Users/thinl/OneDrive/Documents')
       
# install.packages("TH.data") # This data is on german breast cancer study group 2: 686 observations
# install.packages("Ecdat") # data on the unemploment duration "number of observations : 3343
# install.packages("survival")
# install.packages("survminer")
# install.packages("reshape2")

#Survival package
library("survival")
library("survminer")
library("reshape2")

# Moreinformation on data 
# help(UnempDur, package="Ecdat")
# help(GBSG2, package="TH.data")

#Load the data 
data(UnempDur, package="Ecdat")
data(GBSG2, package="TH.data")
dta1<-GBSG2
dta2<-UnempDur


#Summary statistics
summary(GBSG2)
summary(UnempDur)

#GBSG2 

#count the censored and uncensored observation 
num_censored1<-table(dta1$cens)
num_censored1 # cens=1 means person died during the event


 #barplot
barplot(num_censored1)

# Create Surv-Object
sobj <- Surv(GBSG2$time,GBSG2$cens)
sobj

#o/p eg:1814  2018   712  1807   772   448  2172+ 2161+ (+ means dataset is censored is 0 the event has not eventuated)

# Look at 10 first elements
sobj[1:11]


# Look at summary
summary(sobj)

# Look at structure
str(sobj)

#UnempDur
# Count censored and uncensored data
cens_employ_ft <- table(UnempDur$censor1)
cens_employ_ft

# Create barplot of censored and uncensored data
barplot(cens_employ_ft)

# Create Surv-Object
sobj2 <- Surv(UnempDur$spell, UnempDur$censor1) # Spell indicate the the length of time  individual was unemployed in 2 week interval
sobj2
# Look at 10 first elements
sobj[1:10]

## Kaplan-Meier estimate GBSG2
km <- survfit(Surv(time,cens) ~1, data = GBSG2)

# Plot of the Kaplan-Meier estimate
ggsurvplot(km)

# Add the risk table to plot
ggsurvplot(km, risk.table = TRUE)

# Add a line showing the median survival time
ggsurvplot(km,risk.table = TRUE,surv.median.line = "hv")


#weibull model GBSG2
Wb<- survreg(Surv(time,cens) ~1, data = GBSG2)

#to compute a time point were 90% patient survive we use:

predict(Wb, type="quantile", p= 1-0.9, newdata=data.frame(1))
#90% of patients survived 384 days

sur<-seq(0.99, 0.01, by= -.01)
t<-predict(Wb, type="quantile", p= 1-sur, newdata=data.frame(1))

wef<- data.frame(time=t, surv=sur, upper=NA, lower=NA, std.err=NA)
wef

ggsurvplot_df(fit=wef, surv.geom=geom_line)

#WEIBULL with covariates
# Weibull model
wbmod <- survreg(Surv(time, cens) ~ horTh + tsize, data = GBSG2)

# Imaginary patients
newdat <- expand.grid(
  horTh = levels(GBSG2$horTh),
  tsize = quantile(GBSG2$tsize, probs = c(0.25, 0.50, 0.75)))
newdat

surv<-seq(0.99, 0.01, by= -.01)
t2<-predict(wbmod, type="quantile", p= 1-surv, newdata=newdat)

surv_wbmod_wide <- cbind(newdat, t2)

# Use melt() to bring the data.frame to long format
surv_wbmod <- melt(surv_wbmod_wide, id.vars = c("horTh", "tsize"), variable.name = "surv_id", value.name = "time")

# Use surv_wbmod$surv_id to add the correct survival probabilities surv
surv_wbmod$surv <- surv[as.numeric(surv_wbmod$surv_id)]

# Add columns upper, lower, std.err, and strata to the data.frame
surv_wbmod[, c("upper", "lower", "std.err", "strata")] <- NA

# Plot the survival curves
ggsurvplot_df(surv_wbmod, surv.geom = geom_line,
              linetype = "horTh", color = "tsize", legend.title = NULL)
surv_wbmod_wide<- cbind(newdat,t2)

remove(list=ls());cat("\14");
#alternative distribution option

#exponential 
expmode <- survreg(Surv(time, cens) ~ horTh + tsize, data = GBSG2, dist= "exponential")

#lognormal
lognormalmode <- survreg(Surv(time, cens) ~ horTh + tsize, data = GBSG2, dist= "lognormal")


#------Example---------#
# Weibull model
wbmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2)

# Log-Normal model
lnmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2, dist = "lognormal")

# Newdata
newdat <- data.frame(horTh = levels(GBSG2$horTh))

# Surv
surv <- seq(.99, .01, by = -.01)

# Survival curve from Weibull model and log-normal model
wbt <- predict(wbmod, type = "quantile", p = 1-surv, newdata = newdat)
lnt <- predict(lnmod, type = "quantile", p = 1-surv, newdata =newdat)

surv_w<- rbind(wbt,lnt)
surv_wide<-cbind(newdat,surv_w)
surv_wide$dist <- rep(c("Weibull", "Log-Normal"), each = nrow(newdat))

# Melt the data.frame into long format.
surv_long <- melt(surv_wide, id.vars = c("horTh", "dist"), variable.name= "surv_id", value.name = "time")

# Add column for the survival probabilities
surv_long$surv <- surv[as.numeric(surv_long$surv_id)]

# Add columns upper, lower, std.err, and strata contianing NA values
surv_long[, c("upper", "lower", "std.err", "strata")] <- NA


ggsurvplot_df(surv_long, surv.geom = geom_line,
              linetype = "horTh", color = "dist", legend.title = NULL)



#------COX model-----------#
cxmod <- coxph(Surv(time, cens) ~ horTh+tsize, data = GBSG2)
coef(cxmod)
#horThyes is -0.3640099 (-ve values means positive  effect )
#covariates
#
newdat <- expand.grid(
  horTh = levels(GBSG2$horTh),
  tsize = quantile(GBSG2$tsize, probs = c(0.25, 0.50, 0.75)))
newdat

cxsf<-survfit(cxmod, data=GBSG2, newdata=newdat, conf.type="none")
surv_cxmod0<-surv_summary(cxsf)
surv_cxmod<- cbind(surv_cxmod0,newdat[as.character(surv_cxmod0$strata),])

ggsurvplot_df(surv_cxmod, linetype = 'horTh', color="tsize", legend.title=NULL,censor=FALSE)


remove(list=ls());cat("\14");
library("survival")
library("survminer")

#survival curves 
# Create time and event data
example<- data.frame(
time = c(5, 6, 2, 4, 4),
event =c(1, 0, 0, 1, 1))

# Compute Kaplan-Meier estimate
km <- survfit(Surv(time,event) ~ 1, data=example)
km

# Take a look at the structure
str(km)

# Create data.frame
data.frame(time = km$time, n.risk = km$n.risk, n.event = km$n.event,
           n.censor = km$n.censor, surv = km$surv)


#graph
ggsurvplot(fit=km, 
           palette = "blue",
           linetye=1,
           surv.median.line="hv",
           risk.table=TRUE,
           cumevents=TRUE,
           cumcensor=TRUE,
           tables.height=0.1)
