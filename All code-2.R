#set up weight first 
nyQ1.gal <- nb2listw(poly2nb(ny2)) # Notice we don't need to use zero.policy=T
summary(nyQ1.gal)
### Create and attach new ny2.df as the default data set.

ny2.df <- data.frame(ny2)
attach(ny2.df)
class(ny2.df)   # Object class is data.frame, which we need
names(ny2.df)
row.names(ny2.df) <- ny2.df$POSTAL

#logit transformation
lrate <- logit(Rate1000)
#OLS model
reg1 <- lm(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH)
summary(reg1)
# Run the Moran & Geary tests on the residuals:

moran.test(reg1$residuals, nyQ1.gal, alternative="two.sided")
geary(reg1$residuals, nyQ1.gal, n=length(lrate), n1=length(lrate)-1,
      S0=length(lrate))
## GeoDa-style Lagrange multiplier tests can be obtained using the

lm.LMtests(reg1, nyQ1.gal, test=c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA"))
#if you want speed up the procese 
W.eig <- eigenw(nyQ1.gal, quiet=NULL)
1/min(W.eig)   
1/max(W.eig)
# FIRST, THE SPATIAL LAG MODEL:
ny.lag.eig <- lagsarlm(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH,
                       data=ny2, nyQ1.gal, method="eigen", quiet=FALSE)
summary(ny.lag.eig)
#check assymerty of the residual 
me1 <- mean(residuals(ny.lag.eig))
me1    
sd1 <- sd(residuals(ny.lag.eig))
sd1   
summary(residuals(ny.lag.eig))  # reasonably symmetric; good

hist(residuals(ny.lag.eig), breaks=seq(-3, 3, .2),
     col=8, probability=T,
     ylab='Density',
     main='Histogram of Residuals(Spatial Lag Model)',
     xlab='Residuals(Spatial Lag Model)')
box()
curve(dnorm(x, mean=me1, sd=sd1), from=-3, to=3, add=T,
      col='red', lwd=2)
# Check residuals for spatial autocorrelation using Q1 matrix:
names(ny.lag.eig)
moran.test(ny.lag.eig$resid, nyQ1.gal, alternative="two.sided")
#raw plot about where there still have spatial pattern
breaks <-seq(-3,3,1)
colbrks <- findInterval(ny.lag.eig$resid, breaks, all.inside=FALSE)

plot(ny2)
plot(ny2, col=colbrks, add=T)
##  THE SPATIAL ERROR MODEL:
ny.err.eig <- errorsarlm(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH,
                         data=ny2, nyQ1.gal, method="eigen", quiet=FALSE)
summary(ny.err.eig)
# As we did with the spatial lag model, let's look at the residuals:

me2 <- mean(residuals(ny.err.eig))
me2    # Effectively zero
sd2 <- sd(residuals(ny.err.eig))
sd2    # 0.582
summary(residuals(ny.err.eig))
hist(residuals(ny.err.eig), breaks=seq(-3, 3, .2),
     col=8, probability=T,
     ylab='Density', main='Histogram of Residuals(Spatial Error Model)',
     xlab='Residuals(Spatial Error Model)')
box()
curve(dnorm(x, mean=me2, sd=sd2), from=-3, to=3, add=T,
      col='red', lwd=2)
# Check residuals for spatial autocorrelation using Q1 matrix:

names(ny.err.eig)
moran.test(ny.err.eig$residuals, nyQ1.gal, alternative="two.sided")
##### NOW THE SPATIAL ERROR REGRESSION MODEL USING THE K-P METHOD
nyKP <- GMerrorsar(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH,
                   data=ny2, nyQ1.gal, verbose=FALSE)
summary(nyKP)
# Can also try estimating the spatial lag model using a
# two-state-least-squares estimator.  This is not an MLE estimator,
# so, as with the GMM estimator, we don't get a likelihood score
# or AIC.  We should anticipate different results than from the
# spatial lag model estimated with MLE.  It too runs very quickly.

ny2SLS <- stsls(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH,
                data=ny2, nyQ1.gal)
summary(ny2SLS)

# There's even an argument in this function call that gives us
# robust estimates that are adjusted for heteroskedasticity.  That is,
# there is a heteroskedasticity correction to the coefficient
# covariances:

ny2SLSR <- stsls(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH,
                 robust=T, data=ny2, nyQ1.gal)
summary(ny2SLSR)
# Let's compare the autocorrelation in the residuals

moran.test(reg1$residuals, nyQ1.gal, alternative="two.sided")
moran.test(ny.lag.eig$residuals, nyQ1.gal, alternative="two.sided")
moran.test(ny.err.eig$residuals, nyQ1.gal, alternative="two.sided")
moran.test(nyKP$residuals, nyQ1.gal, alternative="two.sided")
moran.test(ny2SLS$residuals, nyQ1.gal, alternative="two.sided")
moran.test(ny2SLSR$residuals, nyQ1.gal, alternative="two.sided")
# Let's also map out the residuals across the different models

# First, set up the plot window
windows()
par(mfrow=c(1,2))
par(mar=c(1.2,1.2,1.2,1.2))

# Choose a color palette - we will use this for all the maps
colors <- brewer.pal(5, "YlOrBr")  # Stores our colors in object color

# Now classify and map the Lag model
windows()
color.cat.lag<-classIntervals(ny.lag.eig$residuals, n=5, style="quantile", dataPrecision=2)
collag <- findColours(color.cat.lag, colors)
plot(ny2, col=collag)
title('Map of Lag Model Residuals')
legend('topleft', legend=c(names(attr(collag, 'table'))), fill=c(attr(collag, 'palette')), 
       title='Regression Residuals')
box()

# Plot Error model
color.cat.err<-classIntervals(ny.err.eig$residuals, n=5, style="quantile", dataPrecision=2)
colerr <- findColours(color.cat.err, colors)
plot(ny2, col=colerr)
title('Map of Error Model Residuals')
legend('topleft', legend=c(names(attr(colerr, 'table'))), fill=c(attr(colerr, 'palette')), 
       title='Regression Residuals')
box()
#how to interpret lag model 
effects.SAR <- impacts(ny.lag.eig, listw=nyQ1.gal)
effects.SAR


############################
#
##### HIGHER ORDER SPATIAL MODELS
#
############################

##### SARAR MODEL

ny.lag.error <- sacsarlm(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH,
                         data=ny2, nyQ1.gal, method = "eigen", quiet=FALSE)
summary(ny.lag.error)

# As we did with the spatial error regression model and the
# spatial lag regression model, let's look at the residuals:

me3 <- mean(residuals(ny.lag.error))
me3    # Effectively zero
sd3 <- sd(residuals(ny.lag.error))
sd3    # 0.5284
summary(residuals(ny.lag.error))
hist(residuals(ny.lag.error), breaks=seq(-3, 3, .2),
     col=8, probability=T,
     ylab='Density', main='Histogram of Residuals(SARAR Model)',
     xlab='Residuals(SARAR Model)')
box()
curve(dnorm(x, mean=me3, sd=sd3), from=-3, to=3, add=T,
      col='red', lwd=2)

plot(ny.lag.error$resid); abline(h=0)

# Check residuals for spatial autocorrelation using Q1 matrix:

moran.test(ny.lag.error$residuals, nyQ1.gal,alternative="two.sided")
##### SPATIAL DURBIN MODEL
#
# This is an important model.  Fresh attention has been brought to it
# in the LeSage & Pace (2009) book.  It will be fit using a the spatial
# lag model (as in the first of the estimated models, above), but it 
# includes spatial lags of the independent variables.
#

# LAG THE INDEPENDENT VARIABLES

lwhite <- lag.listw(nyQ1.gal, PctWht)
lhisp <- lag.listw(nyQ1.gal, PctHisp)
lgini <- lag.listw(nyQ1.gal, Gini)
lhsed <- lag.listw(nyQ1.gal, PctHSEd)
lfemhh <- lag.listw(nyQ1.gal, PctFemHH)

# Add the lagged variables to ny.df:

detach(ny2.df)
lag.cols <- data.frame(lwhite, lhisp, lgini, lhsed, lfemhh)
nylag.df <- cbind(ny2.df, lag.cols)
head(nylag.df)
attach(nylag.df)

# Run the Spatial Durbin Model

ny.SDM <- lagsarlm(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH +
                     lwhite + lhisp + lgini,
                   data=nylag.df, nyQ1.gal, method="eigen", quiet=FALSE)
summary(ny.SDM)

# An alternative way to run this is using 'type="mixed", but that
# automatically includes all the independent variables.

ny.SDM2 <- lagsarlm(lrate ~ PctWht + PctHisp + Gini + PctHSEd + PctFemHH ,
                    data=nylag.df, nyQ1.gal, method="eigen", type="mixed",
                    quiet=FALSE)
summary(ny.SDM2) # it gives slightly different results.
#check residual 
me4 <- mean(residuals(ny.SDM))
me4    # Effectively zero
sd4 <- sd(residuals(ny.SDM))
sd4    # 0.566
summary(residuals(ny.SDM))  # reasonably symmetric; good

hist(residuals(ny.SDM), breaks=seq(-3.2, 3.2, .2),
     col=8, probability=T,
     ylab='Density',
     main='Histogram of Residuals(Spatial Durbin Model)',
     xlab='Residuals(Spatial Durbin Model)')
box()
curve(dnorm(x, mean=me4, sd=sd4), from=-3.2, to=3.2, add=T,
      col='red', lwd=2)
# Check residuals for spatial autocorrelation using Q1 matrix:
names(ny.SDM)
moran.test(ny.SDM$residuals, nyQ1.gal, alternative="two.sided")
# Let's also map out the residuals across the different models

# First, set up the plot window
windows()
par(mfrow=c(2,2))
par(mar=c(1.2,1.2,1.2,1.2))

# Choose a color palette - we will use this for all the maps
colors <- brewer.pal(5, "YlOrBr")  # Stores our colors in object color

# Now classify and map the Lag model
color.cat.lag<-classIntervals(ny.lag.eig$residuals, n=5, style="quantile", dataPrecision=2)
collag <- findColours(color.cat.lag, colors)
plot(ny2, col=collag)
title('Map of Lag Model Residuals')
legend('topleft', legend=c(names(attr(collag, 'table'))), fill=c(attr(collag, 'palette')), 
       title='Regression Residuals')
box()

# Plot Error model
color.cat.err<-classIntervals(ny.err.eig$residuals, n=5, style="quantile", dataPrecision=2)
colerr <- findColours(color.cat.err, colors)
plot(ny2, col=colerr)
title('Map of Error Model Residuals')
legend('topleft', legend=c(names(attr(colerr, 'table'))), fill=c(attr(colerr, 'palette')), 
       title='Regression Residuals')
box()

# Plot SARAR model
color.cat.sarar<-classIntervals(ny.lag.error$residuals, n=5, style="quantile", dataPrecision=2)
colsarar <- findColours(color.cat.sarar, colors)
plot(ny2, col=colsarar)
title('Map of SARAR Model Residuals')
legend('topleft', legend=c(names(attr(colsarar, 'table'))), fill=c(attr(colsarar, 'palette')), 
       title='Regression Residuals')
box()

# Plot Durbin model
color.cat.sdm<-classIntervals(ny.SDM$residuals, n=5, style="quantile", dataPrecision=2)
colsdm <- findColours(color.cat.sdm, colors)
plot(ny2, col=colsdm)
title('Map of SDM Model Residuals')
legend('topleft', legend=c(names(attr(colsdm, 'table'))), fill=c(attr(colsdm, 'palette')), 
       title='Regression Residuals')
box()




#attach dataset first 
attach(ny1)
2*sqrt(length(CASES))         # 33.52

hist(CASES, nclass=33, border=2) # Looks left skewed, but really we have 3 outliers

mean(CASES); median(CASES)

#QQ plot
qqPlot(CASES, distribution="norm",
       xlab='', main='Quantile Comparison Plot CASES',
       envelope=.95, las=0, pch=NA, lwd=2, col="red",
       line="quartiles")
par(new=TRUE)
qqPlot(CASES, distribution="norm", envelope=FALSE, 
       pch=1, cex=1, col="black")
par(new=FALSE)
#QQ plot for OLS residural 
qqPlot(residuals(lm1), distribution="norm",
       xlab='', main='Quantile Comparison Plot reg1 esiduals',
       envelope=.95, las=0, pch=NA, lwd=2, col="red",
       line="quartiles")
par(new=TRUE)
qqPlot(residuals(lm1), distribution="norm", envelope=FALSE,
       pch=1, cex=1, col="black")
par(new=FALSE)
#heteroskedasticity test 
plot(fitted(lm1), residuals(lm1), xlab="Fitted y", ylab= "Residuals",
     main="Plot of Residuals against Fitted y")
abline(h=0)
#moran's test 
ny.q<-nb2listw(poly2nb(ny))
moran.test(CASES, ny.q)   
moran.test(lm1$resid, ny.q) 

# Let's look at the Moran's Plot
par(mfrow=c(1,2))
moran.plot(CASES, ny.q, xlab=NULL, labels=F, ylab=NULL, pch=1)
moran.plot(lm1$resid, ny.q, xlab=NULL, labels=F, ylab=NULL, pch=1)
#####################################################################
###
############ GEOGRAPHICALLY WEIGHTED REGRESSION MODELS
###
#####################################################################

# The first step in GWR is to comput the optimal bandwidth.  There are 
# two methods for this: CV and AICc minimization.  We also have the option
# to choose an adaptive bandwidth or a fixed one.  Finally, we can specify
# the mathematical form we'd like out kernel to take.

# The options for kernels in GWmodel package are: gaussian, exponential,
# bisquare, tricube, boxcar.  Most researchers agree it doesn't matter too
# much, but the most commonly used are gaussian, exponential and bisquare.

# In the code below, let's see how specifying different kernels changes the
# bandwidth estimate.

# First, we create two with fixed kernels
bw.fbi<-bw.gwr(CASES~PCTAGE65P+PCTOWNHOME+PEXPOSURE, 
               data=ny.spdf, approach="CV", kernel="bisquare", adaptive=F)
print(bw.fbi) # 167164

bw.fgau<-bw.gwr(CASES~PCTAGE65P+PCTOWNHOME+PEXPOSURE, 
                data=ny, approach="CV", kernel="gaussian", adaptive=F)
print(bw.fgau) # 167128 really not that different from the bisquare

# Let's run the first one again, but using AICc as the selection
bw.fbi.aic<-bw.gwr(CASES~PCTAGE65P+PCTOWNHOME+PEXPOSURE, 
                   data=ny, approach="AICc", kernel="bisquare", adaptive=F)
print(bw.fbi.aic) # 167171 again, not that different from the first one

# Here we create two with adaptive kernels
bw.abi<-bw.gwr(CASES~PCTAGE65P+PCTOWNHOME+PEXPOSURE, 
               data=ny, approach="CV", kernel="bisquare", adaptive=T)
print(bw.abi) # 126 this is the number of neighbors to include, not the distance

bw.agau<-bw.gwr(CASES~PCTAGE65P+PCTOWNHOME+PEXPOSURE, 
                data=ny, approach="CV", kernel="gaussian", adaptive=T)
print(bw.agau) # 279 this is quite different from the bisquare
#####
##### FIRST, JUST THE BASIC MODEL:
#####

gwr.fbi<-gwr.basic(CASES~PCTAGE65P+PCTOWNHOME+PEXPOSURE, 
                   data=ny, bw=bw.fbi, kernel="bisquare", adaptive=F, F123.test=T)
print(gwr.fbi)
gwr.abi<-gwr.basic(CASES~PCTAGE65P+PCTOWNHOME+PEXPOSURE, 
                   data=ny, bw=bw.abi, kernel="bisquare", adaptive=T, F123.test=T)
gwr.t.adj <- gwr.t.adjust(gwr.abi)
print(gwr.abi)
#MC test
mctest<-gwr.montecarlo(CASES~PCTAGE65P+PCTOWNHOME+PEXPOSURE, 
                       data=ny, nsims=999, kernel="bisquare",adaptive=T, bw.abi)
mctest
# Let's  map the parameter values and the t-values side by side
names(gwr.abi$SDF)

mycol.1<-brewer.pal(6, "Greens")
mycol.2<-brewer.pal(4, "RdBu")

plot.exp<-spplot(gwr.abi$SDF, "PEXPOSURE", col.regions=mycol.1, cuts=5, main="Parameter Values for Exposure")
plot.expt<-spplot(gwr.t.adj$SDF, "PEXPOSURE_p", col.regions=mycol.2, cuts=5, 
                  at = c(0, 0.025, 0.05, 0.1, 1.0000001), main="p-values for Exposure")
windows()
plot(plot.exp, split=c(1,1,2,1), more=T)
plot(plot.expt, split=c(2,1,2,1))

plot.hown<-spplot(gwr.abi$SDF, "PCTOWNHOME", col.regions=mycol.1, cuts=5, main="Parameter Values for Homeownership")
plot.hownt<-spplot(gwr.t.adj$SDF, "PCTOWNHOME_p", col.regions=mycol.2, cuts=5, 
                   at = c(0, 0.025, 0.05, 0.1, 1.0000001), main="p-values for Homeownership")

plot(plot.hown, split=c(1,1,2,1), more=T)
plot(plot.hownt, split=c(2,1,2,1))

plot.age<-spplot(gwr.abi$SDF, "PCTAGE65P", col.regions=mycol.1, cuts=5, main="Parameter Values for Prop 65+")
plot.aget<-spplot(gwr.t.adj$SDF, "PCTAGE65P_p", col.regions=mycol.2, cuts=5, 
                  at = c(0, 0.025, 0.05, 0.1, 1.0000001), main="p-values for Prop 65+")

plot(plot.age, split=c(1,1,2,1), more=T)
plot(plot.aget, split=c(2,1,2,1))
#####
##### NOW LET'S LOOK FOR LOCAL COLLINEARITIES:
#####

gwr.coll<-gwr.collin.diagno(CASES~PCTAGE65P+PCTOWNHOME+PEXPOSURE, 
                            data=ny, bw=bw.abi, kernel="bisquare", adaptive=T)
names(gwr.coll$SDF)

# Once we've run the code, we can look at the VIF and local correlations.  
windows()
par(mfrow=c(3,1))
hist(gwr.coll$SDF$PCTAGE65P_VIF, freq=T, main="VIF Prop. Age 65+", xlab="VIF", xlim=c(1,2), ylim=c(0,200))
abline(v=c(5,10), col=c("red","darkred"),lty=3, lwd=5)
hist(gwr.coll$SDF$PCTOWNHOME_VIF, freq=T, main="VIF Prop. Homeownership", xlab="VIF", xlim=c(1,2), ylim=c(0,200))
abline(v=c(5,10), col=c("red","darkred"),lty=3, lwd=5)
hist(gwr.coll$SDF$PEXPOSURE_VIF, freq=T, main="VIF Exposure to TCE Sites", xlab="VIF", xlim=c(1,2), ylim=c(0,200))
abline(v=c(5,10), col=c("red","darkred"),lty=3, lwd=5)

# Note that I have graphed a vertical line at x=5 and 10, these are traditionally cutoff
# values which we need to pay attention to.  We don't have a problem with our data, so these
# lines don't even show up.

windows()
par(mfrow=c(3,1))
hist(abs(gwr.coll$SDF$Corr_PCTAGE65P.PEXPOSURE), probability=F, main="EXPOSURE-AGE65+", xlab="Corr Coeff |r|", xlim=c(0,1), ylim=c(0,100))
abline(v=.8, col="red",lty=3, lwd=3)
hist(abs(gwr.coll$SDF$Corr_PCTOWNHOME.PEXPOSURE), probability=F, main="EXPOSURE-OWNHOME", xlab="Corr Coeff |r|", xlim=c(0,1), ylim=c(0,100))
abline(v=.8, col="red",lty=3, lwd=3)
hist(abs(gwr.coll$SDF$Corr_PCTAGE65P.PCTOWNHOME), probability=F, main="OWNHOME-AGE65+", xlab="Corr Coeff |r|", xlim=c(0,1), ylim=c(0,100))
abline(v=.8, col="red",lty=3, lwd=3)

# Looking at these, we don't seem to have any problems with local collinearities, either.
# Fotheringham suggests looking at correlations of >0.8, which I've marked here with a 
# straight red line.

#####
##### NOW LET'S RUN THE ROBUST MODEL:
#####

# If we have outliers, we can run the robust regression.  The code is below.  However,
# this can take a VERY long time, so I'd recommend not running it until you have some time
# to dedicate to waiting around for the results

gwr.rob<-gwr.robust(CASES~PCTAGE65P+PCTOWNHOME+PEXPOSURE, 
                    data=ny, bw=bw.abi, kernel="bisquare", adaptive=T, F123.test=T)
print(gwr.rob)

gwr.rob.t.adj <- gwr.t.adjust(gwr.rob)

moran.mc(gwr.rob$SDF$residual, ny.q, nsim=99)
#mixed model 
gwr.mix<-gwr.mixed(CASES~PCTAGE65P+PCTOWNHOME+PEXPOSURE,
                   fixed.vars=c("PCTAGE65P", "PCTOWNHOME"),
                   data=ny.spdf, bw=bw.abi, kernel="bisquare", adaptive=T, 
                   diagnostic=T)
print(gwr.mix)

