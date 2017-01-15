################################################################################
# Updated version of the R code for the analysis in:
#
#   "Temporal variation in heat-mortality associations: a multicountry study"
#   Antonio Gasparrini and collaborators
#   Environmental Health Perspectives - 2015
#   http://www.ag-myresearch.com/2015_gasparrini_ehp.html
#
# Update: 15 January 2017
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/2015_gasparrini_EHP_Rcodedata
################################################################################

################################################################################
# MULTIVARIATE META-ANALYSIS OF THE REDUCED COEF
################################################################################

# CREATE AVERAGE TEMPERATURE AND RANGE AS META-PREDICTORS
avgtmean <- sapply(dlist,function(x) mean(x$tmean,na.rm=T))
rangetmean <- sapply(dlist,function(x) diff(range(x$tmean,na.rm=T)))

################################################################################
# META-ANALYSIS
# NB: COUNTRY EFFECT IS NOT INCLUDED IN THIS EXAMPLE

# RUN THE MODELS FOR AVERAGE EFFECT
mv <- mvmeta(coef~avgtmean+rangetmean,vcov,data=cities,control=list(showiter=T))
summary(mv)
# NB: IN THIS EXAMPLE THE MV-META MODEL IS CLEARLY OVERPARAMETERIZED

# RUN THE MODELS FOR 1993 AND 2006
mv1 <- mvmeta(coef1~avgtmean+rangetmean,vcov1,cities,control=list(showiter=T))
summary(mv1)
mv2 <- mvmeta(coef2~avgtmean+rangetmean,vcov2,cities,control=list(showiter=T))
summary(mv2)
# NOTE THAT THERE ARE ISSUES WITH ESTIMATION OF THE BETWEEN-LOCATION 
# (CO)VARIANCE MATRIX. THESE ARE DUE TO THE SMALL SAMPLE USED HERE

# RUN THE MODELS FOR INTERACTION
mvint <- mvmeta(coefint~avgtmean+rangetmean,vcovint,cities,
  control=list(showiter=T))
summary(mvint)

################################################################################
# PREDICT THE POOLED COEFFICIENTS

datanew <- data.frame(avgtmean=mean(tapply(avgtmean,cities$city,mean)),
  rangetmean=mean(tapply(rangetmean,cities$city,mean)))

mvpred <- predict(mv,datanew,vcov=T,format="list")
mvpred1 <- predict(mv1,datanew,vcov=T,format="list")
mvpred2 <- predict(mv2,datanew,vcov=T,format="list")
mvpredint <- predict(mvint,datanew,vcov=T,format="list")

################################################################################
# BLUPS

blup <- blup(mv,vcov=T)
blup1 <- blup(mv1,vcov=T)
blup2 <- blup(mv2,vcov=T)

################################################################################
# RE-CENTERING

# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
set.seed(13041975)
tmeancountry <- rowMeans(sapply(dlist,function(x) quantile(jitter(x$tmean),
  predper/100,na.rm=T)))

# DEFINE INDICATOR FOR CENTERING PERCENTILE FROM AVERAGE ASSOCIATION
bvar <- onebasis(tmeancountry,fun="bs",degree=2,
  knots=tmeancountry[paste(varper,".0%",sep="")])
cenindcountry <- which.min(bvar%*%mvpred$fit)

# DEFINE CENTERING PERCENTILE FOR EACH COUNTRY AND CITY
cenpercountry <- pmin(pmax(predper[cenindcountry],10),90)

################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS

# OBTAIN THE CENTERED PREDICTIONS
bvar <- onebasis(tmeancountry,fun="bs",degree=2,
  knots=tmeancountry[paste(varper,".0%",sep="")])
cen <- tmeancountry[paste(cenpercountry,".0%",sep="")]

cp <- crosspred(bvar,coef=mvpred$fit,vcov=mvpred$vcov,model.link="log",
  at=tmeancountry,cen=cen)
cp1 <- crosspred(bvar,coef=mvpred1$fit,vcov=mvpred1$vcov,model.link="log",
  at=tmeancountry,cen=cen)
cp2 <- crosspred(bvar,coef=mvpred2$fit,vcov=mvpred2$vcov,model.link="log",
  at=tmeancountry,cen=cen)
cpint <- crosspred(bvar,coef=mvpredint$fit,vcov=mvpredint$vcov,model.link="log",
  at=tmeancountry,cen=cen)

#
