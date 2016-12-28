################################################################################
# Updated version of the R code for the analysis in:
#
#   "Temporal variation in heat-mortality associations: a multicountry study"
#   Antonio Gasparrini and collaborators
#   Environmental Health Perspectives - 2015
#   http://www.ag-myresearch.com/ehp2015.html
#
# Update: 14 March 2016
# For any problem with this code, please contact antonio.gasparrini@lshtm.ac.uk
# Please refer to the original code for any copyright issue
#
#  See www.ag-myresearch.com for future updates
################################################################################

################################################################################
# LAG-RESPONSE ASSOCIATIONS (USING CITY-SPECIFIC CENTERING)
################################################################################

################################################################################
# CREATE THE OBJECTS TO STORE THE RESULTS
# NOW LAG-RESPONSE ASSOCIATION AT THE 99TH PERCENTILE VS THE MMP

# OVERALL CUMULATIVE
coeflag <- coeflag1 <- coeflag2 <- matrix(NA,nrow(cities),4,
  dimnames=list(cities$city))
vcovlag <- vcovlag1 <- vcovlag2 <- vector("list",nrow(cities))
names(vcovlag) <- names(vcovlag1) <- names(vcovlag2) <- cities$city

################################################################################
# RUN THE MODEL FOR EACH CITY (WITH RE-CENTERED VALUES)

# LOOP
for(i in seq(nrow(cities))) {

  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data <- dlist[[i]]
  
  # DEFINE THE CROSS-BASIS (NOW CENTERED AT MINIMUM MORTALITY TEMPERATURE)
  argvar <- list(fun="bs",degree=2,knots=quantile(data$tmean,
    varper/100,na.rm=T))
  cb <- crossbasis(data$tmean,lag=lag,argvar=argvar,
    arglag=list(knots=logknots(lag,lagnk)),group=data$season)
  #summary(cb)
  
  # RUN THE MODELS WITHOUT INTERACTION
  model <- glm(formula,family=quasipoisson,data=data,na.action="na.exclude")
 
  # DEFINE INTERACTION TERMS, STANDARDIZED
  # -INT1: CENTERED AT 1ST OF AUGUST 1993
  # -INT1: CENTERED AT 1ST OF AUGUST 2006
  datenum <- as.numeric(data$date)
  day1 <- as.numeric(as.Date("1993-08-01"))
  day2 <- as.numeric(as.Date("2006-08-01"))
  int1 <- ((datenum-day1)/(day2-day1))*cb
  int2 <- ((datenum-day2)/(day2-day1))*cb
  
  # RUN THE MODELS WITH INTERACTION
  model1 <- glm(formula1,family=quasipoisson,data=data)
  model2 <- glm(formula2,family=quasipoisson,data=data)
  
  # PREDICTIONS AND REDUCTION TO LAG-RESPONSE AT 99TH
  # - AVERAGE, 1993, 2006, INTERACTION TERMS
  # NB: CENTERING NEEDED HERE, AS IT CHANGES COEF-VCOV
  cen <- quantile(data$tmean,cenpercountry/100,na.rm=T)
  redlag <- crossreduce(cb,model,"var",value=quantile(data$tmean,0.99,na.rm=T),
    cen=cen)
  coeflag[i,] <- coef(redlag)
  vcovlag[[i]] <- vcov(redlag)

  redlag1 <- crossreduce(cb,model1,"var",value=quantile(data$tmean,0.99,na.rm=T),
    cen=cen)
  coeflag1[i,] <- coef(redlag1)
  vcovlag1[[i]] <- vcov(redlag1)

  redlag2 <- crossreduce(cb,model2,"var",value=quantile(data$tmean,0.99,na.rm=T),
    cen=cen)
  coeflag2[i,] <- coef(redlag2)
  vcovlag2[[i]] <- vcov(redlag2)
}

################################################################################
# META-ANALYSIS

# RUN THE MODELS FOR AVERAGE EFFECT
mvlag <- mvmeta(coeflag~avgtmean+rangetmean,vcovlag,data=cities,
  control=list(showiter=T))
summary(mvlag)

# RUN THE MODELS FOR 1993 AND 2006
mvlag1 <- mvmeta(coeflag1~avgtmean+rangetmean,vcovlag1,cities,
  control=list(showiter=T))
summary(mvlag1)

mvlag2 <- mvmeta(coeflag2~avgtmean+rangetmean,vcovlag2,cities,
  control=list(showiter=T))
summary(mvlag2)

################################################################################
# PREDICT THE POOLED COEFFICIENTS

mvpredlag <- predict(mvlag,datanew,vcov=T,format="list")
mvpredlag1 <- predict(mvlag1,datanew,vcov=T,format="list")
mvpredlag2 <- predict(mvlag2,datanew,vcov=T,format="list")

################################################################################
# PREDICT THE POOLED LAG-RESPONSE ASSOCIATIONS

# OBTAIN THE PREDICTIONS
blag <- do.call(onebasis,c(list(x=seq(0,lag)),attr(cb,"arglag")))

cplag <- crosspred(blag,coef=mvpredlag$fit,vcov=mvpredlag$vcov,model.link="log",
  at=0:100/10)
cplag1 <- crosspred(blag,coef=mvpredlag1$fit,vcov=mvpredlag1$vcov,
  model.link="log",at=0:100/10)
cplag2 <- crosspred(blag,coef=mvpredlag2$fit,vcov=mvpredlag2$vcov,
  model.link="log",at=0:100/10)

#

