################################################################################
# Updated version of the R code for the analysis in:
#
#   "Temporal variation in heat-mortality associations: a multicountry study"
#   Antonio Gasparrini and collaborators
#   Environmental Health Perspectives - 2015
#   http://www.ag-myresearch.com/2015_gasparrini_ehp.html
#
# Update: 05 December 2017
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/2015_gasparrini_EHP_Rcodedata
################################################################################

################################################################################
# FIRST-STAGE ANALYSIS: RUN THE MODEL IN EACH CITY, REDUCE AND SAVE
################################################################################

################################################################################
# CREATE THE OBJECTS TO STORE THE RESULTS FOR 
# THE REDUCED OVERALL CUMULATIVE EXPOSURE-RESPONSE ONLY

coef <- coef1 <- coef2 <- coefint <- matrix(NA,nrow(cities),4,
  dimnames=list(cities$city))
vcov <- vcov1 <- vcov2 <- vcovint <- vector("list",nrow(cities))
names(vcov) <- names(vcov1) <- names(vcov2) <- names(vcovint) <- cities$city

################################################################################
# RUN THE MODEL FOR EACH CITY

# LOOP
for(i in seq(nrow(cities))) {

  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data <- dlist[[i]]
  
  # DEFINE THE CROSS-BASIS
  argvar <- list(fun="bs",degree=2,knots=quantile(data$tmean,
    varper/100,na.rm=T))
  cb <- crossbasis(data$tmean,lag=lag,argvar=argvar,
    arglag=list(knots=logknots(lag,lagnk)),group=data$year)
  #summary(cb)
  
  # RUN THE MODELS WITHOUT INTERACTION
  model <- glm(formula,family=quasipoisson,data=data,na.action="na.exclude")
 
  # DEFINE INTERACTION TERMS, STANDARDIZED
  # -INT1: CENTERED AT 1ST OF AUGUST 1993
  # -INT2: CENTERED AT 1ST OF AUGUST 2006
  datenum <- as.numeric(data$date)
  day1 <- as.numeric(as.Date("1993-08-01"))
  day2 <- as.numeric(as.Date("2006-08-01"))
  int1 <- ((datenum-day1)/(day2-day1))*cb
  int2 <- ((datenum-day2)/(day2-day1))*cb
  
  # RUN THE MODELS WITH INTERACTION
  model1 <- glm(formula1,family=quasipoisson,data=data)
  model2 <- glm(formula2,family=quasipoisson,data=data)
  
  # PREDICTION AND REDUCTION TO OVERALL CUMULATIVE EXPOSURE-RESPONSE
  # - AVERAGE, 1993, 2006, INTERACTION TERMS
  # NB: CENTERING NOT NEEDED AT THIS STAGE, AS IT DOES NOT CHANGE COEF-VCOV
  cen <- mean(data$tmean,na.rm=T)
  red <- crossreduce(cb,model,cen=cen)
  coef[i,] <- coef(red)
  vcov[[i]] <- vcov(red)

  red1 <- crossreduce(cb,model1,cen=cen)
  coef1[i,] <- coef(red1)
  vcov1[[i]] <- vcov(red1)

  red2 <- crossreduce(cb,model2,cen=cen)
  coef2[i,] <- coef(red2)
  vcov2[[i]] <- vcov(red2)

  redint <- crossreduce(int1,model1,cen=cen)
  coefint[i,] <- coef(redint)
  vcovint[[i]] <- vcov(redint)
}

#

