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
# MAIN RESULTS
################################################################################

################################################################################
# DESCRIPTIVE (TABLE 1 IN THE MANUSCRIPT)

# NUMBER OF LOCATIONS
length(dlist)

# TOTAL NUMBER OF DEATHS
sum(sapply(dlist,function(x) sum(x$death,na.rm=T)))

# TEMPERATURE DISTRIBUTION IN THE FIRST AND LAST PERIOD
rowMeans(sapply(dlist,function(x) summary(subset(x,year<1999.5)$tmean)))
rowMeans(sapply(dlist,function(x) summary(subset(x,year>1999.5)$tmean)))

################################################################################
# EFFECTS BY COUNTRY (TABLE 2 IN THE MANUSCRIPT)
# NB: NOT IDENTICAL TO MANUSCRIPT, AS BASED ON UK ONLY

# FUNCTION FOR MULTIVARIATE WALD TEST (FIRST VERSION, BASED ON COEF-VCOV)
fwald1 <- function(b1,V1,b2=NULL,V2=NULL) {
  invVp <- if(is.null(b2)) solve(V1) else solve(V1+V2)
  b <- if(is.null(b2)) b1 else b1-b2
  stat <- t(b)%*%invVp%*%(b)
  df <- length(b1)
  pchisq(stat,df,lower.tail=FALSE)
}

# MINIMUM MORTALITY PERCENTILE (MMP)
cenpercountry

# PERIOD
(period <- range(dlist[[1]]$year))

# RR AT 90TH AND 99TH VS MMP (WITH 95%CI)
# AVERAGE, 1993 AND 2006
cp$allRRfit[predper==90];cp$allRRlow[predper==90];cp$allRRhigh[predper==90]
cp1$allRRfit[predper==90];cp1$allRRlow[predper==90];cp1$allRRhigh[predper==90]
cp2$allRRfit[predper==90];cp2$allRRlow[predper==90];cp2$allRRhigh[predper==90]

cp$allRRfit[predper==99];cp$allRRlow[predper==99];cp$allRRhigh[predper==99]
cp1$allRRfit[predper==99];cp1$allRRlow[predper==99];cp1$allRRhigh[predper==99]
cp2$allRRfit[predper==99];cp2$allRRlow[predper==99];cp2$allRRhigh[predper==99]

# MULTIVARIATE TESTS FOR A NULL INTERACTION (P-VALUE)
fwald1(coef(cpint),vcov(cpint))
fwald1(coef(cp1),vcov(cp1),coef(cp2),vcov(cp2))

################################################################################
# META-REGRESSION MODEL (TABLE 3 IN THE MANUSCRIPT)
# NB: NOT IDENTICAL TO MANUSCRIPT, AS BASED ON UK ONLY

# FUNCTION FOR MULTIVARIATE WALD TEST (SECOND VERSION, BASED MODEL AND VAR)
# NB: SAME TEST OF VERSION 1, ONLY DIFFERENT WAY TO EXTRACT COEF-VCOV
fwald2 <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}

# FULL MODEL: TEST FOR PREDICTORS
fwald2(mv,"avgtmean")
fwald2(mv,"rangetmean")

# FULL MODEL: Q TEST
qtest(mv)

# FULL MODEL: I-SQUARE
summary(mv)
# OR
Q <- qtest(mv)
(Q$Q[1]-Q$df[1])/Q$Q[1]*100

#
