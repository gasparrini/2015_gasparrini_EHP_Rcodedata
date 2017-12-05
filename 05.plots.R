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
# PLOTS
################################################################################

################################################################################
# MAP

library(maps) ; library(mapdata)
#library(mapproj)

#library(maptools)

summary(avgtmean)
#hist(avgtmean)
levels <- cut(avgtmean,c(-100,7:15*2+1,100),labels=c('<15',8:15*2,'>31'))
funcol <- colorRampPalette(c("yellow","orange","red"))
col <- funcol(length(levels(levels)))[levels]
size <- (c(65,70,26,26,27,29,70)/70)[cities$country]

pdf("figure1.pdf",width=6,height=5)

map("worldHires",xlim=c(-22,12),ylim=c(47,62),mar=c(0,0,0,0),myborder=0.1)
map("worldHires","UK:Great Britain",add=T,fill=T,col=grey(0.9))
symbols(cities$long,cities$lat,circles=rep(0.5,10),inches=F,bg=col,add=T)
map.scale(-15,47,ratio=F,cex=0.8)
text(-12.5,62,"United\nKingdom",font=2,cex=1.3)

legend(-19,53,paste0(c('<15',8:15*2,'>31'),"C"),xjust=0.5,yjust=0.5,pt.cex=1.5,
  pch=21,pt.bg=funcol(length(levels(levels))),bty="n",cex=0.8,
  title="Summer average\ntemperature")

dev.off()

################################################################################
# BY CITY: AVERAGE

pdf("cityavg.pdf",width=9,height=9)
layout(matrix(seq(3*3),nrow=3,byrow=T))
par(mex=0.8,mgp=c(2.5,1,0),las=1)

for(i in seq(dlist)) {
  bvar <- onebasis(dlist[[i]]$tmean,fun="bs",degree=2,
    knots=quantile(dlist[[i]]$tmean,varper/100,na.rm=T))
  cen <- quantile(dlist[[i]]$tmean,cenpercountry/100,na.rm=T)
  pred <- crosspred(bvar,coef=blup[[i]]$blup,vcov=blup[[i]]$vcov,
    model.link="log",cen=cen)
  plot(pred,ylab="RR",xlab="Summer temperature",ylim=c(0.8,2),
    lwd=2,col=2,main=cities$cityname[i])
}

dev.off()

################################################################################
# BY CITY: 1993-2006 (FIGURE S2 IN SUPPLEMENTAL MATERIAL)

pdf("figureS2.pdf",width=9,height=9)
layout(matrix(seq(3*3),nrow=3,byrow=T))
par(mex=0.8,mgp=c(2.5,1,0),las=1)

for(i in seq(dlist)) {
  bvar <- onebasis(dlist[[i]]$tmean,fun="bs",degree=2,
    knots=quantile(dlist[[i]]$tmean,varper/100,na.rm=T))
  cen <- quantile(dlist[[i]]$tmean,cenpercountry/100,na.rm=T)
  pred1 <- crosspred(bvar,coef=blup1[[i]]$blup,vcov=blup1[[i]]$vcov,
    model.link="log",cen=cen)
  pred2 <- crosspred(bvar,coef=blup2[[i]]$blup,vcov=blup2[[i]]$vcov,
    model.link="log",cen=cen)
  plot(pred1,ylab="RR",xlab="Summer temperature ",ylim=c(0.8,2),
    lwd=2,col=3,ci.arg=list(density=20,col=3))
  lines(pred2,ci="area",lwd=2,col=4,ci.arg=list(density=20,angle=-45,col=4))
  title(cities$cityname[i])
  legend("top",as.character(period),col=3:4,lwd=2,bg="white",
    cex=0.6,ncol=2,inset=0.05)
}

dev.off()

################################################################################
# POOLED: AVERAGE (FIGURE 2 AND S3 IN SUPPLEMENTAL MATERIAL)

indlab <- predper%in%c(0,1,10,50,90,99,100)

pdf("figure2.pdf",width=4,height=4)
layout(1)
par(mex=0.8,mgp=c(2.5,1,0),las=1)

plot(cp,ylab="RR",xlab="Summer temperature percentile",axes=F,
  ylim=c(0.8,1.8),lwd=2,col=2,main="UK")
axis(1,at=tmeancountry[indlab],labels=predper[indlab],cex.axis=0.9)
axis(2,cex.axis=0.9)
abline(v=c(tmeancountry[cenindcountry],tmeancountry[c("90.0%","99.0%")]),
  lty=c(3,2,2))
mtext(paste(period,collapse="-"),cex=0.7,line=0)

dev.off()

pdf("figureS4.pdf",width=4,height=4)
layout(1)
par(mex=0.8,mgp=c(2.5,1,0),las=1)

plot(cplag,ylab="RR",xlab="Lag",ylim=c(0.95,1.15),lwd=2,col=2,main="UK",
  cex.axis=0.9)
mtext(paste(period,collapse="-"),cex=0.7,line=0)

dev.off()

################################################################################
# POOLED: 1993-2006 (FIGURE 2 AND 3 IN THE MANUSCRIPT)

pdf("figure3.pdf",width=4,height=4)
layout(1)
par(mex=0.8,mgp=c(2.5,1,0),las=1)

plot(cp1,ylab="RR",xlab="Summer temperature percentile",axes=F,
  ylim=c(0.8,1.8),lwd=2,col=3,ci.arg=list(density=20,col=3),main="UK")
lines(cp2,ci="area",lwd=2,col=4,ci.arg=list(density=20,angle=-45,col=4))
axis(1,at=tmeancountry[indlab],labels=predper[indlab],cex.axis=0.9)
axis(2,cex.axis=0.9)
abline(v=c(tmeancountry[cenindcountry],tmeancountry[c("90.0%","99.0%")]),
  lty=c(3,2,2))
legend("top",as.character(period),col=3:4,lwd=2,bg="white",cex=0.8,ncol=2)

dev.off()

pdf("figure4.pdf",width=4,height=4)
layout(1)
par(mex=0.8,mgp=c(2.5,1,0),las=1)

plot(cplag1,ylab="RR",xlab="Lag",ylim=c(0.95,1.15),lwd=2,col=3,
  ci.arg=list(density=20,col=3),main="UK",cex.axis=0.9)
lines(cplag2,ci="area",lwd=2,col=4,ci.arg=list(density=20,angle=-45,col=4))
legend("top",as.character(period),col=3:4,lwd=2,bg="white",cex=0.8,ncol=2)

dev.off()

################################################################################
# POOLED: ONLY INTERACTION (FIGURE S3 IN SUPPLEMENTAL MATERIAL)

pdf("figureS3.pdf",width=4,height=4)
layout(1)
par(mex=0.8,mgp=c(2.5,1,0),las=1)

plot(cpint,ylab="RR",xlab="Summer temperature percentile",axes=F,
  ylim=c(0.8,1.8),lwd=2,col=2,main="UK")
axis(1,at=tmeancountry[indlab],labels=predper[indlab],cex.axis=0.9)
axis(2,cex.axis=0.9)
abline(v=c(tmeancountry[cenindcountry[i]],tmeancountry[c("90.0%","99.0%")]),
  lty=c(3,2,2))
mtext(paste(period,collapse="-"),cex=0.7,line=0)

dev.off()

#
