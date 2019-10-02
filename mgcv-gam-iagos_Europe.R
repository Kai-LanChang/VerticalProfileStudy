#This example is made for the trend analysis with IAGOS data above Europe
#Kai-Lan Chang (1st Oct 2019)

####################################################
##read data
####################################################
#All the output will be saved at following directory
dir="C:/Users/kchang/Downloads/ssgam/"
library(mgcv)
library(dichromat)
library(viridis)
library(fields)
col=colorschemes$DarkRedtoBlue.18
dcol=c('blue 4',col,'red 4')

mm=read.table("C://Users/kchang/Downloads/ssgam/iagos/o3_profiles_Europe_IAGOS_forTrendsCalculation.txt",header=F, sep=",",skip=9)
colnames(mm)=c("YY", "MM", "DD", "HH_start", "HH_end", "AIRPORT",
    "h100", "h150", "h200", "h250", "h300", "h350", "h400", "h450", "h500",
    "h550", "h600", "h650", "h700", "h750", "h800", "h850", "h900", "h950", "h1000")

cycle=aggregate(.~MM, mm, function(x) mean(x,na.rm=T), na.action=NULL)
for (j in 1:nrow(mm)){
	for (i in 1:12){
	if (mm[j,2]==i) {mm[j,7:ncol(mm)]=mm[j,7:ncol(mm)]-cycle[i,7:ncol(cycle)]}
}}

mm=mm[mm$YY<=2016,]
mm=aggregate(.~DD+MM+YY, mm, function(x) mean(x,na.rm=T), na.action=NULL)
date=cbind(MM=rep(seq(1,12),23), YY=rep(seq(1994,2016), each=12), ind=seq(1,12*23))
dt=merge(date, mm, by=c("YY","MM"), all=T)
dt=dt[,-ncol(dt)]

cm=cor(dt[,11:25], use="complete.obs")
library(corrplot)
png(paste(dir,"o3_iagos_eu_vert_cor.png", sep=""), height=10, width=10, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
corrplot(cm, type="upper", tl.col="black", col=dcol)  
corrplot.mixed(cm, tl.col="black", lower.col=dcol, upper.col=dcol) 
corrplot.mixed(cm, tl.col="black", lower.col=2, upper.col=dcol) 
dev.off()
##########################################

mm=read.table("C://Users/kchang/Downloads/ssgam/iagos/o3_profiles_Europe_IAGOS_forTrendsCalculation.txt",header=F, sep=",",skip=9)
colnames(mm)=c("YY", "MM", "DD", "HH_start", "HH_end", "AIRPORT",
    "h100", "h150", "h200", "h250", "h300", "h350", "h400", "h450", "h500",
    "h550", "h600", "h650", "h700", "h750", "h800", "h850", "h900", "h950", "h1000")


mm=mm[mm$YY<=2017,]
ny=24
mm=aggregate(.~DD+MM+YY, mm, function(x) mean(x,na.rm=T), na.action=NULL)
date=cbind(MM=rep(seq(1,12),ny), YY=rep(seq(1994,2017), each=12), ind=seq(1,12*ny))
dt=merge(date, mm, by=c("YY","MM"), all=T)
dt=dt[,-ncol(dt)]

dt1=cbind(dt[,1:3], o3=dt$h250, alt=0.25)
dt2=cbind(dt[,1:3], o3=dt$h300, alt=0.30)
dt3=cbind(dt[,1:3], o3=dt$h350, alt=0.35)
dt4=cbind(dt[,1:3], o3=dt$h400, alt=0.40)
dt5=cbind(dt[,1:3], o3=dt$h450, alt=0.45)
dt6=cbind(dt[,1:3], o3=dt$h500, alt=0.50)
dt7=cbind(dt[,1:3], o3=dt$h550, alt=0.55)
dt8=cbind(dt[,1:3], o3=dt$h600, alt=0.60)
dt9=cbind(dt[,1:3], o3=dt$h650, alt=0.65)
dt10=cbind(dt[,1:3], o3=dt$h700, alt=0.70)
dt11=cbind(dt[,1:3], o3=dt$h750, alt=0.75)
dt12=cbind(dt[,1:3], o3=dt$h800, alt=0.80)
dt13=cbind(dt[,1:3], o3=dt$h850, alt=0.85)
dt14=cbind(dt[,1:3], o3=dt$h900, alt=0.90)
dt15=cbind(dt[,1:3], o3=dt$h950, alt=0.95)
dt=rbind(dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9,dt10,dt11,dt12,dt13,dt14,dt15)


####################################################
##model formulation
####################################################
dt=aggregate(o3~alt+YY+MM+ind, data=dt, FUN=function(x) {mean(x,na.rm=T)})
mod2=gam(o3 ~ s(MM,alt,bs="ds",k=120) +s(YY,alt,bs="ds",k=300), na.action=na.omit, data=dt)
png(paste(dir,"o3_iagos_eu_gam_check.png", sep=""), height=10, width=10, units="in", res=500)
gam.check(mod2)
dev.off()

png(paste(dir,"o3_iagos_eu_obs_fit.png", sep=""), height=7, width=7, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
plot(mod2$mod$o3, mod2$fitted.values, xlab="Observations", ylab="Fitted values",
    ylim=c(0,160), xlim=c(0,160), cex.lab=1.3)
abline(lm(mod2$fitted.values~mod2$model$o3), lwd=2) 
#lo=loess(mod2$fitted.values~mod2$model$o3)
#j=order(mod2$model$o3)
#lines(mod2$model$o3[j],lo$fitted[j],col="red",lwd=2)
abline(a=0,b=1,lty=2)
dev.off()

png(paste(dir,"o3_iagos_eu_res_fit.png", sep=""), height=7, width=7, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
rr=(residuals(mod2)-mean(residuals(mod2)))/sd(residuals(mod2))
plot(mod2$fitted.values,rr, xlab="Fitted values", ylab="Residuals", cex.lab=1.3)
abline(lm(rr~mod2$fitted.values), lwd=2)#col='red'
dev.off()

png(paste(dir,"o3_iagos_eu_res_month.png", sep=""), height=7, width=7, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
rr=(residuals(mod2)-mean(residuals(mod2)))/sd(residuals(mod2))
boxplot(rr~mod2$model$MM, xlab="Month", ylab="Residuals", ylim=c(-8,8), cex.lab=1.3)
dev.off()

png(paste(dir,"o3_iagos_eu_res_year.png", sep=""), height=7, width=7, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
rr=(residuals(mod2)-mean(residuals(mod2)))/sd(residuals(mod2))
boxplot(rr~mod2$model$YY, xlab="Year", ylab="Residuals", ylim=c(-8,8), xaxt="n", cex.lab=1.3)
axis(1,at=seq(1,ny,by=5),labels=seq(1994,2017, by=5),cex=1.5)
dev.off()

png(paste(dir,"o3_iagos_eu_res_pressure.png", sep=""), height=7, width=7, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
rr=(residuals(mod2)-mean(residuals(mod2)))/sd(residuals(mod2))
alt=mod2$model$alt*1000
boxplot(rr~alt, xlab="Pressure [hPa]", ylab="Residuals", ylim=c(-8,8), cex.lab=1.3)
dev.off()



aa=plot(mod2, select=2,scheme=2,n2=150,too.far=0.2)
surf=matrix(aa[[2]]$fit,150) 
surf=surf[,ncol(surf):1]
LL=aa[[2]]$x
LL=seq(1994,2017, length=150)
ll=aa[[2]]$y*1000
png(paste(dir,"o3_iagos_eu_trend.png", sep=""), height=10, width=14, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
image.plot(LL,ll,surf, xlab="Year", ylab="Pressure [hPa]", horizontal=F, col=dcol,
    yaxt="n", breaks=seq(-25, 25, length.out=21),ylim=c(250,950), cex.lab=1.5)
contour(LL,ll,surf,drawlabels=T, add=TRUE)
axis(2,at=seq(250,950,by=50),labels=rev(seq(250,950,by=50)))
mtext('[ppb]', side=3, line=0, at=2018, cex=1.5)
dev.off()

pr=aa
prm=matrix(pr[[2]]$fit,150)
prm=prm[,seq(1,150,length=15)] 
pry=seq(1994,2017,length=150)
surf[surf>=25]=25
png(paste(dir,"o3_iagos_eu_tstrend.png", sep=""), height=10, width=14, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
par(fig=c(0,10,4,10)/10)
image.plot(LL,ll,surf, xlab="", ylab="Pressure [hPa]", horizontal=F, col=dcol,
    yaxt="n", breaks=seq(-25, 25, length.out=21),ylim=c(250,950), cex.lab=1.5)

par(fig=c(0,10,4,10)/10)
par(new=T)
contour(LL,ll,surf,drawlabels=T, add=TRUE)
axis(2,at=seq(250,950,by=50),labels=rev(seq(250,950,by=50)))
mtext('[ppb]', side=3, line=0, at=2017.8, cex=1.5)

par(fig=c(0,10,1,4.5)/10)
par(new=T)
matplot(pry, prm, xlab="Year", ylab="Anomalies [ppb]", cex.lab=1.5, xlim=c(1995,2016), col=cividis(ncol(prm)), type="l", lwd=3)
image.plot(legend.only=T, nlevel=ncol(prm), col=rev(cividis(ncol(prm))),
    zlim=range(prm),
    lab.breaks=c(NA,950,NA,850,NA,750,NA,650,NA,550,NA,450,NA,350,NA,250))
par(fig=c(0,10,1,4.5)/10)
par(new=T)
mtext('[hPa]', side=3, line=1, at=2017.6, cex=1.5)
dev.off()


aa=plot(mod2, select=1,scheme=2,n2=150,too.far=0.2)
surf=matrix(aa[[1]]$fit,150) 
surf=surf[,ncol(surf):1]
LL=aa[[1]]$x
ll=aa[[1]]$y*1000
png(paste(dir,"o3_iagos_eu_month.png", sep=""), height=10, width=14, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
image.plot(LL,ll,surf, xlab="Month", ylab="Pressure [hPa]", horizontal=F, col=dcol,
    yaxt="n", breaks=seq(-40, 40, length.out=21), ylim=c(250,950), cex.lab=1.5)
contour(LL,ll,surf,drawlabels=T, add=TRUE)
axis(2,at=seq(250,950,by=50),labels=rev(seq(250,950,by=50)))
mtext('[ppb]', side=3, line=0, at=12.4, cex=1.5)
dev.off()


pr=aa
prm=matrix(pr[[1]]$fit,150)
prm=prm[,seq(1,150,length=15)]  
pry=pr[[1]]$x
png(paste(dir,"o3_iagos_eu_tsmonth.png", sep=""), height=10, width=14, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
par(fig=c(0,10,4,10)/10)
image.plot(LL,ll,surf, xlab="", ylab="Pressure [hPa]", horizontal=F, col=dcol,
    yaxt="n", xaxt="n",
	breaks=seq(-40, 40, length.out=21), ylim=c(250,950), cex.lab=1.5)
par(fig=c(0,10,4,10)/10)
par(new=T)
contour(LL,ll,surf,drawlabels=T, add=TRUE)
axis(1,at=seq(1,12,by=1),labels=seq(1,12,by=1))
axis(2,at=seq(250,950,by=50),labels=rev(seq(250,950,by=50)))
mtext('[ppb]', side=3, line=0, at=12.4, cex=1.5)

par(fig=c(0,10,1,4.5)/10)
par(new=T)
matplot(pry, prm, xlab="Year", ylab="Anomalies [ppb]", cex.lab=1.5,
    xlim=c(1.35,11.65), col=cividis(ncol(prm)), type="l", lwd=3)
axis(1,at=seq(1,12,by=1),labels=seq(1,12,by=1))
image.plot(legend.only=T, nlevel=ncol(prm), col=rev(cividis(ncol(prm))),
    zlim=range(prm),
    lab.breaks=c(NA,950,NA,850,NA,750,NA,650,NA,550,NA,450,NA,350,NA,250))
par(fig=c(0,10,1,4.5)/10)
par(new=T)
mtext('[hPa]', side=3, line=1, at=12.45, cex=1.5)
dev.off()


aa=plot(mod2, select=1,scheme=2,n2=150,too.far=0.2)
surf=matrix(aa[[1]]$fit,150) +unlist(summary(mod2)$p.coeff[1])
surf=surf[,ncol(surf):1]
LL=aa[[1]]$x
ll=aa[[1]]$y*1000
png(paste(dir,"o3_iagos_eu_season.png", sep=""), height=10, width=14, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
image.plot(LL,ll,surf, xlab="Month", ylab="Pressure [hPa]", horizontal=F,
    col=viridis(20), yaxt="n", xaxt="n",
	breaks=seq(20, 100, length.out=21), ylim=c(250,950), cex.lab=1.5)
contour(LL,ll,surf,drawlabels=T, add=TRUE)
axis(1,at=seq(1,12,by=1),labels=seq(1,12,by=1))
axis(2,at=seq(250,950,by=50),labels=rev(seq(250,950,by=50)))
mtext('[ppb]', side=3, line=0, at=12.4, cex=1.5)
dev.off()


pr=aa
prm=matrix(pr[[1]]$fit,150) +unlist(summary(mod2)$p.coeff[1])
prm=prm[,seq(1,150,length=15)]
pry=pr[[1]]$x
png(paste(dir,"o3_iagos_eu_tsseason.png", sep=""), height=10, width=14, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
par(fig=c(0,10,4,10)/10)
image.plot(LL,ll,surf, xlab="", ylab="Pressure [hPa]", horizontal=F,
    col=viridis(20), yaxt="n", xaxt="n",
	breaks=seq(20, 100, length.out=21), ylim=c(250,950), cex.lab=1.5)
par(fig=c(0,10,4,10)/10)
par(new=T)
contour(LL,ll,surf,drawlabels=T, add=TRUE)
axis(1,at=seq(1,12,by=1),labels=seq(1,12,by=1))
axis(2,at=seq(250,950,by=50),labels=rev(seq(250,950,by=50)))
mtext('[ppb]', side=3, line=0, at=12.4, cex=1.5)

par(fig=c(0,10,1,4.5)/10)
par(new=T)
matplot(pry, prm, xlab="Year", ylab="climatologies [ppb]", cex.lab=1.5,
    xlim=c(1.35,11.65), col=cividis(ncol(prm)), type="l", lwd=3)
axis(1,at=seq(1,12,by=1),labels=seq(1,12,by=1))
image.plot(legend.only=T, nlevel=ncol(prm), col=rev(cividis(ncol(prm))),
    zlim=range(prm),
    lab.breaks=c(NA,950,NA,850,NA,750,NA,650,NA,550,NA,450,NA,350,NA,250))
par(fig=c(0,10,1,4.5)/10)
par(new=T)
mtext('[hPa]', side=3, line=1, at=12.45, cex=1.5)
dev.off()
