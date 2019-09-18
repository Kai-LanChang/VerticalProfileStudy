#The ozonesonde data can be found in NOAA ESRL Global Monitoring Division (GMD) database
#ftp://aftp.cmdl.noaa.gov/data/ozwv/Ozonesonde/. Preparation work includes installation of
#all R packages by using install.packages("."), and changes working directory to the folder
#which the data have been saved. All the output will be saved at following directory "dir".
#This example is made for the trend analysis at Trinidad Head, CA (1998-2018).

setwd("C:/Users/kchang/Downloads/ozonesonde-gmd/Trinidad Head, California/")
dir="C:/Users/kchang/Downloads/ssgam/"
####################################################
##read data
####################################################
library(viridis)
library(mgcv)
library(fields)
library(corrplot)
library(dichromat)
col=colorschemes$DarkRedtoBlue.18
col=c('blue 4',col,'red 4')
11
filelist = list.files(pattern = NULL)
datalist = lapply(filelist, function(x)read.table(x, header=F,sep="", skip=29, fill=T)) 

for(i in seq_along(filelist)){
	nn=data.frame(do.call('rbind', strsplit(as.character(filelist[[i]]),'_',fixed=TRUE)))[1,2]
	datalist[[i]]$Year = rep(nn, nrow(datalist[[i]]))
} 
for(i in seq_along(filelist)){
	nn=data.frame(do.call('rbind', strsplit(as.character(filelist[[i]]),'_',fixed=TRUE)))[1,3]
	datalist[[i]]$Month = rep(nn, nrow(datalist[[i]]))
} 
for(i in seq_along(filelist)){
	nn=data.frame(do.call('rbind', strsplit(as.character(filelist[[i]]),'_',fixed=TRUE)))[1,4]
	datalist[[i]]$Day = rep(nn, nrow(datalist[[i]]))
} 

dt = do.call("rbind", datalist) 
colnames(dt)=c('Level', 'Press', 'Alt', 'Pottp', 'Temp', 'FtempV', 'Hum', 
	'Ozone_mpa', 'Ozone_ppbv', 'Ozone_atm', 'Ptemp', 'O3DN', 'O3Res', 'O3Uncert',
	'Year','Month','Day')
dt$Year=as.numeric(as.character(dt$Year))
dt$Month=as.numeric(as.character(dt$Month))
dt$Day=as.numeric(as.character(dt$Day))
dt=dt[dt$Alt<=15,]
dt[dt[,2]>=9990,2]=NA
dt[dt[,9]>=90,9]=NA
dt[,9]=dt[,9]*1000
dt=dt[dt$Year>=1998,]
dt=dt[dt$Year<=2018,]

pr=dt[,c(2,3)]
pr$Alt=round(pr$Alt,1)
pr=aggregate(Press~Alt, data=pr, FUN=function(x) {mean(x,na.rm=T)})

####################################################
##model formulation
####################################################

date=cbind(Year=rep(seq(1998,2018), each=12), Month=rep(seq(1,12),21), ind=seq(1,12*21))
dt=aggregate(Ozone_ppbv~Alt+Year+Month, data=dt, FUN=function(x) {mean(x,na.rm=T)})

mod1=gam(Ozone_ppbv ~ s(Month,Alt,bs="ds",k=120) +s(Year,Alt,bs="ds",k=300), na.action=na.omit, data=dt)


aa=plot(mod1, select=2,scheme=2,n2=150,too.far=0.2)
surf=matrix(aa[[2]]$fit,150) 
LL=aa[[2]]$x
LL=seq(1998,2018, length=150)
ll=aa[[2]]$y
ll=rev(approx(x=pr$Alt,y=pr$Press,xout=ll)$y)
png(paste(dir,"o3_trinidad_0_15km_tstrend_pr.png", sep=""), height=10, width=14, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
par(fig=c(0,10,4,10)/10)
image.plot(LL,ll,surf, xlab="", ylab="Pressure [hPa]", horizontal=F, col=col, 
	breaks=seq(-40, 40, length.out=21),ylim=c(130,1000), cex.lab=1.5, yaxt="n")

par(fig=c(0,10,4,10)/10)
par(new=T)
contour(LL,ll,surf,drawlabels=T, add=TRUE)
axis(2,at=seq(150,1000,by=50),labels=rev(seq(150,1000,by=50)))
mtext('[ppb]', side=3, line=0, at=2019.5, cex=1.5)

par(fig=c(0,10,1,4.5)/10)
par(new=T)
matplot(LL, surf[,ncol(surf):1], xlab="Year", ylab="Anomaly [ppb]", cex.lab=1.5,
    xlim=c(1999,2017), col=rainbow(ncol(surf)), type="l")
image.plot(legend.only=T, nlevel=ncol(surf), col=rev(rainbow(ncol(surf))), zlim=range(surf),
    axis.args=list(at=seq(min(surf),max(surf),length=10),
    labels=c(NA,rev(seq(150,950,length=9)))))
par(fig=c(0,10,1,4.5)/10)
par(new=T)
mtext('[hPa]', side=3, line=1, at=2019.5, cex=1.5)
dev.off()


aa=plot(mod1, select=1,scheme=2,n2=150,too.far=0.2)
surf=matrix(aa[[1]]$fit,150) +unlist(summary(mod1)$p.coeff[1]) 
LL=aa[[1]]$x
ll=aa[[1]]$y
ll=rev(approx(x=pr$Alt,y=pr$Press,xout=ll)$y)
png(paste(dir,"o3_trinidad_0_15km_tsseason_pr.png", sep=""), height=10, width=14, units="in", res=500)
par(mar=c(3.5, 3.5, 2, 1), mgp=c(2.4, 0.8, 0), las=1)
par(fig=c(0,10,4,10)/10)
image.plot(LL,ll,surf, xlab="", ylab="Pressure [hPa]", horizontal=F, col=viridis(20), xaxt="n",
	breaks=seq(20, 180, length.out=21), ylim=c(130,1000), cex.lab=1.5, yaxt="n")
par(fig=c(0,10,4,10)/10)
par(new=T)
contour(LL,ll,surf,drawlabels=T, add=TRUE)
axis(1,at=seq(1,12,by=1),labels=seq(1,12,by=1))
axis(2,at=seq(150,1000,by=50),labels=rev(seq(150,1000,by=50)))
mtext('[ppb]', side=3, line=0, at=12.4, cex=1.5)

par(fig=c(0,10,1,4.5)/10)
par(new=T)
matplot(LL, surf[,ncol(surf):1], xlab="Year", ylab="Climatology [ppb]", cex.lab=1.5,
    xlim=c(1.35,11.65), col=rainbow(ncol(surf)), type="l")
axis(1,at=seq(1,12,by=1),labels=seq(1,12,by=1))
image.plot(legend.only=T, nlevel=ncol(surf), col=rev(rainbow(ncol(surf))), zlim=range(surf),
    axis.args=list(at=seq(min(surf),max(surf),length=10),
    labels=c(NA,rev(seq(150,950,length=9)))))
par(fig=c(0,10,1,4.5)/10)
par(new=T)
mtext('[hPa]', side=3, line=1, at=12.4, cex=1.5)
dev.off()

