# This was originally a copy of Jim's 'fiaAllRegionsAnalysis.r' script. Now, the first few hundred lines are my own more or less garbled testbed for snippets that later got copied out to their own scripts. Keeping, because many chunks (mine or Jim's) could be useful later, but don't expect this to actually run as a script. 


rm(list=ls())

library(ade4)
library(fpc)

setwd('/Work/Research/Macrosystems/Clarkmodel')

source('./Scripts/fn/clarkRcppFunctions.r')
source('./Scripts/fn/demFunctions.r')
source('./Scripts/fn/clarkFunctions.R')
source('./Scripts/fn/FIAfunctions.r')


load('./Data/allData_demPlots.Rdata')
regSections <- c('SW','NE','SE','MW')
nsecs <- length(regSections)

# demographic data: n by nsnb (size-species)(nbreak, nspec)
# reproDatAll,reproEstAll,phiEstAll,ytotDataAll,growMeanDatAll,growMeanEstAll,
# growSdDatAll,survDatAll,survEstAll
# breaks,increments are in breaksTot, dxTot



# savefolder <- 'easternUS/'
# 
# 
# savefolder <- '/md2/clark/pdeOutput/'


# maplon <- c(-96.5,-67)       #location by site
# maplat <- range(lonLatAll[,2] )
# topo <- subETOPO5(maplon,maplat)


# xmap <- c(-100,-84.5,-65)    #define four sections lat/lon
# ymap <- c(25,37,49)



# tmp        <- matrix( unlist( strsplit(rownames(bioDatAll),'_ereg_') ),ncol=2,byrow=T)
# plotEcoreg <- tmp[,2]
# #tmp        <- matrix( unlist( strsplit(tmp[,1],'-') ),ncol=2,byrow=T)
# #plotReg    <- tmp[,1]
# plotName   <- tmp[,2]
# 
# ecoRegs <- sort(unique(plotEcoreg))
# 
# 
# domainNames <- read.table('./extra/fiaProvinces.txt',header=T)
# 
# wname <- match(plotEcoreg,domainNames[,'code'])
# ename <- as.character( domainNames[wname,'short'] )
# ename[is.na(ename)] <- plotEcoreg[is.na(ename)]
# 
# domains <- sort(unique(ename))
# domainIndex <- match(ename,domains)
# 
# nname <- nchar(ename)
# 
# 
# ncol   <- 100
# colF   <- colorRampPalette(c('wheat','orange','brown','black','black'))

######## indices

ii <- matrix( unlist( strsplit(colnames(ytotDataAll),'-') ),ncol=2,byrow=T)
bk <- as.numeric(ii[,2])

specs  <- unique(ii[,1])
breaks <- unique(bk)
sdex   <- match(ii[,1],specs)   #species index columns
kdex   <- match(bk,breaks)      #size species columns

nspec <- length(specs)
nb    <- length(breaks)



####### demography: growMeanDatAll,reproEstAll,phiEstAll,
####### growSdDatAll,survDatAll,survEstAll

#presTreeNCAll, presTreeAll,p50 = presPopAll

# dim(ytotDataAll) = Plots X (spp*size)

nall <- nrow(ytotDataAll)          #no. plots
nsnb <- ncol(ytotDataAll)          # Sp X # Size class
jj   <- as.vector( matrix( c(1:nall),nall,nsnb) )      #plot index
ss   <- as.vector( matrix( sdex,nall,nsnb,byrow=T) )   #species index


# Works: Sum over all size classes
yall <- byFunctionRcpp(as.vector(ytotDataAll),jj,ss,
                       matrix(0,nall,nspec),matrix(0,nall,nspec),
                       MEAN=T)



# Try: Aggregate size classes
breaks.new = c(0,10,50,100)
nb.new     = length(breaks.new)
kdex.new   = findInterval( kdex, breaks.new )



y.i = ytotDataAll[1,]



spXsize.i = byFunctionRcpp(y.i, sdex, kdex.new,
                       matrix(0, nspec, nb.new), MEAN=T)


# Check byFunctionRcpp is doing what I think
#   spXsize.i = byFunctionRcpp(y.i, sdex, kdex.new,
#                          matrix(0, nspec, nb.new), MEAN=T)
#   
#   spXsize.i2 = matrix(0, nspec, nb.new)
#   tots       = matrix(0, nspec, nb.new)
#   for(i in 1:length(y.i)) {
#     spXsize.i2[ sdex[i], kdex.new[i] ] = spXsize.i2[ sdex[i], kdex.new[i] ] + y.i[i]
#     tots[ sdex[i], kdex.new[i] ] = tots[ sdex[i], kdex.new[i] ] + 1
#   }
#   spXsize.i2 = spXsize.i2 / tots # Convert to mean
#   
#   spXsize.i2 - spXsize.i






#  growth rates for demogrphy plots (idem) where there are trees

ndem <- length(idem)

sindex <- as.vector( matrix( sdex,ndem, nsnb,byrow=T))
jindex <- as.vector( matrix( 1:ndem,ndem,nsnb) )

# growDemPlots[growDemPlots <= 0] <- NA

# Same as before, this just aggregates over size classes
growMeanEst_plXsp <- byFunctionRcpp( as.vector(growMeanEstAll),jindex,sindex,
                       matrix(0,ndem,nspec),matrix(0,ndem,nspec),MEAN=T)
                       
                       
tmp2 <- byFunctionRcpp( as.vector(growDemPlots)*0 + 1,jindex,sindex,
                       matrix(0,ndem,nspec),matrix(0,ndem,nspec),MEAN=T)

meanGrowth <- tmp1/tmp2


#y0 <- yall
#y0[yall > 0]  <- 0
#y0[yall == 0] <- 1
#colnames(y0) <- specs

mapDemogr(nn=ytotDataAll,demName='reproEst',
          specI = 'liquStyr',
          zrange=c(0,.1),PRES=T)

mapDemogr(nn=ytotDataAll,demName='fecnEst',
          specI = 'liquStyr',
          zrange=c(0,.1),PRES=T)

mapDemogr(nn=ytotDataAll,demName='growEst',
          specI = 'liriTuli', sizeRange = c(20,200),
          zrange=c(0,.4),PRES=T)

mapDemogr(nn=ytotDataAll,demName='growDat',
          specI = 'liriTuli', sizeRange = c(20,200),
          zrange=c(0,.4),PRES=T)

mapDemogr(nn=ytotDataAll,demName='survEst',
          specI = 'querRubr', sizeRange = c(20,200),
          zrange=c(0,.1),PRES=T)

mapDemogr(nn=ytotDataAll,demName='survDat',
          specI = 'acerRubr', sizeRange = c(20,200),
          zrange=c(0,1),PRES=T)

mapDemogr(nn=ytotDataAll,demName='survEst',
          specI = 'acerRubr', sizeRange = c(20,200),
          zrange=c(0,.06),PRES=T)

mapDemogr(nn=ytotDataAll,demName='reproEst',
          specI = 'liriTuli',
          zrange=c(0,.1),PRES=T)


specj <- rep('acerRubr',4)
demj  <- c('fecnEst','reproEst','growEst','survEst')
zj    <- c(.008,.1,.2,.04)
pj    <- c(.4,.4,.5,.5)
mapx  <- maplon
mapy  <- maplat
insetPos <- NULL

specj <- rep('querAlba',4)
demj  <- c('fecnEst','reproEst','growEst','survEst')
zj    <- c(.008,.1,.2,.04)
pj    <- c(.3,.3,.5,.5)
mapx  <- maplon
mapy  <- maplat
insetPos <- NULL


specj <- rep('liquStyr',4)
demj  <- c('fecnEst','reproEst','growEst','survEst')
zj    <- c(.008,.1,.2,.04)
pj    <- c(.3,.3,.5,.5)
mapx  <- c(-95,-74)
mapy  <- c(26,40)
insetPos <- NULL


# compare demographic rates for a species
specj <- rep('pinuTaed',4)
demj  <- c('fecnEst','reproEst','growEst','survEst')
zj    <- c(.01,.1,.2,.04)
pj    <- c(.4,.4,.5,.5)



# compare demographic rates for a species
specj <- rep('querPrin',4)
demj  <- c('fecnEst','reproEst','growEst','survEst')
zj    <- c(.005,.1,.2,.04)
pj    <- c(.3,.3,.5,.5)
mapx  <- c(-94,-66)
mapy  <- c(32,45)
insetPos <- NULL


graphics.off()
par(mfrow=c(2,2),mar=c(.1,.1,.1,.1))

for(j in 1:4){

   mapDemogr(nn=ytotDataAll,demName=demj[j],
             specI = specj[j], presProb=pj[j],
             zrange=c(0,zj[j]),mapscale=20,mapx=mapx,mapy=mapy,PRES=T,insetPos=insetPos)
   
   tt <- paste(specj[j],demj[j])
   title(tt,adj=0,line=1)
}

fname <- paste('demMap',specj[1],demj[1],'.pdf',sep='_')
dev.copy2pdf(file=outFile(savefolder,fname ))


# compare mortality rates across size classes

specj <- rep('acerRubr',4)
demj  <- c('survEst','survEst','survEst','survEst')
size <- rbind(c(0,10),c(10,20),c(20,30),c(30,1000))
zj    <- c(.06,.06,.06,.06)     # upper range for values
pj    <- c(.5,.5,.5,.5)         # Pr present

graphics.off()
par(mfrow=c(2,2),mar=c(.1,.1,.1,.1))

for(j in 1:4){
  
  mapDemogr(nn=ytotDataAll,demName=demj[j],sizeRange=size[j,],
            specI = specj[j], presProb=pj[j],
            zrange=c(0,zj[j]),mapscale=10,PRES=T)
  
  tt <- paste(specj[j],demj[j])
  title(tt,adj=0)
}

fname <- paste('demMap',specj[1],demj[1],'.pdf',sep='_')
dev.copy2pdf(file=outFile(savefolder,fname ))


# compare growth rates across size classes

specj <- rep('pinuTaed',4)
demj  <- c('growEst','growEst','growEst','growEst')
size <- rbind(c(0,10),c(10,20),c(20,30),c(30,1000))
zj    <- c(.4,.4,.4,.4)
pj    <- c(.5,.5,.5,.5)

graphics.off()
par(mfrow=c(2,2),mar=c(.1,.1,.1,.1))

for(j in 1:4){
  
  mapDemogr(nn=ytotDataAll,demName=demj[j],sizeRange=size[j,],
            specI = specj[j], presProb=pj[j],
            zrange=c(0,zj[j]),mapscale=6,PRES=T)
  
  tt <- paste(specj[j],demj[j])
  title(tt,adj=0)
}

fname <- paste('demMap',specj[1],demj[1],'.pdf',sep='_')
dev.copy2pdf(file=outFile(savefolder,fname ))


#oak recruitment

specj <- c('querVirg','querAlba','querPrin','querRubr')
demj  <- c('reproEst','reproEst','reproEst','reproEst')
zj    <- c(.1,.1,.1,.1)
pj    <- c(.3,.3,.3,.3)

par(mfrow=c(2,2),mar=c(1,1,1,1))

for(j in 1:4){
  
  mapDemogr(nn=ytotDataAll,demName=demj[j],
            specI = specj[j], presProb=pj[j],
            zrange=c(0,zj[j]),mapscale=10,PRES=T)
  
  tt <- paste(specj[j],demj[j])
  title(tt,adj=0)
}

fname <- paste('demMap',specj[1],demj[1],'.pdf',sep='_')
dev.copy2pdf(file=outFile(savefolder,fname ))



#oak growth

specj <- c('querAlba','querAlba','querRubr','querRubr')
demj  <- c('growDat','growEst','growDat','growEst')
zj    <- c(.3,.3,.3,.3)
pj    <- c(.5,.5,.5,.5)

par(mfrow=c(2,2),mar=c(1,1,1,1))

for(j in 1:4){
  
  mapDemogr(nn=ytotDataAll,demName=demj[j],
            specI = specj[j], presProb=pj[j],
            zrange=c(0,zj[j]),mapscale=8,PRES=T)
  
  tt <- paste(specj[j],demj[j])
  title(tt,adj=0)
}

fname <- paste('demMap',specj[1],demj[1],'.pdf',sep='_')
dev.copy2pdf(file=outFile(savefolder,fname ))



#pinuTaed growth

specj <- c('pinuTaed','pinuTaed',
           'pinuEchi','pinuEchi',
           'pinuPalu','pinuPalu')
demj  <- c('survEst','growEst',
           'survEst','growEst',
           'survEst','growEst')
zj    <- c(.04,.4,.04,.4,.04,.4)
pj    <- c(.5,.5,.5,.5,.5,.5)

par(mfrow=c(3,3),mar=c(1,1,1,1))

for(j in 1:6){
  
  mapDemogr(nn=ytotDataAll,demName=demj[j],
            specI = specj[j], presProb=pj[j],
            zrange=c(0,zj[j]),mapscale=8,PRES=T,
            mapx=c(-96,-75),mapy=c(28,39))
  
  tt <- paste(specj[j],demj[j])
  title(tt,adj=0)
}

fname <- paste('demMap',specj[1],demj[1],'.pdf',sep='_')
dev.copy2pdf(file=outFile(savefolder,fname ))

###########which stage limiting?

specj <- c('querAlba','querAlba','querAlba','querAlba','querAlba')
demj  <- c('fecnEst','reproDat','growEst','survEst')
size <- rbind(c(0,20),c(0,20),c(0,20),c(0,20))
zj    <- c(.03,.03,.3,.05)
pj    <- c(.3,.3,.3,.3)

par(mfrow=c(2,2),mar=c(1,1,1,1))

for(j in 1:4){
  
  mapDemogr(nn=ytotDataAll,demName=demj[j],sizeRange=size[j,],
            specI = specj[j], presProb=pj[j],
            zrange=c(0,zj[j]),mapscale=8,PRES=T)
  
  tt <- paste(specj[j],demj[j])
  title(tt,adj=0)
}

fname <- paste('demMap',specj[1],demj[1],'.pdf',sep='_')
dev.copy2pdf(file=outFile(savefolder,fname ))

###################### KL divergence on density #######################

#predG2Y <- gammaMuA*dxMat*sampArea
#predG2Y <- gamPredA*dxTot[,kdex]*sampArea


jj   <- as.vector(matrix( c(1:nall),nall,nsnb))
ss   <- as.vector(matrix( sdex,nall,nsnb,byrow=T))


yall <- byFunctionRcpp(as.vector(ytotDataAll),jj,ss,
                       matrix(0,nall,nspec),matrix(0,nall,nspec),
                       MEAN=T)
gall  <- byFunctionRcpp(as.vector(gamPredA),jj,ss,
                        matrix(0,nall,nspec),matrix(0,nall,nspec),
                        MEAN=T)

klMat <- matrix(NA,nall,nspec)
colnames(klMat) <- specs

maxClass <- 40               #no. of size classes for comparison
agClass  <- round(seq(1,maxClass,length=nb),0)


wj <- as.vector( matrix(c(1:nall),nall,nb) )
wk <- as.vector( matrix(agClass,nall,nb,byrow=T) )
mm <- matrix(0,nall,nb)

for(s in 1:nspec){
  
  ws <- which(sdex == s)
  
  ys <- as.vector(ytotDataAll[,ws])
  gs <- as.vector(gamPredA[,ws])
  
  gsums <- byFunctionRcpp(as.vector(gamPredA[,ws]),wj,wk,
                          mm,mm,MEAN=T)
  ysums <- byFunctionRcpp(as.vector(ytotDataAll[,ws]),wj,wk,
                          mm,mm,MEAN=T)
  
  tmp <- KLdivergence(gsums,ysums,NORM=F)
  
  klMat[,s] <- tmp$KL
  klMat[tmp$zero,s] <- NA
  print(s)
}

klMat[klMat == 0] <- NA

#if(min(klMat,na.rm=T) < 0)klMat <- klMat + -min(klMat,na.rm=T) + .001 

#klMat[!is.finite(klMat)] <- 0


graphics.off()

klPlot <- rowSums(klMat,na.rm=T)
if(min(klPlot,na.rm=T) < 0)klPlot <- klPlot + -min(klPlot) + .001   #positive values for mapping


xx <- range(topo$x)
yy <- range(topo$y)
xy <- c(xx[2] - diff(xx)/10, yy[2] - diff(yy)/10)
 

#ii <- list(eco = ecoReg4,type = c(rep(2,nplotDem),rep(1,nplotFIA)))
mkl <- byIndex(klPlot,domainIndex,mean)
skl <- byIndex(klPlot,domainIndex,sd)
#colnames(mkl) <- colnames(skl) <- c('FIA','LTFD')

#scin <- sort(unique(ecoReg4))
#scol <- match(rownames(mkl),scin)

okl <- order(mkl)

eregFit <- domains[okl]

colF   <- colorRampPalette(c('wheat','tan','brown','black'))  #order by fit
colq   <- colF(length(domains))

scol2   <- colq[match(domains,eregFit)]


#if(!REMOTE){
plotfile <- outFile(savefolder,'KLbyEcoReg.pdf')
#  plotstart(plotfile,REMOTE)

#enames <- ename[match(rownames(mkl),ecoDomains)]
par(mfrow=c(1,1),bty='n')
myBoxPlot(mu = mkl,
          limits = list(cbind(mkl - skl, mkl + skl),
                        cbind(mkl - 1.96*skl, mkl + 1.96*skl) ),
          myorder=okl,boxcol=colq,
          xlab=' ',ylab='K-L divergence',
          ytic=seq(-2,4,by=1),plabels=eregFit,add=F,widthScale=.5)
#    myBoxPlot(mkl[,2],(cbind(mkl[,2] - skl[,2], mkl[,2] + skl[,2]) ),
#          myorder=o,boxcol='red',add=T)
title('+/- 1 sd')
abline(h=0,lty=2)


plotfile <- outFile(savefolder,'KLByY.pdf')

graphics.off()

if(min(klPlot,na.rm=T) < 0)klPlot <- klPlot + -min(klPlot,na.rm=T) + .001 
klPlot[is.na(klPlot)] <- 0

qmax <- quantile(klPlot,.95,na.rm=T)

ssk  <- seq(0,qmax+1,length=20)
klPlot[klPlot > qmax] <- qmax
scol <- findInterval(klPlot,ssk)
scol <- colq[scol]

regMap(topo$x,topo$y,topo$z,IMAGE=F,mapscale=5)
#regMap(lonLatAll[,1],lonLatAll[,2],z=klPlot,ADD=T,bg=scol,fg=scol)

values2contour(x=lonLatAll[,1],y=lonLatAll[,2],
               z=klPlot,nx=50,ny=50,col=scol,
               zlevs=ssk,lwd=1,add=T,lty=2,fill=T)

mapMask(lonLatAll[,1],lonLatAll[,2],dx=.1,dy=.16,whiteOutDist=.15,col='white')

map('state',add=T,col='white',lwd=6,interior=F)
map('state',add=T,col='grey',lwd=3,interior=T)


title('KL divergence size-species-space-time')

rect(maplon[1] + diff(maplon)*.65,
     maplat[1] + diff(maplat)*.002,
     maplon[1] + diff(maplon)*.99,
     maplat[1] + diff(maplat)*.2,col='white',border='white')

colorLegend(maplon[1] + diff(maplon)*c(.92,.97),
            maplat[1] + diff(maplat)*c(.02,.15),
            scale=c(1:20),cols=rev(colq),bg='white',
            labside='left',endLabels=c('worst','best'))


dev.copy2pdf(file=outFile(savefolder,'KLmap1.pdf'))


##############species richness (of dominant species)

y0 <- treeAll*0
y0[treeAll == 0] <- 1
y0[treeAll > 0] <- 0

xmap <- c(-100,-84.5,-65)    #define four sections lat/lon
ymap <- c(25,37,49)

specsByReg <- vector( mode = 'list',length=4)
names(specsByReg) <- regSections

for(j in 1:2){
  mlon <- xmap[c(j,j+1)]
  for(k in 1:2){
    mlat <- ymap[c(k,k+1)]
    if(j == 1 & k == 1)ss <- 'SW'
    if(j == 1 & k == 2)ss <- 'MW'
    if(j == 2 & k == 1)ss <- 'SE'
    if(j == 2 & k == 2)ss <- 'NE'
    tmp <- specsInLonLat(mlon,mlat,lonLatAll[,1],lonLatAll[,2],1 - y0,minPlot=30)
    specsByReg[names(specsByReg) == ss] <- list( names(tmp) )
  }
}


richness <- rowSums(1 - y0)

qmax <- quantile(richness,.95,na.rm=T)

ssk  <- seq(0,qmax+1,length=20)
richness[richness > qmax] <- qmax
scol <- findInterval(richness,ssk)
scol <- colq[scol]

regMap(topo$x,topo$y,topo$z,IMAGE=F,mapscale=5)

values2contour(x=lonLatAll[,1],y=lonLatAll[,2],
               z=richness,nx=70,ny=70,col=scol,
               zlevs=ssk,lwd=1,add=T,lty=2,fill=T)

mapMask(lonLatAll[,1],lonLatAll[,2],dx=.1,dy=.16,whiteOutDist=.15,col='white')

map('state',add=T,col='white',lwd=6,interior=F)
map('state',add=T,col='grey',lwd=3,interior=T)


title('richness of abundant species')

rect(maplon[1] + diff(maplon)*.65,
     maplat[1] + diff(maplat)*.002,
     maplon[1] + diff(maplon)*.99,
     maplat[1] + diff(maplat)*.2,col='white',border='white')

colorLegend(maplon[1] + diff(maplon)*c(.92,.97),
            maplat[1] + diff(maplat)*c(.02,.15),
            scale=c(1:20),cols=rev(colq),bg='white',
            labside='left',endLabels=c('low','high'))


dev.copy2pdf(file=outFile(savefolder,'speciesRichness.pdf'))


#######  biomass: bioMuAllA,bioVrAllA,bioMuPredA,bioVrPredA #######################

graphics.off()

par(mfrow=c(2,2),bty='n')

bioMu <- rowSums(bioMuAllA)
bioSd <- sqrt(rowSums(bioVrAllA))
bioDat <- rowSums(bioDatAll,na.rm=T)

hist(bioDat,nclass=100,xlab='Aboveground biomass (T/ha)',
     ylim=c(0,.01),main='',probability=T,col='grey',border='grey')
tmp <- hist(bioMu,nclass=100,plot=F)
lines(tmp$mids,tmp$density,type='s',lwd=6,col='white')
lines(tmp$mids,tmp$density,type='s',lwd=3,col='black')

dev.copy2pdf(file=outFile(savefolder,'bmassHist.pdf'))


byFunctionRcpp(bioDat,domainIndex,rep(1,nall),
               matrix(length(domains),ncol=1),matrix(length(domains),ncol=1),MEAN=F)
               


rr <- "Coastal_Plain"
add.scatter(hist2Add(bioDat[ename == rr],xlim=c(0,800)),posi='bottomright',
            ratio=.15,inset=c(.02,.1),bg.col='wheat')

rr <- "E_Broadleaf_Cont"
add.scatter(hist2Add(bioDat[ename == rr],xlim=c(0,800)),posi='topright',
            ratio=.15,inset=c(.02,.01),bg.col='wheat')


zlevs <- c(seq(0,300,by=10),10000)
regMap(topo$x,topo$y,topo$z,IMAGE=F,mapscale=10,lineCol='grey')

values2contour(lonLatAll[,1],lonLatAll[,2],
               bioMu,nx=60,ny=60,col=colF(length(zlevs)),
               zlevs=zlevs,lwd=2,add=T,drawlabels=F,fill=T)

mapMask(lonLatAll[,1],lonLatAll[,2],dx=.1,dy=.1,whiteOutDist=.3,col='white')

map('state',add=T,col='white',lwd=6,interior=F)
map('state',add=T,col='grey',lwd=3,interior=T)

add.scatter(hist2Add(bioMu,xlim=c(0,400)),posi='bottomright',
            ratio=.12,inset=c(.02,.1),bg.col='grey')
              

zlevs <- c(seq(0,100,by=1),1000)

regMap(topo$x,topo$y,topo$z,IMAGE=F,mapscale=10,lineCol='grey')

values2contour(lonLatAll[,1],lonLatAll[,2],
               bioSd,nx=60,ny=60,col=colF(length(zlevs)),
               zlevs=zlevs,lwd=2,add=T,drawlabels=F,fill=T)

mapMask(lonLatAll[,1],lonLatAll[,2],dx=.1,dy=.1,whiteOutDist=.4,col='white')

map('state',add=T,col='white',lwd=6,interior=F)
map('state',add=T,col='grey',lwd=3,interior=T)

add.scatter(hist2Add(bioSd,xlim=c(0,120)),posi='bottomright',
            ratio=.12,inset=c(.02,.1),bg.col='grey')

dev.copy2pdf(file=outFile(savefolder,'tonsPerHaEst.pdf'))






              
              

######################climate maps



graphics.off()

par(mfrow=c(2,2),mar=c(.5,.5,.5,.5),bty='n')

#surplus
xr <- maplon[1] + diff(maplon)*c(.68,.92)
yr <- maplat[1] + diff(maplat)*c(.03,.27)
xx <- maplon[1] + diff(maplon)*c(.85,.9)
yy <- maplat[1] + diff(maplat)*c(.07,.2)

climateMap(maplon,maplat,lonLatAll[,1],lonLatAll[,2],topo,climVec=climAll[,'therm'],
           zlevs=seq(min(climAll[,'therm'])+50,2000,by=100),
           colorRamp=c('wheat','green','darkgreen','black'),nx=100,ny=100,
           xbox=cbind(xr,yr),xscale=cbind(xx,yy),
           labside='left',legside='bottomright',endLabels=c('200','2000'),
           mapscale=10)
title('a) Hydrothermal Surplus',adj=0,line=1)

#dev.copy2pdf(file='thermMap.pdf')


#deficit
climateMap(maplon,maplat,lonLatAll[,1],lonLatAll[,2],topo,
           climVec=climAll[,'deficit'],
           zlevs=seq(min(climAll[,'deficit'])+100,4000,by=200),
           colorRamp=c('wheat','orange','brown','red','black'),nx=100,ny=100,
           xbox=cbind(xr,yr),xscale=cbind(xx,yy),
           labside='left',legside='bottomright',endLabels=c('100','4000'),
           mapscale=10)
#dev.copy2pdf(file='deficitMap.pdf')
title('b) Hydrothermal Deficit',adj=0,line=1)

#winter temp
climateMap(maplon,maplat,lonLatAll[,1],lonLatAll[,2],topo,climVec=climAll[,'temp'],
           zlevs=seq(-10,30,by=2),
           colorRamp=rev(c('azure','darkblue')),nx=100,ny=100,
           xbox=cbind(xr,yr),xscale=cbind(xx,yy),
           labside='left',legside='bottomright',endLabels=c('-10','20'),
           mapscale=10)
#dev.copy2pdf(file='winterTempMap.pdf')
title('c) Winter temperature',adj=0,line=1)

ncol   <- 30
colF   <- colorRampPalette(c('red','green','blue'))
cols   <- colF(ncol)
colh   <- colF(3)

colm <- cols[findInterval(climAll[,'moisture'],seq(-1.1,1.1,length=ncol))]

regMap(topo$x,topo$y,topo$z,IMAGE=F,mapscale=10)
points(lonLatAll[,1],lonLatAll[,2],col=colm,cex=.15)

legend('bottomright',c('xeric','intermediate','mesic'),
       text.col=c('red','green','blue'),bty='n')

title('d) Moisture classes',adj=0,line=1)
dev.copy2pdf(file='climateMap.pdf')

tp <- cbind(climAll[,'temp'],climAll[,'prec'])
tphull <- chull( tp )
tphull <- tp[tphull,]


td <- cbind(climAll[,'therm'],climAll[,'deficit'])         #climAll is climate by location
ww <- which(is.finite(td[,1]) & is.finite(td[,2]))
tdhull <- chull( td[ww,] )
tdhull <- td[ww,][tdhull,]

specs <- colnames(presTreeNCAll)
nspec <- length(specs)

presenceMat <- numeric(0)

tpx <- seq(min(tphull[,1]),max(tphull[,1]),length=30)
tpy <- seq(min(tphull[,2]),max(tphull[,2]),length=30)
tpGrid <- as.matrix( expand.grid(tpx,tpy) )

tpx <- seq(min(tdhull[,1]),max(tdhull[,1]),length=50)
tpy <- seq(min(tdhull[,2]),max(tdhull[,2]),length=50)
tdGrid <- as.matrix( expand.grid(tpx,tpy) )


for(s in 1:nspec){
  
  print(s)
  
  graphics.off()
  
  plotfile <- outFile( savefolder,paste('climateSpace_',specs[s],'.pdf',sep='') )
  
#  plotstart(plotfile,REMOTE=F)
  
  par(mfrow=c(2,2),bty='n')
  plot(0,0,xlim=c(-20,25),ylim=c(200,2800),xlab='Temperature C',
       ylab='Precipitation (mm)',cex=.01)
  
  title(specs[s])
  
  mt <- which(presTreeNCAll[,s] > .5)
  if(length(mt) < 5)next
  symbols(tp[mt,1],tp[mt,2],squares=mt*0+.2,add=T,fg='orange',bg='orange',inches=F)
  
  ps <- points2grid(xx=tp[mt,1],yy=tp[mt,2],grid=tpGrid)$gridFraction
  ns <- 'treeNC'
  
  mt <- which(presTreeAll[,s] > .5)
  symbols(tp[mt,1],tp[mt,2],squares=mt*0+.2,add=T,fg='tan',bg='tan',inches=F)
  
  mm <- points2grid(xx=tp[mt,1],yy=tp[mt,2],grid=tpGrid)$gridFraction
  ps <- c(ps, mm)
  ns <- c(ns,'tree')
  
  mt <- which(presPopAll[,s] > .5)
  symbols(tp[mt,1],tp[mt,2],squares=mt*0+.2,add=T,fg='brown',bg='brown',inches=F)
  
  mm <- points2grid(xx=tp[mt,1],yy=tp[mt,2],grid=tpGrid)$gridFraction
  ps <- c(ps, mm)
  ns <- c(ns,'popn')
  
  mt <- which(treeAll[,s] > 0)
  symbols(tp[mt,1],tp[mt,2],squares=mt*0+.2,add=T,fg='red',bg='red',inches=F)
  
  mm <- points2grid(xx=tp[mt,1],yy=tp[mt,2],grid=tpGrid)$gridFraction
  ps <- c(ps, mm)
  ns <- c(ns,'obsTree')
  
  polygon(tphull[,1],tphull[,2],border='white',lwd=5)
  polygon(tphull[,1],tphull[,2],border='grey',lwd=3,lty=2)
  
  legend('topleft',c('w/o competition','trees','population','observed'),bty='n',
         text.col=c('orange','tan','brown','red') )
  
  
  plot(0,0,xlim=c(-20,25),ylim=c(200,2800),xlab='Temperature C',
       ylab=' ',cex=.01)
  mt <- which(presTreeNCAll[,s] > .5)
  symbols(tp[mt,1],tp[mt,2],squares=mt*0+.2,add=T,fg='orange',bg='orange',inches=F)
  
  mt <- which(presFecAll[,s] > .5)
  symbols(tp[mt,1],tp[mt,2],squares=mt*0+.2,add=T,fg='lightblue',bg='lightblue',inches=F)
  
  mm <- points2grid(xx=tp[mt,1],yy=tp[mt,2],grid=tpGrid)$gridFraction
  ps <- c(ps,mm)
  ns <- c(ns,'fecn')
  
  mt <- which(presReprAll[,s] > .5)
  symbols(tp[mt,1],tp[mt,2],squares=mt*0+.2,add=T,fg='darkblue',bg='darkblue',inches=F)
  
  mm <- points2grid(xx=tp[mt,1],yy=tp[mt,2],grid=tpGrid)$gridFraction
  ps <- c(ps,mm)
  ns <- c(ns,'repro')
  
  mt <- which(reproAll[,s] > 0)
  symbols(tp[mt,1],tp[mt,2],squares=mt*0+.2,add=T,fg='red',bg='red',inches=F)
  
  mm <- points2grid(xx=tp[mt,1],yy=tp[mt,2],grid=tpGrid)$gridFraction
  ps <- c(ps,mm)
  ns <- c(ns,'obsRepro')
  
  polygon(tphull[,1],tphull[,2],border='white',lwd=5)
  polygon(tphull[,1],tphull[,2],border='grey',lwd=3,lty=2)
  
  legend('topleft',c('fecundity','recruitment','observed'),bty='n',
         text.col=c('lightblue','darkblue','red') )
  

  
   ############################################
  plot(0,0,xlim=c(0,2500),ylim=c(0,4000),xlab='Surplus',
       ylab='Deficit',cex=.01)
  
#  polygon(tdhull[,1],tdhull[,2],col='azure',border='grey',lwd=5)
  
  mt <- which(presTreeNCAll[,s] > .5)
  if(length(mt) < 5)next
  symbols(td[mt,1],td[mt,2],squares=mt*0+10,add=T,fg='orange',bg='orange',inches=F)
  
  ds <- points2grid(xx=td[mt,1],yy=td[mt,2],grid=tdGrid)$gridFraction
  
  mt <- which(presTreeAll[,s] > .5)
  symbols(td[mt,1],td[mt,2],squares=mt*0+10,add=T,fg='tan',bg='tan',inches=F)
  
  ds <- c(ds,points2grid(xx=td[mt,1],yy=td[mt,2],grid=tdGrid)$gridFraction)
  
  mt <- which(presPopAll[,s] > .5)
  symbols(td[mt,1],td[mt,2],squares=mt*0+10,add=T,fg='brown',bg='brown',inches=F)
  
  ds <- c(ds,points2grid(xx=td[mt,1],yy=td[mt,2],grid=tdGrid)$gridFraction)
  
  mt <- which(treeAll[,s] > 0)
  symbols(td[mt,1],td[mt,2],squares=mt*0+10,add=T,fg='red',bg='red',inches=F)
  
  ds <- c(ds,points2grid(xx=td[mt,1],yy=td[mt,2],grid=tdGrid)$gridFraction)
  
  polygon(tdhull[,1],tdhull[,2],border='white',lwd=5)
  polygon(tdhull[,1],tdhull[,2],border='grey',lwd=3,lty=2)
  
  legend('topright',c('w/o competition','trees','population','observed'),bty='n',
         text.col=c('orange','tan','brown','red') )
  
  plot(0,0,xlim=c(0,3000),ylim=c(0,4000),xlab='Surplus',
       ylab=' ',cex=.01)
  
  mt <- which(presTreeNCAll[,s] > .5)
  symbols(td[mt,1],td[mt,2],squares=mt*0+10,add=T,fg='orange',bg='orange',inches=F)
  
  mt <- which(presFecAll[,s] > .5)
  symbols(td[mt,1],td[mt,2],squares=mt*0+10,add=T,fg='lightblue',bg='lightblue',inches=F)
  
  ds <- c(ds,points2grid(xx=td[mt,1],yy=td[mt,2],grid=tdGrid)$gridFraction)
  
  mt <- which(presReprAll[,s] > .5)
  symbols(td[mt,1],td[mt,2],squares=mt*0+10,add=T,fg='blue',bg='blue',inches=F)
  
  ds <- c(ds,points2grid(xx=td[mt,1],yy=td[mt,2],grid=tdGrid)$gridFraction)
  
  mt <- which(reproAll[,s] > 0)
  symbols(td[mt,1],td[mt,2],squares=mt*0+10,add=T,fg='red',bg='red',inches=F)
  
  ds <- c(ds,points2grid(xx=td[mt,1],yy=td[mt,2],grid=tdGrid)$gridFraction)
  
  polygon(tdhull[,1],tdhull[,2],border='white',lwd=5)
  polygon(tdhull[,1],tdhull[,2],border='grey',lwd=3,lty=2)
  
  legend('topright',c('fecundity','recruitment','observed'),bty='n',
         text.col=c('lightblue','darkblue','red') )
  
  names(ps) <- paste('tempPrec',ns,sep='-')
  names(ds) <- paste('surpDef',ns,sep='-')
  
  presenceMat <- rowBind(presenceMat,round(c(ps,ds),4),specs[s]) 
  
  print(s)
 # print(presenceMat[s,])
  
  dev.copy2pdf(file=plotfile )
}

presenceMat <- presenceMat[rownames(presenceMat) != 'other',]
write.table(presenceMat,file=paste(savefolder,'presenceMat.txt',sep=''),quote=F)

dev.off()

par(mfrow=c(1,1))
v1 <- 'tempPrec'
tmp <- presenceMat[,grep(v1,colnames(presenceMat))]
myorder <- order(tmp[,1],decreasing=T)

par(bty='n')
myBarPlot(tmp[myorder,],ylim=c(0,4),widthFactor=.5,
          barLabels=rownames(tmp)[myorder],textFactor=.7,stack=F)


tmp2 <- tmp/tmp[,1]
myorder <- order(tmp2[,2],decreasing=T)

par(bty='n')
tmp <- myBarPlot(tmp2[myorder,-1],ylim=c(0,7),widthFactor=.5,
                 xlab='Decreasing climate space',ylab='Fraction of range',
          barLabels=rownames(tmp2)[myorder],textFactor=.7,stack=F)
text(rep(nspec,length(tmp)),tmp-.3,
         c('individual','population','obs trees','fecundity','reproduction','obs repro'),
         col=c(1:length(tmp)),pos=2) 
title('fraction of temp/prec space')

plotfile <- outFile( savefolder,paste('tempPrecSpace.pdf',sep='') )
dev.copy2pdf(file=plotfile )


par(mfrow=c(1,1))
v1 <- 'surpDef'
tmp <- presenceMat[,grep(v1,colnames(presenceMat))]

tmp2 <- tmp/tmp[,1]
myorder <- order(tmp2[,2],decreasing=T)

par(bty='n')
tmp <- myBarPlot(tmp2[myorder,-1],ylim=c(0,7),widthFactor=.5,
                 xlab='Decreasing climate space',ylab='Fraction of range',
                 barLabels=rownames(tmp2)[myorder],textFactor=.7,stack=F)
text(rep(nspec,length(tmp)),tmp-.3,
     c('individual','population','obs trees','fecundity','reproduction','obs repro'),
     col=c(1:length(tmp)),pos=2) 
title('fraction of surplus/defict space')

plotfile <- outFile( savefolder,paste('surpDeficitSpace.pdf',sep='') )
dev.copy2pdf(file=plotfile )

tmp <- points2grid(xx=td[mt,1],yy=td[mt,2],grid=tdGrid)

xnames <- xnamesAll[[1]]
for(j in 2:nsecs){
  xnames <- unique( c(xnames,xnamesAll[[j]]) )
}


##############sensitivity

chainSens <- numeric(0)
vnames <- c("therm","temp-DecJanFebMar","deficit","plotBA","xeric","mesic")
mvars  <- c('betaChain','betaMortChain','betaFecChain','betaGamChain')

betaMats <- numeric(0)

for(j in 1:nsecs){         #regions
  
  ww <- which(names(specsByReg) == regSections[[j]])
  specKeep <- specsByReg[[ww]]
  specKeep <- specKeep[specKeep != 'other']
  
  chainReg  <- numeric(0)

  for(m in 1:length(mvars)){       #demography
    
  betaMu    <- numeric(0)
    
  parChain <- chainsAll[[j]][[mvars[m]]]
  specj    <- matrix( unlist( strsplit(colnames(parChain),'_') ),ncol=2,byrow=T)[,2]
  sall     <- unique(specj)
  
  sall <-  sall[sall %in% specKeep]
  
  xrangej <- xrangeAll[[j]]
  jmat    <- matrix(NA,nrow(parChain),length(sall)*length(vnames))
  cjn     <- as.vector(outer(vnames,sall,paste,sep='_'))
  colnames(jmat) <- cjn
  
  ks <- 0
  
  for(s in 1:length(sall)){
    
    pchain <- parChain[,grep(sall[s],colnames(parChain))]
    cc     <- colnames(pchain)
    colnames(pchain) <- matrix( unlist(strsplit(cc,'_')),ncol=2,byrow=T)[,1]
    
    betaMu <- cbind(betaMu,colMeans(pchain))
    
    for(k in 1:length(vnames)){
      
      ks <- ks + 1
      
      if(!vnames[k] %in% colnames(pchain))next
      
      jmat[,ks] <- par2sens(pchain,vname=vnames[k],xrangej,reverseSign)
    }
  }
  colnames(betaMu) <- sall
  chainReg <- append(chainReg,list(jmat))
  betaMats <- append(betaMats,list(betaMu))
  names(betaMats)[[length(betaMats)]] <- paste(regSections[j],mvars[m],sep='_')
  }
  names(chainReg) <- mvars
  chainSens <- append(chainSens,list(chainReg))
}

names(chainSens) <- regSections


chains <- sens <- character(0)
for(r in 1:nsecs){
  chains <- c(chains,names(chainsAll[[1]]) )
  sens <- c(sens,names(chainSens[[1]]) )
}
  
chains <- sort(unique(chains))
chains <- unlist( strsplit(chains,'Chain') )
sens   <- sort(unique(sens))
sens   <- unlist( strsplit(sens,'Chain') )

omitSpecs <- c("acerBarb","acerNegu","betuNigr","betuPopu","carpCaro","caryAqua",
               "magnGran","ostrVirg","other","pinuBank","platOcci","querFalc",
               "querMari","querMich","querNigr","querStel","sassAlbi")


################## Sensitivity


###############NE - SE

# xeric vs mesic growth

dev.off()

# beta - growth
# betaFec - fecundity
# betaMort - mortality


includeSections <- c('NE','SW')

includeSections <- c('SE','SW')

includeSections <- c('NE','SE')
           
includeSections <- c('NE','MW')

includeSections <- c('MW','SW')

#includeSections <- c('transectEast','transectWest')
#includeSections <- c('transectNorth','transectSouth')
           
# beta - growth
           # beta - growth          

resp <- c('beta','beta','beta')        
vars <- c('deficit','xeric','mesic')
nv   <- length(vars)
ht   <- c(1,.0005,.0025)     # ht of plots
xlo  <- c(-.00005,-.1,-.05)  # lo x
xhi  <- c(0,0,.1)     # hi x


resp <- c('beta','beta','beta')
vars <- c('deficit','therm','temp')
nv   <- length(vars)
ht   <- c(1,1,2)     # ht of plots
xlo  <- c(-.00005,-.00001,0)  # lo x
xhi  <- c(0,.00005,.01)     # hi x



resp <- c('beta','beta')
vars <- c('xeric','mesic')
nv   <- length(vars)
ht   <- c(.0005,.0025)     # ht of plots
xlo  <- c(-.2,-.05)  # lo x
xhi  <- c(.1,.1)     # hi x


rr <- c(unique(resp),unique(vars))
plotfile <- outFile( savefolder,
                     paste('postSens_',includeSections[1],includeSections[1]
                           ,'_',paste0(rr,collapse='-') ,'.pdf',sep='') )

par(mfrow=c(1,nv),bty='n')

ord <- NULL

kchains <- getRegionChains(includeSections,chainSens)

colF   <- colorRampPalette(c('blue','brown'))
colorcode   <- colF(length(kchains))


for(k in 1:nv){
  
  xtic <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=length(kchains),ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtic,textSize=1,COLCODE=colorcode,
                           reverseSign=NULL,omitSpecs=omitSpecs,FILL=T   )
  
  legend('topright',legend=tmp$percPos,text.col=colorcode,cex=1.0,
         border='white',box.col='white',xjust=1)
  legend('topleft',legend=tmp$percNeg,text.col=colorcode,cex=1.0,
         border='white',box.col='white')
  
  if(k == 1){
    ord <- tmp$specOrd
    legend('bottomleft',names(kchains),text.col=tmp$repSpecs,box.col='white',cex=1.4)
  }
}

dev.copy2pdf(file=plotfile )






resp <- c('beta','beta')
vars <- c('temp-DecJanFebMar','therm')
nv   <- length(vars)
ht   <- c(.9,2)
xlo  <- c(-.002,-.00001)
xhi  <- c(.01,.0001)


rr <- c(unique(resp),unique(vars))
plotfile <- outFile( savefolder,
                     paste('postSens_',includeSections[1],includeSections[2],
                           '_',paste0(rr,collapse='-') ,'.pdf',sep='') )

par(mfrow=c(1,nv),bty='n')

ord <- NULL

kchains <- getRegionChains(includeSections,chainSens)

for(k in 1:nv){
  
  xtic <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=length(kchains),ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtic,textSize=.7,COLCODE=colorcode,
                           reverseSign=NULL,omitSpecs=omitSpecs,FILL=T   )
  
  legend('topright',legend=tmp$percPos,text.col=colorcode,cex=.8,
         border='white',box.col='white',xjust=1)
  legend('topleft',legend=tmp$percNeg,text.col=colorcode,cex=.8,
         border='white',box.col='white')
  
  if(k == 1){
    ord <- tmp$specOrd
    legend('bottomleft',names(kchains),text.col=tmp$repSpecs,box.col='white',cex=.7)
  }
}

dev.copy2pdf(file=plotfile )



#############################################33

includeSections <- c('SW')
kchains <- getRegionChains(includeSections,chainSens)

resp <- c('beta','beta')
vars <- c('therm','plotBA')
nv   <- length(vars)
ht   <- c(.9,2)
xlo  <- c(-.00002,-.001)
xhi  <- c(.00005,.0005)


rr <- c(unique(resp),unique(vars))
plotfile <- outFile( savefolder,
                     paste('postSens_',includeSections[1],includeSections[2],
                           '_',paste0(rr,collapse='-') ,'.pdf',sep='') )

par(mfrow=c(1,nv),bty='n')

ord <- NULL

kchains <- getRegionChains(includeSections,chainSens)

for(k in 1:nv){
  
  xtic <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=length(kchains),ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtic,textSize=.7,COLCODE=colorcode,
                           reverseSign=NULL,omitSpecs=omitSpecs,FILL=T   )
  
  legend('topright',legend=tmp$percPos,text.col=colorcode,cex=.8,
         border='white',box.col='white',xjust=1)
  legend('topleft',legend=tmp$percNeg,text.col=colorcode,cex=.8,
         border='white',box.col='white')
  
  if(k == 1){
    ord <- tmp$specOrd
    legend('bottomleft',names(kchains),text.col=tmp$repSpecs,box.col='white',cex=.7)
  }
}

dev.copy2pdf(file=plotfile )


#includeSections <- c('transectEast','transectWest')
#includeSections <- c('transectNorth','transectSouth')

resp <- c('betaFec','betaFec')
vars <- c('temp-DecJanFebMar','therm')
nv   <- length(vars)
ht   <- c(.1,.1)
xlo  <- c(0,0)
xhi  <- c(.2,.001)

rr <- c(unique(resp),unique(vars))
plotfile <- outFile( savefolder,
                     paste('postSens_',includeSections[1],includeSections[2],
                           '_',paste0(rr,collapse='-') ,'.pdf',sep='') )


par(mfrow=c(1,nv),bty='n')

ord <- NULL

kchains <- getRegionChains(includeSections,chainSens)

for(k in 1:nv){
  
  xtic <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=length(kchains),ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtic,textSize=.7,COLCODE=colorcode,
                           reverseSign=NULL,omitSpecs=omitSpecs,FILL=T   )
  
  legend('topright',legend=tmp$percPos,text.col=colorcode,cex=.8,
         border='white',box.col='white',xjust=1)
  legend('topleft',legend=tmp$percNeg,text.col=colorcode,cex=.8,
         border='white',box.col='white')
  
  if(k == 1){
    ord <- tmp$specOrd
    legend('bottomleft',names(kchains),text.col=tmp$repSpecs,box.col='white',cex=.7)
  }
}
dev.copy2pdf(file=plotfile )



#includeSections <- c('transectEast','transectWest')
#includeSections <- c('transectNorth','transectSouth')

resp <- c('betaMort','betaMort','betaMort')
vars <- c('deficit','xeric','mesic')
nv   <- length(vars)
ht   <- c(1,.0001,.001,1)
xlo  <- c(-.0001,-.2,-.2)
xhi  <- c(.0001,.6,.2)

#includeSections <- c('NE','SE')

rr <- c(unique(resp),unique(vars))
plotfile <- outFile( savefolder,paste('postSens_NES', paste0(rr,collapse='-') ,'.pdf',sep='') )

par(mfrow=c(1,nv),bty='n')

ord <- NULL

kchains <- getRegionChains(includeSections,chainSens)

for(k in 1:nv){
  
  xtic <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=length(kchains),ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtic,textSize=.8,COLCODE=colorcode,
                           reverseSign=NULL,omitSpecs=omitSpecs,FILL=T   )
  
  legend('topright',legend=tmp$percPos,text.col=colorcode,cex=.8,
         border='white',box.col='white',xjust=1)
  legend('topleft',legend=tmp$percNeg,text.col=colorcode,cex=.8,
         border='white',box.col='white')
  
  if(k == 1){
    ord <- tmp$specOrd
    legend('bottomleft',names(kchains),text.col=tmp$repSpecs,box.col='white',cex=.7)
  }
}
dev.copy2pdf(file=plotfile )






######################### end NE-SE--all regions


kchains <- getRegionChains(regSections,chainSens)


# xeric vs mesic mortality

resp <- c('betaMort','betaMort','betaMort')
vars <- c('deficit','xeric','mesic')
nv   <- length(vars)
ht   <- c(.9,.01,.8)
xlo  <- c(-.0001,-.5 ,-.5)
xhi  <- c(.0001,1,.5)

rr <- c(unique(resp),unique(vars))
plotfile <- outFile( savefolder,paste('postSens_', paste0(rr,collapse='-') ,'.pdf',sep=''))

par(mfrow=c(1,nv))

ord <- NULL

for(k in 1:nv){
  
  xtic <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=length(kchains),ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtic,textSize=textSize,
                           reverseSign=NULL,omitSpecs=omitSpecs,FILL=T   )
  ord <- tmp$specOrd
}
legend('topright',regSections,text.col=tmp$repSpecs,box.col='white')
dev.copy2pdf(file=plotfile )




resp <- c('betaGam','beta','betaMort','betaFec')
vars <- c('temp-DecJanFebMar','temp-DecJanFebMar','temp-DecJanFebMar','temp-DecJanFebMar')
nv   <- length(vars)
ht   <- c(.07,.001,.9,1)
xlo  <- c(0, -.001,-.020,0)
xhi  <- c(2,.01,  .01,.4)

rr <- c(unique(resp),unique(vars))
plotfile <- outFile( savefolder,paste('postSens_', paste0(rr,collapse='-') ,'.pdf',sep=''))

par(mfrow=c(1,nv))

ord <- NULL
for(k in 1:nv){
  
  xtic <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=length(kchains),ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtic,textSize=textSize,
                           reverseSign=NULL,omitSpecs=omitSpecs,FILL=T)
  ord <- tmp$specOrd
}
legend('topright',regSections,text.col=tmp$repSpecs,box.col='white')
dev.copy2pdf(file=plotfile )


resp <- c('beta','betaMort','betaFec')
vars <- c('plotBA','plotBA','plotBA')
nv   <- length(vars)
ht   <- c(.9,1,1.8)
xlo  <- c(-.002,-.004,-.1)
xhi  <- c(  0,.004,  0)
plotfile <- outFile( savefolder,paste('postSens_',paste0(resp,vars,collapse='-'),'.pdf',sep=''))

par(mfrow=c(1,nv))

ord <- NULL
for(k in 1:nv){
  
  
  xtic <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=length(kchains),ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtic,textSize=textSize,
                           reverseSign=NULL,omitSpecs=omitSpecs,FILL=T   )
  ord <- tmp$specOrd
  if(k == nv)legend('topleft',regSections,text.col=tmp$repSpecs,box.col='white')
}
dev.copy2pdf(file=plotfile )


resp <- c('beta','beta')
vars <- c('temp','plotBA')
nv   <- length(vars)
ht   <- c(.9,1)
xlo  <- c(-.002,-.002)
xhi  <- c(  .01,.0005)
plotfile <- outFile( savefolder,paste('postSens_',paste0(resp,vars,collapse='-'),'.pdf',sep=''))

par(mfrow=c(1,nv))

ord <- NULL
for(k in 1:nv){
  
  
  xtic <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=length(kchains),ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtic,textSize=textSize,
                           reverseSign=NULL,omitSpecs=omitSpecs,FILL=T   )
  ord <- tmp$specOrd
  if(k == nv)legend('topleft',regSections,text.col=tmp$repSpecs,box.col='white')
}
dev.copy2pdf(file=plotfile )




cnames <- c('therm','deficit')
cnames <- c('deficit','temp-DecJanFebMar')
cnames <- c('therm','temp-DecJanFebMar')

cnames <- c('deficit','xeric')


cnames <- c('xeric','mesic')

cnames <- c('plotBA','deficit')

cnames <- c('plotBA','therm')
bnames <- 'betaChain'


includeSections <- c('MW','NE','SW','SE')
nsecs <- length(includeSections)

kchains <- getRegionChains(includeSections,chainSens)

par(mfrow=c(2,2),mar=c(5,4,2,2))

colF   <- colorRampPalette(c('orange','brown','black','blue','darkgreen'))
cseq   <- colF(length(specs))


xlim <- ylim <- numeric(0)

for(k in 1:nsecs){
  
  ccc <- chainSens[names(chainSens) == includeSections[k]]
  
  ck <- ccc[[1]][[bnames]]
  
  ch1 <- ck[,grep(cnames[1],colnames(ck))]
  ch2 <- ck[,grep(cnames[2],colnames(ck))]
  
  if(is.na(max(ch1)) | is.na(max(ch2)))next
  
  xlim <- range(c(xlim,ch1))
  ylim <- range(c(ylim,ch2))
}

for(k in 1:nsecs){
  
  
  ccc <- chainSens[names(chainSens) == includeSections[k]]
  
  ck <- ccc[[1]][[bnames]]
  
  
#  ck <- chainSens[[k]][[bnames]]
  
  ch1 <- ck[,grep(cnames[1],colnames(ck))]
  ch2 <- ck[,grep(cnames[2],colnames(ck))]
  
  if(is.na(max(ch1)) | is.na(max(ch2)))next
  
  tmp <- biVarMoments(ch1[,1],ch2[,1],wt=1,pr = .95)
  
  sknames <- matrix( unlist( strsplit(colnames(ch1),'_') ),ncol=2,byrow=T)[,2]
  
  plot(tmp$ellipse[,1],tmp$ellipse[,2],xlim=xlim,ylim=ylim,type='l',
                  xlab=cnames[1],ylab=cnames[2])
  abline(v=0,lty=2,lwd=2,col='grey')
  abline(h=0,lty=2,lwd=2,col='grey')
  for(s in 1:length(sknames)){
    if(sknames[s] %in% omitSpecs)next
    wcol <- match(sknames[s],specs)
 
    tmp <- biVarMoments(ch1[,s],ch2[,s],wt=1,pr = .95)
    lines(tmp$ellipse[,1],tmp$ellipse[,2],col=cseq[wcol])
#    text(max(tmp$ellipse[,1]),mean(tmp$ellipse[,2]),
#         sknames[s],col=cseq[wcol],pos=4,cex=.8)
  }
  title(includeSections[k])
}

plotfile <- outFile( savefolder,
                     paste('ellipse_',paste0(c(bnames,cnames),collapse='_') ,'.pdf',sep=''))
dev.copy2pdf(file=plotfile )



reverseSign <- c('deficit','xeric','plotBA')


################## GROWTH RESPONSES

# xeric vs mesic growth

dev.off()

resp <- c('beta','beta','beta')
vars <- c('deficit','xeric','mesic')
nv   <- length(vars)
ht   <- c(.8,.0005,.0025)
xlo  <- c(-.0001,-.15,0)
xhi  <- c(0,.05,.1)

rr <- c(unique(resp),unique(vars))
plotfile <- outFile( savefolder,paste('post_', paste0(rr,collapse='-') ,'.pdf',sep=''))

par(mfrow=c(1,nv),bty='n')

ord <- NULL


for(k in 1:nv){
  
  xtic <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=nsecs,ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtic,textSize=textSize,
                           reverseSign=NULL,omitSpecs=omitSpecs,FILL=T   )
  ord <- tmp$specOrd
}
legend('topright',regSections,text.col=tmp$repSpecs,box.col='white')
dev.copy2pdf(file=plotfile )


################SENSITIVITY IS FULL EFFECT, SO NO INTERACTIONS

resp <- c('beta','beta','beta')
vars <- c('deficit','xericXdeficit','mesicXdeficit')
nv   <- length(vars)
ht   <- c(.4,.01,.01)
xlo  <- c(-.0001,-.1,-.1)
xhi  <- c( 0, .1,.1)

rr <- c(unique(resp),unique(vars))
plotfile <- outFile( savefolder,paste('post_', paste0(rr,collapse='-') ,'.pdf',sep=''))

par(mfrow=c(1,nv))

ord <- NULL
for(k in 1:nv){
  
  xtic <- c(xlo[k],xhi[k] )
            
            tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                                     allChains=kchains,chainLength=nsecs,ORD=ord,
                                     vnames= vnames,MAINEFFECT=F, htFraction=ht[k],
                                     xtic=xtic,textSize=textSize,
                                     reverseSign=NULL,omitSpecs=omitSpecs,FILL=T   )
            ord <- rev(tmp$specOrd)
}
legend('topright',regSections,text.col=tmp$repSpecs,box.col='white')
dev.copy2pdf(file=plotfile )



################## MORT RESPONSES

# xeric vs mesic growth

resp <- c('betaMort','betaMort','betaMort')
vars <- c('deficit','xeric','mesic')
nv   <- length(vars)
ht   <- c(.8,.00001,.005)
xlo  <- c(-.0002,-.5 ,-.5)
xhi  <- c(.0001,1,.5)

rr <- c(unique(resp),unique(vars))
plotfile <- outFile( savefolder,paste('post_', paste0(rr,collapse='-') ,'.pdf',sep=''))

par(mfrow=c(1,nv))

ord <- NULL

for(k in 1:nv){
  
  xtk <- c( xlo[k],xhi[k] )
            
            tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                                     allChains=kchains,chainLength=nsecs,ORD=ord,
                                     vnames= vnames, htFraction=ht[k],
                                     xtic=xtk,textSize=textSize,
                                     reverseSign=reverseSign,omitSpecs=omitSpecs,FILL=T   )
            ord <- rev(tmp$specOrd)
}
legend('topright',regSections,text.col=tmp$repSpecs,box.col='white')
dev.copy2pdf(file=plotfile )

#######################winterTEMP

resp <- c('betaGam','beta','betaMort','betaFec')
vars <- c('temp-DecJanFebMar','temp-DecJanFebMar','temp-DecJanFebMar','temp-DecJanFebMar')
nv   <- length(vars)
ht   <- c(.01,.0000001,.000001,.0000000001)
xlo  <- c(-.5, 0,-.01,0)
xhi  <- c(2,.01,  .01,.2)

rr <- c(unique(resp),unique(vars))
plotfile <- outFile( savefolder,paste('post_', paste0(rr,collapse='-') ,'.pdf',sep=''))

par(mfrow=c(1,nv))

ord <- NULL
for(k in 1:nv){
  

  xtk <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=nsecs,ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtk,textSize=textSize,
                           reverseSign=reverseSign,omitSpecs=omitSpecs,FILL=T   )
  ord <- tmp$specOrd
}
legend('topright',regSections,text.col=tmp$repSpecs,box.col='white')
dev.copy2pdf(file=plotfile )



#######################plotBA

resp <- c('beta','betaMort','betaFec')
vars <- c('plotBA','plotBA','plotBA')
nv   <- length(vars)
ht   <- c(.8,.6,.1)
xlo  <- c(-.002,-.006,-.20)
xhi  <- c(  0,  .006,  0)
plotfile <- outFile( savefolder,paste('post_',paste0(resp,vars,collapse='-'),'.pdf',sep=''))

par(mfrow=c(1,nv))

ord <- NULL
for(k in 1:nv){
  
  xtk <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=nsecs,ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtk,textSize=textSize,
                           reverseSign=NULL,omitSpecs=omitSpecs,FILL=T   )
  ord <- tmp$specOrd
  if(k == 1)legend('topleft',regSections,text.col=tmp$repSpecs,box.col='white')
}
dev.copy2pdf(file=plotfile )

####################### deficit

resp <- c('beta','betaMort','betaFec')
vars <- c('deficit','deficit','deficit')
nv   <- length(vars)
ht   <- c(.9,.9,.01)
xlo  <- c(-.0001,-.0001,-.0040)
xhi  <- c(  0,  .0001,  0)
plotfile <- outFile( savefolder,paste('post_',paste0(resp,vars,collapse='-'),'.pdf',sep=''))

par(mfrow=c(1,nv))

ord <- NULL
for(k in 1:nv){
  
  if(k == 1)legend('topleft',regSections,text.col=tmp$repSpecs,bty='n',box.col='white')
  xtk <- c( xlo[k],xhi[k] )
  
  tmp <-   plotDensColumns(var2plot=vars[k],chain2plot=resp[k],
                           allChains=kchains,chainLength=nsecs,ORD=ord,
                           vnames= vnames, htFraction=ht[k],
                           xtic=xtk,textSize=textSize,
                           reverseSign=NULL,omitSpecs=omitSpecs,FILL=T   )
  ord <- tmp$specOrd
}

dev.copy2pdf(file=plotfile )

############################## maps

f0 <- reproDatAll*0
f0[reproDatAll == 0] <- 1
f0[reproDatAll > 0] <- 0

temp <- climAll[,'temp']
prec <- climAll[,'prec']

outfolder <- savefolder

 lims <- rbind(maplon,maplat)

for(s in 1:nspec){
  
  if( specs[s] %in% omitSpecs )next
  
  graphics.off()
  
  par(mfrow=c(2,2),bty='n')
  
  plotRate(s,demName='grow',demog=growMeanDatAll,
           temp = temp, prec = prec,
           maplon=maplon,maplat=maplat,LEGEND=F)
  
  plotRate(s,demName='surv',demog=survDatAll,
           temp = temp, prec = prec,
           maplon=maplon,maplat=maplat,LEGEND=F)
  
  dev.copy2pdf(file=outFile(savefolder,paste('rates_',specs[s],'.pdf',sep='')))
  
}

for(s in 1:nspec){

  graphics.off()
  
  data(countyMapEnv)
  data(stateMapEnv)
  
  varList <- list(i50no_competition = presTreeNCAll, 
                  i50competition = presTreeAll,
                  p50 = presPopAll)
  varCols <- c('orange','brown','blue')
  varFill <- c(T,T,T)
  varLineLwd  <- c(3,1,1)
  varPoints   <- c(T,T,T)
  zlevs       <- c(.5,.5,.5)
  inLegend    <- c(T,T,T)
  
  data(countyMapEnv)
  data(stateMapEnv)
  
  tmp <- presenceMap(spec=specs[s],varList=varList,
                     varCols=varCols,varFill=varFill,
                     TOPOMAP=F,varLineLwd=varLineLwd,varPoints=varPoints,
                     xlim=lims[1,], ylim=lims[2,],scaleSym=.1,
                     noY=y0[,s],topoList=topo,
                     lon=lonLatAll[,1],lat=lonLatAll[,2],
                     inLegend=inLegend,legendCorner='bottomright',
                     zlevs=zlevs,little=F,REMOTE=REMOTE,PRESENCE=T,mapscale=3.5)
  
 # sTree <- tmp$i50no_competition
 # sPop  <- tmp$p50
  
  graphics.off()
  data(countyMapEnv)
  data(stateMapEnv)
  
  varList <- list(p5_fecundity = presFecAll,p5_recruitment = presReprAll)
  varCols <- c('tan','brown')
  varFill <- c(T,T)
  varLineLwd  <- c(3,1)
  varPoints   <- c(T,T)
  zlevs       <- c(.4,.4)
  inLegend    <- c(T,T)
  
  data(countyMapEnv)
  data(stateMapEnv)
  
  tmp <- presenceMap(spec=specs[s],varList=varList,varCols=varCols,varFill=varFill,
                     TOPOMAP=F,varLineLwd=varLineLwd,varPoints=varPoints,
                     xlim=lims[1,], ylim=lims[2,],scaleSym=.1,
                     lon=lonLatAll[,1],lat=lonLatAll[,2],
                     noY=f0[,s],inLegend=inLegend,zlevs=zlevs,
                     little=F,legendCorner='bottomright',
                     REMOTE=REMOTE,PRESENCE=T,mapscale=3.5)
  
  #  ww <- which(sTree$z > .5)
 # sFec  <- tmp$p5_fecundity
  
 # matTree[,s] <- as.vector(sTree$z)
 # matPop[,s] <- as.vector(sPop$z)
 # matRec[,s] <- as.vector(sFec$z)
  
}

#######################scores


plotfile <- outFile(savefolder,'scoreMap.pdf')

score <- scoreAll

ncol   <- 50
colF   <- colorRampPalette(c('wheat','brown','black'))
colq   <- mapColors(ncol)

colC <- colorRampPalette(c('blue','darkgreen','brown','tan','orange','red','magenta'))
colq <- colC(ncol)

graphics.off()
par(mfrow=c(2,2),mar=c(.5,.5,.5,.5))

scaleCoords <- rbind(maplon[1] + diff(maplon)*c(.9,.97),
                     maplat[1] + diff(maplat)*c(.05,.25))

pvals  <- c(.1,.9)
qscore <- quantile(score,pvals)

scaleCoords <- rbind(maplon[1] + diff(maplon)*c(.9,.97),
                     maplat[1] + diff(maplat)*c(.05,.25))

for(j in 1:ncol(score)){
  
  pj <- colnames(score)[j]
  
  ssj <- score[,j]
  
  qscore <- quantile(ssj,pvals)
  
  ssj[ssj < (qscore[1] + .01)] <- qscore[1] + .01
  ssj[ssj > qscore[2]] <- qscore[2] - 1e-5
  
  ss <- qscore
  
  sc <- seq(ss[1],ss[2],length=ncol)
  sc[sc < ss[1]] <- sc[1]
  sc[sc > ss[2]] <- sc[2]
  
  data(countyMapEnv)
  data(stateMapEnv)
  

  
  regMap(x=topo$x,y=topo$y,z=topo$z,IMAGE=F,xl=c(-96,-67),
         yl=maplat,mapscale=20)
  
 
  zlevs <- round(seq(ss[1],ss[2],length=ncol),5)

  
  fi <- findInterval(zlevs,sc)
  fi[fi == 1] <- 2
  
  colz  <- colq[fi]
  
  values2contour(xx=lonLatAll[,1],yy=lonLatAll[,2],
                 z=ssj,nx=50,ny=50,col=colz,
                 zlevs=zlevs,lwd=4,add=T,fill=T)
  
  wz <- which(topo$z > -1,arr.ind=T)
 # mapMask(xx=topo$x[wz[,1]],yy=topo$y[wz[,2]],dx=.4,dy=.4,whiteOutDist=.2,col='white')
  
  mapMask(lonLatAll[,1],lonLatAll[,2],dx=.3,dy=.3,whiteOutDist=.4,col='white')
  
  map(add=T,col='white',lwd=8)
  map(add=T,col='black',lwd=2)
  
  text(scaleCoords[1,2],scaleCoords[2,2] + 2,pj,pos=2)
  
  colorLegend(scaleCoords[1,],scaleCoords[2,],
                        scale=c(1:ncol),cols=colq,
                        labside='left',endLabels=round(ss,1))
}


dev.copy2pdf(file=plotfile)

################### SPECIES GROUPS


mvars  <- c('gsigmaChain','sigmaChain','fsigmaChain','msigmaChain')
mnames <- c('SDM','growth','fecundity','mortality')

omitFromClus <- c('other')

covAll <- cov(treeAll)
nall   <- nrow(covAll)

tmp <- treeAll
tmp[treeAll == 0] <- NA

corAll <- cor(tmp,use='pairwise.complete.obs')
corAll[is.na(corAll)] <- 0                      #only where both spp occur

corAll <- cor(treeAll)
           
clusStats <- numeric(0)
cname <- character(0)

for(j in 1:nsecs){         #regions
  
  par(mfrow=c(1,5), mar = c(5,2,1,8))
  
  jm <- regSections[j]
  plotfile <- outFile(savefolder,paste('cluster_',jm,'.pdf',sep='') )
  
  varVec <- sampSize[match(jm,rownames(sampSize)),]
  
  specm    <- specsAll[[j]]
  specm    <- specm[!specm %in% omitFromClus]
  mspec    <- length(specm)
  
  pm <- length(xnamesAll[[m]])
  
  mm <- match(specm,rownames(covAll))
  acov <- conditionalMVN(rep(0,nall),rep(0,nall),covAll,mm)$vr
  
  colCode <- NULL
  
 # acor <- cov2cor(acov)
  acor <- corAll[mm,mm]
  dist <- cov2Dist(covAll[mm,mm])
  
#  dist <- adis
#  dist <- 1 - acor
  
  tmp <- clusterPlot( dcor = acor ,method='complete',main='Data',
                      cex=.4,ncluster=6, colCode=NULL)
  colCode <- tmp$colCode
  ciMain <- tmp$clusterIndex
  cim    <- numeric(0)
  
  for(m in 1:length(mnames)){       #correlations
    
    if(mnames[m] == 'SDM'){
      vname <- 'gsigmaChain'
      bname <- 'betaGamChain'
      bcov  <- xgamCov[[j]]
    }
    if(mnames[m] == 'growth'){
      vname <- 'sigmaChain'
      bname <- 'betaChain'
      bcov  <- xgroCov[[j]]
    }
    if(mnames[m] == 'fecundity'){
      vname <- 'fsigmaChain'
      bname <- 'betaFecChain'
      bcov  <- xfecCov[[j]]
    }
    if(mnames[m] == 'mortality'){
      vname <- 'msigmaChain'
      bname <- 'betaMortChain'
      bcov  <- xmorCov[[j]]
    }
    
    nvar <- varVec[mnames[m]]
    
    wm <- match(paste(regSections[j],bname,sep='_'),names(betaMats))
    
    bm <- betaMats[[wm]][,specm]
    
    bvar <- t(bm)%*%bcov%*%bm
    
    dist <- cov2Dist(bvar)
    rownames(dist) <- colnames(dist) <- colnames(bm)
    
    acor <- cov2cor(bvar)
    specCols <- colCode[match(rownames(bvar),names(colCode))]
  #  dist <- cov2Dist(acov)
    
    #  dist <- adis
  #  dist <- 1 - acor
    
    
    tmp <- clusterPlot( dcor = acor ,method='complete',
                        main='',cex=.3,ncluster=6, 
                        colCode=specCols)
    cim <- append(cim,list(tmp$clusterIndex))
    
    tmp1 <- cluster.stats(clustering=ciMain,alt.clustering=cim[[m]],compareonly=T)
    
    
    cname <- c(cname,paste(regSections[j],vname,sep='-'))
    
    tmp2 <- compareCorMatrix(corAll[mm,mm],acor,
                             varVec['nplot'],nvar,alpha=.05)  #NEED REAL N1, N2
    percDiff <- tmp2$fraction*100
    
    clusStats <- rbind(clusStats,c(tmp1$corrected.rand,tmp1$vi,
                                   percDiff,tmp2$meanDis,tmp2$meanDisST))
    
    
    mp <- signif(tmp1$corrected.rand*100,2)
    md <- signif(tmp2$meanDis,2)
    mmain <- paste(mnames[m],' ',mp,'%',sep='')
    title(main=mmain)
    
    
    
    #  colCode <- tmp$colCode
  
    
  }
  dev.copy2pdf(file=plotfile )
}

rownames(clusStats) <- cname
colnames(clusStats) <- c('corrected Rand','VI','% different','mean Dist','standard Dist')

signif(clusStats,3)

# high means aggrement
    


