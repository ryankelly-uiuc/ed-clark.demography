library(maps)


getTraits <- function(STANDARDIZE=F){
  
  tmpN <- read.table('../allocationmodel/datafiles/traitN.txt',header=T)
  tmpP <- read.table('../allocationmodel/datafiles/traitP.txt',header=T)
  tSLA <- read.table('../allocationmodel/datafiles/traitSLA.txt',header=T)
  
  tmp <- read.table('../allocationmodel/datafiles/treeCodesDuke.txt',header=T)
  
  nindex <- match(tmpN[,'species'],tmp[,'code'])
  pindex <- match(tmpP[,'species'],tmp[,'code'])
  sindex <- match(tSLA[,'species'],tmp[,'code'])
  
  tmpN <- tmpN[is.finite(nindex),]
  tmpP <- tmpP[is.finite(pindex),]
  tSLA <- tSLA[is.finite(sindex),]
  
  nindex <- match(tmpN[,'species'],tmp[,'code'])
  pindex <- match(tmpP[,'species'],tmp[,'code'])
  sindex <- match(tSLA[,'species'],tmp[,'code'])
  
  moreCols <- matrix(NA,nrow(tmp),3)
  colnames(moreCols) <- c('leafN','leafP','SLA')
  moreCols[nindex,1] <- tmpN[,'mu']
  moreCols[pindex,2] <- tmpP[,'mu']
  moreCols[sindex,3] <- tSLA[,'mu']
  
  tmp <- cbind(tmp,moreCols)
  
  
  quanT <- c('kgM3Wood','gm1000Seed','maxHt','leafN','leafP','SLA')
  qualT <- c('monoecious','leaves')
  
  nq <- length(quanT)
  
  genera <- sort(unique(tmp[,'genus']))
  gindex <- match(tmp[,'genus'],genera)
  genTable <- numeric(0)
  
  fam      <- sort(unique(tmp[,'family']))
  findex   <- match(tmp[,'family'],fam)
  famTable <- numeric(0)
  
  for(j in 1:nq){
    
    tj <- tmp[,quanT[j]]
    gj <- byIndex(tj,gindex,mean,na.rm=T,coerce=T)
    gj[gj == 0] <- NA
    
    tg <- gj[gindex]                 #value for genus
    tj[is.na(tj)] <- tg[is.na(tj)]
    tmp[,quanT[j]] <- signif(tj,3)
    genTable <- cbind(genTable,signif(gj,3))
    
    tj <- tmp[,quanT[j]]
    gj <- byIndex(tj,findex,mean,na.rm=T,coerce=T)
    gj[gj == 0] <- NA
    
    tg <- gj[findex]
    tj[is.na(tj)] <- tg[is.na(tj)]
    tmp[,quanT[j]] <- signif(tj,3)
    famTable <- cbind(famTable,signif(gj,3))
    
  }
  
  rownames(genTable) <- genera
  colnames(genTable) <- quanT
  
  rownames(famTable) <- fam
  colnames(famTable) <- quanT
  
  cmu <- colMeans(tmp[,quanT],na.rm=T)
  ww  <- which(is.na(tmp[,quanT]),arr.ind=T)
  tmp[,quanT][ww] <- cmu[ww[,2]]
  
  tmp[,'gm1000Seed'] <- log10(tmp[,'gm1000Seed']/1000)
  tmp[,'kgM3Wood']   <- tmp[,'kgM3Wood']/1000
  
  colnames(tmp)[colnames(tmp) == 'kgM3Wood'] <- 'gmPerCm'
  colnames(tmp)[colnames(tmp) == 'gm1000Seed'] <- 'gmPerSeed'
  
  genTable[,'gm1000Seed'] <- log10(genTable[,'gm1000Seed']/1000)
  genTable[,'kgM3Wood']   <- genTable[,'kgM3Wood']/1000
  
  colnames(genTable)[colnames(genTable) == 'kgM3Wood'] <- 'gmPerCm'
  colnames(genTable)[colnames(genTable) == 'gm1000Seed'] <- 'gmPerSeed'
  
  famTable[,'gm1000Seed'] <- log10(famTable[,'gm1000Seed']/1000)
  famTable[,'kgM3Wood']   <- famTable[,'kgM3Wood']/1000
  
  colnames(famTable)[colnames(famTable) == 'kgM3Wood'] <- 'gmPerCm'
  colnames(famTable)[colnames(famTable) == 'gm1000Seed'] <- 'gmPerSeed'
  
  #########################
  
  traitNames <- c('gmPerSeed','gmPerCm','maxHt','leafN','leafP','SLA')
  lt     <- tmp[,'leaves']
  lTot   <- sort(unique(lt))
  leaf   <- match(lt,lTot)
  
  lmat <- matrix(0,length(leaf),max(leaf))
  lmat[cbind( c(1:length(leaf)),leaf) ] <- 1
  colnames(lmat) <- lTot
  other <- rep(0,length(lTot))
  other[1] <- 1
  lmat <- rbind(lmat,other)
  
  ln <- matrix( unlist(strsplit(colnames(lmat),'_')), 2,4)
  for(j in 1:4)colnames(lmat)[j] <- paste(ln[1,j],ln[2,j],sep='')
  
  lmat <- lmat[,colnames(lmat) != 'needledeciduous']    # only a few species
  
  tnames <- c(traitNames,colnames(lmat))
 # M      <- length(tnames)
  
  traits <- tmp[,traitNames]
  other <- colMeans(traits)
  traits <- rbind(traits,other)
  traits <- cbind(traits,lmat)
  rownames(traits) <- c(as.character(tmp[,'code']),'other')
  
  tmu <- colMeans(traits[,traitNames],na.rm=T)
  tsd <- apply(traits[,traitNames],2,sd,na.rm=T)
  
  if(STANDARDIZE){
    
    tm1 <- matrix(tmu,nrow(tmp),length(traitNames),byrow=T)
    ts1 <- matrix(tsd,nrow(tmp),length(traitNames),byrow=T)
    tmp[,traitNames] <- (tmp[,traitNames] - tm1)/ts1
    
    tm2 <- matrix(tmu,nrow(traits),length(traitNames),byrow=T)
    ts2 <- matrix(tsd,nrow(traits),length(traitNames),byrow=T)
    traits[,traitNames] <- (traits[,traitNames] - tm2)/ts2
  }
  
  list( traitTable = tmp, traitsByGenus = genTable, traitsByFam = famTable,
        traitMatrix = traits,
        traitMeans = tmu, traitSds = tsd)
}

traitLabel <- function(tname){
  
  if(tname == 'gmPerSeed')tname <- 'seed mass'
  if(tname == 'gmPerCm')  tname <- 'wood density'
  if(tname == 'maxHt')    tname <- 'maximum height'
  if(tname == 'leafN')    tname <- 'leaf N'
  if(tname == 'leafP')    tname <- 'leaf P'
  tname
}

getTraitValues <- function(ww,betaMu=NULL,sigMu=NULL,covX){
  
  tmp <- getTraits(STANDARDIZE=T)$traitMatrix
  traitMat <- t( tmp[colnames(ww),] )
  
  nn   <- nrow(ww)
  traitMu  <- rowMeans(traitMat)
  traitDev <- traitMat - matrix(traitMu,M,S)
  TD   <- traitDev%*%t(traitDev)/S
  tbar <- 1/nn*traitMat%*%t(ww)%*%matrix(1,nn)   #weighted trait norm
  
  tbarMat <- matrix(tbar,M,S)
  tdev <- traitMat - tbarMat
  WTD  <- matrix(1,M)%*%matrix(1,ncol=nn)%*%ww
  WTD  <- (tdev*WTD)%*%t(tdev)/nn
  
  sig <- cov(ww,use='pairwise.complete.obs')
  sig[is.na(sig)] <- 0
  TV  <- traitMat%*%sig%*%t(traitMat)
  
  sumG <- rbind(t(tbar),diag(TD),diag(WTD),diag(TV))
  rownames(sumG) <- c('WTN','TD','WTD','TV')
  
  XTV <- RTV <- pXV <- numeric(0)
  
  if(!is.null(betaMu)){
    ETV    <- traitMat%*%t(betaMu)%*%covX%*%betaMu%*%t(traitMat)
    tsigma <- traitMat%*%sigMu%*%t(traitMat)
    XTV    <- diag(ETV)
    RTV    <- diag(tsigma)
   # sumG   <- rbind(sumG,XTV,RTV)
    RTV <- tsigma
    PTV  <- XTV/(XTV + diag(RTV))
    sumG <- rbind(sumG,PTV)
  }
  
  list(table = sumG, TD = TD, WTD = WTD, TV = TV,
       ETV = ETV)
}

setupPrior <- function(){
  
  if('hydro' %in% posPrior & hydroVariable == 'factor'){
    posPrior <- posPrior[!posPrior == 'hydro']
    posPrior <- c(posPrior,'mesic')
    negPrior <- c(negPrior,'xeric')
  }
  if('hydro' %in% negPrior & hydroVariable == 'factor'){
    negPrior <- negPrior[!negPrior == 'hydro']
    posPrior <- c(posPrior,'xeric')
  }
  list(posPrior = posPrior, negPrior = negPrior)
}


setupHydroFactor <- function(){
  
  if('hydro' %in% xnames){
    hvars <- c('xeric','mesic')
    xtmp   <- xnames[xnames != 'hydro']
    xnames <- c(xtmp,hvars)
    
    xtmp   <- fnames[fnames != 'hydro']
    fnames <- c(xtmp,hvars)
    
    xtmp   <- mnames[mnames != 'hydro']
    mnames <- c(xtmp,hvars)
    
    xtmp   <- gnames[gnames != 'hydro']
    gnames <- c(xtmp,hvars)
  }
  
  list(xnames = xnames, fnames = fnames, mnames = mnames, gnames = gnames)
}


regMap <- function(x,y,z=NULL,fg=rep(1,length(x)),bg=rep(1,length(x)),
                   xl=range(x),yl=range(y),lineCol='black',axes=T,lwd=1,
                   zminmax=NULL,scaleSym=1,IMAGE=F,ADD=F,mapscale=1,county=F,
                   stateVec=c('north carolina','south carolina','tennessee',
                              'georgia','virginia')){
  
  # if IMAGE z is a matrix with length(x) rows, length(y) columns)
  # if !IMAGE, xx, yy, and z are all same length
  require(maps)
  require(maptools)
  
  if(!ADD){
    mapSetup(xl,yl,scale=mapscale)
  }
  
  colSea   <- colorRampPalette(c('black','darkblue','blue','lightblue','white'))
  
  cols <- c(colSea(5),terrain.colors(100))
  
  zr <- max(z)
  bb <- c( -10000,seq(-300,-5,length=5), seq(-2,160,length=50), seq(170,zr+1,length=50))
  
  if(!is.null(zminmax)){
    bb <- seq(zminmax[1],zminmax[2],length=106)
  }
  
  if(IMAGE){
    image(x,y,z,col=cols,breaks=bb,add=ADD,xlab=' ',ylab=' ',xlim=xl,ylim=yl) 
    ADD <- T
  }
  
#  map('county','north carolina',boundary=F,col='grey',add=T)

  if(!is.null(stateVec) & county)map('county',stateVec,boundary=T,col='grey',lwd=lwd,
                                     add=ADD,xlim=range(x),ylim=range(y))
  map('state',add=ADD,col=lineCol,lwd=lwd,xlim=xl,ylim=yl)
#  map('county',stateVec,interior=F,col='black',add=T)
  if(axes)map.axes()
  
  if(IMAGE)return()
  
  if(is.null(z)) points(x,y,col=bg)
  if(!IMAGE & !is.null(z) & length(x) == length(y))
    symbols(x,y,circles=z*scaleSym,fg=fg,bg=bg,inches=F,add=ADD)
}



specEllipse <- function(chainName,v1,v2,xlim=NULL,ylim=NULL,PLOT=F){
  
  require(cluster)
  
  ww <- which(chains == chainName)
  xx <- chainList[[ww]]
  
  cc <- colnames(xx)
  cc[grep('X',cc)] <- 'BOGUS'   #do not include terms with interactions
  
  g1 <- grep(v1,cc)
  g2 <- grep(v2,cc)
  
  if(is.null(xlim))xlim <- range(xx[,g1])
  if(is.null(ylim))ylim <- range(xx[,g2])
  
  if(PLOT){
    plot(0,0,xlim=xlim,ylim=ylim,cex=.01,xlab=v1,ylab=v2)
    title(chainName)
  
    for(j in 1:nspec){
      
      x1 <- intersect(g1,grep(specs[j],colnames(xx)))
      x2 <- intersect(g2,grep(specs[j],colnames(xx)))
      
      x1 <- xx[,x1]
      x2 <- xx[,x2]
      
      mu <- c(mean(x1),mean(x2))
      vr <- cov(cbind(x1,x2))
      
      tmp  <- list(loc = mu, cov = vr, d2 = qchisq(.9,1) )
      tmp  <- predict.ellipsoid(tmp)
      lines(tmp[,1],tmp[,2],col=j,lwd=2)
      text(mu[1],mu[2],specs[j],col=j)
    } 
  }
}

abundanceMap <- function(z,spec,nx=100,ny=100,zlevs=c(0,.1),title=' ',topoList=topo,
                         lon=plotLon,lat=plotLat,little=F,REMOTE=F,nHoldOut=0,mapscale=3){
  
  cols <- colorRampPalette(c('wheat','orange','brown','maroon'))
  
  colors <- cols(length(zlevs))
  
  ss <- match(spec,specs)
  
  zr <- quantile(z[,ss],.95)
  
  plotfile <- outFile(savefolder,paste(spec,'_',title,'.pdf',sep=''))
#  plotstart(plotfile,REMOTE)
  
  regMap(topoList$x,topoList$y,topoList$z,mapscale=mapscale,IMAGE=F)
  
  
  keep <- c(1:nplot)
  if(nHoldOut > 0)keep <- keep[-holdOutPlots]
  
  llon <- lon[keep]
  llat <- lat[keep]
  z    <- z[keep,ss]
  
  ww <- which(topoList$z < 0,arr.ind=T)  #ocean
  
  if(length(ww) > 0){
    
    llon <- c(llon,topoList$x[ww[,1]])
    llat <- c(llat,topoList$y[ww[,2]])
    z   <- c(z,rep(-.1,nrow(ww)))
  }
  
  xr    <- range(topoList$x,na.rm=T)
  yr    <- range(topoList$y,na.rm=T)
  d     <- max(c(diff(xr),diff(yr)),na.rm=T)
  
  for(k in 1:length(zlevs)){
    if(zlevs[k] > zr)break
#    z1 <- zlevs[k-1]
    values2contour(llon,llat,z,nx=nx,ny=ny,col=colors[k],
                   lwd=2,zlevs=c(zlevs[k],max(zlevs)),add=T,fill=T)
  }
  
  if(little){
    tmp <- distributionMapLittle(spec)
    plot(tmp,add=T,lwd=3,lty=2)
  }
  
  mapMask(lon,lat,dx=.2,dy=.2,whiteOutDist=.3,col='white')
  
  
  map('state',add=T,col='white',lwd=6,interior=F)
  map('state',add=T,col='grey',lwd=3,interior=T)
  
  add.scatter(hist2Add(z),posi='bottomright',ratio=.15,inset=c(.02,.07),bg.col='grey')
  
  
  title(paste(title,spec,sep='-'))
  
  dev.copy2pdf(file=plotfile)
  
#  plotend(plotfile,REMOTE)
}

diffDem <- function(dem1=1-prAbsReg,dem2=1-prAbsent){         #difference in demographic Probs
  
  demDif <- dem1 - dem2
  wpos   <- which(demDif > 0)
}


plotDemArrows <- function(lonlat = demLonLat){
  
  if(is.list(lonlat))   lonlat <- as.matrix(lonlat)
  if(!is.matrix(lonlat))lonlat <- matrix(lonlat,1)
  
  points(lonlat[,1],lonlat[,2],lwd=10,col='white',pch=7)
  points(lonlat[,1],lonlat[,2],lwd=5,col='brown',pch=7)
  
#  arrows(lonlat[,1] - .3,lonlat[,2] - .3, lonlat[,1],lonlat[,2],
#         length=.2,angle=20,lwd=8,col='white')
#  arrows(lonlat[,1] - .3,lonlat[,2] - .3, lonlat[,1],lonlat[,2],
#         length=.2,angle=20,lwd=4)
  
}


presenceMap <- function(spec,varList,
                        varCols=rep(1,length(varList)),
                        varFill=rep(F,length(varList)),
                        varLineLwd=rep(1,length(varList)),
                        varPoints = rep(F,length(varList)),
                        scaleSym=1/15,
                        xlim=NULL, ylim=NULL,
                        inLegend = rep(F,length(varList)),legendCorner='bottomright',
                        noY=y0,TOPOMAP=F,topoList=topo,
                        lon=plotLon,lat=plotLat,
                        zlevs = rep(.5,length(varList)),little=F,REMOTE=F,
                        drawlabels=T, PRESENCE=F,mapscale=2,nHoldOut=0){   
  require(mapdata)
 # require(MBA)
  require(sp)

  nplot    <- length(lon)
  
  plotfile <- outFile(outfolder,paste(spec,'_',names(varList)[1],
                                      'presence.pdf',sep=''))
#  plotstart(plotfile,REMOTE)
  
  nsurf <- length(varList)
  
  surface <- vector(mode='list',length=nsurf) 
  names(surface) <- names(varList)
  
  if(nsurf > 1 & length(varCols) == 1)    varCols <- rep(varCols,nsurf)
  if(nsurf > 1 & length(varLineLwd) == 1) varLineLwd <- rep(varLineLwd,nsurf)
  
  if(is.null(xlim))xlim <- range(topoList$x)
  if(is.null(ylim))ylim <- range(topoList$y)
  
  print(xlim)
  print(ylim)
  
  wx <- which(topoList$x >= xlim[1] & topoList$x <= xlim[2])
  wy <- which(topoList$y >= ylim[1] & topoList$y <= ylim[2])
  
  xx <- topoList$x[wx]
  yy <- topoList$y[wy]
  zz <- topoList$z[wx,wy]
  
#  cols <- c('wheat','orange','brown')
  
  ss <- match(spec,specs)
  
  if(is.na(ss))stop('species not in data')
  
  graphics.off()
  
  IM <- F
  if(TOPOMAP)IM <- T
  
  regMap(xx,yy,zz,IMAGE=IM,mapscale=mapscale)

  keep <- c(1:nplot)
  if(nHoldOut > 0)keep <- keep[-holdOutPlots]
  
  wk <- which(lon >= xlim[1] & lon <= xlim[2] &
              lat >= ylim[1] & lat <= ylim[2])
  
  keep <- keep[keep %in% wk]
  
  vList <- varList
  
  for(m in 1:nsurf)vList[[m]] <- varList[[m]][keep,ss]
  
  z0 <- noY[keep]
  llon <- lon[keep]
  llat <- lat[keep]
  
  ww <- which(zz < 0,arr.ind=T)  #ocean
  
  if(length(ww) > 0){
    
    llon <- c(llon,xx[ww[,1]])
    llat <- c(llat,yy[ww[,2]])
    for(m in 1:nsurf)vList[[m]] <- c(vList[[m]],rep(-.1,nrow(ww)))
  }
  
  xr <- range(xx,na.rm=T)
  yr <- range(yy,na.rm=T)
 # d  <- max(c(diff(xr),diff(yr)),na.rm=T)
 # q  <- d/30
  
  
  ww <- which(noY == 1)
 # fg <- bg <- 'grey'
 # regMap(lon[ww],lat[ww],noY[ww],fg=fg, bg=bg, 
 #        scaleSym=scaleSym/10,ADD=T,stateVec=NULL)
  
  for(m in 1:nsurf){
    
    surface[[m]]  <- values2contour(xx=llon,yy=llat,z=vList[[m]],
                                    nx=130,ny=130,
                                    col=varCols[m],
                                    lwd=varLineLwd[m],
                                    zlevs=c(zlevs[m],1.1),add=T,fill=varFill[m],
                                    drawlabels=drawlabels)
  }
  

  
  ww <- which(  varList[[1]][,ss] > zlevs[1] 
              & varList[[2]][,ss] < zlevs[2])
  fg <- bg <- varCols[1]
  regMap(lon[ww],lat[ww],ww*0+1,fg=fg, bg=bg, 
         scaleSym=scaleSym/1.6,ADD=T,stateVec=NULL)

  
  if(nsurf > 2){
    for(kk in 2:(nsurf-1)){
      ww <- which(varList[[kk]][,ss] > zlevs[kk] & 
                  varList[[kk+1]][,ss] < zlevs[kk])
      regMap(lon[ww],lat[ww],ww*0+1,fg=varCols[kk], bg=varCols[kk], 
             scaleSym=scaleSym*kk*.1,ADD=T,stateVec=NULL)
    }
  }
  ww <- which(varList[[1]][,ss] < zlevs[1])
  fg <- bg <- 'white'   #completely absent
  regMap(lon[ww],lat[ww],ww*0+1,fg=fg, bg=bg, 
         scaleSym=scaleSym/1.7,ADD=T,stateVec=NULL)
  
  if(PRESENCE){
    ww <- which(noY == 0)
    if(length(ww) > 0){
      fg <- bg <- 'red'
      regMap(lon[ww],lat[ww],noY[ww] + 1,
             fg=varCols[m], bg=varCols[m], xl=xlim,yl=ylim,
             scaleSym=scaleSym*.9,ADD=T,stateVec=NULL)
      regMap(lon[ww],lat[ww],noY[ww] + 1,fg=fg, bg=bg, xl=xlim,yl=ylim,
             scaleSym=scaleSym*.7,ADD=T,stateVec=NULL)
    }
  }
    
  if(little){
    tmp <- distributionMapLittle(spec)
    plot(tmp,add=T,lwd=9,border='white')
    plot(tmp,add=T,lwd=7,border='grey')
    plot(tmp,add=T,lwd=4,lty=2,border='black')
  }
  
  title(spec)
  
  
  map('world','canada',fill=T,add=T,boundary=F,col='white')
  
  lnames <- replaceString(names(varList),now='_',new=' ')
    
  if(length(lnames[inLegend]) > 0)legend(legendCorner,lnames[inLegend],
                          col='white',box.col='white',text.col=varCols[inLegend],cex=.7)
  
  
 # if(region == 'SE')plotDemArrows(demLonLat)
 # if(region == 'NE')plotDemArrows(HFLonLat)
  
  dev.copy2pdf(file=plotfile)
  
#  plotend(plotfile,REMOTE)
  
  invisible( surface )
}

nearest <- function(x1,x2){   #closest point in vector x2 for each element of vector x1
  
  n1 <- length(x1)
  n2 <- length(x2)
  
  tmp <- abs( matrix(x1,n1,n2) - matrix(x2,n1,n2,byrow=T) )
  ww  <- apply(tmp,1,which.min)
  mm  <- apply(tmp,1,min)
  
  list(i2 = ww, m2 = mm)
}


compareEst <- function(variable,xlim=NULL,ylim=NULL){
  
  p1 <- postSummary[intersect(grep('beta_',rownames(postSummary)),grep(variable,rownames(postSummary))),]
  p2mu <- betaInd[grep(variable,rownames(betaInd)),]
  p2se <- betaIndSE[grep(variable,rownames(betaInd)),]
  
  if(is.null(xlim))xlim <- range(p1[,3:4])
  if(is.null(ylim))ylim <- range(c(p2mu-2*p2se,p2mu+2*p2se),na.rm=T)
  
  plot(p1[-nspec,1],p2mu[-nspec],xlim=xlim,ylim=ylim,xlab='This analysis',ylab='Traditional')
  
  for(s in 1:(nspec-1)){
    ws <- grep(specs[s],rownames(p1))
    lines(p1[ws,3:4],c(p2mu[s],p2mu[s]))
    lines(c(p1[ws,1],p1[ws,1]),c(p2mu[s]-1.96*p2se[2],p2mu[s]+1.96*p2se[2]))
  }
  abline(0,1,lty=2)
  title(variable)
  
}

makeRateMaps <- function(maplon=c(-98,-68),maplat=c(25,50),dx=1,dy=1,xbin=20,MAP=F){

  xgrid <- seq(maplon[1],maplon[2],by=dx)
  ygrid <- seq(maplat[1],maplat[2],by=dy)

  nx <- length(xgrid)
  ny <- length(ygrid)

  dcol <- c('dia1','dia2')
  
  allData <- numeric(0)

  for(k in 1:nspec){

    dxdt <- dldt <- dldx <- ij <- mask <- numeric(0)

    for(i in 1:(nx-1)){                                 #longitude

       dmat <- cmat <- jj <- numeric(0)

       for(j in 1:(ny-1)){                              #latitude

           wk <- which(data[,'lon'] >= xgrid[i] &
                       data[,'lon'] <  xgrid[i+1] &
                       data[,'lat'] >= ygrid[j] &
                       data[,'lat'] <  ygrid[j+1] &
                       data[,'sp'] ==  specs[k])
           if(length(wk) == 0)next

           dk <- data[wk,]

          tmp  <- histByYr(dk[,dcol],yr=c(1,2),xbin=xbin,dens=F,dt=1,weight=NULL)
          dj   <- matrix(tmp$dist[,1],nrow=1)
          cj   <- matrix(tmp$change[,1],nrow=1)
          colnames(dj) <- colnames(cj) <- tmp$size
          dmat <- appendMatrix(dmat,dj)
          cmat <- appendMatrix(cmat,cj)
  
          nij  <- nrow(dmat)
          jj   <- c(jj,j)
          rownames(dmat)[nij] <- rownames(cmat)[nij] <- paste('grid',i,j,sep='_')
       }

       if(length(jj) < 2)next

       dmat <- dmat[,order(as.numeric(colnames(dmat)))]
       cmat <- cmat[,order(as.numeric(colnames(cmat)))]

       ival <- rep(i,nrow(dmat))
       jval <- c(1:nrow(dmat))
       kval <- rep(k,nrow(dmat))
       ijk  <- cbind(kval,jval,ival)
       colnames(ijk) <- c('spec','lat','lon')

       allData <- appendMatrix(allData,cbind(ijk,dmat))

       dl <- cmat[-nrow(cmat),]             #dl/dt
       if(!is.matrix(dl))dl <- matrix(dl,nrow=1)

       dlx <- apply(dmat,2,diff,na.rm=T)/diff(jj)/dy
       if(!is.matrix(dlx))dlx <- matrix(dlx,nrow=1)

       ddd <- dl/dlx
       if(!is.matrix(ddd))ddd <- matrix(ddd,nrow=1)

       dxdt <- appendMatrix(dxdt,ddd)
       dldt <- appendMatrix(dldt,dl)
       dldx <- appendMatrix(dldx,dlx)

       dxdt <- dxdt[,order(as.numeric(colnames(dxdt)))]
       dldt <- dldt[,order(as.numeric(colnames(dldt)))]
       dldx <- dldx[,order(as.numeric(colnames(dldx)))]

       ii <- jj*0 + i
       ij <- rbind(ij,cbind(ii,jj)[-length(jj),])
    }

    rss <- which(rowSums(dxdt,na.rm=T) == 0)
    if(length(rss) > 0){
      dxdt <- dxdt[-rss,]
      dldt <- dldt[-rss,]
      dldx <- dldx[-rss,]
    }

    print(k)
    print(specs[k])
    
    if(!MAP)next


    xyi <- matrix( as.numeric( unlist( strsplit(rownames(dxdt),'_') )),ncol=3,byrow=T)[,2:3]

    mm <- matrix(0,nx,ny)
    mm[xyi] <- 1
    mask <- which(mm == 0,arr.ind=T)
    rk   <- which(mm == 1,arr.ind=T)

    xrange <- range(xgrid[rk[,1]])
    yrange <- range(ygrid[rk[,2]])

    if(!is.finite(xrange[1]) | !is.finite(xrange[2]) | 
       !is.finite(yrange[1]) | !is.finite(yrange[2]) ) next

    graphics.off()
    par(mfcol=c(3,3),bty='n',mar=c(1,1,2,1))

   mat <- dldt

   var <- c('dl/dt','dl/dx','dx/dt')

   for(ll in 1:3){

    if(ll == 2)mat <- dldx
    if(ll == 3)mat <- dxdt

    for(kk in 1:3){     #diameter classes

      xymat <- matrix(0,nx,ny)
      xymat[xyi] <- mat[,kk]
      xymat[is.na(xymat)] <- 0

      map('state',interior=F,col='grey',xlim=xrange, ylim=yrange)
      image(xgrid,ygrid,xymat,col=terrain.colors(100),add=T)
      contour(xgrid,ygrid,xymat,add=T,levels=c(-5,0,5))
      contour(xgrid,ygrid,xymat,add=T,levels=c(0),lwd=2)
      symbols(xgrid[mask[,1]],ygrid[mask[,2]],squares=rep(dx,nrow(mask)),
                 inches=F,fg='white',bg='white',add=T)
      map('state',boundary=F,col='grey',add=T)
      map('state',interior=F,add=T)
      title(paste(specs[k],colnames(dxdt)[kk],sep='-'))
    }
    text(xrange[2]-1,yrange[1]+.5,var[ll],pos=3)
  }

   dev.copy2pdf(file=paste('rateMap_',specs[k],'.pdf',sep=''))

  }

  allData
}

getH <- function(b,ss){

  tmat <- matrix(c(1,xmean['temp'],xmean['prec']),3,r)
  pmat <- matrix(c(1,xmean['prec'],xmean['temp']),3,r)

#  if(CENTER | CENSTAND){
#    tmat <- tmat*0
#    pmat <- pmat*0
#  }

  delfT <- apply(b[c(2,4,6),]*tmat,2,sum)
  delfP <- apply(b[c(3,5,6),]*pmat,2,sum)
  delf  <- cbind(delfT,delfP)

  t(delf)%*%invMat(ss)%*%delf
}

makeX1 <- function(temp,prec){

  x    <- cbind(temp,prec)
  xmean <- apply(x,2,mean)
  xsd   <- apply(x,2,sd)

  if(CENTER)  xx   <- (x - matrix(xmean,nrow(x),ncol(x),byrow=T))
  if(CENSTAND)xx   <- (x - matrix(xmean,nrow(x),ncol(x),byrow=T))/matrix(xsd,nrow(x),ncol(x),byrow=T) 

  x <- xx

  tempXtemp <- x[,1]^2
  precXprec <- x[,2]^2
  tempXprec <- x[,1]*x[,2]
  n   <- nrow(x)
  int <- rep(1,n)
  cbind(int,x,tempXtemp,precXprec,tempXprec)

}

yStats <- function(y){

  ytab <- rbind(apply(y,2,mean,na.rm=T),apply(y,2,range,na.rm=T))
  rownames(ytab) <- c('mean','min','max')

  yi   <- y
  yi[y == 0] <- NA

  yi  <- yi*0 + 1
  plots  <- apply(yi,2,sum,na.rm=T)
  signif(rbind(ytab,plots),3)
}
 


getSpecs <- function(SPECS=NULL, minA = 1, minI = 5){

  specs <- sp.tab[,'code']
  specs <- specs[!specs %in% notree]
  y <- as.matrix(ba.mat[,specs])

  ytab <- yStats(y)

  wy <- which(ytab['max',] > minA & ytab['plots',] > minI)

  if(!is.null(SPECS))wy <- intersect(wy,which(colnames(y) %in% SPECS))

  y[,wy]
}

fitLinear <- function(){

  yy <- y
  yy[yy == 0] <- NA
  yy <- apply(yy,2,var,na.rm=T)/3

  linearModsA <- matrix(NA,r,q+4)
  rownames(linearModsA) <- specs
  colnames(linearModsA) <- c('n','r^2','sigma','mspe',xnames)
  linearModsB <- linearModsA

  linearStandB <- linearModsB

  xt  <- seq(min(x[,'temp']),max(x[,'temp']),length=30)
  xt2 <- xt^2
  xp  <- seq(min(x[,'prec']),max(x[,'prec']),length=30)
  xp2 <- xp^2

  newTemp <- data.frame(x1 = xt, x2 = xp*0, x12 = xt2, x22 = xp2*0, x3 = xt*0)
  newPrec <- data.frame(x1 = xt*0, x2 = xp, x12 = xt2*0, x22 = xp2, x3 = xt*0)

  tempPredMu <- precPredMu <- matrix(NA,r,length(xt))
  tempLo <- tempHi <- precLo <- precHi <- tempSE <- precSE <- tempPredMu

  rtemp <- range(x[,'temp'])
  rprec <- range(x[,'prec'])

  tempmin <- tempmax <- precmin <- precmax <- rep(NA,r)

  fitSummary <- numeric(0)

  sigmaAll <- matrix(0,r,2)
  colnames(sigmaAll) <- c('sigma','mspe')

  predXtemp <- predXprec <- numeric(0)
  
  for(j in 1:r){

    wj   <- which(y[,j] > 0)    #fit non-zeros
 
    x1  <- x[wj,'temp']
    x2  <- x[wj,'prec']
    x12 <- x[wj,'temp']^2
    x22 <- x[wj,'prec']^2
    x3  <- x[wj,'temp']*x[wj,'prec']

    trange <- range(x1)
    prange <- range(x2)
    tempmin[j] <- trange[1]
    tempmax[j] <- trange[2]
    precmin[j] <- prange[1]
    precmax[j] <- prange[2]

    lfit <- lm(y[wj,j] ~ x1 + x2 + x12 + x22 + x3)

    tmp <- summary(lfit)

    pred <- predict.lm(lfit,interval='prediction')
    mspe <- mean( ((pred[,3] - pred[,1])/1.96)^2 )

    fitSummary <- append(fitSummary,list(tmp))
    names(fitSummary)[[j]] <- specs[j]
 
    linearModsB[j,]  <- c(length(wj),tmp$r.squared,tmp$sigma,mspe,tmp$coefficients[,1])

    xx1  <- x[,'temp']   #fit all data 
    xx2  <- x[,'prec']
    xx12 <- x[,'temp']^2
    xx22 <- x[,'prec']^2
    xx3  <- x[,'temp']*x[,'prec']

    lfita <- lm(y[,j] ~ xx1 + xx2 + xx12 + xx22 + xx3)
    tmpa <- summary(lfita)

    pred <- predict.lm(lfita,interval='prediction')
    mspe <- mean( ((pred[,3] - pred[,1])/1.96)^2 )
    sigmaAll[j,] <- c(tmpa$sigma,mspe)

    xs1 <- (x1 - rtemp[1])/(rtemp[2] - rtemp[1])   #fit standardized to range of x
    xs2 <- (x2 - rprec[1])/(rprec[2] - rprec[1])
    xs12 <- xs1^2
    xs22 <- xs2^2
    xs3  <- xs1*xs2

    tmp    <- lm(y[wj,j] ~ xs1 + xs2 + xs12 + xs22 + xs3)
    lstand <- summary(tmp)
    pred <- predict.lm(tmp,interval='prediction')
    mspe <- mean( ((pred[,3] - pred[,1])/1.96)^2 )

    linearStandB[j,] <- c(length(wj),lstand$r.squared,lstand$sigma,mspe,lstand$coefficients[,1])


    xtt         <- mean(x1)
    xpp         <- mean(x2)
    xtp         <- newTemp
    xtp[,'x2']  <- xpp
    xtp[,'x22'] <- xpp^2
    xtp[,'x3']  <- newTemp[,'x1']*xpp
    xpp <- newPrec
    xpp[,'x1']  <- xtt
    xpp[,'x12'] <- xtt^2
    xpp[,'x3']  <- newPrec[,'x2']*xtt
  
    te <- predict.lm(lfit,xtp, interval='confidence')
    pe <- predict.lm(lfit,xpp, interval='confidence')

    te[xtp[,'x1'] < trange[1],] <- 0
    te[xtp[,'x1'] > trange[2],] <- 0

    pe[xpp[,'x2'] < prange[1],] <- 0
    pe[xpp[,'x2'] > prange[2],] <- 0

    predXtemp <- append(predXtemp,list(xtp))
    predXprec <- append(predXprec,list(xpp))

    tempPredMu[j,] <- te[,1]
    tempLo[j,]     <- te[,2]
    tempHi[j,]     <- te[,3]
    precPredMu[j,] <- pe[,1]
    precLo[j,]     <- pe[,2]
    precHi[j,]     <- pe[,3]

    te <- predict.lm(lfit,newTemp, interval='prediction')
    pe <- predict.lm(lfit,newPrec, interval='prediction')

    tempSE[j,] <- (te[,'upr'] - te[,'lwr'])/1.96
    precSE[j,] <- (pe[,'upr'] - pe[,'lwr'])/1.96

    lfit <- lm(zg[wj,j] ~ x1 + x2 + x12 + x22 + x3)
    tmp <- summary(lfit)
    pred <- predict.lm(lfit,interval='prediction')
    mspe <- mean( ((pred[,3] - pred[,1])/1.96)^2 )
    linearModsA[j,]  <- c(length(wj),tmp$r.squared,tmp$sigma,mspe,tmp$coefficients[,1])
  }
  
  tempSE[is.na(tempSE)] <- 10
  precSE[is.na(precSE)] <- 100

  linearModsA[,-1]  <- signif(linearModsA[,-1],3)
  linearModsB[,-1]  <- signif(linearModsB[,-1],3)
  linearStandB[,-1] <- signif(linearStandB[,-1],3)

  predTempBA <- matrix(NA,r,30)
  rownames(predTempBA) <- specs
  predPrecBA <- predPrecLo <- predPrecHi <- predTempLo <- predTempHi <- predTempBA

  nsim <- 1000
  
  for(k in 1:30){
    
    tm <- pm <- matrix(0,nsim,r)
    wt <- which(tempmin > xt[k] | tempmax < xt[k])
    wp <- which(precmin > xp[k] | precmax < xp[k])

    for(m in 1:nsim){
      tk <- rnorm(r,tempPredMu[,k],tempSE[,k])
      tk[tk < 0] <- 0
      tk[wt] <- 0
      tm[m,] <- cumsum(tk)
      tk <- rnorm(r,precPredMu[,k],precSE[,k])
      tk[tk < 0] <- 0
      tk[wp] <- 0
      pm[m,] <- cumsum(tk)
    }

    t1 <- apply(tm,2,quantile,c(.5,.025,.975))
    predTempBA[,k] <- t1[1,]
    predTempLo[,k] <- t1[2,]
    predTempHi[,k] <- t1[3,]
    t1 <- apply(pm,2,quantile,c(.5,.025,.975))
    predPrecBA[,k] <- t1[1,]
    predPrecLo[,k] <- t1[2,]
    predPrecHi[,k] <- t1[3,]
  }

  cols <- mapColors(r)


  xtt <- xt*xscale['temp'] + xadd['temp']
  xpp <- xp*xscale['prec'] + xadd['prec']
  par(mfrow=c(2,2),bty='n')

  BAplot(xtt,predTempBA,'Temperature (C)','Mean basal area',xlimit=c(-10,15))
  BAplot(xtt,predTempHi,'Temperature (C)','95th quantile',xlimit=c(-10,15))
  BAplot(xpp,predPrecBA,'Precipitation (mm)','Mean basal area',xlimit=c(600,1600))
  BAplot(xpp,predPrecHi,'Precipitation (mm)','95th quantile',xlimit=c(600,1600))

 # dev.print(device=postscript,file='gradientPred.ps',width=6,horizontal=F)

  list(precPredMu = precPredMu, precLo = precLo, precHi = precHi,
       tempPredMu = tempPredMu, tempLo = tempLo, tempHi = tempHi,
       fitSummary = fitSummary, linearModsA = linearModsA, linearModsB = linearModsB, 
       sigmaAll = sigmaAll, predXtemp = predXtemp, predXprec = predXprec)
}

regPriors <- function(){
  
  beta <- matrix(0,nix,nspec)
  colnames(beta) <- specs
  rownames(beta) <- xnames
  loB <- beta - 1000
  hiB <- beta + 1000
  loB[xnames %in% posPrior,] <- 0
  hiB[xnames %in% negPrior,] <- 0
#  loB[xnames == 'intercept',] <- -5
#  hiB[xnames == 'intercept',] <- 5
  if('size2' %in% xnames)hiB[xnames == 'size2',] <- 0
  
  ww <- grep('temp',xnames)
  w2 <- grep('temp2',xnames)
  ww <- ww[!ww %in% w2]
  
  if(length(ww) > 0){
    loB[ww[1],] <- 0
    w2 <- grep('temp2',xnames)
    if(length(w2) > 0)hiB[w2,] <- 0
  }
  ww <- grep('prec',xnames)
  w2 <- grep('prec2',xnames)
  ww <- ww[!ww %in% w2]
  if(length(ww) > 0){
    loB[ww[1],] <- 0
    w2 <- grep('prec2',xnames)
    if(length(w2) > 0)hiB[w2,] <- 0
  }
  
  betaFec <- matrix(0,length(fnames),nspec)
  colnames(betaFec) <- specs
  rownames(betaFec) <- fnames
  loFec <- betaFec*0 - 1000
  hiFec <- betaFec*0 + 1000
  loFec[fnames %in% posPrior,] <- 0
  hiFec[fnames %in% negPrior,] <- 0
 # loFec[fnames == 'intercept',] <- -5
 # hiFec[fnames == 'intercept',] <- 5
  
  ww <- grep('temp',fnames)
  w2 <- grep('temp2',fnames)
  ww <- ww[!ww %in% w2]
  if(length(ww) > 0){
    loFec[ww[1],] <- 0
    w2 <- grep('temp2',fnames)
    if(length(w2) > 0)hiFec[w2,] <- 0
  }
  ww <- grep('prec',fnames)
  w2 <- grep('prec2',fnames)
  ww <- ww[!ww %in% w2]
  if(length(ww) > 0){
    loFec[ww[1],] <- 0
    w2 <- grep('prec2',fnames)
    if(length(w2) > 0)hiFec[w2,] <- 0
  }
  
  rmort <- logit( c(.0001,.1) )  #range of possible mort rates
  
  betaMort <- solve(crossprod(xall[,mnames]))%*%crossprod(xall[,mnames],
                    logit(make2d( thetaMat+.01 )))
  rownames(betaMort) <- mnames
  colnames(betaMort) <- specs
  
  loMmat <- betaMort*0 - 100
  hiMmat <- betaMort*0 + 100
  hiMmat[mnames %in% posPrior,] <- 0
  loMmat[mnames %in% c('plotBA','deficit'),] <- 0  #effects reversed for mortality
  hiMmat[mnames %in% c('plotBA','deficit'),] <- 1000
  
#  hiMmat[mnames %in% negPrior,] <- 0  
#  loMmat[mnames %in% posPrior,] <- 0
  
#  loMmat['intercept',] <- rmort[1]
#  hiMmat['intercept',] <- rmort[2]
  
  ww <- grep('temp',mnames)
  w2 <- grep('temp2',mnames)
  ww <- ww[!ww %in% w2]
  
  if(length(ww) > 0){
    hiMmat[ww,] <- 0
    if(length(w2) > 0)loMmat[w2,] <- 0
  }
  
  ww <- grep('prec',mnames)
  w2 <- grep('prec2',mnames)
  ww <- ww[!ww %in% w2]
  if(length(ww) > 0){
    loMmat[ww[1],] <- 0
    w2 <- grep('prec2',mnames)
    if(length(w2) > 0)hiMmat[w2,] <- 0
  }
  
  betaMort[betaMort < loMmat] <- loMmat[betaMort < loMmat]
  betaMort[betaMort > hiMmat] <- hiMmat[betaMort > hiMmat]
  
  loM <- loMmat
  hiM <- hiMmat
  
  betaGam <- matrix(0,length(gnames),nspec)
  colnames(betaGam) <- specs
  rownames(betaGam) <- gnames
  
  loGam <- betaGam*0 - 10000
  hiGam <- betaGam*0 + 10000
  
  loGam[gnames %in% posPrior,] <- 0
  hiGam[gnames %in% negPrior,] <- 0
  
  ww <- grep('temp',gnames)
  w2 <- grep('temp2',gnames)
  ww <- ww[!ww %in% w2]
  if(length(ww) > 0){
 #   loGam[ww[1],] <- 0
    w2 <- grep('temp2',gnames)
    if(length(w2) > 0)hiGam[w2,] <- 0
  }
  ww <- grep('prec',gnames)
  w2 <- grep('prec2',gnames)
  ww <- ww[!ww %in% w2]
  if(length(ww) > 0){
 #   loGam[ww[1],] <- 0
    w2 <- grep('prec2',gnames)
    if(length(w2) > 0)hiGam[w2,] <- 0
  }
  
  list(beta = beta, loB = loB, hiB = hiB,
       betaFec = betaFec, loFec = loFec, hiFec = hiFec,
       betaMort = betaMort, loM = loM, hiM = hiM,
       betaGam = betaGam, loGam = loGam, hiGam = hiGam)
}


  

BAplot <- function(predx,predy,xlabel=' ',ylabel=' ',xlimit){    #predy is spec by climate value

  cols <- mapColors(r)
  
  plot(predx,predy[1,],type='l',ylim=c(0,100),xlim=xlimit,xlab=xlabel,ylab=ylabel)
  for(j in 2:r){
     polygon(c(predx,rev(predx)),c(predy[j,],rev(predy[j-1,])),col=cols[j],border=cols[j])
  }

}

plotFromUnivariate <- function(j,vname){

  print( fitSummary[[j]] )

  par(mfrow=c(2,1),bty='n')

  minx <- min(x[y[,jplot] > 0,vname])
  maxx <- max(x[y[,jplot] > 0,vname])

  if(vname == 'temp'){
    xj <- predXtemp[[j]]
    pm <- tempprecPredMu[j,]
    pl <- tempLo[j,]
    ph <- tempHi[j,]
  }

  if(vname == 'prec'){
    xj <- predXprec[[j]]
    pm <- precPredMu[j,]
    pl <- precLo[j,]
    ph <- precHi[j,]
  }

  plot(x[,vname]*xsd[vname] + xadd[vname],y[,j],xlab='Precipitation (mm)',
        ylab='Basal area (m2/ha)',cex=.4,
        xlim=c(minx*xsd[vname] + xadd[vname],maxx*xsd[vname] + xadd[vname]))

  wk <- which(xj[,'x2'] > minx & xj[,'x2'] < maxx)
  lines(xj[wk,'x2']*xsd[vname] + xadd[vname],pm[wk],lwd=2)
  lines(xj[wk,'x2']*xsd[vname] + xadd[vname],pl[wk],lty=2)
  lines(xj[wk,'x2']*xsd[vname] + xadd[vname],ph[wk],lty=2)

  dev.print(device=postscript,file=paste('linearFit',specs[j],'.ps',sep=''),width=6,horizontal=F)
}

######################################################3
predY <- function(nsim=100,specs2Pred=specs,DATA=F,PRESENCE=F){  #if PRESENCE, assign zero where absent
  
  ws    <- match(specs2Pred,specs)
  rs    <- length(ws)

  gg    <- sample(burnin:ng,nsim,replace=T) 

  npred <- 20

 # tempbyspec <- apply(matrix(x[,'temp'],n,r)*y/matrix(apply(y,2,sum),n,r,byrow=T),2,sum)
 # precbyspec <- apply(matrix(x[,'prec'],n,r)*y/matrix(apply(y,2,sum),n,r,byrow=T),2,sum)

  tempgrad <- seq(min(x[,'temp']),max(x[,'temp']),length=npred)
  precgrad <- seq(min(x[,'prec']),max(x[,'prec']),length=npred)

  xtemp <- xprec <- matrix(0,npred,q)

  xtemp[,1] <- xprec[,1] <- 1
  xprec[,xnames == 'prec']        <- precgrad
  xprec[,xnames == 'precXprec']   <- precgrad^2
  xtemp[,xnames == 'temp']        <- tempgrad
  xtemp[,xnames == 'tempXtemp']   <- tempgrad^2

  ytemp <- yprec <- matrix(0,nsim,npred*rs)

  if(PRESENCE){
     temp0 <- x[,'temp']*y
     prec0 <- x[,'prec']*y
     temp0[temp0 == 0] <- NA
     prec0[prec0 == 0] <- NA
     tmin <- apply((temp0*0 + 1)*x[,'temp'],2,min,na.rm=T)
     tmax <- apply((temp0*0 + 1)*x[,'temp'],2,max,na.rm=T)
     pmin <- apply((prec0*0 + 1)*x[,'prec'],2,min,na.rm=T)
     pmax <- apply((prec0*0 + 1)*x[,'prec'],2,max,na.rm=T)
  }

  for(gk in 1:nsim){

     bg <- matrix(bgibbs[gg[gk],],q,r)
     sg <- matrix(sgibbs[gg[gk],],r,r)

     ag <- matrix(agibbs[gg[gk],],q,r)
     mg <- matrix(mgibbs[gg[gk],],r,r)

     mut <- xtemp%*%bg
     mup <- xprec%*%bg

     wt <- myrmvnorm(npred,mut,sg)
     wp <- myrmvnorm(npred,mup,sg)

     aut <- xtemp%*%ag
     aup <- xprec%*%ag

     at <- myrmvnorm(npred,aut,mg)
     ap <- myrmvnorm(npred,aup,mg)

     probt <- 1 - exp(-exp(at)*areaRareSpec)
     probp <- 1 - exp(-exp(ap)*areaRareSpec)

 #    wrt <- wrp <- matrix(0,npred,r)
 #    for(j in 1:rs){
 #       wrt[,j] <- conditionalMVNcdf(0,wt,mut,sigma=sg,ws[j])
 #       wrp[,j] <- conditionalMVNcdf(0,wp,mup,sigma=sg,ws[j])
 #    }

  #   pgg   <- matrix(pgibbs[gg[gk],],npred,r,byrow=T)
  #   probt <- 1 - (pgg + (1 - pgg)*wrt)
  #   probp <- 1 - (pgg + (1 - pgg)*wrp)

     wtt <- matrix(rbinom(npred*r,1,probt),npred,r)
     wpp <- matrix(rbinom(npred*r,1,probp),npred,r)

     yt <- (wt*wtt)[,ws]
     yp <- (wp*wpp)[,ws]

     yt[yt < 0] <- 0
     yp[yp < 0] <- 0

     ytemp[gk,] <- as.vector(yt)
     yprec[gk,] <- as.vector(yp)
  }

  mu <- matrix(apply(ytemp,2,mean),npred,rs)
  ci <- apply(ytemp,2,quantile,c(.025,.975))
  lo <- matrix(ci[1,],npred,rs)
  hi <- matrix(ci[2,],npred,rs)

  if(PRESENCE){
    tempMat <- matrix(tempgrad,npred,rs)
    wt      <- which(tempMat < matrix(tmin,npred,rs,byrow=T) | tempMat > matrix(tmax,npred,rs,byrow=T))
    mu[wt] <- 0
    hi[wt] <- 0
  }

  plot(tempgrad,mu[,1],type='l',ylim=c(0,10),xlab='Temperature',ylab='m2/ha')
  for(j in 1:rs){
    lines(tempgrad,mu[,j],type='l',lwd=2,col=j)
    lines(tempgrad,lo[,j],lty=2,col=j)
    lines(tempgrad,hi[,j],lty=2,col=j)
    if(DATA)points(x[,'temp'],y[,ws[j]],col=j)
  }

  tmean <- apply(mu,1,cumsum)
  thi   <- apply(hi,1,cumsum)

  mu <- matrix(apply(yprec,2,mean),npred,rs)
  ci <- apply(yprec,2,quantile,c(.025,.975))
  lo <- matrix(ci[1,],npred,rs)
  hi <- matrix(ci[2,],npred,rs)

  if(PRESENCE){
    precMat <- matrix(precgrad,npred,rs)
    wt      <- which(precMat < matrix(pmin,npred,rs,byrow=T) | precMat > matrix(pmax,npred,rs,byrow=T))
    mu[wt] <- 0
    hi[wt] <- 0
  }

  plot(precgrad,mu[,1],type='l',ylim=c(0,10),xlab='Precipitation',ylab=' ')
  for(j in 1:rs){
    lines(precgrad,mu[,j],type='l',lwd=2,col=j)
    lines(precgrad,lo[,j],lty=2,col=j)
    lines(precgrad,hi[,j],lty=2,col=j)
    if(DATA)points(x[,'prec'],y[,ws[j]],col=j)
  }

  legend('topright',specs2Pred,text.col=c(1:rs),bty='n')

  pmean <- apply(mu,1,cumsum)
  phi   <- apply(hi,1,cumsum)

  if(length(specs2Pred) > 10){
    par(mfrow=c(2,2),bty='n')
    BAplot(tempgrad + xmean['temp'],tmean,xlimit=c(-10,15))
    BAplot(tempgrad + xmean['temp'],thi,xlimit=c(-10,15))
    BAplot(precgrad + xmean['prec'],pmean,xlimit=c(600,1600))
    BAplot(precgrad + xmean['prec'],phi,xlimit=c(600,1600))
  }
 
}
   
#######################################

DIPold <- function(meanTP=c(0,0),nsim=10,scoreScale=c(-100,10)){

  par(mfrow=c(2,2),bty='n')

  tiny <- 1e-10

  ywt <- y
  ywt[y > 0] <- 1
  ywt <- apply(ywt,1,sum)

  priorMu <- 0
  priorV  <- 100000

  wmean <- w1/ntot
  wvr   <- w2/ntot - wmean^2
  wvr[wvr < tiny] <- tiny
  wsd   <- sqrt(wvr)

  scoreMat <- matrix(NA,length(main),1)

  scores <- matrix(NA,n,6)
  colnames(scores) <- as.vector(outer(c('mu','var','score'),c('temp','prec'),paste,sep='-'))

  xmean <- xvar <- xscore <- matrix(NA,n,nsim)

  gg   <- sample(burnin:ng,nsim,replace=T)

  for(j in 1:length(main)){

    for(gk in 1:nsim){

     A <- matrix(bgibbs[gg[gk],],q,r)
     sigma <- matrix(sgibbs[gg[gk],],r,r)
     sinv  <- invMat(sigma)

     qq <- match(main[j],xnames)
     qj <- match(get(main[j]),xnames)
     qi <- numeric(0)
     B  <- matrix(A[qq,],n,r,byrow=T)

     if(length(qj) > 0){  #interact with q
       for(jj in 1:length(qj)){
         B <- B + matrix(x[,qj[jj]]/x[,qq],n,1)%*%A[qj[jj],]
       }
     }
     qn <- c(1:q)[-c(qq,qj)]  # do not interact with q
 
     term1 <- matrix(0,n,r)
     for(jj in 1:length(qn))term1 <- term1 + matrix(x[,qn[jj]],n,1)%*%matrix(A[qn[jj],],1,r) 

     term2 <- matrix(x[,qq],n,1)%*%matrix(A[qq,],1,r)
     if(length(qj) > 0){
        for(jj in 1:length(qj)){
            term2 <- term2 + matrix(x[,qj[jj]]/x[,qq],n,1)%*%matrix(A[qj[jj],],1,r) 
        }
     }

     mu <- term1 + term2

     term2 <- matrix(x[,qq],n,r)*B

     for(i in 1:n){

        wimu <- rnorm(r,wmean[i,],wsd[i,])
     #   wimu <- wmean[i,]

        bi <- matrix(B[i,],1,r)
        Vi <- invMat(bi%*%sinv%*%t(bi) + 1/priorV)
        vi <- bi%*%sinv%*%(wimu - term1[i,]) + priorMu/priorV
        score     <- -(x[i,qq] - vi*Vi)^2/Vi - log(Vi)
        xmean[i,gk] <- vi*Vi
        xvar[i,gk]  <- Vi
        xscore[i,gk] <- score
     }
   }

   xmu   <- apply(xmean,1,mean)
   xvr  <- apply(xvar,1,mean)
   score <- apply(xscore,1,mean)

   if(j == 1)jcol <- c(1:3)
   if(j == 2)jcol <- c(4:6)
   scores[,jcol] <- cbind(xmu,xvr,score)

      if(j == length(main))xlabel <- 'Observed'
      plotObsPred(x[,qq],xmu,sqrt(xvr),xlabel=main[j],ylim=range(x[,qq]))
      abline(0,1,lty=2)
      abline(h=priorMu,lty=2)
      text(min(x[,qq]),max(xmu),main[j],pos=4)

      ss <- score
      ss[ss > scoreScale[2]] <- scoreScale[2]
      ss[ss < scoreScale[1]] <- scoreScale[1]
   #   hist(ss,breaks=seq(scoreScale[1],scoreScale[2],length=30),xlim=scoreScale,main=' ',xlab='Score')
   #   text(scoreScale[1],2,signif(mean(ss),2),pos=4)

      ylimit <- quantile(score,c(.05,.95))

      plotObsPred(x[,qq],score,nbin=20,xlabel=main[j],ylabel='score',ylim=ylimit)

   #   plot(x[,qq],score,xlab=main[j],ylab='score',ylim=ylimit)
      hh <- hist(x[,qq],nclass=100,plot=F)
      hy <- hh$density*diff(ylimit)/max(hh$density)*.2
      lines(hh$mids,hy + ylimit[1],type='s')

      hh <- hist(rep(x[,qq],ywt),nclass=100,plot=F)
      lines(hh$mids,hh$density*diff(ylimit)/max(hh$density)*.2+ylimit[1],type='s',col='orange')
      
      scoreMat[j,] <- mean(ss)
   }
  dev.print(device=postscript,file='dip.ps',width=6,horizontal=F)

  scores
}

#################################

DIP <- function(meanTP=c(0,0),scoreScale=c(-100,10)){

  par(mfrow=c(2,2),bty='n')

  tiny <- 1e-10

  ywt <- y
  ywt[y > 0] <- 1
  ywt <- apply(ywt,1,sum)

 # wmean <- w1/ntot
  wvr   <- w2/ntot - wmean^2
  wvr[wvr < tiny] <- tiny
  wsd   <- sqrt(wvr)

  scoreMat <- matrix(NA,length(main),1)

  xmean <- cbind( apply(xtemp,2,mean), apply(xprec,2,mean) )
  xvar  <- cbind( apply(xtemp,2,var), apply(xprec,2,var) )

  colnames(xmean) <- colnames(xvar) <- c('temp','prec')

  wlo   <- which(xvar < tiny)
  xvar[wlo] <- tiny

  score <- getScoreNorm(x[,c('temp','prec')],xmean,xvar)

  score[wlo] <- NA

  scores <- cbind(xmean[,'temp'],xvar[,'temp'],score[,'temp'],xmean[,'prec'],xvar[,'prec'],score[,'prec'])
  colnames(scores) <- as.vector(outer(c('mu','var','score'),c('temp','prec'),paste,sep='-'))

  plotObsPred(x[,'temp'],xmean[,'temp'],sqrt(xvar[,'temp']),
              xlabel='temp',ylim=range(x[,'temp']))
  abline(0,1,lty=2)
  abline(h=priorXmu,lty=2)
  text(min(x[,'temp']),max(x[,'temp']),'temp',pos=4)

  ylimit <- quantile(score[,'temp'],c(.1,1),na.rm=T)

  plotObsPred(x[,'temp'],score[,'temp'],nPerBin=30,xlabel='temp',
              ylabel='score',ylim=ylimit)

      hh <- hist(x[,'temp'],nclass=100,plot=F)
      hy <- hh$density*diff(ylimit)/max(hh$density)*.2
      lines(hh$mids,hy + ylimit[1],type='s')

      hh <- hist(rep(x[,'temp'],ywt),nclass=100,plot=F)
      lines(hh$mids,hh$density*diff(ylimit)/max(hh$density)*.2+ylimit[1],type='s',col='orange')

  plotObsPred(x[,'prec'],xmean[,'prec'],sqrt(xvar[,'prec']),
              xlabel='prec',ylim=range(x[,'prec']))
  abline(0,1,lty=2)
  abline(h=priorXmu,lty=2)
  text(min(x[,'prec']),max(x[,'prec']),'prec',pos=4)

  ylimit <- quantile(score[,'prec'],c(.1,1),na.rm=T)

  plotObsPred(x[,'prec'],score[,'prec'],nPerBin=45,xlabel='prec',
              ylabel='score',ylim=ylimit)

      hh <- hist(x[,'prec'],nclass=100,plot=F)
      hy <- hh$density*diff(ylimit)/max(hh$density)*.2
      lines(hh$mids,hy + ylimit[1],type='s')

      hh <- hist(rep(x[,'prec'],ywt),nclass=100,plot=F)
      lines(hh$mids,hh$density*diff(ylimit)/max(hh$density)*.2+ylimit[1],type='s',col='orange')


  dev.print(device=postscript,file='dip.ps',width=6,horizontal=F)


  scores
}






predStates <- function(){

  par(mfrow=c(2,2),bty='n')

  ytab <- yStats(y) 

  ymean <- y1/ntot
  ysd   <- sqrt(y2/ntot - ymean^2)

 # plotObsPred(y,ymean,ysd,xlabel='Observed BA',ylim=range(y))
  plotObsPred(y,ymean,nbin=20,xlabel='Observed BA',ylim=range(y))
  abline(0,1,lty=2)
  title('all species')

  rmse <- rep(0,r)

  plot(y[,1],ymean[,1],xlim=c(0,14),ylim=c(0,12),xlab='Observed',ylab='Predicted',cex=.1)
  for(j in 1:r){
    yj <- y[,j]
    wj <- which(yj > 0)
    rmse[j] <- sqrt( sum( (yj[wj] - ymean[wj,j])^2) /length(wj) )
    points(y[,j],ymean[,j],col=j,cex=.1)
    
  }
  abline(0,1,lty=2)
  title('individual species')

  ytab <- rbind(ytab,rmse)

  wmean <- w1/ntot
  wsd   <- sqrt(w2/ntot - wmean^2)

  plotObsPred(ymean,wmean,nbin=20,xlabel='predicted y',ylabel='estimated w',ylim=range(y))
  abline(0,1,lty=2)

  rmean <- r1/ntot

  plotObsPred(wmean[y == 0],rmean[y == 0],nbin=20,xlabel='w',ylabel='r|y = 0')
  abline(h=0,lty=2)

  dev.print(device=postscript,file='states.ps',width=6,horizontal=F)

  ytab

}

updateB <- function(bnow,lo,hi,yy,sigma,smallv=NULL,bigv=NULL){  # update b's for mvnorm
  
  if(is.null(bigv))  bigv   <- invMat(crossprod(x))      #multivariate
  if(is.null(smallv))smallv <- crossprod(x,yy)
  mu     <- bigv%*%smallv
  vaa    <- kronecker(sigma,bigv)   
  b      <- matrix(myrmvnorm(1,as.vector(mu),sigma=vaa),q,r,byrow=F)
 # b <- matrix(tnorm.mvt(as.vector(bnow),as.vector(mu),vaa,lo,hi,times=1),q,r) #if prior is truncated (lo, hi)
 
  b
}

updateXold <- function(){
  
  xsd <- propXTP
  xsd[sample(length(propXTP),20)] <- .1

  propp <- tnorm(n,loPrec,hiPrec,xg[,'prec'],xsd[,'prec'])
  propt <- tnorm(n,loTemp,hiTemp,xg[,'temp'],xsd[,'temp'])

  propp <- propp - mean(propp)
  propt <- propt - mean(propt)

  xnew  <- xg

  xnew[,'temp'] <- propt
  xnew[,'prec'] <- propp
  xnew[,'tempXtemp'] <- propt^2
  xnew[,'precXprec'] <- propp^2
  xnew[,'tempXprec'] <- propt*propp

  munow <- xg%*%bg
  munew <- xnew%*%bg

  pnow <- dmvnormZeroMean(w - munow,sg) + dmvnormZeroMean(zg - xg%*%ag,amat) +
          dnorm(xg[,'temp'],priorXmu,sqrt(priorXvr),log=T) +
          dnorm(xg[,'prec'],priorXmu,sqrt(priorXvr),log=T) 

  pnew <- dmvnormZeroMean(w - munew,sg) + dmvnormZeroMean(zg - xnew%*%ag,amat) +
          dnorm(xnew[,'temp'],priorXmu,sqrt(priorXvr),log=T) +
          dnorm(xnew[,'prec'],priorXmu,sqrt(priorXvr),log=T) 

  aa <- exp(pnew - pnow)
  z  <- runif(n,0,1)

  accept <- rep(0,n)
  accept[z < aa] <- 1

  xg[z < aa,] <- xnew[z < aa,]
  xg[,'temp'] <- xg[,'temp'] - mean(xg[,'temp'])
  xg[,'prec'] <- xg[,'prec'] - mean(xg[,'prec'])
  xg[,'tempXtemp'] <- xg[,'temp']^2
  xg[,'precXprec'] <- xg[,'prec']^2
  xg[,'tempXprec'] <- xg[,'temp']*xg[,'prec']

  list(xg = xg, accept = accept)
}

proposeX <- function(){

  xsd <- propXTP
  xsd[sample(length(propXTP),20)] <- .1

  propp <- tnorm(n,loPrec,hiPrec,xg[,'prec'],xsd[,'prec'])
  propt <- tnorm(n,loTemp,hiTemp,xg[,'temp'],xsd[,'temp'])

  propp <- propp - mean(propp)
  propt <- propt - mean(propt)

  xnew  <- xg

  xnew[,'temp'] <- propt
  xnew[,'prec'] <- propp
  xnew[,'tempXtemp'] <- propt^2
  xnew[,'precXprec'] <- propp^2
  xnew[,'tempXprec'] <- propt*propp
  xnew
}

updateXnew <- function(){
  
  xnew <- proposeX()

  pnow <- dmvnormZeroMean(w - xg%*%bg,sg) + dmvnormZeroMean(zg - xg%*%ag,amat) +
          dnorm(xg[,'temp'],priorXmu,sqrt(priorXvr),log=T) +
          dnorm(xg[,'prec'],priorXmu,sqrt(priorXvr),log=T) 

  pnew <- dmvnormZeroMean(w - xnew%*%bg,sg) + dmvnormZeroMean(zg - xnew%*%ag,amat) +
          dnorm(xnew[,'temp'],priorXmu,sqrt(priorXvr),log=T) +
          dnorm(xnew[,'prec'],priorXmu,sqrt(priorXvr),log=T) 

  aa <- exp(pnew - pnow)
  z  <- runif(n,0,1)

  accept <- rep(0,n)
  accept[z < aa] <- 1

  xg[z < aa,] <- xnew[z < aa,]
  xg[,'temp'] <- xg[,'temp'] - mean(xg[,'temp'])
  xg[,'prec'] <- xg[,'prec'] - mean(xg[,'prec'])
  xg[,'tempXtemp'] <- xg[,'temp']^2
  xg[,'precXprec'] <- xg[,'prec']^2
  xg[,'tempXprec'] <- xg[,'temp']*xg[,'prec']

  list(xg = xg, accept = accept)
}


updateZ <- function(LIKE){  #log intensity

  tiny <- 1e-15

  sprop <- rexp(n*r,50)

  propz <- matrix(rnorm(n*r,zg,sprop),n,r)

  psnow <- exp(-exp(zg)*area)  
  psnew <- exp(-exp(propz)*area)  

  wr <- w*0
  mu <- x%*%bg

  for(j in 1:r){
     wr[,j] <- conditionalMVNcdf(0,w,mu,sigma=sg,j)
  }

  psnow <- psnow/(psnow + (1 - psnow)*wr)  # Pr that rr = 1 (absent)
  psnew <- psnew/(psnew + (1 - psnew)*wr)

  psnow[psnow < tiny] <- tiny
  psnow[psnow > (1 - tiny)] <- 1 - tiny
  psnew[psnew < tiny] <- tiny
  psnew[psnew > (1 - tiny)] <- 1 - tiny

  pnow1 <- dmvnormZeroMean(zg - x%*%ag,amat) 
  if(LIKE == 'pois') pnow2 <- dpois(ntreeAll,area*exp(zg),log=T) #eqn A.5
  if(LIKE == 'binom')pnow2 <- dbinom(zpresent,1,1 - exp(-exp(zg)*area))

  pnow4 <- rr*log(psnow) + (1 - rr)*log(1 - psnow)  #g_is in Appendix
  pnow4[y > 0] <- 0

  pnew1 <- dmvnormZeroMean(propz - x%*%ag,amat)     #proposed z
  if(LIKE == 'pois') pnew2 <- dpois(ntreeAll,area*exp(propz),log=T)
  if(LIKE == 'binom')pnew2 <- dbinom(zpresent,1,1 - exp(-exp(propz)*area))
  pnew4 <- rr*log(psnew) + (1 - rr)*log(1 - psnew)
  pnew4[y > 0] <- 0

  pnow5 <- rowSums(pnow2 + pnow4)
  pnew5 <- rowSums(pnew2 + pnew4)

  pnow <- pnow1 + pnow5
  pnew <- pnew1 + pnew5

  a <- exp(pnew - pnow)
  z <- runif(n,0,1)
  zg[z < a,] <- propz[z < a,]
  zg
}

updateDensity <- function(z,x,y,rr,area,sg,alpha,amat){  #log intensity for Poisson sampling
  
  #area  - vector of sample areas
  #x     - predictors n by p
  #y     - response n by r
  #z     - current log density
  #sg    - r by r covariance matrix for abundance
  #alpha - p by r coefficient matrix
  #amat  - r by r covariance for log density
  #rr    - n by r absence indicators
  
  tiny <- 1e-15
  
#  area  <- tarea + sarea
  
  r <- ncol(y)
  n <- nrow(y)
  p <- ncol(x)
  
  sprop <- rexp(n*r,50)
  
  propz <- matrix(rnorm(n*r,z,sprop),n,r)
  
  psnow <- exp(-exp(z)*area)  
  psnew <- exp(-exp(propz)*area)  
  
  wr <- w*0
  mu <- x%*%bg
  
  for(j in 1:r){
    wr[,j] <- conditionalMVNcdf(0,w,mu,sigma=sg,j)
  }
  
  psnow <- psnow/(psnow + (1 - psnow)*wr)  # Pr that r = 1 (absent)
  psnew <- psnew/(psnew + (1 - psnew)*wr)
  
  psnow[psnow < tiny] <- tiny
  psnow[psnow > (1 - tiny)] <- 1 - tiny
  psnew[psnew < tiny] <- tiny
  psnew[psnew > (1 - tiny)] <- 1 - tiny
  
  pnow1  <- dmvnormZeroMean(z - x%*%alpha,amat) 
  pnow2  <- dpois(y,area*exp(z),log=T) #eqn A.5
  pnow4  <- rr*log(psnow) + (1 - rr)*log(1 - psnow)  #g_is in Appendix
  pnow4[y > 0] <- 0
  
  pnew1  <- dmvnormZeroMean(propz - x%*%alpha,amat)  #proposed z
  pnew2  <- dpois(y,area*exp(propz),log=T)
  pnew4  <- rr*log(psnew) + (1 - rr)*log(1 - psnew)
  pnew4[y > 0] <- 0
  
  pnow5 <- rowSums(pnow2 + pnow4)
  pnew5 <- rowSums(pnew2 + pnew4)
  
  pnow <- pnow1 + pnow5
  pnew <- pnew1 + pnew5
  
  a <- exp(pnew - pnow)
  zz <- runif(n,0,1)
  z[zz < a,] <- propz[zz < a,]
  z
}

  

updateW <- function(){  #w is drawn from a MVN, values > 0 set equal to y

  mu <- x%*%bg
  ww <- w

  tiny <- 1e-10

  for(k in 1:r){

     lok <- rep(-Inf,n)
     hik <- rep(Inf,n)

     absent  <- which(rr[,k] == 1)
     present <- which(rr[,k] == 0)

     hik[absent]  <- 0
     lok[present] <- 0

    testv <- try(chol(sg[-k,-k]),T)
    if(inherits(testv,'try-error')){
      return( list(mu = numeric(0), vr = numeric(0)) )
    }

    sin <- chol2inv(testv)
    p1  <- sg[-k,k]%*%sin

    muk <- mu[,k] + p1%*%t(ww[,-k] - mu[,-k])
    vrk <- sg[k,k] - p1%*%sg[k,-k]
    if(vrk < tiny)vrk <- tiny

     ww[,k]  <- tnorm(n,lok,hik,muk,sqrt(vrk))
  }

  ww[rr == 0] <- y[rr == 0]
  ww
}


updateP <- function(){   #Pr species is absent

  u1 <- colSums(rr)
  u2 <- n - u1
  rbeta(r,u1 + 1,u2 + 1)
}

updateR <- function(muw){
  
  tiny <- 1e-6
  
  wr <- w*0              # Pr(w < 0)
  for(j in 1:r)wr[,j] <- conditionalMVNcdf(0,w,muw,sigma=sg,j)
  
  wr[wr < tiny] <- tiny
  
  rr[y > 0] <- 0
  
#  theta <- pz/(pz + (1 - pz)*wr)  # g_is in Appendix
  
  theta <- wr/(wr + (1 - wr)*pz)  # g_is in Appendix
  
  w0 <- which(y == 0)
  rr[w0] <- rbinom(n*r,1,theta)[w0]
  
  rr
}

updateR1 <- function(muw){
  
  tiny <- 1e-6

  wr <- w*0              # Pr(w < 0)
  for(j in 1:r)wr[,j] <- conditionalMVNcdf(0,w,muw,sigma=sg,j)
  
  wr[wr < tiny] <- tiny

  rr[y > 0] <- 0

  theta <- pz/(pz + (1 - pz)*wr)  # g_is in Appendix
    
  w0 <- which(y == 0)
  rr[w0] <- rbinom(n*r,1,theta)[w0]

  rr
}

rhoMLE <- function(rho){  #for optimization

  tiny <- 1e-10
  ee <- exp(-rho*D)
  ll <- -(1 - S)*rho*D + S*log(1 - ee + tiny)
  sum(ll)
}

fiaTimes <- function(xx){  #fia matrix has columns 'date_t1','date_t2','plt1'
  
  pnames <- sort(as.character(unique(xx[,'plt1'])))
  
  t11 <- as.Date(xx[,'date_t1'])
  t1 <- as.numeric(format(t11,'%Y'))
  
  t22 <- as.Date(xx[,'date_t2'])
  t2 <-as.numeric(format(t22,'%Y'))
  
  t11 <- aggregate(t1,by=list(xx[,'plt1']),min,na.rm=T)
  t1 <- t11[[2]][ match(t11[[1]],pnames) ]
  
  t11 <- aggregate(t2,by=list(xx[,'plt1']),min,na.rm=T)
  t2 <- t11[[2]][ match(t11[[1]],pnames) ]
  dt <- t2 - t1
  
  lon <- aggregate(xx[,'lon'],by=list(xx[,'plt1']),min,na.rm=T)
  lon <- lon[[2]][ match(lon[[1]],pnames) ]
  
  lat <- aggregate(xx[,'lat'],by=list(xx[,'plt1']),min,na.rm=T)
  lat <- lat[[2]][ match(lat[[1]],pnames) ]
  
  yy <- cbind(t1,t2,dt,lon,lat)
  rownames(yy) <- pnames
  yy
}



groByDiamSpec <- function(){           #not used
  
  #demography plots
  diam <- apply(demDiam,1,mean,na.rm=T)              #mean diameter by tree
  tmp  <- matrix( unlist(strsplit(rownames(demDiam),'-')),ncol=2,byrow=T)
  tmp[!tmp[,2] %in% specs,2] <- 'other'
  si   <- match(tmp[,2],specs)
  ji   <- match(tmp[,1],plotnames)
  ti   <- as.numeric(matrix( unlist( strsplit(colnames(demDiam),'diam') ),ncol=2,byrow=2)[,2])
  tmat <- matrix(ti,nrow(demDiam),ncol(demDiam),byrow=T)
  t1   <- apply(tmat*(demDiam*0 + 1),1,min,na.rm=T)
  t2   <- apply(tmat*(demDiam*0 + 1),1,max,na.rm=T)
  
  wdead  <- deathyrDem[,'died'] + 1
  wd     <- which(is.finite(wdead))
  
  t2[wd] <- years[wdead[wd]]
  
  dmat  <- interpRows(demDiam,startIndex=match(t1,years),endIndex=match(t2,years),
                      INCREASING=T,minVal=0,maxVal=Inf,
                      defaultValue=NULL,tinySlope=.001)
  
  
  
  
  #fia
  fcode <- as.character(treeCodes[match(fiaDiam[,'spec'],treeCodes[,'fiaCode']),'code'])
  fcode[!fcode %in% specs] <- 'other'

  diam <- c( diam,rowMeans(fiaDiam[,c('diam.dia1','diam.dia2')],na.rm=T) )
  si  <- c( si,match(fcode,specs) )
  ji  <- c( ji,match(fiaDiam[,'plot'],plotnames) )
  dinc <- c(demDinc,fiaDiam[,'dinc'])
  
  tt1 <- c(t1,fiaDiam[,'yr.t1'])
  tt2 <- c(t2,fiaDiam[,'yr.t2'])
  
  dincAll <- numeric(0)
  
  wdead <- deathyrDem[,2] + 1
  wd    <- which(is.finite(wdead))
  
  fcode <- as.character(treeCodes[match(fiaDiam[,'spec'],treeCodes[,'fiaCode']),'code'])
  fcode[!fcode %in% specs] <- 'other'
  
  tmp  <- matrix( unlist(strsplit(rownames(demDiam),'-')),ncol=2,byrow=T)
  tmp[!tmp[,2] %in% specs,2] <- 'other'
  si   <- match(tmp[,2],specs)
  ji   <- match(tmp[,1],plotnames)
  
  smat <- matrix(si,nrow(demDiam),nt)
  jmat <- matrix(ji,nrow(demDiam),nt)
  
  tmat  <- matrix(1:nt,nrow(demDiam),nt,byrow=T)
  t1    <- apply(tmat*(demDiam*0 + 1),1,min,na.rm=T)
  t2    <- apply(tmat*(demDiam*0 + 1),1,max,na.rm=T)
  
  lmat  <- demDiam*0
  lmat[cbind(1:nrow(demDiam),t1)] <- 1
  lmat[cbind(wd,wdead[wd] + 1)] <- -1
  
  lmat[is.na(lmat)] <- 0
  lmat  <- t(apply(lmat,1,cumsum))     #live
  omat  <- lmat*0
  omat[cbind(wd,wdead[wd] + 1)] <- 1   #dead
  
  www <- which(lmat == -1 & omat == 1,arr.ind=T)
  wwz <- www
  wwz[,2] <- wwz[,2] + 1
  lmat[www] <- 1
  omat[www] <- 0
  omat[wwz] <- 1
  
  www <- which(lmat == -1 & omat == 0,arr.ind=T)
  lmat[www] <- 0
  
  dmat  <- interpRows(demDiam,startIndex=t1,endIndex=t2,
                      INCREASING=T,minVal=0,maxVal=Inf,
                      defaultValue=NULL,tinySlope=.001)
  imat <- diam2dinc(dmat,minInc=0,firstTime=t1)
  imat <- cbind(rep(NA,nrow(imat)),imat)
  imat <- imat[,-ncol(imat)]
  colnames(imat) <- years
  
  ww   <- which(is.finite(imat[,-nt]),arr.ind=T)
  
  t1 <- tmat[ww]
  t2 <- t1 + 1
  jj <- jmat[ww]
  ss <- smat[ww]
  dd <- dmat[ww]
  ii <- imat[ww]
  
  ll <- lmat[ww]
  oo <- omat[ww]
  
  ifia <- fiaDiam[,'dinc']
  ifia[ifia < .00001] <- .00001
  
  
  dd <- c( dd,rowMeans(fiaDiam[,c('diam.dia1','diam.dia2')],na.rm=T) )
  ss  <- c( ss,match(fcode,specs) )
  jj  <- c( jj,match(fiaDiam[,'plot'],plotnames) )
  ii <- c(ii,ifia)
  
  t1 <- c(t1,match(fiaDiam[,'yr.t1'],years))
  t2 <- c(t2,match(fiaDiam[,'yr.t2'],years))
  
  kmat <- matrix(0,nt,nbreak)
  
  gmat <- livemat <- deadmat <- y*0
  
  for(j in 1:nplot){
    
    wj <- which(jtdex[,'j'] == j)
    db <- breakMat[wj[1],1:nbreak]
    
    wji <- which(ji == j)
    
    diall <- cbind(rep(j,length(wji)),si[wji],tt1[wji],tt2[wji],diam[wji],dinc[wji])
    dincAll <- rbind(dincAll,diall)
    
    for(s in 1:nspec){
      
      ww <- which(jj == j & ss == s)
      if(length(ww) == 0)next
      
      kk <- findInterval(dd[ww],db) + 1
      kk[kk > (nbreak-1)] <- nbreak - 1
      
      inc <- c(ii[ww],ii[ww])
      ttt <- c(t1[ww],t2[ww])
      
      tmp <- byFunctionRcpp(inc,ttt,kk,kmat*0,kmat*0,MEAN=T)
      
      gmat[wj,sdex == s] <- tmp[jtdex[wj,'t'],]
      
      if(j <= nplotDem){
        
        live <- byFunctionRcpp(ll[ww],ttt,kk,kmat*0,kmat*0,MEAN=F)
        dead <- byFunctionRcpp(oo[ww],ttt,kk,kmat*0,kmat*0,MEAN=F)
        livemat[wj,sdex == s] <- live[jtdex[wj,'t'],]
        deadmat[wj,sdex == s] <- dead[jtdex[wj,'t'],]
      }
      
    }
  }
  
  colnames(dincAll) <- c('j','spec','t1','t2','diam','dinc')
  
  deadmat <- rbind(deadmat[jtdex[,'j'] <= nplotDem,],deadmatFIA)
  livemat <- rbind(livemat[jtdex[,'j'] <= nplotDem,],livematFIA)
  
  list(groMat = gmat, liveMat = livemat, deadMat = deadmat, dincAll = dincAll)
}



groByDiamSpecOld <- function(){
  
  diam <- apply(demDiam,1,mean,na.rm=T)
  tmp  <- matrix( unlist(strsplit(rownames(demDiam),'-')),ncol=2,byrow=T)
  tmp[!tmp[,2] %in% specs,2] <- 'other'
  si   <- match(tmp[,2],specs)
  ji   <- match(tmp[,1],plotnames)
  ti   <- as.numeric(matrix( unlist( strsplit(colnames(demDiam),'diam') ),ncol=2,byrow=2)[,2])
  tmat <- matrix(ti,nrow(demDiam),ncol(demDiam),byrow=T)
  t1   <- apply(tmat*(demDiam*0 + 1),1,min,na.rm=T)
  t2   <- apply(tmat*(demDiam*0 + 1),1,max,na.rm=T)
  
  fcode <- as.character(treeCodes[match(fiaDiam[,'spec'],treeCodes[,'fiaCode']),'code'])
  fcode[!fcode %in% specs] <- 'other'
  
  
  diam <- c( diam,rowMeans(fiaDiam[,c('diam.dia1','diam.dia2')],na.rm=T) )
  si  <- c( si,match(fcode,specs) )
  ji  <- c( ji,match(fiaDiam[,'plot'],plotnames) )
  dinc <- c(demDinc,fiaDiam[,'dinc'])
  
  tt1 <- c(t1,fiaDiam[,'yr.t1'])
  tt2 <- c(t2,fiaDiam[,'yr.t2'])
  
  dincAll <- numeric(0)
  
  wdead <- deathyrDem[,2] + 1
  wd    <- which(is.finite(wdead))

  fcode <- as.character(treeCodes[match(fiaDiam[,'spec'],treeCodes[,'fiaCode']),'code'])
  fcode[!fcode %in% specs] <- 'other'

  tmp  <- matrix( unlist(strsplit(rownames(demDiam),'-')),ncol=2,byrow=T)
  tmp[!tmp[,2] %in% specs,2] <- 'other'
  si   <- match(tmp[,2],specs)
  ji   <- match(tmp[,1],plotnames)

  smat <- matrix(si,nrow(demDiam),nt)
  jmat <- matrix(ji,nrow(demDiam),nt)

  tmat  <- matrix(1:nt,nrow(demDiam),nt,byrow=T)
  t1    <- apply(tmat*(demDiam*0 + 1),1,min,na.rm=T)
  t2    <- apply(tmat*(demDiam*0 + 1),1,max,na.rm=T)
  
  lmat  <- demDiam*0
  lmat[cbind(1:nrow(demDiam),t1)] <- 1
  lmat[cbind(wd,wdead[wd] + 1)] <- -1
  
  lmat[is.na(lmat)] <- 0
  lmat  <- t(apply(lmat,1,cumsum))     #live
  omat  <- lmat*0
  omat[cbind(wd,wdead[wd] + 1)] <- 1   #dead
  
  www <- which(lmat == -1 & omat == 1,arr.ind=T)
  wwz <- www
  wwz[,2] <- wwz[,2] + 1
  lmat[www] <- 1
  omat[www] <- 0
  omat[wwz] <- 1
  
  www <- which(lmat == -1 & omat == 0,arr.ind=T)
  lmat[www] <- 0
  
  dmat  <- interpRows(demDiam,startIndex=t1,endIndex=t2,
                      INCREASING=T,minVal=0,maxVal=Inf,
                      defaultValue=NULL,tinySlope=.001)
  imat <- diam2dinc(dmat,minInc=0,firstTime=t1)
  imat <- cbind(rep(NA,nrow(imat)),imat)
  imat <- imat[,-ncol(imat)]
  colnames(imat) <- years
  
  ww   <- which(is.finite(imat[,-nt]),arr.ind=T)
  
  t1 <- tmat[ww]
  t2 <- t1 + 1
  jj <- jmat[ww]
  ss <- smat[ww]
  dd <- dmat[ww]
  ii <- imat[ww]
  
  ll <- lmat[ww]
  oo <- omat[ww]
  
  ifia <- fiaDiam[,'dinc']
  ifia[ifia < .00001] <- .00001
  
  
  dd <- c( dd,rowMeans(fiaDiam[,c('diam.dia1','diam.dia2')],na.rm=T) )
  ss  <- c( ss,match(fcode,specs) )
  jj  <- c( jj,match(fiaDiam[,'plot'],plotnames) )
  ii <- c(ii,ifia)
  
  t1 <- c(t1,match(fiaDiam[,'yr.t1'],years))
  t2 <- c(t2,match(fiaDiam[,'yr.t2'],years))
  
  kmat <- matrix(0,nt,nbreak)
  
  gmat <- livemat <- deadmat <- y*0
  
  for(j in 1:nplot){
    
    wj <- which(jtdex[,'j'] == j)
    db <- breakMat[wj[1],1:nbreak]
    
    wji <- which(ji == j)
    
    diall <- cbind(rep(j,length(wji)),si[wji],tt1[wji],tt2[wji],diam[wji],dinc[wji])
    dincAll <- rbind(dincAll,diall)
    
    for(s in 1:nspec){
      
      ww <- which(jj == j & ss == s)
      if(length(ww) == 0)next
      
      kk <- findInterval(dd[ww],db) + 1
      kk[kk > (nbreak-1)] <- nbreak - 1
      
      inc <- c(ii[ww],ii[ww])
      ttt <- c(t1[ww],t2[ww])
      
      tmp <- byFunctionRcpp(inc,ttt,kk,kmat*0,kmat*0,MEAN=T)
      
      gmat[wj,sdex == s] <- tmp[jtdex[wj,'t'],]
      
      if(j <= nplotDem){
        
        live <- byFunctionRcpp(ll[ww],ttt,kk,kmat*0,kmat*0,MEAN=F)
        dead <- byFunctionRcpp(oo[ww],ttt,kk,kmat*0,kmat*0,MEAN=F)
        livemat[wj,sdex == s] <- live[jtdex[wj,'t'],]
        deadmat[wj,sdex == s] <- dead[jtdex[wj,'t'],]
      }
      
    }
  }
  
  colnames(dincAll) <- c('j','spec','t1','t2','diam','dinc')
  
  deadmat <- rbind(deadmat[jtdex[,'j'] <= nplotDem,],deadmatFIA)
  livemat <- rbind(livemat[jtdex[,'j'] <= nplotDem,],livematFIA)
  
  list(groMat = gmat, liveMat = livemat, deadMat = deadmat, dincAll = dincAll)
}


groByDiamSpecOld <- function(){
  
  diam <- apply(demDiam,1,mean,na.rm=T)
  tmp  <- matrix( unlist(strsplit(rownames(demDiam),'-')),ncol=2,byrow=T)
  tmp[!tmp[,2] %in% specs,2] <- 'other'
  si   <- match(tmp[,2],specs)
  ji   <- match(tmp[,1],plotnames)
  ti   <- as.numeric(matrix( unlist( strsplit(colnames(demDiam),'diam') ),ncol=2,byrow=2)[,2])
  tmat <- matrix(ti,nrow(demDiam),ncol(demDiam),byrow=T)
  t1   <- apply(tmat*(demDiam*0 + 1),1,min,na.rm=T)
  t2   <- apply(tmat*(demDiam*0 + 1),1,max,na.rm=T)
  
  fcode <- as.character(treeCodes[match(fiaDiam[,'spec'],treeCodes[,'fiaCode']),'code'])
  fcode[!fcode %in% specs] <- 'other'
  
  
  diam <- c( diam,rowMeans(fiaDiam[,c('diam.dia1','diam.dia2')],na.rm=T) )
  si  <- c( si,match(fcode,specs) )
  ji  <- c( ji,match(fiaDiam[,'plot'],plotnames) )
  dinc <- c(demDinc,fiaDiam[,'dinc'])
  
  t1 <- c(t1,fiaDiam[,'yr.t1'])
  t2 <- c(t2,fiaDiam[,'yr.t2'])
  
  gmat <- matrix(NA,ntnp,nbreak*nspec)
  amat <- matrix(0,nspec,nbreak)
  
  dincAll <- numeric(0)
  
  for(j in 1:nplot){
    
    wj <- which(ji == j)
    gj <- which(jtdex[,'j'] == j)
    dx <- breakMat[gj[1],1:nbreak]
    di <- findInterval(diam[wj],dx)
    di[di == 0] <- 1
    
    tmp <- byFunction(dinc[wj],si[wj],di,amat*0,mean)
    tmp <- as.vector(t(tmp))
    
    gmat[gj,] <- matrix(tmp,length(gj),nbreak*nspec,byrow=T)
    
    diall <- cbind(rep(j,length(wj)),si[wj],t1[wj],t2[wj],diam[wj],dinc[wj])
    dincAll <- rbind(dincAll,diall)
    
  }
  
  colnames(dincAll) <- c('j','spec','t1','t2','diam','dinc')
  
  list(groMat = gmat, dincAll = dincAll)
}

clim4FIA <- function(lat,lon,plVec,yrVec,nplot=max(plVec),
                     monTemp=NULL,monPrec=NULL,monThermIndex=c(1:12),
                     start=1969,end=2012,GINDEX=F,ONEYRVEC=F,LST=T){
  
  # if ONEYRVEC there is one yrVec that applies to all plots
  # if !ONEYRVEC there is a vector containing a year for each lat,lon
  
  tmp <- inClimate(lat=lat,lon=lon,start=start,end=end,LST=LST,PADLAST=T)
  temp <- tmp$temp
  prec <- tmp$prec
  
  end  <- max( round( as.numeric(colnames(temp),0 ) ) )
  npp <- max(plVec)
  nyy <- end - start + 1
  yiIndex <- start:end
  
  gindex <- numeric(0)
  
  monthNames <- getMonthNames()
  
  regClim <- matrix(NA,length(plVec),2)
  colnames(regClim) <- c('temp','prec')
  
  if(GINDEX){
    yvec   <- rep(yiIndex,each=12)
    mvec   <- rep(1:12,nyy)
    ccols  <- as.character(round( yvec + mvec*.01,2 ))
 
    wshort <- which(nchar(ccols) == 6)
    if(length(wshort) > 0)ccols[wshort] <- paste(ccols[wshort],'0',sep='')
    ccols  <- ccols[ccols %in% colnames(temp)]
    tmp    <- monthlyPHr(yi=yvec,mi=mvec,
                         tempMatrix=temp[,ccols],precMatrix=prec[,ccols],lat=lat)
    gindex <- tmp$degHrPosYr
    gindex <- gindex[,yvec - start + 1]
    dindex <- tmp$degHrNegYr
    dindex <- dindex[,yvec - start + 1]
    regClim <- matrix(NA,length(plVec),4)
    colnames(regClim) <- c('temp','prec','therm','deficit')
  }
  
  
  lastClimYr <- max( floor( as.numeric(colnames(temp)) ) )
  
  monTempIndex <- monPrecIndex <- c(1:12)
  
  if(!is.null(monTemp)){
    moVec <- numeric(0)
    s1 <- seq(1,nchar(monTemp),by=3)
    s2 <- s1 + 2
    monTempIndex <- match( substring(monTemp,s1 , s2 ), monthNames)
  }
  if(!is.null(monPrec)){
    moVec <- numeric(0)
    s1 <- seq(1,nchar(monPrec),by=3)
    s2 <- s1 + 2
    monPrecIndex <- match( substring(monPrec,s1 , s2 ), monthNames)
  }
  
  cCols <- matrix( unlist( strsplit(colnames(temp),'[.]') ), ncol=2,byrow=T)
  yrCol <- as.numeric(cCols[,1])
  moCol <- as.numeric(cCols[,2])
  
  for(j in 1:npp){
    
    wj <- which(plVec == j)   #yrs for plot j
    nw <- length(wj)
    if(nw == 0)next
    
    if(ONEYRVEC) yjj <- yrVec
    if(!ONEYRVEC)yjj <- yrVec[wj]
    
    if(nw == 1){       #only one row for plot
      
      yjj[yjj > end] <- end
      
      wt <- which( yrCol %in% yjj & moCol %in% monTempIndex )
      wp <- which( yrCol %in% yjj & moCol %in% monPrecIndex )
      wg <- which( yrCol %in% yjj & moCol %in% c(1:12) )
      
      tj <- temp[j,wt]
      pj <- prec[j,wp]
      
      meanTemp <- byIndex(tj,yrCol[wt],mean,na.rm=T)
      meanPrec <- byIndex(pj,yrCol[wp],sum,na.rm=T)
      
      regClim[wj,'temp']    <- mean( meanTemp, na.rm=T)
  #    regClim[wj[k],'temp'] <- meanTemp[length(meanTemp)]
      
      regClim[wj,'prec'] <- mean( meanPrec, na.rm=T)
  #    regClim[wj[k],'prec'] <- meanPrec[length(meanPrec)]
      
      if(GINDEX){
  #      wg <- which( yiIndex %in% yjj  )
  #      gj <- gindex[j,wg]
        
        gj <- gindex[j,wt]
        gj <- byIndex(gj,yrCol[wp],mean,na.rm=T)
        regClim[wj,'therm']    <- mean(gj)
        
        dj <- dindex[j,wt]
        dj <- byIndex(dj,yrCol[wp],mean,na.rm=T)
        regClim[wj,'deficit']    <- mean(dj)
        
  #      regClim[wj[k],'therm'] <- gj[length(gj)]
      }
      
      next
    }
      
    
    if(nw > 1){
    
      for(k in 2:nw){
        
        kk <- wj[k-1]:wj[k]
        yj <- yrVec[kk]
        yj[yj >= lastClimYr] <- lastClimYr
        
        wt <- which( yrCol %in% yj & moCol %in% monTempIndex )
        wp <- which( yrCol %in% yj & moCol %in% monPrecIndex )
        wg <- which( yrCol %in% yj & moCol %in% c(1:12) )
        
        tj <- temp[j,wt]
        pj <- prec[j,wp]
        
        meanTemp <- byIndex(tj,yrCol[wt],mean,na.rm=T)
        meanPrec <- byIndex(pj,yrCol[wp],sum,na.rm=T)
        
        regClim[kk,'temp']    <- mean( meanTemp, na.rm=T)
        regClim[wj[k],'temp'] <- meanTemp[length(meanTemp)]
        
        regClim[kk,'prec'] <- mean( meanPrec, na.rm=T)
        regClim[wj[k],'prec'] <- meanPrec[length(meanPrec)]
        
        if(GINDEX){
     #     yj[yj > end] <- end
     #     wg <- which( yiIndex %in% yj  )
     #     gj <- gindex[wj[k],wg]
          
          gj <- gindex[j,wg]
          gj <- byIndex(gj,yrCol[wg],mean,na.rm=T)
          regClim[kk,'therm']    <- gj
          regClim[wj[k],'therm'] <- gj[length(gj)]
          
          dj <- dindex[j,wg]
          dj <- byIndex(dj,yrCol[wg],mean,na.rm=T)
          regClim[kk,'deficit']    <- dj
          regClim[wj[k],'deficit'] <- dj[length(dj)]
        }
      }
    } #end if nw > 1
      
  }
  
  if(!is.null(monTemp)){
    colnames(regClim)[colnames(regClim) == 'temp'] <- paste('temp',monTemp,sep='-')
  }
  if(!is.null(monPrec)){
    colnames(regClim)[colnames(regClim) == 'prec'] <- paste('prec',monPrec,sep='-')
  }
               
  list(regClim = regClim, gindex = gindex)
}
  
  
loadFIA <- function(maplon=NULL,maplat=NULL,mapPolygon=NULL,centers=NULL){
  
  #maplon, maplat each have two values for boundaries
  #mapPolygon has two rows, first row is xcorners, 2nd row is ycorners
  
  if(is.null(mapPolygon)){
    xcorners <- maplon[c(1,2,2,1)]
    ycorners <- maplat[c(1,1,2,2)]
    mapPolygon <- rbind(xcorners,ycorners)
  }
  
  require(sp)
  
  load('tree/all.rdata')
  
  # tmpData <- read.csv('tree/remeasured.csv',header=T,na.strings="NULL")
  
  tmp <- point.in.polygon(remeas.tre[,'lon'],remeas.tre[,'lat'],
                           mapPolygon[1,],mapPolygon[2,])
  wm <- which(tmp > 0)
    
  
  if(!is.null(maplon)){
    wm  <- which(remeas.tre[,'lon'] > maplon[1] & 
                 remeas.tre[,'lon'] < maplon[2] &
                 remeas.tre[,'lat'] > maplat[1] & 
                 remeas.tre[,'lat'] < maplat[2])
  }
  if(!is.null(centers)){
    #   near <- nn2(tmp[,c('lon','lat')],centers,
    #                searchtype='radius',radius=buffer,k=ceiling(buffer*10000))
    near <- nn2(remeas.tre[,c('lon','lat')],centers,k=ceiling(buffer*5000))
    www  <- which(near[[2]] <= buffer,arr.ind=T)
    wm   <- sort(unique(near[[1]][near[[2]] <= buffer]))
  }
  
  slope  <- (remeas.tre[wm,'slope_t1'] + remeas.tre[wm,'slope_t2'])/2
  aspect <- (remeas.tre[wm,'aspect_t1'] + remeas.tre[wm,'aspect_t2'])/2
  
  data <- remeas.tre[wm,c('plt1','eco_sub','lon','lat','date_t1','date_t2','sp',
                          'dia1','dia2','physclcd_t2','stat1','stat2')]
  data <- cbind(data,slope,aspect)
  
  ww <- which( data[,c('dia1','dia2')] == 0 ,arr.ind=T)
  
  if(length(ww) > 0)data[,c('dia1','dia2')][ww] <- NA
  data
}


inPlotFIA <- function(maplon=NULL,maplat=NULL,mapPolygon=NULL,
                      numSpec=NULL,minMaxBA=NULL,nplots=NULL,
                      centers=NULL,buffer=NULL){

  # requires either (maplon,maplat) for bounds or centers (lat,lon) and buffer
  # centers  - (lat, lon) for current plots
  # buffer   - all FIA plots within buffer degrees from centers
  
  require(RANN)
  require(maps)
  require(maptools)
  
  fnames <- as.vector( t(outer(specs,1:nbreak,paste,sep='-')) )
  
  if(is.null(mapPolygon)){
    xcorners <- maplon[c(1,2,2,1)]
    ycorners <- maplat[c(1,1,2,2)]
    mapPolygon <- rbind(xcorners,ycorners)
  }

  treeCodes <- read.table('../allocationmodel/datafiles/treeCodesDuke.txt',header=T)
  
  data <- loadFIA(mapPolygon=mapPolygon,centers)
  
  fspecs <- fiaSpecs[!fiaSpecs == 'other']
  
  tdata <- fiaTimes(data)
  
  ww <- which(tdata[,'dt'] < 1)
  if(length(ww) > 0){
    wd <- which(data[,'plt1'] %in% rownames(tdata)[ww])
    data <- data[-wd,]
    tdata <- tdata[-ww,]
  }
  
  wws <- which(data[,'sp'] %in% fiaSpecs)
  dd  <- data[wws,c('dia1','dia2')]
  
  t1 <- as.Date(data[wws,'date_t1'])
  t1 <- as.numeric(format(t1,'%Y'))
  t2 <- as.Date(data[wws,'date_t2'])
  t2 <- as.numeric(format(t2,'%Y'))
  
  tmp <- paste(as.character(data[wws,'plt1']),data[wws,'sp'],sep='-')
  dinc <- (dd[,2] - dd[,1])/(t2 - t1)
#  rownames(dd) <- tmp
  diamAll <- data.frame(plot=data[wws,'plt1'], spec = data[wws,'sp'], 
                   yr = cbind(t1,t2), diam = as.matrix(dd), dinc = dinc)
  
  
  recruits <- which(is.na(data[,'dia1']) & is.finite(data[,'dia2']) &
                      data[,'dia2'] < minTree)  #new on sapling plots
  
  pnames <- rownames(tdata)
  pindex <- match(pnames,data[,'plt1'])
  eco    <- data[pindex,'eco_sub']
  lon    <- data[pindex,'lon']
  lat    <- data[pindex,'lat']
  phys   <- data[pindex,'physclcd_t2']
  slope  <- data[pindex,'slope']
  aspect <- data[pindex,'aspect']
                 
  plotData <- data.frame( list(cbind(tdata,lon,lat,phys,slope,aspect),eco = eco) )
  
  if(!REMOTE){
    data(countyMapEnv)
    data(stateMapEnv)
    graphics.off()
    mapscale <- min(c(diff(maplon),diff(maplat)))/2
    mapSetup(maplon,maplat,scale=mapscale)
  
    map('state',lwd=1,xlim=maplon,ylim=maplat)
    points(plotData[,'lon'],plotData[,'lat'],cex=.2)
    polygon(mapPolygon[1,],mapPolygon[2,],lwd=3,border='grey')
    
    print('number of FIA plots')
    print(nrow(plotData))
  }
  
  data[,c('date_t1','date_t2')] <- tdata[match(data[,'plt1'],rownames(tdata)),c('t1','t2')]
  
  dxm <- getBreaks(dxSeq,minSap,maxTree,tdata[,'dt'],fia=T)
  dxx <- dxm$dx
  bxx <- dxm$br
                 
  np <- nrow(tdata)
  
  yy <- jtdexFIA <- numeric(0)
  
  liveMat <- deadMat <- deadInt <- numeric(0)
  growMu <-   growVr <- numeric(0)
  
  rvec <- rep(0,length(specs))
  names(rvec) <- specs
  rfia <- rother <- numeric(0)
  
  growth <- survive <- numeric(0)
  
  ws <- match(data[,'sp'],treeCodes[,'fiaCode'])
  dspec <- as.character(treeCodes[ws,'code'])
  dspec[!dspec %in% specs] <- 'other'
  data[,'sp'] <- dspec
  
  for(j in 1:np){
    
    bj <- bxx[j,]
    wj <- which(data[,'plt1'] == pnames[j])
    wr <- recruits[recruits %in% wj]
    rvec <- rvec*0
    
    #recruits, only second time
    if(length(wr) > 0){
      ij <- data[wr,'sp']
      tmp <- table(ij)
      names(tmp)[!names(tmp) %in% specs] <- 'other'
      rt  <- tmp[names(tmp) %in% specs]
      rvec[ match(names(rt),specs) ] <- rt
    }
    rfia <- rbind(rfia,rvec*0,rvec)
    
    #first time
    f1  <- wj[is.finite(data[wj,'dia1']) & data[wj,'stat1'] == 1]
    w1  <- findInterval(data[f1,'dia1'],bj)
    yy1 <- dat2Tab(dvec=data[f1,'dia1']*0+1,ispec=data[f1,'sp'],jj=w1,fnames,fiaSpecs,nbreak,sum)
    
    #second time
    f2  <- wj[is.finite(data[wj,'dia2']) & data[wj,'stat2'] == 1]
    w2  <- findInterval(data[f2,'dia2'],bj)
    yy2 <- dat2Tab(data[f2,'dia2']*0+1,data[f2,'sp'],w2,fnames,fiaSpecs,nbreak,sum)
    
    tj <- data[wj[1],c('date_t1','date_t2')]
    jj <- rep(j,2)
    jtdexFIA <- rbind(jtdexFIA,cbind(jj,match(tj,years)))
    
    yj <- rbind(yy1,yy2)
    rownames(yj) <- paste(pnames[j],tj,sep='-')
    yy <- rbind(yy,yj)
    
    #survive: 1st & 2nd
    wsurv <- intersect(f1,f2)
    w3    <- findInterval(data[wsurv,'dia1'],bj)
    surv  <- dat2Tab(data[wsurv,'dia1']*0+1,data[wsurv,'sp'],w3,fnames,fiaSpecs,nbreak,sum)
    surv  <- rbind(surv,yy2)
    rownames(surv) <- paste(pnames[j],tj,sep='-')
    
    survive <- rbind(survive,surv)
    
    # growth: 1st & 2nd, increase bin
    wgrow <- wj[is.finite(data[wj,'dia1']) & data[wj,'stat1'] == 1 &
                is.finite(data[wj,'dia2']) & data[wj,'stat2'] == 1 &
                data[wj,'dia2'] > data[wj,'dia1']]

    grow  <- dat2Tab(dvec=data[wj[wgrow],'dia1']*0+1,data[wj[wgrow],'sp'],
                     w1[wgrow],fnames,fiaSpecs,nbreak,sum)
    grow  <- rbind(grow,grow)
    rownames(grow) <- paste(pnames[j],tj,sep='-')
    growth <- rbind(growth,grow)
    
    w3   <- findInterval(data[wsurv,'dia1'],bj)
    incr <- (data[wsurv,'dia2'] - data[wsurv,'dia1'])/
            (data[wsurv,'date_t2'] - data[wsurv,'date_t1'])
    
    gmu <- dat2Tab(dvec=incr,data[wsurv,'sp'],w3,fnames,fiaSpecs,nbreak,mean)
    gvr <- dat2Tab(dvec=incr,data[wsurv,'sp'],w3,fnames,fiaSpecs,nbreak,var)
    
    gMu <- rbind(gmu,gmu)
    rownames(gMu) <- paste(pnames[j],tj,sep='-')

    gVr <- rbind(gvr,gvr)
    rownames(gVr) <- paste(pnames[j],tj,sep='-')
    
    growMu <- rbind(growMu,gMu)
    growVr <- rbind(growVr,gVr)
    
    #death: 1st, not 2nd
    wl <- which(is.finite(data[wj,'dia1']) & data[wj,'dia1'] > 0 & 
                is.finite(data[wj,'dia2']) & data[wj,'dia2'] > 0)
    wd <- which(is.finite(data[wj,'dia1']) & data[wj,'dia1'] > 0 & 
               (!is.finite(data[wj,'dia2'] | data[wj,'dia2'] == 0)))
    lage  <- data[wj[wl],'date_t2'] - data[wj[wl],'t1']
    dage  <- data[wj[wd],'date_t2'] - data[wj[wd],'t1']
    lspec <- data[wj[wl],'sp']
    dspec <- data[wj[wd],'sp']
    
    ls <- dat2Tab(data[wj[wl],'dia1']*0+1,data[wj[wl],'sp'],w1[wl],fnames,fiaSpecs,nbreak,sum)
    ds <- ls*0
    
    if(length(wd) > 0){
      ds <- dat2Tab(data[wj[wd],'dia1']*0+1,data[wj[wd],'sp'],w1[wd],fnames,fiaSpecs,nbreak,sum)
    }
    
    liveMat <- rbind(liveMat,ls)
    deadMat <- rbind(deadMat,ds)
    deadInt <- c(deadInt,lage[1])
    
  }
  
  recruit <- rfia
  colnames(recruit) <- specs
  
  colnames(jtdexFIA) <- c('j','t')
  
  survive[survive > yy] <- yy[survive > yy]
  
  rm(remeas.tre)
  
  growMu[yy == 0] <- growVr[yy == 0] <- NA
  growVr[ (growMu != 0 & growVr == 0) | (growMu > 0 & is.na(growVr)) ] <- 1
  
  
  list(plotData = plotData, y = yy, recruit = recruit, jtdex = jtdexFIA, diamAll = diamAll,
       breaks = bxx, dxMat = dxx, liveMat = liveMat, deadMat = deadMat, deadInt = deadInt,
       survive = survive, growth = growth, growMu = growMu, growVr = growVr)
}
  

dat2Tab <- function(dvec,ispec,jj,fnames,fiaSpecs,nbreak,FUN){
  
  FUN <- match.fun(FUN)
  
  dtmp <- matrix(0,length(specs),nbreak)
  
 # ispec[!ispec %in% fiaSpecs] <- 'other'
  wspec <- match(ispec,specs)
  
  wj <- which(is.finite(jj))
  dvec <- dvec[wj]
  wspec <- wspec[wj]
  jj   <- jj[wj]
  
  xx <- byFunction(dvec,wspec,jj,dtmp,FUN)
  xx <- as.vector(t(xx))
  names(xx) <- fnames
  xx
}

  

getBreaks <- function(dxSeq,minD,maxD,meanInc,fia=F){    
  
  #requires reference dxSeq from main program, which has length nbreak-1
  
  nr <- length(meanInc)
  nc <- nbreak-1
  
  imat <- matrix(dxSeq,nr,nc,byrow=T) * matrix(meanInc,nr,nc)
  imat <- cbind(rep(minD,nr),imat)
  br   <- t(apply(imat,1,cumsum))
  br[,nbreak] <- maxD
  
  if(fia){
    fb <- 2.54*5
    ww <- apply( (fb - br)^2 ,1,which.min)
    br[cbind(1:nr,ww)] <- fb
  }
    
  dj <- t(apply(br,1,diff))
  dj[,ncol(dj)] <- dj[,ncol(dj)-1]
  dj <- cbind(dj,dj[,ncol(dj)])
  colnames(dj) <- colnames(br) <- c(1:nbreak)
  
  list(breaks = br, dx = dj)
}







inPlotDem <- function(minDiam=0,SPECS=NULL,years=c(1992:2013)){  
  #input demography plots


  dataPath <- '../allocationmodel/datafiles/'
  treePath  <- paste(dataPath,'trees',sep='')
  hydroPath <- paste(dataPath,'hydro/',sep='')

  pNames <- character(0)
  nt     <- length(years)

  sapl <- tree <- baTab <- plotTab <- sizeBySpec <- 
    survBySpec <- growBySpec <- sizeSpecData <- numeric(0)
  
  growMuSpec <-growVrSpec <- growNoSpec <- numeric(0)
  
  deathage <- censage <- numeric(0)
  utm <- numeric(0)
  startYr <- numeric(0)
  
  pData <- getPlotNames(dataPath,treeOnly=T,regions=regions)

  plotnames <- pData$plotnames
  plotData  <- pData$data
  
  alt <- rep(NA,nrow(pData$data))
  
  ww <- which(is.finite(pData$data[,'lon']))
  
  if(length(ww) > 0){
    
    xy <- expand.grid(topo$x,topo$y)
    wx <- nn2(xy,pData$data[ww,c('lon','lat')],k=1)$nn.idx
    
    wr <- match(xy[wx,1],topo$x)
    wc <- match(xy[wx,2],topo$y)
    alt[ww] <- topo$z[cbind(wr,wc)]
    plotData[!is.finite(plotData[,'elev']),'elev'] <- alt[!is.finite(plotData[,'elev'])]
  }

  tfiles <- list.files(treePath,pattern='_active.txt')   #in tree directory
  backup <- grep("~",tfiles)
  if(length(backup) > 0)tfiles <- tfiles[-backup]
  tplots <- filename2plot(tfiles)

  wj <- which(tplots %in% plotnames)
  tplots <- tplots[wj]
  tfiles <- tfiles[wj]

  nfile <- length(tplots)
 # wide <- sqrt(sampleArea*20000)

  kp <- jp <- 0
  
  treePlotArea   <- read.table(paste(dataPath,'treePlotAreaDuke.txt',sep=''),header=T)
  tmp    <- colnames(treePlotArea)[-1]
  treeYr <- as.numeric( matrix( unlist(strsplit(tmp,'X')),ncol=2,byrow=T)[,2] )
  
  tmp <- as.numeric(matrix( unlist( strsplit( colnames(treePlotArea)[-1],'X' ) ),
                            ncol=2,byrow=T)[,2] )
  
  if(min(tmp) > years[1])years <- c(min(tmp):max(years))
  
  sampAreaYr <- numeric(0)
  
  censusYears <- deathyr <- numeric(0)
  plotSummary <- numeric(0)
  dxMat <- breakMat <- numeric(0)
  jtdex <- numeric(0)
  
  recruit <- numeric(0)
  diamAll <- dincAll <- numeric(0)
  
  
  for(j in 1:nfile){

    jplot <- unlist(strsplit(tfiles[j],'_'))[1:2]
    jplot <- paste(jplot[1],jplot[2],sep='_')

  #  print(j)
    print(jplot)

    ff     <- tfiles[j]
    fname  <- paste(treePath,ff,sep='/')
    ball   <- read.table(fname,header=T,sep="\t",fill=T)
    ball[,'species'] <- as.character( ball[,'species'] )
    
    wu <- grep('UNKN',ball[,'species'])
    if(length(wu) > 0)ball <- ball[-wu,]
    
    if(exists('rename')){
      for(kk in 1:nrow(rename)){
        wk <- which(ball[,'species'] == rename[kk,1])
        if(length(wk) > 0)ball[wk,'species'] <- rep( as.character( rename[kk,2] ), 
                                                     length(wk) )
      }
    }

    dtmp    <- dateColumn(b=ball,vname='diam',yrvec=years)
    diamMat <- dtmp$x 
    dyr     <- dtmp$yr
    
    yrStats <- c(range(dyr),length(dyr))
    
    #has diam measurements
    tmp  <- apply(diamMat,1,min,na.rm=T)
    ww   <- which(is.finite(tmp))
    ball <- ball[ww,]
    diamMat <- diamMat[ww,]

    #first/last by censoring/death
    first   <- apply(ball[,c('growinyr','censinyr')],1,min,na.rm=T)
    fmin    <- min(first,na.rm=T)
    last    <- ball[,'deathyr']
    if('censoryr' %in% colnames(ball)){
      last    <- apply(ball[,c('deathyr','censoryr')],1,min,na.rm=T)
    }
    died    <- ball[,'deathyr']

    #censoring errors & last measurements for survivors
    tm <- matrix(years,nrow(diamMat),ncol(diamMat),byrow=T)
    tm <- (diamMat*0 + 1)*tm
    fd <- apply(tm,1,min,na.rm=T)
    first[!is.finite(first)] <- fd[!is.finite(first)]
    ld <- apply(tm,1,max,na.rm=T)
    last[!is.finite(last)] <- ld[!is.finite(last)]
    
    first[first < years[1]] <- years[1]
    last[last > max(years)] <- max(years)
    
    #censored last sample
    mm <- max(ld,na.rm=T)

    rr <- range( (last - first),na.rm=T)
    if(rr[2] < 1)next
    
    #ingrowth too late
    wf <- which(is.finite(ball[,'growinyr']) & ball[,'growinyr'] > fd)
    if(length(wf) > 0)ball[wf,'growinyr'] <- first[wf] <- fd[wf]
    
    #measurements taken after death
    ww <- which(ld > died)
    if(length(ww) > 0){
      ld[is.finite(died) & ld > died]   <- died[is.finite(died) & ld > died]
      last[is.finite(died) & ld > died] <- ld[is.finite(died) & ld > died]
    }

    
    #new recruits by specs/year
    rmat <- matrix(0,length(dyr),length(SPECS))   
    colnames(rmat) <- SPECS
    rownames(rmat) <- paste(jplot,dyr,sep='-')
    
    rf <- byIndex(ball[,'growinyr']*0 + 1,
                  list(yr=ball[,'growinyr'],spec=ball[,'species']),sum,na.rm=T)
    
    ww <- which(is.finite(ball[,'growinyr']) & ball[,'growinyr'] %in% dyr
                & ball[,'species'] %in% SPECS)
    sw <- match(ball[ww,'species'],SPECS)
    xm <- matrix(0,length(dyr),length(SPECS))
    rf <- byFunctionRcpp(ball[ww,'growinyr']*0 + 1,match(ball[ww,'growinyr'],dyr),sw,
                         xm,xm,MEAN=F)
    rownames(rf) <- dyr
    colnames(rf) <- SPECS
    
    if(length(rf) == 0)rm <- rmat
    if(length(rf) > 0){
      rn <- rownames(rf)
      r1 <- rf[,colnames(rf) %in% SPECS]
      if(!is.matrix(r1)){
        cn <- names(r1)
        r1 <- matrix(r1,1)
        rownames(r1) <- rn
        colnames(r1) <- cn
      }
      rmat[match(rownames(r1),dyr),match(colnames(r1),SPECS)] <- r1
      wrr    <- which(!colnames(rf) %in% SPECS)
      rother <- rowSums(matrix(rf[,wrr],nrow(rf)))
      other <- rep(0,length(dyr))
      other[match(rownames(rf),dyr)] <- rother
      rmat <- cbind(rmat,other)
    }
    
    recruit <- appendMatrix(recruit,rmat)
    
    

    #more than 0 yr
    ww <- which( (ld - fd) > 0 & (last - first) > 0)  
    if(length(ww) == 0)next
    ball    <- ball[ww,]
    diamMat <- diamMat[ww,]
    first   <- first[ww]
    last    <- last[ww]
    fd      <- fd[ww]
    ld      <- ld[ww]
    died    <- died[ww]
    
    #mean diameter increment by tree
    diamIncr <- diamMat[cbind( 1:length(ww),match(ld,colnames(diamMat)) )] - 
                diamMat[cbind( 1:length(ww),match(fd,colnames(diamMat)))]
    diamIncr <- diamIncr/(ld - fd)
    
    
    
    startYr <- c(startYr,min(first,na.rm=T))
    names( startYr )[length(startYr)] <- jplot
    

    diamMax <- diamMat
    if(ncol(diamMat) > 1)diamMax <- apply(diamMat,1,max,na.rm=T)
    
    sp <- as.character(ball[,'species'])
    sp[!sp %in% SPECS] <- 'other'
    
    jjj <- rep(jplot,length(sp))
    dd  <- diamMat
    colnames(dd) <- paste('diam',colnames(dd),sep='')
    rownames(dd) <- paste(jjj,sp,sep='-')
    
    diamAll <- appendMatrix(diamAll,dd)   #tree by yr
    dincAll <- c(dincAll,diamIncr)   #mean increment for tree
    cens    <- ld
    last    <- apply(cbind(died,cens),1,min,na.rm=T)

    diamMat <- interpRows(x=diamMat,startIndex=match(first,years),
                          endIndex=match(last,years),INCREASING=F,minVal=0)
    
    if('UTMx' %in% colnames(ball)){
      ha    <- xy2area(ball[,c('UTMx','UTMy')])
      areaj <- ha$area/10000
      hullj <- ha$hull
      
      ux <- mean(ball[,'UTMx'],na.rm=T)
      uy <- mean(ball[,'UTMy'],na.rm=T)
    }
    if(!'UTMx' %in% colnames(ball)){
      ha    <- areaj <- hullj <- NA
      ux <- uy <- NA
    }

    k <- 0
    
    plotArea <- treePlotArea[treePlotArea[,'plot'] == jplot,-1]
    plotArea <- unlist( plotArea[match(years,treeYr)] )
    names(plotArea) <- paste(jplot,years,sep='-')

    # below and above minimum diameter
    ws   <- which(diamMax < minDiam)
    wt   <- which(diamMax >= minDiam)
    
    baj  <- as.matrix( by(baPerHa(diamMax[wt],max(plotArea,na.rm=T)),sp[wt],sum,na.rm=T) ) 
    sap  <- as.matrix(table(sp[ws]) )
    tre  <- as.matrix(table(sp[wt]) )
    
    slope  <- runif(1,0,.01)
    aspect <- runif(1,0,360)
    
    ele <- alt[j] 
    
    wpd <- which( plotData[,'plotname'] == jplot )
    lonlat <- unlist(  pData$data[wpd,c('lon','lat')] )
    

    
    if('elev' %in% colnames(ball)){
      ele  <- round(mean(ball[,'elev'],na.rm=T),0)
      wa   <- which(is.finite(ball[,'UTMx']) & is.finite(ball[,'UTMy']) & is.finite(ball[,'elev']))
      if(length(wa) > 0){
        tmp    <- localSlope(ball[wa,'UTMx'],ball[wa,'UTMy'],ball[wa,'elev'])
        slope  <- round(tmp$grade,3)
        aspect <- round(tmp$aspect,1)
      }
    }
    
    
    censusYears <- append(censusYears,list(dyr))      # census years
    meanInc     <- mean(diff(dyr))
    tmpBreak    <- getBreaks(dxSeq,minSap,maxTree,meanInc) 
    breakj      <- tmpBreak$breaks
    dxj         <- tmpBreak$dx
  
    sapl    <- appendMatrix(sapl,t(sap))
    tree    <- appendMatrix(tree,t(tre))
    baTab   <- appendMatrix(baTab,t(baj))
    ha      <- max(plotArea,na.rm=T)
    ba_ha   <- round(sum(baj,na.rm=T),2)
    stems   <- sum(sap,na.rm=T)+sum(tre,na.rm=T)
    plotTab <- rbind(plotTab,c(j,k,ux,uy,lonlat,ele,slope,aspect,ha,
                               ba_ha,stems) )
                     
    names(yrStats) <- c('first','last','no.')
    pj <- c( plotData[wpd,c('plotname','state','contact','lon','lat')] ) 
    pj <- as.data.frame(c(pj,ha = ha,ba_ha = ba_ha,stems = stems,yrStats))
    plotSummary <- rbind(plotSummary, pj )
                                        
                                        
    dxMat    <- appendMatrix(dxMat,dxj)
    breakMat <- appendMatrix(breakMat,breakj)
    rownames(dxMat)[nrow(dxMat)] <- rownames(breakMat)[nrow(dxMat)] <- jplot
    deathyr  <- rbind(deathyr,cbind(cens,died))
    
    rname <- tplots[j]
    rownames(tree)[nrow(tree)] <- rownames(baTab)[nrow(baTab)] <- 
                                  rownames(plotTab)[nrow(plotTab)] <-  rname
    
    if(ncol(diamMat) > 1){
      
      diamSamp <- diamMat[wt,colnames(diamMat) %in% dyr]   #only sample years
      
      tmp <- sizeBySpecDistribution( dmat=diamSamp,sss=sp[wt],breaks=breakj,
                                     yr=dyr,die=died[wt])
      atmp <- tmp$all
      stmp <- tmp$surv
      groMu <- tmp$groMu
      groVr <- tmp$groVr
      groNo <- tmp$groNo
      rownames(atmp) <- rownames(stmp) <- 
        rownames(groMu) <- rownames(groVr) <- rownames(groNo) <- 
        paste(jplot,rownames(atmp),sep='-')
      
    #  gtmp[nrow(gtmp),]  <- gtmp[nrow(gtmp)-1,]   #last yr gets previous yr
      groMu[nrow(atmp),] <- groMu[nrow(atmp)-1,]   #last yr gets previous yr
      groVr[nrow(atmp),] <- groVr[nrow(atmp)-1,]   #last yr gets previous yr
      groNo[nrow(atmp),] <- groNo[nrow(atmp)-1,]   #last yr gets previous yr
       
      sizeBySpec <- appendMatrix(sizeBySpec,atmp,fill=0)
      survBySpec <- appendMatrix(survBySpec,stmp,fill=0)
  #    growBySpec <- appendMatrix(growBySpec,gtmp,fill=0)
      growMuSpec <- appendMatrix(growMuSpec,groMu,fill=0)
      growVrSpec <- appendMatrix(growVrSpec,groVr,fill=0)
      growNoSpec <- appendMatrix(growNoSpec,groNo,fill=0)
      
      
      jy <- cbind( rep(match(jplot,plotnames),length(years)),
                   rep(k,length(years)),rep(jp,length(years)),years )
      jy <- jy[years %in% dyr,]
      sizeSpecData <- rbind(sizeSpecData,jy)
      utj <- c(ux,uy)
      
      utm <- rbind(utm,utj )
      pNames <- c(pNames,jplot)
      
      sampAreaYr <- c(sampAreaYr,plotArea[years %in% dyr])
      
    }
    
  }###################################end plot loop

  colnames(plotTab) <- c('j','k','UTMx','UTMy','lon','lat',
                         'elev','slope','aspect','ha','ba/ha','stems')
  sapl[is.na(sapl)] <- 0
  tree[is.na(tree)] <- 0
  baTab[is.na(baTab)] <- 0
  
  phys <- pData$data[match(rownames(plotTab),pData$data[,'plotname']),'phys']
  plotTab <- cbind(plotTab,phys)

  colnames(sizeSpecData) <- c('plot','subplot','j','year')
  
  recruit[is.na(recruit)] <- 0
  recruit <- recruit[rownames(recruit) %in% rownames(sizeBySpec),]
  
  sizeBySpec[is.na(sizeBySpec)] <- 0
  
  o <- orderYnames(colnames(sizeBySpec))$o
  
  sizeBySpec <- sizeBySpec[,o]
  survBySpec <- survBySpec[,o]
#  growBySpec <- growBySpec[,o]
  growMuSpec <- growMuSpec[,o]
  growVrSpec <- growVrSpec[,o]
  growNoSpec <- growNoSpec[,o]

  list(plotData = as.matrix(plotTab), sapStems = sapl, treeStems = tree, recruit = recruit,
       basalArea = baTab, plotnames = pNames, diamAll = diamAll,
       dincAll = dincAll,deathyr = deathyr,
       sizeBySpec = sizeBySpec, survBySpec = survBySpec,
       growMuSpec = growMuSpec, growVrSpec = growVrSpec, growNoSpec = growNoSpec,
       sizeSpecData = sizeSpecData, utm = utm,
       sampAreaYr = sampAreaYr, censusYears = censusYears, breaks = breakMat, 
       dxMat = dxMat,years = years, plotSummary = plotSummary)
}


mortFunction <- function(b){
  
  theta <- exp(xmat[ij,]%*%b)/(1 + exp(xmat[ij,]%*%b))
  pr <- theta*0
  pr <- i[ij]*yy[ij]*log(1 - theta) +
        (1 - i[ij])*( (yy[ij] - 1)*log(1 - theta) + log(theta) )
  
  -sum(pr)
}



orderYnames <- function(yn){
  
  cs <- matrix( unlist( strsplit(yn,'-') ),length(yn),2,byrow=T)
  ot <- which(cs[,1] == 'other')
  cs[ot,1] <- paste('zzz',cs[ot,1],sep='')
  ss <- as.numeric(cs[,2])
  bs <- cs[,1]
  ac <- sort(unique(cs[,1]))
  ac[ac == 'zzzother'] <- 'other'
  
  list(o = order(bs,ss), specs = ac)
}

diam2biomass <- function(allomFile=NULL){
  
  if(is.null(allomFile))allomFile <- '../allocationmodel/datafiles/allometricCoeffs.txt'
  
  bb  <- breakMat
  bb[bb > 200] <- 200
  bb[,kdex == nbreak] <- bb[,kdex == (nbreak-1)] + 1
  
#  mids <- (bb[,kdex != 1] + bb[,kdex != nbreak])/2
#  bb[,kdex != nbreak] <- mids
  
  specmat <- matrix(specs[sdex],ntnp,nsnb,byrow=T)
  
  tmp <- allomConvert(bb,specmat,allomFile=allomFile,
                      int='bmassInt',slope='bmassSlope',
                      defaultSpec='other')
  tmp
}

bmassPerPlot <- function(gmu,gsd=NULL,BYPLOT=F){  #tonnes per ha
  
  # if BYPLOT, then by plot
  # if!BYPLOT, then by plot-yr
  
  bvar <- numeric(0)
  
  tmp   <- gmu*biomass*dxMat/1000
  bvar  <- tmp*biomass/1000/sampArea
  
  bmu  <- byFunctionRcpp( as.vector(tmp), as.vector(wideJTmat),as.vector(wideSmat),
                          matSpec*0,matSpec*0,MEAN=F)
  bvar <- byFunctionRcpp( as.vector(bvar), as.vector(wideJTmat),as.vector(wideSmat),
                          matSpec*0,matSpec*0,MEAN=F)
  
  if(BYPLOT){
    ji <- rep(jtdex[,'j'],nspec)
    si <- as.vector( matrix( c(1:nspec), ntnp,nspec,byrow=T) )
    
    bmu <- byFunctionRcpp( as.vector(bmu), ji  ,si ,
                           matrix(0,nplot,nspec),matrix(0,nplot,nspec),MEAN=T)   
    bvar <- byFunctionRcpp( as.vector(bvar), ji,si,
                            matrix(0,nplot,nspec),matrix(0,nplot,nspec),MEAN=T)
  }

  list(bmu = bmu, bvar = bvar)
}

getClimateData <- function(x = 'prec',months=NULL,whichplot,whichsub=NULL,years,
                           climatePath='climateData/',outFile=NULL,ANNUAL=T){
  
  monthnames <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  
  jfile   <- paste(x,'_',whichplot,sep='')
  jfolder <- paste(climatePath,'monthly_',x,'/',sep='')
  lf      <- list.files(jfolder,pattern=jfile)
  
  if(!is.null(whichsub)){
    kfile <- paste(jfile,whichsub,sep='-')   #is there a separate file for subplot?
    lk     <- list.files(jfolder,pattern=kfile)
    if(length(lk) > 0)lf <- lk
  }
  
  backup <- grep('~',lf)
  if(length(backup) > 0)lf <- lf[-backup]
  
  if(length(lf) == 0){
    warning( paste('no local climate:',jfile,'for variable',x,sep=' ') )
    return()
  }
  
  ff <- paste(jfolder,lf,sep='')
  if(!is.null(outFile))write.table(ff,file=outFile,append=F)
  
  print(ff)
  
  #  zz  <- file(ff,'r')
  tmp <- read.table(ff,header=T)
  
  #  close(zz)
  
  #  tmp    <- read.table(ff,header=F)
  colnames(tmp) <- c('yr',monthnames)
  
  mm <- c(2:13)
  if(!is.null(months)){
    mm     <- numeric(0)
    for(m in 1:12)if(length(grep(monthnames[m],months)) > 0)mm <- c(mm,m)
    mm <- mm + 1
  }
  meank <- rowMeans(tmp[,mm],na.rm=T)        #mean temp or pdsi
  if(x == 'prec')meank <- meank*length(mm)   #total rainfall
  
  ytmp   <- tmp[,'yr']
  ymatch <- match(years,ytmp)
  wm     <- which(is.na(ymatch))
  if(length(wm) > 0){
    ymatch[wm] <- length(ytmp)
  }
  
  yy     <- cbind(years,meank[ymatch])
  colnames(yy) <- c('yr',paste(x,months,sep='_'))
  
  if(!ANNUAL){
    yy <- tmp[ymatch,mm]
    rownames(yy) <- years
  }
  
  list(climvar = as.matrix(yy), monthVec = mm - 1)
}

biomassPlots <- function(gmu,gsd=NULL,trimBiomass=Inf,mapscale=1.8){
  
  require(ade4)
  
  graphics.off()
  tmp <- bmassPerPlot(gmu,gsd) 
  bmass <- tmp$bmu
  bvar  <- tmp$bvar
  
  bmassTot <- rowSums(bmass,na.rm=T)
  bmassVar <- rowSums(bvar,na.rm=T)
  bmassTot[bmassTot > trimBiomass] <- trimBiomass
  
  if(nHoldOut > 0){
    whold <- which(jtdex[,'j'] %in% holdOutPlots)
    hist(bmassTot[-whold],nclass=100)
  }
  dev.copy2pdf(file=outFile(outfolder,'plotBmassHist.pdf'))
  
  graphics.off()
  
  bmassByPlot <- by(bmassTot,jtdex[,'j'],mean)
  bVarByPlot  <- by(bmassVar,jtdex[,'j'],mean)
  
  tmp <- bmassTot
  tmp[dtVec == 0] <- NA
  corr <- cor(tmp[-1],tmp[-ntnp],use='complete.obs')
  
  tmp1 <- (bmassTot[-1] - bmassTot[-ntnp])/dtVec[-1]
  jt   <- jtdex[dtVec > 0,'j']
  tmp  <- tmp1[is.finite(tmp1)]
  prodByPlot <- by(tmp,jt,mean)
  
  tmp2     <- (bmassVar[-1] + bmassVar[-ntnp] - 
                 2*corr*sqrt(bmassVar[-1]*bmassVar[-ntnp]))/dtVec[-1]/dtVec[-1]
  tmp2    <- tmp2[is.finite(tmp2)]
  prodVar <- by(tmp2,jt,mean)
  
  bioOut <- cbind(bmassByPlot,bVarByPlot,prodByPlot,prodVar)
  colnames(bioOut) <-c('T/ha','bVar','T/ha/yr','pVar')
  
  ll <- cbind(plotLon,plotLat)
  colnames(ll) <- c('lon','lat')
  bioOut <- cbind(ll,bioOut)
  rownames(bioOut) <- NULL
  
  out <- data.frame(plotnames,ecoReg,ecoReg4,bioOut)
  
  save(out,file=outFile(savefolder,'biomassByPlot.Rdata'))
  
  np <- nplot - nHoldOut
  
  scaleSym <- .0003
  regMap(topo$x,topo$y,topo$z,IMAGE=T,mapscale=mapscale,lineCol='grey')
  regMap(plotLon[-holdOutPlots],plotLat[-holdOutPlots],bmassByPlot[-holdOutPlots],
         lineCol='grey',fg=rep('black',np),scaleSym=scaleSym,IMAGE=F,ADD=T)
  
  xr <- maplon[1] + diff(maplon)*c(.77,1)
  yr <- maplat[1] + diff(maplat)*c(0,.15)
 # xx <- maplon[1] + diff(maplon)*c(.92,.97)
 # yy <- maplat[1] + diff(maplat)*c(.02,.15)
  
  rect(xr[1],yr[1],xr[2],yr[2],col='white',border='white')
  
  bvals <- c(50,250)
  xpt   <- rep(xr[1] + diff(xr)*.1,2)
  ypt   <- yr[1] + c(.25,.75)*diff(yr)
  
  symbols(xpt,ypt,circles=scaleSym*bvals,bg='black',inches=F,add=T)
  text(xpt,ypt,bvals,pos=4)
  ylabel <- expression( "T ha"^{-1} )
  text( xr[1] + diff(xr)*.6,mean(ypt),ylabel)
  title('Biomass')
  
  
  dev.copy2pdf(file=outFile(outfolder,'tonsPerHa.pdf'))
  

  ncol   <- 100
  colF   <- colorRampPalette(c('white','wheat','orange','brown','black'))

  xx <- bmassByPlot[-holdOutPlots]
  zlevs <- c(0,seq(100,250,by=10),10000)
  
  regMap(topo$x,topo$y,topo$z,IMAGE=F,mapscale=mapscale,lineCol='grey')
  
  values2contour(plotLon[-holdOutPlots],plotLat[-holdOutPlots],
                 xx,nx=50,ny=50,col=colF(length(zlevs)),
                 zlevs=zlevs,lwd=2,add=T,drawlabels=F,fill=T)
  
  mapMask(plotLon,plotLat,dx=.1,dy=.1,whiteOutDist=.25,col='white')
  
  map('state',add=T,col='white',lwd=6,interior=F)
  map('state',add=T,col='grey',lwd=3,interior=T)
  
  add.scatter(hist2Add(xx,xlim=c(0,400)),posi='bottomright',ratio=.15,inset=c(.02,.07),bg.col='grey')
  
  dev.copy2pdf(file=outFile(outfolder,'tonsPerHaContour.pdf'))
  
  #std error
  
  regMap(topo$x,topo$y,topo$z,IMAGE=T,mapscale=mapscale,lineCol='grey')
  
  regMap(plotLon[-holdOutPlots],plotLat[-holdOutPlots],
         sqrt(bVarByPlot[-holdOutPlots]),lineCol='grey',fg=rep('black',np),
         scaleSym=scaleSym,IMAGE=F,ADD=T)
  
  rect(xr[1],yr[1],xr[2],yr[2],col='white',border='white')
  
  bvals <- c(50,250)
  xpt   <- rep(xr[1] + diff(xr)*.1,2)
  ypt   <- yr[1] + c(.25,.75)*diff(yr)
  
  symbols(xpt,ypt,circles=scaleSym*bvals,bg='black',inches=F,add=T)
  text(xpt,ypt,bvals,pos=4)
  ylabel <- expression( "T ha"^{-1} )
  text( xr[1] + diff(xr)*.6,mean(ypt),ylabel)
  title('Biomass standard error')
 # text( xr[1] + diff(xr)*.6,ypt[1],"SD")
 # values2contour(plotLon[-holdOutPlots],plotLat[-holdOutPlots],
 #                sqrt(bVarByPlot[-holdOutPlots]),nx=50,ny=50,col='white',
 #                zlevs=c(0,.1,.5,2,4,8),lwd=4,add=T,drawlabels=F)
 # values2contour(plotLon[-holdOutPlots],plotLat[-holdOutPlots],
 #                sqrt(bVarByPlot[-holdOutPlots]),nx=50,ny=50,col='blue',
 #                zlevs=c(0,.1,.5,2,4,8),lwd=2,add=T,drawlabels=T)
  
  dev.copy2pdf(file=outFile(outfolder,'tonsPerHaSE.pdf'))
  
  
  xx <- sqrt(bVarByPlot[-holdOutPlots])/bmassByPlot[-holdOutPlots]
  xx[xx > 1] <- 1
  zlevs <- c(0,seq(.2,.8,length=100))
  
  regMap(topo$x,topo$y,topo$z,IMAGE=F,mapscale=mapscale,lineCol='grey')
  
  values2contour(plotLon[-holdOutPlots],plotLat[-holdOutPlots],
                 xx,nx=50,ny=50,col=colF(length(zlevs)),
                 zlevs=zlevs,lwd=2,add=T,drawlabels=F,fill=T)
  
  mapMask(plotLon,plotLat,dx=.1,dy=.1,whiteOutDist=.25,col='white')
  
  map('state',add=T,col='white',lwd=6,interior=F)
  map('state',add=T,col='grey',lwd=3,interior=T)
  
  add.scatter(hist2Add(xx,xlim=c(0,1)),
              posi='bottomright',ratio=.15,inset=c(.02,.07),bg.col='grey')

  dev.copy2pdf(file=outFile(outfolder,'tonsPerHaCvContour.pdf'))
  
  
  
  xx <- sqrt(bVarByPlot[-holdOutPlots])
  zlevs <- c(0,seq(10,100,by=10),300)
  
  regMap(topo$x,topo$y,topo$z,IMAGE=F,mapscale=mapscale,lineCol='grey')
  
  values2contour(plotLon[-holdOutPlots],plotLat[-holdOutPlots],
                 xx,nx=50,ny=50,col=colF(length(zlevs)),
                 zlevs=zlevs,lwd=2,add=T,drawlabels=F,fill=T)
  
  mapMask(plotLon,plotLat,dx=.1,dy=.1,whiteOutDist=.25,col='white')
  
  map('state',add=T,col='white',lwd=6,interior=F)
  map('state',add=T,col='grey',lwd=3,interior=T)
  
  add.scatter(hist2Add(xx,xlim=c(0,100)),
              posi='bottomright',ratio=.15,inset=c(.02,.07),bg.col='grey')
  
  dev.copy2pdf(file=outFile(outfolder,'tonsPerHaSdContour.pdf'))
  
  
  
  
  
  scaleSym=.007
  regMap(topo$x,topo$y,topo$z,IMAGE=T,mapscale=mapscale,lineCol='grey')
  wp <- which(prodByPlot > 0)
  wp <- wp[!wp %in% holdOutPlots]
  regMap(plotLon[wp],plotLat[wp],lineCol='grey',
         prodByPlot[wp],fg='blue',bg='black',scaleSym=scaleSym,IMAGE=F,ADD=T)
  wn <- which(prodByPlot < 0)
  wn <- wn[!wn %in% holdOutPlots]
  regMap(plotLon[wn],plotLat[wn],lineCol='grey',
         -prodByPlot[wn],fg='red',bg='red',scaleSym=scaleSym,IMAGE=F,ADD=T)
  
  
  rect(xr[1],yr[1],xr[2],yr[2],col='white',border='white')
  
  bvals <- c(5,15)
  xpt   <- rep(xr[1] + diff(xr)*.1,2)
  ypt   <- yr[1] + c(.25,.75)*diff(yr)
  
  symbols(xpt,ypt,circles=scaleSym*bvals,bg='black',inches=F,add=T)
  text(xpt,ypt,bvals,pos=4)
  ylabel <- expression( "T ha"^{-1}*yr^{-1} )
  text( xr[1] + diff(xr)*.6,mean(ypt),ylabel)
  title('Productivity')
  
  dev.copy2pdf(file=outFile(outfolder,'tonsPerHaPerYr.pdf'))
  
  scaleSym <- .004
  regMap(topo$x,topo$y,topo$z,IMAGE=T,mapscale=mapscale,lineCol='grey')
  regMap(plotLon[wp],plotLat[wp],lineCol='grey',
         sqrt(prodVar[wp]),fg=rep('blue',np),bg=rep('black',np),
         scaleSym=scaleSym,IMAGE=F,ADD=T)
  regMap(plotLon[wn],plotLat[wn],lineCol='grey',
         sqrt(prodVar[wn]),fg=rep('red',np),bg=rep('red',np),
         scaleSym=scaleSym,IMAGE=F,ADD=T)
 # values2contour(plotLon[-holdOutPlots],plotLat[-holdOutPlots],
 #                sqrt(prodVar[-holdOutPlots]),nx=50,ny=50,col='white',
 #                zlevs=c(2,4,8,16),lwd=4,add=T,drawlabels=F)
 # values2contour(plotLon[-holdOutPlots],plotLat[-holdOutPlots],
 #                sqrt(prodVar[-holdOutPlots]),nx=50,ny=50,col='blue',
 #                zlevs=c(2,4,8,16),lwd=2,add=T,drawlabels=T)
  
  rect(xr[1],yr[1],xr[2],yr[2],col='white',border='white')
  bvals <- c(5,15)
  xpt   <- rep(xr[1] + diff(xr)*.1,2)
  ypt   <- yr[1] + c(.25,.75)*diff(yr)
  
  symbols(xpt,ypt,circles=scaleSym*bvals,bg='black',inches=F,add=T)
  text(xpt,ypt,bvals,pos=4)
  ylabel <- expression( "T ha"^{-1}*yr^{-1} )
  text( xr[1] + diff(xr)*.6,mean(ypt),ylabel)
#  text( xr[1] + diff(xr)*.6,ypt[1],"SD")
  title('Productivity standard error')
  dev.copy2pdf(file=outFile(outfolder,'tonsPerHaPerYrSE.pdf'))
  
  
  
  xx <- prodByPlot[wp]
  xx[xx > 11] <- 11
  zlevs <- seq(0,11,length=100)
  
  regMap(topo$x,topo$y,topo$z,IMAGE=F,mapscale=mapscale,lineCol='grey')
  
  values2contour(plotLon[wp],plotLat[wp],
                 xx,nx=50,ny=50,col=colF(length(zlevs)),
                 zlevs=zlevs,lwd=2,add=T,drawlabels=F,fill=T)
  
  mapMask(plotLon,plotLat,dx=.1,dy=.1,whiteOutDist=.25,col='white')
  
  map('state',add=T,col='white',lwd=6,interior=F)
  map('state',add=T,col='grey',lwd=3,interior=T)
  
  add.scatter(hist2Add(xx,xlim=c(0,10)),
              posi='bottomright',ratio=.15,inset=c(.02,.07),bg.col='grey')
  
  dev.copy2pdf(file=outFile(outfolder,'tonsPerHaYrContour.pdf'))
  
  
  xx <- sqrt(prodVar[wp])
  xx[xx > 11] <- 11
  zlevs <- seq(0,11,length=100)
  
  regMap(topo$x,topo$y,topo$z,IMAGE=F,mapscale=mapscale,lineCol='grey')
  
  values2contour(plotLon[wp],plotLat[wp],
                 xx,nx=50,ny=50,col=colF(length(zlevs)),
                 zlevs=zlevs,lwd=2,add=T,drawlabels=F,fill=T)
  
  mapMask(plotLon,plotLat,dx=.1,dy=.1,whiteOutDist=.25,col='white')
  
  map('state',add=T,col='white',lwd=6,interior=F)
  map('state',add=T,col='grey',lwd=3,interior=T)
  
  add.scatter(hist2Add(xx,xlim=c(0,10)),
              posi='bottomright',ratio=.15,inset=c(.02,.07),bg.col='grey')
  
  dev.copy2pdf(file=outFile(outfolder,'tonsPerHaYrSeContour.pdf'))
  
  
  
   list(biomass = bmassByPlot, prod = prodByPlot, 
       biomassVar = bVarByPlot, prodVar = prodVar)
}


gatherOther <- function(xx,keepCol){  #shift rare spp into other
  
  scol <- which( !c(1:ncol(xx)) %in% keepCol )
  ocol <- grep('other',colnames(xx))
  
  noto <- keepCol[!keepCol %in% ocol]
  
  xnew <- xx[,noto]
  
  oindex <- matrix( 1:nbreak,nrow(xx),length(scol),byrow=T)
  rindex <- matrix( 1:nrow(xx),nrow(xx),length(scol))
  mat    <- matrix(0,nrow(xx),nbreak)
  
  tmp  <- byFunction(as.vector(xx[,scol]),i=as.vector(rindex),j=as.vector(oindex),mat)
  omat <- xx[,ocol] + tmp
  
  xnew <- cbind(xnew,omat)
  sp   <- unique(matrix( unlist( strsplit(colnames(xnew),'-') ), ncol=2,byrow=T)[,1])
  list(x = xnew, specs = sp)
}



getDemX <- function(vnames,xplot,xyear,CENTER=F,STANDARD=F,jdex,hvar='wetness'){
  
  # matrix demData includes column 'phys' for hvar = 'factor'

  pall <- unique(xplot)
  np   <- length(pall)
                                                                                
  hfile <- paste('../allocationmodel/datafiles/hydro/',hvar,'/',sep='')
  p <- length(vnames)

  x <- matrix(1,length(xplot),p)
  colnames(x) <- vnames
  
  monTemp <- monPrec <- NULL
                                                                                
  wtt <- grep('temp',vnames)
  if(length(wtt) > 0){
    wmo <- grep('-',vnames[wtt[1]])
    if(length(wmo) > 0){
       tjk <- unlist( strsplit(vnames[wtt[1]],'-') )
       yjk <- tjk[1]
       monTemp <- tjk[2]
    }
  }
  wtt <- grep('prec',vnames)
  if(length(wtt) > 0){
    wmo <- grep('-',vnames[wtt[1]])
    if(length(wmo) > 0){
      tjk <- unlist( strsplit(vnames[wtt[1]],'-') )
      yjk <- tjk[1]
      monPrec <- tjk[2]
    }
  }
  
#  centSize <- (breaks - mean(breaks))/sd(breaks)
#  centBA   <- (baLarge - mean(baLarge))/sd(baLarge)
                                                                                
  cvars <- unique( c(grep('temp',vnames),grep('prec',vnames),
                     grep('therm',vnames),grep('pdsi',vnames),
                     grep('deficit',vnames)) )
  GINDEX <- F
  if('therm' %in% vnames | 'deficit' %in% vnames)GINDEX <- T
                                                                                
  tmp <- clim4FIA(lat=demLonLat[pall,'lat'],lon=demLonLat[pall,'lon'],
                plVec=jdex[,'j'],yrVec=years[ jdex[,'t'] ],
                monTemp=monTemp,monPrec=monPrec,GINDEX=GINDEX,
                start=firstClimYr,end=lastClimYr )
  
  climData <- tmp$regClim
  rownames(climData) <- paste(pall[jdex[,'j']],years[jdex[,'t']],sep='-')
                                                                                
  
  for(k in 1:p){

    for(j in 1:np){

      wj   <- which(xplot == pall[j])
      latj <- demLonLat[j,'lat']

      if(k %in% cvars){
        
        months <- NULL
        yjk    <- vnames[k]
        
        wmo <- grep('-',vnames[k])
        if(length(wmo) > 0){
          tjk <- unlist( strsplit(vnames[k],'-') )
          yjk <- tjk[1]
          months <- tjk[2]
        }
        
        if(yjk == 'temp2')yjk <- 'temp'
        
        if(yjk != 'therm'){
          climk <- getClimateData(x=yjk,months=months,whichplot=pall[j],
                              years=years[jdex[wj,'t']],
                              climatePath='../allocationmodel/climateData/')$climvar
          
          if(!is.null(climk))x[wj,k] <- climk[,2]
          if(is.null(climk)){   #no local climate file
            if(yjk == 'temp')x[wj,k] <- climData[wj,1]
            if(yjk == 'prec')x[wj,k] <- climData[wj,2]
          }
          wna <- which(is.na(climk[,2]))
          if(length(wna) > 0){
            if(yjk == 'temp')x[wj[wna],k] <- climData[wj[wna],1]
            if(yjk == 'prec')x[wj[wna],k] <- climData[wj[wna],2]
          }
          
        }
        if(yjk == 'therm' | yjk == 'deficit'){
          temp <- getClimateData(x='temp',months=months,whichplot=pall[j],
                                years=years[jdex[wj,'t']],
                                climatePath='../allocationmodel/climateData/',
                                 ANNUAL=F)$climvar
 
          prec <- getClimateData(x='prec',months=months,whichplot=pall[j],
                                 years=years[jdex[wj,'t']],
                                 climatePath='../allocationmodel/climateData/',
                                 ANNUAL=F)$climvar
          if(is.null(temp))x[wj,k] <- climData[wj,3]
          if(!is.null(temp)){
            tmp <-  monthlyPHr(yi=rep(1,12),mi=c(1:12),tempMatrix=temp,precMatrix=prec,
                               lat=rep(latj,nrow(temp)))
            
            if(yjk == 'therm')  x[wj,k] <- tmp$degHrPosYr
            if(yjk == 'deficit')x[wj,k] <-tmp$degHrNegYr
          }
          
        }
          
      }
      
      if(vnames[k] == 'xeric'){
        wd <- which(rownames(demData) == pall[j])
        hh <- demData[wd,'phys']
        x[wj,k] <- 0
        if(hh == 1)x[wj,k] <- 1
      }
      if(vnames[k] == 'mesic'){
        wd <- which(rownames(demData) == pall[j])
        hh <- demData[wd,'phys']
        x[wj,k] <- 0
        if(hh == 3)x[wj,k] <- 1
      }

      if(vnames[k] %in% c('phys','hydro')){

        utm    <- matrix(demUTM[j,],length(wj),ncol=2,byrow=T)
        region <- matrix( unlist(strsplit(pall[j],'_')),ncol=2,byrow=2)[,1]
        jfile  <- paste(hfile,region,'.asc',sep='')

        x[wj,k] <- nearestHydro(utm[,1],utm[,2],hydroFile=jfile) #neg. values for wetness
      }
    }
  }
  
 # linearXpred <- c('temp','prec','hydro')              #inputs that have priors, linear
 # downScale   <- c('phys','elev','slope','aspect')     #downscale prior variables
                                                                                
  
 if(hvar == 'tcinew')x[,'hydro'] <- log(x[,'hydro'])
  
  tmp <- inClimate(lon=demLonLat[,'lon'],lat=demLonLat[,'lat'],
                   start=firstClimYr,end=lastClimYr,PADLAST=T)
  rTemp  <- tmp$temp
  rPrec  <- tmp$prec
  
  yri <- trunc( as.numeric(colnames(rTemp)), 0)
  yi <- yri - yri[1] + 1
  mi <- rep(1:12,length(yi)/12)
  tmp <- monthlyPHr(yi=yi,mi=mi,tempMatrix=rTemp,precMatrix=rPrec,
                       lat=demLonLat[,'lat'])
  rTherm <- tmp$degreeHrPos
  rDef <- tmp$degreeHrNeg
  
  lastClimYr <- max( floor( as.numeric(colnames(rTemp)) ) )
  
  wplot <- match(xplot,rownames(demData))
  
  regClim <- matrix(NA,length(wplot),4)
  colnames(regClim) <- c('temp','prec','therm','deficit')
  
  for(j in 1:length(xyear)){
    
    yj <- grep(xyear[j],colnames(rTemp))
    if(length(yj) == 0)yj <- grep(lastClimYr,colnames(rTemp))
    regClim[j,'temp'] <- mean( rTemp[wplot[j],yj] ,na.rm=T)
    regClim[j,'prec'] <- sum( rPrec[wplot[j],yj] ,na.rm=T)
    regClim[j,'therm'] <- sum( rTherm[wplot[j],yj] ,na.rm=T)
    regClim[j,'deficit'] <- sum( rDef[wplot[j],yj] ,na.rm=T)
  }
  
  
  xpriorX <-  xpriorY <- bhat <- bsig <- numeric(0)
  
  if(!is.null(linearXpred)){
    
    downScale <- downScale[downScale != 'aspect']
    if(!'u1' %in% downScale)downScale <- c(downScale,c('u1','u2'))

    pX <- length(downScale)
    xpriorX <- matrix(1,nrow(x),pX+1)
    colnames(xpriorX) <- c('intercept',downScale)
    
  
    xpriorY <- x[,linearXpred]  #contains local variables
 #   xpriorY[,c('temp','prec')] <- xpriorY[,c('temp','prec')] - regClim  #local anomalies
    
    xpriorX[,'phys']  <- demData[wplot,'phys']
    if('elev' %in% downScale)xpriorX[,'elev']  <- demData[wplot,'elev']
    slope  <- demData[wplot,'slope']
    
    asp   <- degrees2radians( demData[wplot,'aspect'] )
    u1    <- slope*sin(asp)
    u2    <- slope*cos(asp)
    xpriorX[,'slope'] <- slope
    xpriorX[,c('u1','u2')] <- cbind(u1,u2)
  
    xpriorX <- cbind(xpriorX,regClim)
    
    XX <- crossprod(xpriorX)
    XY <- crossprod(xpriorX,xpriorY)
    bhat <- solve(XX)%*%XY

    bsig   <- crossprod(xpriorY - xpriorX%*%bhat)/nrow(xpriorX)
  }
 
 #   xpriorY[,c('temp','prec')] <- xpriorY[,c('temp','prec')] + regClim
      
  if(CENTER){
       means <- colMeans(x)
       x <- x - matrix(means,nrow(x),p,byrow=T)
  }
  if(STANDARD){
       sds <- apply(x,2,sd)
       x <- x/matrix(sds,nrow(x),p,byrow=T)
  }

  x[,'intercept'] <- 1                 # nt*nplot by p
  
  breakTmp <- demBreaks
  breakTmp[,nbreak] <- breakTmp[,(nbreak-1)] + demDx[,(nbreak-2)]
  

  list( x = x,xpriorX = xpriorX, xpriorY = xpriorY, regClim = regClim, 
        betax = bhat, bsig = bsig)
}


initRhoPars <- function(){

  rtmp <- rho
  rtmp[rtmp == 0] <- .05
  rtmp[rtmp == 1] <- .95

  yy <- logit(rtmp)

  xx <- crossprod(xmatRho)
  bb <- invMat(xx)%*%crossprod(xmatRho,yy)
  ss <- yy - xmatRho%*%bb
  ss <- crossprod(ss)/nplot

  list(beta = bb, sigma = ss)
}


    
plotSizeDist <- function(count,values,weight,jt,s2plot,p2plot=1,xmax=NULL,ymax=NULL,
                         ADD=F,LINES=T,POINTS=F,colseq=NULL,smoothBW = 1,lty='l'){
  
  #count  - ntnp matrix in counts in bin
  #values - values at each break
  #wt     - the wt for each count (e.g., 1/area/binwidth)
  #nt     - number of years to plot
  
    
    cc <- count
    xx <- cc*weight
    
    if(is.character(s2plot))wr <- which(specs == s2plot)
    
    if(is.numeric(s2plot))wr <- s2plot
    
    wc <- grep(specs[wr],colnames(xx))
      
    wp <- which(jt[,'j'] %in% p2plot)

    xi <- xx[wp,wc]
    ji <- jt[wp,'j']
    ti <- years[jt[wp,'t']]
    nt <- length(years)
                
    fraction <- xi/matrix(rowSums(xi,na.rm=T),length(wp),length(wc))
    
    breakVals <- values[wp,wc]
    
    if(length(colseq) == 1)colseq <- rep(colseq,nt)
    if(is.null(colseq))    colseq <- mapColors(nt)
    
    if(is.null(ymax))ymax <- max(xi,na.rm=T)
    if(is.null(xmax))xmax <- max(breakVals)
  

    if(!ADD) plot(breaks,breaks*0-10,xlim=c(0,xmax),ylim=c(0,ymax),cex=.1)

    for(t in 1:nt){
      
        ww <- which(ti == years[t])
        if(length(ww) == 0)next
        
        xt <- breakVals[ww,]
        xq <- which(xt < xmax)
      #  xq <- 1:max(xq[,2])
        xt <- as.vector(xt[xq])
        
        
        wt <- as.vector(xi[ww,][xq])
        if(max(wt,na.rm=T) == 0)next
        
        wt <- wt/sum(wt,na.rm=T)

        fr <- as.vector(fraction[ww,][xq])
        tmp <- density(xt,weights=wt,bw=smoothBW)
        if(LINES){
          if(lty == 'l')lines(tmp$x,tmp$y,col=colseq[t],lwd=2)
          if(lty == 's')lines(xt,fr/sum(fr,na.rm=T),col=colseq[t],type='s',lwd=2)
        }
        if(POINTS) points(jitter(xt,factor=.2),fr/sum(fr,na.rm=T),cex=.8,col=colseq[t])
    }
    
    invisible(xi)
    
}



  

propGamma <- function(gam=gam){

  tiny <- 1e-4

  gam[is.na(gam)] <- 0

  tseq <- c(nt:1)
  tseq <- sample(tseq)

  for(t in tseq){

    ptiny <- 1
    if(t < nt)ptiny <- tiny

    wt   <- which(jtdex[,'t'] == t)

    k0t0 <- k1t0 <- k1t2 <- k2t2 <- gam[wt,adex[,2]]*0

    k1t1 <- gam[wt,adex[,2]]
    k0t1 <- gam[wt,adex[,1]]
    k2t1 <- gam[wt,adex[,3]]

    if(t > 1){
      w0   <- wt-1
      k0t0 <- gam[w0,adex[,1]]
      k1t0 <- gam[w0,adex[,2]]
    }

    if(t < nt){
      w2   <- wt+1
      k1t2 <- gam[w2,adex[,2]]
      k2t2 <- gam[w2,adex[,3]]
    } 

    k0t0[,ksdex[,1]] <- fecundity[wt,]  #first stage
    k1t0[,ksdex[,1]] <- 0  

    k2t1[,ksdex[,nbreak]] <- 0
    k2t2[,ksdex[,nbreak]] <- 0

    ht <- k0t0 - k0t1 + k1t0 - k2t1
    lt <- k1t2 - k0t1 + k2t2 - k2t1

    lt[lt < 0]  <- 0
    ht[ht < lt] <- lt[ht < lt] + ptiny

    gt <- gam[wt,]
    gt[gt < lt] <- lt[gt < lt]
    gt[gt > ht] <- ht[gt > ht]

    gam[wt,] <- tnorm(length(gt),lt,ht,gt,.1)
  }

  gam
}
  

  


sizeBySpecDistribution <- function(dmat,sss,breaks,yr,die){   #b is a plot data.frame

    break1 <- breaks
    break1[1] <- 0
    nbreak <- length(breaks)

    dmat <- row2Mat(dmat)

 #   sampleYears <- dateColumn(bb,'diam',yrvec=years)$yr

    speck  <- sort(unique(sss))
    colk   <- as.vector( t(outer(speck,c(1:nbreak),paste,sep='-')) )
    amk    <- matrix(NA,length(yr),nbreak*length(speck))
    colnames(amk) <- colk
    rownames(amk) <- yr
    
    groMu <- groVr <- groNo <- smk <- amk

    for(jj in 1:length(speck)){
      
      wj    <- which(sss == speck[jj])
      if(length(wj) < 1)next
      
      djmat <- row2Mat( dmat[wj,] )
      
      #size distribution
      jcol <- paste(speck[jj],c(1:nbreak),sep='-')
      
      tmp  <- histByYr( xx=djmat,yr,dens=F,breaks=break1,deathyr=die[wj])
      dtab <- tmp$dist
      stab <- tmp$surv
      g1 <- tmp$groMu
      g2 <- tmp$groVr
      g3 <- tmp$groNo
      amk[,jcol]   <- t(dtab[1:nbreak,])
      smk[,jcol]   <- t(stab[1:nbreak,])
      groMu[,jcol] <- t(g1[1:nbreak,])
      groVr[,jcol] <- t(g2[1:nbreak,])
      groNo[,jcol] <- t(g3[1:nbreak,])
      
    }
    
    list(all = amk, surv = smk, groMu = groMu, groVr = groVr, groNo = groNo)
}


makeZ <- function(minStem=0,minPlot=1,numSpec=NULL,ALLYEARS=F,mustHave=NULL){
  
  # ALLYEARS produces matrix with years inclusive
  # !ALLYEARS includes only census years
  # mustHave - character vector of species that must be included
  
 # nspec <- nspecFIA
  
 # treeCodes <- read.table('../allocationmodel/datafiles/treeCodesDuke.txt',header=T)
 # snames    <- as.character(treeCodes[match(fiaspecs,treeCodes[,'fiaCode']),'code'])
 # snames[is.na(snames)] <- fiaspecs[is.na(snames)]
  
  ns   <- length(fiaSpecs)
  nsnb <- ns*nbreak
  
 # fspecs <- c(fiaSpecs,'other')
  
  yTree <- matrix(0,nplotFIA,nsnb)
  rownames(yTree) <- plotnamesFIA
  cn <- as.vector(t(outer(fiaSpecs,c(1:nbreak),paste,sep='-')))
  colnames(yTree) <- cn
  yTree1 <- yTree
 
  yi <- ysi  <- seq(1,nsnb,by=nbreak)
  for(i in 1:(nbreak-1))ysi  <- cbind(ysi,yi+i)
  
  maxSpec <- rep(0,nspec)
  ecoReg  <- data[match(plots,data[,'plt1']),'eco']

  for(k in 1:ns){
    
    wtree <- which( data[,'sp'] == fiaSpecs[k] & is.finite(data[,'dia1']))
    wtb1  <- cut(data[wtree,'dia1'],breaks,label=F )
    tmp   <- byIndex( wtb1*0 + 1, list(bin = wtb1,plot = data[wtree,'plt1']), sum,na.rm=T )
    yTree[colnames(tmp),ysi[k,as.numeric(rownames(tmp))]] <- t(tmp)
    
    mk <- matrix(t(tmp),ncol(tmp),nrow(tmp))

    wtree <- which( data[,'sp'] == fiaSpecs[k] & is.finite(data[,'dia2']))
    wtb1 <- cut(data[wtree,'dia2'],breaks,label=F )
    tmp  <- byIndex( wtb1*0 + 1, list(bin = wtb1,plot = data[wtree,'plt1']), sum,na.rm=T )
    yTree1[colnames(tmp),ysi[k,as.numeric(rownames(tmp))]] <- t(tmp)
    
    mk <- c(tmp,mk)
    maxSpec[k] <-  max(mk,na.rm=T)

  }
  
  yTree[is.na(yTree)]   <- 0
  yTree1[is.na(yTree1)] <- 0
  
  #################################3333
  allspecs <- fiaSpecs
  nspec    <- 0
  specs    <- character(0)
  if(!is.null(mustHave)){
    specs <- fiaSpecs[snames %in% mustHave]
    nspec <- length(specs)
  }
  if(!is.null(numSpec)){
    ws      <- which(fiaSpecs %in% specs)
    ofs      <- fiaSpecs[order(maxSpec,decreasing=T)]
    fs      <- ofs[!ofs %in% specs]
    ns      <- numSpec - length(ws)
    if(ns > 0)specs   <- sort(c(specs, fs[1:ns]))
  }
  
  ###################################3
 # newNames  <- sort(unique( c('other',as.character( treeCodes[match(specs,treeCodes[,'fiaCode']),'code']) ) ))

  
  yn1 <- paste(plotnamesFIA,t1,sep='-')
  yn2 <- paste(plotnamesFIA,t2,sep='-')
  
  if(!ALLYEARS)rnames <- as.vector( t( cbind(yn1,yn2) ) )
  if(ALLYEARS) rnames <- as.vector(t(outer(plotnamesFIA,years,paste,sep='-')))
    
#  ns <- length(newNames)
  
  yt <- matrix(0,nrow(jtdexFIA),nbreak*ns)
  colnames(yt) <- colnames(yTree)
  
  other <- paste('other',c(1:nbreak),sep='-')
  
  for(j in 1:ns){
    
    print(fiaspecs[j])
    
    wc   <- grep(paste(fiaspecs[j],'-',sep=''),cn)
    
    sj <- snames[j]
    cj <- paste(sj,c(1:nbreak),sep='-')
    
    if(length(wc) == 0)next
    
    if(!fiaspecs[j] %in% specs){
      
      yt[match(yn1,rnames),other] <- yt[match(yn1,rnames),other] + yTree[,wc]
      yt[match(yn2,rnames),other] <- yt[match(yn2,rnames),other] + yTree1[,wc]
      yt[is.na(yt)]  <- 0

      next
    }
    
    ymat <- matrix(NA,length(rnames),nbreak)
    
    colnames(ymat) <- cj
    
    ymat[match(yn1,rnames),] <- yTree[,wc]
    ymat[match(yn2,rnames),] <- yTree1[,wc]
    yt[,cj] <- ymat
    yt[is.na(yt)] <- 0
    
  }
    
  rownames(yt) <- rnames
  
  ns <- length(newNames)
  
  sampAreaFIA <- rep(sampAreaFIAtree,nbreak)
  sampAreaFIA[breaks < minTree] <- sampAreaFIAsap
  sampArea <- matrix(rep(sampAreaFIA,ns),2*nplot,nbreak*ns,byrow=T)
  

  list(yTree1 = yTree, yTree2 = yTree1, ksIndex = ysi,
       yt = yt, specs = newNames, sampArea = sampArea, ecoReg = ecoReg)
}

getFIAspecs <- function(specs){
  
  treeCodes <- read.table('../allocationmodel/datafiles/treeCodesDuke.txt',header=T)
  fiaSpecs  <- as.character(treeCodes[match(specs,treeCodes[,'code']),'fiaCode'])
  fiaSpecs[specs == 'other'] <- 'other'
  if(!'other' %in% fiaSpecs)fiaSpecs <- c(fiaSpecs,'other')
  
  fiaSpecs
}


groupFIAregs <- function(ecr,level=3,noM = T){    
  
  #level 3 remove no trailing letters ('Province')
  #level 2 remove 1 letter
  #level 1 remove 2 letters
  
  if(!level %in% c(1:3))stop('level must be 1:3')
  
  ee <- enom <- as.character(ecr)
  ws <- grep(' ',ee)
  if(length(ws) > 0){             # remove leading spaces
    ee[ws] <- substring(ee[ws],2)
    enom   <- ee
  }
  
  wm <- grep('M',enom)
  if(length(wm) > 0){
    notm <- matrix( unlist( strsplit(enom[wm],'M') ), ncol=2,byrow=T)[,2]
    enom[wm] <- notm
  }
  
  if(level == 3){
    if(noM)return(enom)
    if(!noM){
      enol <- substr(enom,1,3)
      enol[wm] <- paste('M',enol[wm],sep='')
      return(enol)
    }
  }
  
  enol <- substr(enom,1,4)
  if(level == 1)enol <- substr(enom,1,3)
  
  if(noM)return(enol)
  if(!noM){
    enol[wm] <- paste('M',enol[wm],sep='')
    return(enol)
  }
  
}
  
inputLST <- function(){
  
  read.csv('climate/LST_estimates/LST_0829_2014.csv',header=T)
}

calibratePrism2LSF <- function(lat=plotLat, lon=plotLon,start=1969,end=2012){
  
  se <- paste(start,end,sep='_')
  
  rname <- paste('climate/tmax_800m',se,'monthly.RData',sep='_')
  sname <- paste('climate/tmin_800m',se,'monthly.RData',sep='_')
  load(rname)
    #  tmax.dat <- out
  load(sname)
    #  tmin.dat <- out
  prismPlots  <- as.numeric( rownames(tmax.dat) )
  prismLonLat <- tmin.dat[,c('lon','lat')]
  tmax.dat   <- tmax.dat[,-c(1,2)]
  tmin.dat   <- tmin.dat[,-c(1,2)]
    
  tmp       <- inputLST()
  lstPlots  <- tmp[,'plot']
  lstLonLat <- tmp[,c('lon','lat')]
  lst.dat   <- as.matrix(tmp[,-c(1:3)])
  
  months <- getMonthNames()
  
  lname  <- capFirstLetter(colnames(lst.dat))
  tmp    <- matrix( unlist( strsplit(lname,'_') ),ncol=2,byrow=2)
  lmonth <- match(tmp[,1],months)
  lyear  <- as.numeric(tmp[,2])
  
  tmp <- nn2(prismLonLat[,c('lat','lon')],lstLonLat[,c('lat','lon')],k = 1  )$nn.idx

  tmax    <- as.matrix( tmax.dat[tmp,] )
  tmin    <- as.matrix( tmin.dat[tmp,] )
  pLonLat <- as.matrix(prismLonLat[tmp,])
  
  pname  <- matrix( unlist( strsplit(colnames(tmax),'[.]') ),ncol=2,byrow=2)
  pmonth <- as.numeric(pname[,2])
  pyear  <- as.numeric(pname[,1])
                       
  lcols    <- lyear + lmonth/12
  pmatch   <- match(lcols,pyear + pmonth/12)
  moreCols <- which(!is.finite(pmatch))
  lmatch   <- which(is.finite(pmatch))
  pmatch   <- pmatch[lmatch]
  
  tminP <- as.vector(tmin[,pmatch])
  tmaxP <- as.vector(tmax[,pmatch])
  
  xx <- cbind( tminP, tmaxP )  #tmin,tmax
  
  mm <- matrix( pLonLat[,'lon'],nrow(tmin),length(pmatch))   #lon
  lo <- as.vector(mm)
  mm <- matrix( pLonLat[,'lat'],nrow(tmin),length(pmatch))   #lat
  la <- as.vector(mm)
  
  xx <- cbind(xx,lo,la)
  colnames(xx) <- c('tmin','tmax','lon','lat')
  
  meanLo <- mean(lo)
  meanLa <- mean(la)
  meanMin <- mean(tminP)
  meanMax <- mean(tmaxP)
  
  ml <- cbind( (lo - meanLo)*(tminP - meanMin),
               (lo - meanLo)*(tmaxP - meanMax),
               (la - meanLa)*(tminP - meanMin),
               (la - meanLa)*(tmaxP - meanMax))
  colnames(ml) <- c('lonXtmin','lonXtmax','latXtmin','latXtmax')
  
  xx <- cbind(xx,ml)
  
  mm <- matrix( pmonth[pmatch], nrow(tmin),length(pmatch),byrow=T)   #month
  mm <- as.vector(mm)
  tmp <- matrix(0,nrow(xx),12)
  tmp[cbind(c(1:nrow(xx)),mm)] <- 1
  colnames(tmp) <- months
  
  xx <- cbind(xx,tmp)
  
 # lobm <- tmp*(lo - meanLo)
 # colnames(lobm) <- paste('lon',months,sep='_')
 # labm <- tmp*(la - meanLa)
 # colnames(labm) <- paste('lat',months,sep='_')
  
 # lola <- cbind( tmp,lobm,labm )
  
 # xx <- cbind(xx,lola)
  yy <- as.vector(lst.dat[,lmatch])
  
  ww   <- which(is.finite(yy))
  bhat <- solve(crossprod(xx[ww,]))%*%crossprod(xx[ww,],yy[ww])
  SSreg <- sum( (yy[ww] - xx[ww,]%*%bhat)^2 )
  SStot <- sum(yy[ww]^2)
  sig  <- SSreg/(length(ww) - length(bhat))
  rsq  <- 1 - SSreg/SStot
  bvar <- sig*solve(crossprod(xx[ww,]))
  bse  <- sqrt(diag(bvar))
  
  bcoeff <- cbind(signif(bhat,4),signif(bse,4))
  colnames(bcoeff) <- c('estimate','SE')
  
  meanVars <- c(meanLo,meanLa,meanMin,meanMax)
  names(meanVars) <- c('lon','lat','tmin','tmax')

  
  tminP <- as.vector(tmin)
  tmaxP <- as.vector(tmax)
  
  xx <- cbind( tminP,tmaxP )
  
  mm <- matrix( pLonLat[,'lon'],nrow(tmin),length(pmonth))   #lon
  lo <- as.vector(mm)
  mm <- matrix( pLonLat[,'lat'],nrow(tmin),length(pmonth))   #lat
  la <- as.vector(mm)
  
  xx <- cbind(xx,lo,la)
  colnames(xx) <- c('tmin','tmax','lon','lat')
  
  ml <- cbind( (lo - meanLo)*(tminP - meanMin),
               (lo - meanLo)*(tmaxP - meanMax),
               (la - meanLa)*(tminP - meanMin),
               (la - meanLa)*(tmaxP - meanMax))
  colnames(ml) <- c('lonXtmin','lonXtmax','latXtmin','latXtmax')
  
  xx <- cbind(xx,ml)
  
  mm <- matrix( pmonth, nrow(tmin),length(pmonth),byrow=T)
  mm <- as.vector(mm)
  tmp <- matrix(0,nrow(xx),12)
  tmp[cbind(c(1:nrow(xx)),mm)] <- 1
  xx <- cbind(xx,tmp)
  
  yy <- xx%*%bhat
  calibLST <- matrix(yy,nrow(tmin),ncol(tmin))
  
  ww <- which(is.finite(calibLST[,pmatch]) & is.finite(lst.dat[,lmatch]) )
  
  rmspe <- sqrt( mean( ( calibLST[,pmatch][ww] - lst.dat[,lmatch][ww] )^2 ) )
  
  calib <- list(coefficients = bcoeff, sigma = signif(sig,4),
                rsq = round(rsq,3),bVAR = signif(bvar,3),meanVars = meanVars,
                rmspe = rmspe)
  print(calib)
  
  calibLST[,pmatch] <- lst.dat[,lmatch]
  colnames(calibLST) <- colnames(tmin)
  moreL    <- lst.dat[,moreCols]
  colnames(moreL) <- round( lyear[moreCols] + lmonth[moreCols]*.01,2 )
  calibLST <- round(cbind(calibLST,moreL),2)
  
  calibLST <- cbind(lstLonLat,calibLST)
  
  rownames(calibLST) <- rownames(tmin)
  
  save(calibLST,file='climate/LST_estimates/calibLST.Rdata')
  invisible( calib )

}
  

inClimate <- function(lat=plotLat, lon=plotLon,start=1969,end=2012,LST=T,PADLAST=F){
  
  #LST     - use modis LST
  #PADLAST - pad missing values at end with most recent available

  require(RANN)
  
  se <- paste(start,end,sep='_')
  pname <- paste('climate/ppt_800m',se,'monthly.RData',sep='_')
  load(pname)
  # ppt.dat <- out
  precDat    <- ppt.dat[,-c(1:2)]
  precPlots  <- as.numeric( rownames(ppt.dat) )
  precLonLat <- ppt.dat[,c('lon','lat')]
  
  if(!LST){
    rname <- paste('climate/tmax_800m',se,'monthly.RData',sep='_')
    sname <- paste('climate/tmin_800m',se,'monthly.RData',sep='_')
    load(rname)
    #  tmax.dat <- out
    load(sname)
    #  tmin.dat <- out
  
    temp <- tmax.dat
    temp[,-c(1:2)] <- ( tmax.dat[,-c(1:2)] + tmin.dat[,-c(1:2)] )/2
    tempLonLat <- temp[,c('lon','lat')]
    tmp        <- nn2(tempLonLat[,c('lat','lon')], cbind(lat,lon),k = 1  )$nn.idx
    temp <- tempDat <- as.matrix(temp[tmp,-c(1:2)])
  }

  if(LST){
    
    load('climate/LST_estimates/calibLST.Rdata')
    
    tempPlots  <- rownames( calibLST )
    tempLonLat <- calibLST[,c('lon','lat')]
    tempDat    <- calibLST[,-c(1:2)]
    
    tmp  <- nn2(tempLonLat[,c('lat','lon')], cbind(lat,lon),k = 1  )$nn.idx
    temp <- as.matrix(tempDat[tmp,])
  }
  
  wrr <- 1
  
  while(length(wrr) > 0){
    tmp  <- nn2(ppt.dat[,c('lat','lon')], cbind(lat,lon),k = 1  )$nn.idx
    prec <- as.matrix(precDat[tmp,])
    
    ps <- rowSums(prec*0 + 1)
    wrr <- which(is.na(ps))
    
    if(length(wrr) > 0)precDat <- precDat[-tmp[wrr],]
  }

  mcol <- match(colnames(temp),colnames(prec))  #temp has more data (columns) than prec
  ww   <- which(is.na(mcol))
  if(length(ww) > 0){
    tmp <- matrix(NA,nrow(prec),length(ww))
    colnames(tmp) <- colnames(temp)[ww]
    prec <- cbind(prec,tmp)
  }
  
  wshort <- which(nchar(colnames(temp)) == 6)
  if(length(wshort) > 0)colnames(temp)[wshort] <- pasteChar2End(colnames(temp)[wshort],'0')
  wshort <- which(nchar(colnames(prec)) == 6)
  if(length(wshort) > 0)colnames(prec)[wshort] <- pasteChar2End(colnames(prec)[wshort],'0')
  
  if(PADLAST){
    prec <- padClimateData(prec,startYr = 2011)
    temp <- padClimateData(temp,startYr = 2011)
  }
    
  list(temp = temp, prec = prec) 
    
}

padClimateData <- function(data,startYr=2010,endYr){
  
  # pad missing values at end with most recent available
  # colnames - "1969.01" "1969.02" "1969.03" "1969.04" "1969.05"...
  # one row per location
  # startYr - earliest yr to pad
  
  numericNames <- as.numeric(colnames(data))
  yrCol <- round( numericNames,0 )
  moCol <- round( 100*( numericNames - yrCol ), 0 )
  ncol  <- length(yrCol)
  
  if(missing(endYr))endYr <- max( yrCol )
  
  maxYrData <- yrCol[ncol]
  cols4Last <- numeric(0)         
  
  if(moCol[ncol] != 12){               #remaining months of last year
    mm        <- 12 - moCol[ncol]
    moreCols <- matrix(NA,nrow(data),mm) 
    mnames   <- substr(as.character(c(1:12)/100),2,4)[(moCol[ncol]+1):12]
    colnames(moreCols) <- paste(maxYrData,mnames,sep='')
    data <- cbind(data,moreCols)
  }
  if(maxYrData < endYr){              #remaining years
    moreCols <- matrix(NA,nrow(data),12*(endYr - maxYrData))
    mnames   <- substr(as.character(c(1:12)/100),2,4)
    cnames   <- as.vector( t( outer( (maxYrData+1):endYr,mnames,paste,sep='') ) )
    colnames(moreCols) <- cnames
    data <- cbind(data,moreCols)
  }
  wshort <- which(nchar(colnames(data)) == 6)
  if(length(wshort) > 0)colnames(data)[wshort] <- pasteChar2End(colnames(data)[wshort],'0')
    
  wyr   <- max( grep(startYr + 1,colnames(data)) ) 
  wmiss <- which(!is.finite(data),arr.ind=T)
  wmiss <- wmiss[wmiss[,2] > wyr,]
  
  while(length(wmiss) > 0){
    wcols <- sort(unique(wmiss[,2]))
    yrMiss <- colnames(data)[wcols]
    yrB4   <- as.character( as.numeric(yrMiss) - 1 )
    wshort <- which(nchar(yrB4) == 6)
    if(length(wshort) > 0)yrB4[wshort] <- pasteChar2End(yrB4[wshort],'0')
    
    yrB4   <- match(yrB4,colnames(data))
    wfrom  <- yrB4[ match( colnames(data)[wmiss[,2]],yrMiss ) ]
    data[wmiss] <- data[ cbind(wmiss[,1],wfrom) ]
    wmiss <- which(!is.finite(data),arr.ind=T)
    if(length(wmiss) == 2)wmiss <- matrix(wmiss,1,2)
  }
  data
}

size2BA <- function(count,diam,plotarea,kk){  
  
  #input is matrix with species-size
  #returns BA
  #kdex is column index for each size class
  
  nb <- max(kk)
  
  baBySize <- matrix(NA,nrow(count),nb)
  
  for(k in 1:nb){
    
    matk  <- matrix(0,nrow(count),1)
    ck    <- count[,kk == k]
    diamk <- diam[,kk == k]
    pk    <- plotarea[,kk == k]
    w1    <- which(ck > 0,arr.ind=T)
    dd    <- diamk[w1]
    cc    <- ck[w1]
    tmp   <- rep(dd,times=cc)
    pa    <- rep(pk[w1],times=cc)
    loc   <- rep(w1[,1],times=cc)
    bk    <- baPerHa(tmp,pa)
    tmp   <- byIndex(bk,loc,sum)
    matk[as.numeric(rownames(tmp))] <- tmp
    
    baBySize[,k] <-  matk
  }
  
  baBySize
}

getFIAX <- function(vnames,xplot,xyear,jt=jtFIA,plotData=plotDataFIA,
                    regClim=regClimFIA,CENTER=F,STANDARD=F){
  
  pall <- unique(xplot)
  np   <- length(pall)
  
  p <- length(vnames)
  
  x <- matrix(1,length(xplot),p)
  colnames(x) <- vnames
  
 # centSize <- (breaks - mean(breaks))/sd(breaks)
#  centBA   <- (baLarge - mean(baLarge))/sd(baLarge)
  
  if('phys' %in% vnames){
    phys <- plotData[jt[,'j'],'phys']
    x[,'phys'] <- trunc(phys/10,0)
  }
  
  wt <- grep('temp',vnames)
  w2 <- grep('temp2',vnames)
  if(length(wt) > 0)x[,wt[1]] <- regClim[,grep('temp',colnames(regClim))]
  if(length(w2) > 0)x[,w2]    <- regClim[,grep('temp',colnames(regClim))]^2
  
  wt <- grep('prec',vnames)
  w2 <- grep('prec2',vnames)
  if(length(wt) > 0)x[,wt[1]] <- regClim[,grep('prec',colnames(regClim))]
  if(length(w2) > 0)x[,w2]    <- regClim[,grep('prec',colnames(regClim))]^2
  
  wt <- grep('therm',vnames)
  w2 <- grep('therm2',vnames)
  if(length(wt) > 0)x[,wt[1]] <- regClim[,grep('therm',colnames(regClim))]
  if(length(w2) > 0)x[,w2]    <- regClim[,grep('therm',colnames(regClim))]^2
  
  wt <- grep('deficit',vnames)
  w2 <- grep('deficit2',vnames)
  if(length(wt) > 0)x[,wt[1]] <- regClim[,grep('deficit',colnames(regClim))]
  if(length(w2) > 0)x[,w2]    <- regClim[,grep('deficit',colnames(regClim))]^2
  
  
  if('mesic' %in% vnames){
    x[,'xeric'] <- x[,'mesic'] <- 0
    hvec <- plotDataFIA[,'phys']
    hvec[hvec < 20] <- 1
    hvec[hvec %in% c(21:22)] <- 2
    hvec[hvec > 2] <- 3
    plotData[,'phys'] <- hvec
    
    hvec <- hvec[jt[,'j']]
    w1 <- grep('xeric',vnames)
    w2 <- grep('mesic',vnames)
    x[hvec == 1,w1[1]] <- 1
    x[hvec == 3,w2[1]] <- 1

  }
    xpriorXf <- xpriorYf  <- numeric(0)
    
  if(!is.null(linearXpred)){
    
    phys <- plotData[jt[,'j'],'phys']
    # elev <- plotData[jt[,'j'],'elev']
    slope <- plotData[jt[,'j'],'slope']/100
    asp   <- degrees2radians(plotData[jt[,'j'],'aspect'])
    u1    <- slope*sin(asp)
    u2    <- slope*cos(asp)
    intercept <- rep(1,length(u1))
    xpriorXf <- cbind(intercept,phys,slope,u1,u2,regClim[jt[,'j'],])
    xpriorYf <- xpriorXf%*%betax
  }
  
  list( x = x, xpriorX = xpriorXf, xpriorY = xpriorYf, regClim = regClim)
}

updateX <- function(){       #CHANGE TO ANOMALIES
  
  xx <- xall
  
  wt <- which(linearXpred %in% c('temp','prec'))
  wtt <- linearXcols[wt]
  x0 <- xall[,-wtt]
  b0 <- beta[-wtt,]
  b1 <- beta[wtt,] 
  
  v  <- xpriorY[,wt]*0
  if(!is.matrix(v))v <- matrix(v,ncol=1)
  if(!is.matrix(b1)) b1 <- matrix(b1,nrow=1)
  
  bx <- b1%*%sinv
  v[groRegIndexF,]  <- t( bx%*%t(gall[groRegIndexF,] - x0[groRegIndexL,]%*%b0) )
  VI  <- bx%*%t(b1)
  
  #time varying (temp, prec)
  vt <- matrix(0,ntnp,length(wt))
  
  for(t in 1:length(wt)){
    
    tmp <- byFunctionRcpp(v[,t],longSpecJ[1:ntnpnb],longSpecT[1:ntnpnb],
                          matPlotTime*0,matPlotTime*0,MEAN=F) 
    vt[,t] <- tmp[jtdex]
  }
  
  if(length(wt) == 1){
    vt <- vt + xpriorY[1:ntnp,wt]*priorIVX[wt,wt] 
    V  <- 1/( nbreak*VI + priorIVX[wt,wt] )
  }
  if(length(wt) > 1){
    vt <- vt + xpriorY[1:ntnp,wt]%*%priorIVX[wt,wt]
    V  <- solve(nbreak*VI + priorIVX[wt,wt])
  }
  
  mu <- vt%*%V
  
  tmp <- myrmvnorm(nrow(mu),mu,V)
  tmp <- tmp[rep( c(1:ntnp),nbreak),]
  tmp[flong,] <- xpriorY[flong,wt]
  
  tmp[longSpecJ[1:ntnpnb] <= nplotDem,] <-  xpriorY[longSpecJ[1:ntnpnb] <= nplotDem,wt]
  xx[,wtt] <- tmp
  
  
  # constant for plot, only hydro at the moment
  wt <- which(linearXpred == 'hydro')
  wtt <- linearXcols[wt]
  
  x0 <- xall[,-wtt]
  b0 <- beta[-wtt,]
  b1 <- matrix(beta[wtt,],nrow=length(wt))
  
  v  <- matrix(xpriorY[,wt]*0,ncol=length(wt))
  
  bx <- b1%*%sinv
  v[-flong,]  <- t( bx%*%t(gall[-flong,] - x0[-llong,]%*%b0) )
  VI  <- bx%*%t(b1)
  
  vt <- matrix(0,nplot,length(wt))
  
  for(t in 1:length(wt)){
    vt[,t]   <- byIndex(v[,t],list(plot=longSpecJ[1:ntnpnb]),coerce=T,sum,na.rm=T)
  }
  
#  if(length(wt) > 1){
#    vt <- vt + xpriorY[1:nplot,wt]%*%priorIVX[wt,wt]
#    V  <- solve(nbreak*VI + priorIVX[wt,wt])
#  }
  if(length(wt) == 1){
    xprior <- as.vector(unlist(by(xpriorY[1:ntnp,wt],jtdex[,'j'],mean)) )
    vt <- vt + xprior*priorIVX[wt,wt]
    V  <- 1/(nbreak*yrByPlot*VI[1] + priorIVX[wt,wt])
    mu <- as.vector(vt)*V
    tmp <- rnorm(nplot,mu,sqrt(V))
    tmp[1:nplotDem] <- xprior[1:nplotDem]
    tmp <- tmp[jtdex[,'j']]
    tmp <- rep(tmp,nbreak)
  }
  xx[,wtt] <- tmp

  xx
}

whichHydroInt <- function(vnames){
  
  wint <- c( grep('hydroX',vnames), grep('Xhydro',vnames) )
  main <- unlist(strsplit(vnames[wint],'X'))
  main <- main[main != 'hydro']
  
  list(wint = wint, main = match(main,vnames) )
  
}

updateXhydro <- function(){          #in progress
  
  hnow  <- xall[xind[firstTime],'hydro']
  hprop <- tnorm(msample,0,1,hnow,.01)
  
  tmp <- whichHydroInt(xnames)
  wint1 <- tmp$wint
  main1 <- tmp$main
  
  x1 <- xall[xind,]
  x1[,'hydro'] <- hprop[longstageJ]
  if(length(wint1) > 0)x1[,wint1]    <- matrix(x1[,'hydro'],nrow(x1),length(main1))*x1[,main1]
  
  pnow1 <- dmvnormZeroMean(gallG[-xfirstG,] - xall[xind[-xlastG],]%*%beta ,sinv=sinv)
  pnew1 <- dmvnormZeroMean(gallG[-xfirstG,] - x1[-xlastG,]%*%beta ,sinv=sinv)
  
  yy = short2long( logit(thetaG),n1=ntnpG,n2=ntnpnbG,n3=ntnpnsG )[-xlastG,]
  
  pnow2 <- dmvnormZeroMean(yy - xall[xind[-xlastG],mnames]%*%betaMort ,sinv=minv)
  pnew2 <- dmvnormZeroMean(yy - x1[-xlastG,mnames]%*%betaMort ,sinv=minv)
  
  pnow1 <- pnow1 + pnow2
  pnew1 <- pnew1 + pnew2
  
  tmp <- whichHydroInt(fnames)
  wint2 <- tmp$wint
  main2 <- tmp$main
  
  x2 <- fall[wind,][fecRegIndexG,]
  x2[,'hydro'] <- hprop[jtind[fecRegIndexG,'j']]
  if(length(wint2) > 0)x2[,wint2]    <- matrix(x2[,'hydro'],nrow(x2),length(main2))*x2[,main2]
  
  pnow2 <- dmvnormZeroMean(phi[wind,][fecRegIndexG,] - fall[fecRegIndexG,]%*%betaFec ,sinv=finv)
  pnew2 <- dmvnormZeroMean(phi[wind,][fecRegIndexG,] - x2%*%betaFec ,sinv=finv)
  
  tmp <- whichHydroInt(gnames)
  wint3 <- tmp$wint
  main3 <- tmp$main
  
  x3 <- xgall[msamp,]
  x3[,'hydro'] <- hprop
  if(length(wint3) > 0)x3[,wint3]    <- matrix(x3[,'hydro'],nrow(x3),length(main3))*x3[,main3]
  
  pnow3 <- dmvnormZeroMean(wmean - xgall[msamp,]%*%betaGam ,sinv=ginv)
  pnew3 <- dmvnormZeroMean(wmean - x3%*%betaGam ,sinv=ginv)
  
  pr1 <- byIndex(pnew1 - pnow1, longstageJ[xfirstG],sum,coerce=T)
  pr2 <- as.vector( byIndex(pnew2 - pnow2, jtind[fecRegIndexG,'j'],sum,coerce=T) )
  
  a <- exp(pr1 + pr2 + pnew3 - pnow3)
  
  ww <- which(runif(msample,0,1) < a)
  hnow[ww] <- hprop[ww]
  
  x1[,'hydro'] <- hnow[longstageJ]
  if(length(wint1) > 0)x1[,wint1]    <- matrix(x1[,'hydro'],nrow(x1),length(main1))*x1[,main1]
  
  x2[,'hydro'] <- hprop[jtind[fecRegIndexG,'j']]
  if(length(wint2) > 0)x2[,wint2]    <- matrix(x2[,'hydro'],nrow(x2),length(main2))*x2[,main2]
  
  x3[,'hydro'] <- hnow
  if(length(wint3) > 0)x3[,wint3]    <- matrix(x3[,'hydro'],nrow(x3),length(main3))*x3[,main3]
  
  list(x = x1, xfec = x2, xg = x3)
  
}
 
makeX <- function(){

  xmat <- matrix(1,nplot*nt,length(xnames))
  colnames(xmat) <- xnames
  rownames(xmat) <- plots

  phys <- by(data[,'physclcd1'],data[,'plt1'],min,na.rm=T)
  phys <- trunc(phys/10,0)

  xmat[names(phys),'phys'] <- as.vector(phys)

  if('temp' %in% xnames)xmat[as.character(plots),'temp'] <- temp
  if('prec' %in% xnames)xmat[as.character(plots),'prec'] <- prec

  xmat
}


makeM <- function(){              #propagator matrix for beta sampling

  up   <- (gam2 - gam1*(1 - rho*matrix(dt,n,nsnb)))[,k0Index]
  down <- (gam1[,k1Index] - gam2[,k2Index])
  M    <- up/down
  zero <- which(down == 0,arr.ind=T)

  list(M = M, zero = zero)

}

getTrans <- function(tindex,kindex){
  
  groTrans(tindex,kindex) - survTrans(tindex,kindex)
}

groTrans <- function(tindex,kindex,gro){
  
  gro[tindex,kindex]*dtMat[tindex,kindex]/dxMat[tindex,kindex] 
}

survTrans <-  function(tindex,kindex,tmat=thetaMat){
   
    tmat[tindex,kindex]/2*dtMat[tindex,kindex]
}


initGamma <- function(){
  
  gam <- y*0
  
  for(j in 1:nspec){
    
    sumList <- list(plot=sizeSpecData[,'j'],yr=sizeSpecData[,'year'])
    
    tmp <- as.matrix( by(sizeBySpec[,sdex == j],sumList,sum,na.rm=T) )
    tmp <- tmp[1:nrow(tmp),]
    
    fj <- apply(tmp,1,max,na.rm=T)
    fj[fj < 1] <- 1
    
    gam[,sdex == j] <- fj[jtdex[,'j']]
  }
  
  
  gam[tnext,know] <- gam[tnow,klast]*getTrans(tnow,klast) + gam[tnow,know]*(1 - getTrans(tnow,klast))
  
  gam <- jitter(gam)/sampleArea
  
  gam
}


make2d <- function(xx,nr=ntnpnb,nc=nspec){
  
  # xx is nt*nplot, nsnb or 3d array nt*nplot, nbreak, nspec
  # returns nt*nplot*nbreak, nspec
  
  matrix(xx,nr,nc )
}
  

make3d  <- function(xx,n1=ntnp){     
  
  # xx is nt*nplot, nsnb
  
  tmp <- array(xx,c(n1,nbreak,nspec))
  
  dimnames(tmp)[[2]] <- breaks
  dimnames(tmp)[[3]] <- specs
#  if(!is.null(rownames(y)))dimnames(tmp)[[1]] <- rownames(y)[wind]
  
  tmp
  
}

stack3d <- function(xx,margin='SPEC',n1=ntnp,n2=ntnpnb,n3=ntnpns){  
  
  # xx is 3d array nt*nplot, nbreak, nspec
  # output is nt*nplot*nbreak, nspec (margin = 'SPEC') or
  #           nt*nplot*nspec, nbreak (margin = 'STAGE')
  
  if(margin == 'SPEC') tmp <- make2d( xx,nr=n2,nc=nspec )
  
  if(margin == 'STAGE'){
    
    ii <- cbind( rep( c(1:n1),c(nsnb) ), 
                 rep( c(1:nbreak),each=(n3)),
                 rep( rep( c(1:nspec), each=n1,nbreak) ) )

    tmp <- matrix( xx[ii],n3,nbreak )
  }
  tmp
}



predGro <- function(nsim=100,xx=xall,gmat=groMat,np=nplot,jt=jtdex){
  
  gs    <- sample(burnin:ng,size=nsim,replace=T)
  gpred <- gpred2 <- gmat*0
  nall  <- nrow(gmat)
  pjs   <- matrix(0,np,nspec)
  
  ji <- rep(jt[,'j'],(nbreak*nspec))
  si <- rep(c(1:nspec),each=nall*nbreak)
  
  for(j in 1:nsim){
    
    b  <- matrix(chainList$betaChain[gs[j],],p,nspec)
    s  <- matrix(chainList$sigmaChain[gs[j],],nspec,nspec)
    m  <- xx%*%b
    yy <- matrix( myrmvnorm(nrow(xx),m,s), nrow(gmat),ncol(gmat))
    gpred  <- gpred + yy
    gpred2 <- gpred2 + yy^2
    
    gg  <- byFunction(as.vector(yy),ji,si,mat=matrix(0,np,nspec),mean)
    mu  <- byFunction(as.vector(m),ji,si,mat=matrix(0,np,nspec),mean)
    
    pjs <- pjs +  exp( pmvnormCondRcpp(q=0,x=gg,mu,smat=s,whichVar=c(1:nspec)) )
    
  }
  gmu <- gpred/nsim
  gsd <- sqrt(gpred2/nsim - gmu^2)
  
  pjs <- pjs/nsim
  
  list(mu = gmu, sd = gsd, pabsGro = pjs)
}


groThetaInit <- function(){
  
  groMu <- make2d(yGrowMu,nr=ntnpnb,nc=nspec) 
  logitTheta  <- short2long( logit(thetaMat),n1=ntnp,n2=ntnpnb,n3=ntnpns )

  groRange <- matrix(NA,2,nspec)
  colnames(groRange) <- specs
  groInit <- thetaInit <- yGrowMu*0
  
  for(j in 1:nspec){
    
    absIndex <- rep(y0[jtdex[,'j'],j],nbreak)
    hi <- (1 - absIndex)*2
    lo <- -absIndex*2
    
    wm  <- which(hi > 0)
    
    muj <- mean(groMu[wm,j])
    rnj <- range(groMu[wm,j])
    if(rnj[1] < .001)rnj[1] <- .001
    if(rnj[2] > 2)rnj[2] <- 2
    wj <- wm[wm %in% groRegIndexL]
    wy <- wm[wm %in% groRegIndexF]
    xx <- xall[wj,]
    rx <- qr(xx)$rank
    
    if(rx < ncol(xx))mu <- muj
    
    if(rx == ncol(xx)){
      bj <- solve(crossprod(xx))%*%crossprod(xx,groMu[wy,j])
      bj <- tnorm.mvtRcpp(bj, bj, diag(.001,rx), loB[,1], hiB[,1])
      mu <- xall%*%t(bj)
    }
    
    pj <- rep(0,ntnpnb)
    pj     <- tnorm(ntnpnb,lo,hi,muj,.1)
    pj[wm] <- groMu[wm,j]
    groMu[,j] <- pj
    groRange[,j] <- rnj
    groInit[,sdex == j] <- pj
    
    muj <- mean(logitTheta[wm,j])
    xx  <- xall[wj,mnames]
    rx  <- qr(xx)$rank
    
    if(rx < ncol(xx))mu <- muj
    
    if(rx == ncol(xx)){
      bt  <- solve(crossprod(xx))%*%crossprod(xx,logitTheta[wy,j])
      bt <- tnorm.mvtRcpp(bt, bt, diag(.001,rx), loM[,1], hiM[,1])
      muj <- xall[,mnames]%*%t(bt)
    }

    pj     <- tnorm(ntnpnb,logit(.002),0,muj,.2)
    logitTheta[,j] <- pj
    thetaInit[,sdex == j] <- invlogit(pj)
  }
  
  groRange[groRange < .001] <- .001
  loGro <- matrix(rep(groRange[1,],each=nbreak),ntnp,nsnb,byrow=T)
  hiGro <- matrix(rep(groRange[2,],each=nbreak),ntnp,nsnb,byrow=T)
  
  hiGro[whichY] <- yGrowMu[whichY] + 2*yGrowSd[whichY]
  loGro[whichY] <- yGrowMu[whichY] - 2*yGrowSd[whichY]
  loGro[loGro < .001] <- .001
  
  hiGro[hiGro > 2] <- 2
  hiGro[hiGro < loGro] <- loGro[hiGro < loGro]*1.2
  groInit[groInit > hiGro] <- hiGro[groInit > hiGro]

  list(loGro = loGro, hiGro = hiGro,groMat = groInit, 
       logitTheta = logitTheta, thetaMat = thetaMat)
}
  
mortLimits <- function(){
  
 # tmeans <- 1 - colMeans( 1 - (1 - ySurv/y)^(1/dtMat),na.rm=T)
  
  lotheta <- .05*exp(-breaks[wideKmat]) + .005
  hitheta <- .4*exp(-.03*breaks[wideKmat])

  survInit <- ySurv/y
  survInit[!is.finite(survInit)] <- .99
  
  thetaMat <- 1 - survInit^(1/dtMat)
  ww <- which(yGrowNo > 4)
  lotheta[ww] <- thetaMat[ww] - .2
  hitheta[ww] <- thetaMat[ww] + .2
  
  lotheta[lotheta < 1e-3] <- 1e-3
  hitheta[hitheta > .9999] <- .9999
  
  hitheta[hitheta < lotheta] <- lotheta[hitheta < lotheta]*1.000001
  
  list( lotheta = matrix(lotheta,ntnp,nsnb), hitheta = matrix(hitheta,ntnp,nsnb) )
}
  
  
  
recruitLimit <- function(gg=gam){
  
 # tiny <- 1e-6
  yno  <- y0[jtdex[,'j'],]
  fno  <- f0[jtdex[,'j'],]
 # wyno <- which(yno == 1)
  wf   <- which(fno == 0)
  
  rr  <- getRecruit(gg)   # basal area above maturation diameter
  cc <- 1/dtMat[,fstage]/sampArea[,fstage]/(rr + 1)
  
  xfec <- fecMat + 1
  
  sdf <- cc*sqrt(xfec)
  
  hiP <- xfec - 1 + 2*sdf
  loP <- xfec - 1 - 2*sdf
  
  phi <- (hiP + loP)/2
  phi[fecMat > 0] <- cc[fecMat > 0]*fecMat[fecMat > 0]
  
  hiP[hiP > 1000] <- 1000
  
  hiI <- xfec/(10 + dtMat[,fstage]*sampArea[,fstage])
  hiI[hiI > 5] <- 2
  hiI[yno == 1 & fno == 0] <- 5         # recruits by no adults
  loI <- hiI*0
  im  <- hiI*0 + .0000001
  
  list(phi = phi, loP = loP, hiP = hiP, im = im, loI = loI, hiI = hiI)
}
  
  
  
getRecruit <- function(gg = gam){
  
  cc <- rep(c(1:ntnp),nspec*nbreak)
  
  tmp <- gg*dxMat*pi*(fecSizeMat/2)^2/10000    #basal area by tree size
 # tmp <- byFunction(as.vector(tmp),cc, wideS,matSpec*0,FUN=sum)
  tmp <- byFunctionRcpp(as.vector(tmp),cc,as.vector(wideSmat),
                        matSpec*0,matSpec*0,MEAN=F)
  colnames(tmp) <- specs
  tmp
}


updateAbsent <- function(np=msample,ggam=gamG,pGam,jdex,jtiSamp=jtind,
                         xb=xgroBeta[xfast,][-lastTime,],xbno=NULL,
                         gro2d,FJ,FS,mat=matJS,
                         windex=wind,ftimeG=firstTime,whichPlot,gamVarG = gamVar){
  #indicator for absence
  
#  pieSamp <- pie[whichPlot,]  #[pop can grow|trees can grow]   should be low outside pop range
  
  tiny <- 1e-10

  tmp  <- getProbGroAbs(xb=xb,xbno,gro2d,FJ,FS,mat,COMP=T)
  pjs  <- tmp$pjs                                                          #trees cannot grow
  pnon <- tmp$pNon                                             #trees cannot grow w/0 competition
  qjs  <- getProbPopAbs(ggam,pGam,gamVarG,windex,whichPlot,jdex,jtiSamp)   #population cannot grow
  ujs  <- getProbZeroCount(np,pGam,jdex,jtiSamp,mat,windex,ftimeG)                    #[unobs|pop can grow]
         
 # qjs[qjs < (1 - pieSamp)] <- pieSamp[qjs < (1 - pieSamp)]
  
  ww <- which(pjs > qjs)
  
  if(length(ww) > 0){
    switch <- rbinom(1,1,.5)
    if(switch == 0)pjs[ww]   <- qjs[ww]
    if(switch == 1)qjs[ww]   <- pjs[ww]
      
    wn <- which(pnon > pjs)
    pnon[wn] <- pjs[wn]
  }
  
  ww <- which(y0[whichPlot,] == 0)
  
  pjs[ww] <- pnon[ww] <- qjs[ww] <- 0
  
  pieSamp <- (1 - qjs)/(1 - pjs + tiny)
  
  vjs <- ujs*pieSamp + 1 - pieSamp                          #[not observe|trees can grow]
                
  m1 <- vjs*(1 - pjs)
  Mi <- m1/(m1 + pjs)                                       #individuals can grow
  
  m1  <- vjs*(1 - pnon)
  MiN <- m1/(m1 + pnon)
  
  m1 <- ujs*(1 - qjs)
  Mp <- m1/(m1 + qjs + tiny)                              #pop can grow
  
  Mi[Mi < Mp]   <- Mp[Mi < Mp]
  MiN[MiN < Mp] <- Mp[MiN < Mp]
  
  #Probs are [bp=0,bi=0], [bp=1,bi=0], [bp=1,bi=1]; note[bp=0|bi=1] = 0
  prob <- cbind( as.vector(Mp), as.vector(Mi - Mp), as.vector(1 - Mi) )
  tmp  <- myrmultinom(1,prob)
  
  bii <- bpp <- pjs*0
  
  bii[tmp[,3] == 1] <- 1
  bpp[tmp[,2] == 1] <- 1
  bpp[tmp[,3] == 1] <- 1
  bii[ww] <- 0
  bii[knownTreeGrow[whichPlot,] == 1]  <- 1
  bii[knownTreeGrow[whichPlot,] == -1] <- 0
  bpp[ww] <- 0
  bpp[knownAbsent[whichPlot,] == 1]    <- 1

  bpp[bii == 1] <- 1
  
  prob <- cbind( as.vector(Mp), as.vector(MiN - Mp), as.vector(1 - MiN) )
  tmp  <- myrmultinom(1,prob)
  
  biN <- pjs*0
  
  biN[tmp[,3] == 1] <- 1
  biN[ww] <- 0
  biN[bii == 0] <- 0
  bpp[biN == 1] <- 1
  
#  hi  <- qjs*0 + .95
#  lo  <- as.vector(1 - qjs)
#  lo[lo > hi] <- hi[lo > hi]*.99
#  tmp <- rtrunc(nn=length(bii),lo=lo,hi=hi,
#                p1=as.vector(1 + (1 - bpp)*(1 - bii)),
#                p2=as.vector(10 + bpp*(1 - bii)),'beta')
#  pieSamp <- matrix(tmp,nrow(bii),ncol(bii))
  
  
  list(bAbsS = bpp, bAbsG = bii, bAbsGnoComp = biN, 
       pie = pieSamp, pjs = pjs, qjs = qjs, ujs = ujs)
}

getProbGroAbs <- function(xb=xgroBeta[xfast,][-lastTime,],xbno=NULL,gro2d,  #trees cannot grow
                          FJ=fastJ,FS=fastS,mat=matJS,COMP=F){  
  
  # evaluated for a single stage, indexed by xfast
  
  pNon <- numeric(0)
  
#  gg <- ( make2d(gro,nr=ntnpnbG,nc=nspec)[xfast,] )[-lastTime,]
  
  MEAN <- T
  
  mub <- byFunctionRcpp(as.vector(xb),FJ,FS,mat*0,mat*0 ,MEAN=T)
  gg  <- byFunctionRcpp(x=as.vector(gro2d),FJ,FS,mat*0,mat*0 ,MEAN=T)
  
  pjs <- exp( pmvnormCondRcpp(q=0,xx=gg,mub,smat=sigma,whichVar=c(1:nspec)) )     #missing due to growth
  
  if(COMP){

 #   xc  <- (xallG[xfast,pc]%*%beta[pc,])[-lastTime,]   #no competition
    muc <- byFunctionRcpp(as.vector(xbno),FJ,FS,
                          summat=mat*0,totmat=mat*0 ,MEAN=T)
    pNon <- exp( pmvnormCondRcpp(q=0,xx=gg,muc,smat=sigma,whichVar=c(1:nspec)) )     #missing w/o comp
  }
  
  list(mu = mub, gg = gg, pjs = pjs, pNon = pNon)
}

getProbPopAbsOld2 <- function(xg=xgamBeta,wmeangam=wmean){
  
  qjs <- pmvnormCondRcpp(q=0,xx=wmeangam,mu=xg,smat=gsigma,whichVar=c(1:nspec))  
  
  qjs
}

getProbPopAbs <- function(ggam,pGam,gamVar,windex,whichPlot,jdex,jtiSamp){
  
  nn  <- length(whichPlot)
  ji  <- rep(jdex,nsnb)
  si  <- rep(c(1:nspec),each=length(jdex)*nbreak)
  
  tmp <- pnorm(0,pGam,sqrt(gamVar))
  qjs <- byFunctionRcpp(as.vector(tmp),ji,si,matrix(0,nn,nspec),matrix(0,nn,nspec),MEAN=T)
  qjs
}


getProbPopAbs1 <- function(ggam,pGam,gamVar,windex,whichPlot,jdex,jtiSamp,n1,n2,n3){
  
  nn     <- length(whichPlot)
  ji  <- rep(jdex,nsnb)
  si  <- rep(c(1:nspec),each=length(jdex)*nbreak)
  
  tmp <- pnorm(0,pGam,sqrt(gamVar))
  qjs <- byFunctionRcpp(as.vector(tmp),ji,si,matrix(0,nn,nspec),matrix(0,nn,nspec),MEAN=T)
  qjs
}

getProbPopAbsOld2 <- function(ggam,pGam,windex,whichPlot,jtiSamp,n1,n2,n3){
  
  nn    <- length(whichPlot)
  ggnow <- short2long(pGam,'STAGE',n1,n2,n3)    #check n2, n3

  oinv <- invertAR1(nbreak,psi,ovar)
  ev   <- eigen(oinv, only.values = TRUE)$values
  
 # lmat <- lambda[jtdex[windex,'j'],]*tmat <- dtMat[windex,kdex==1]
  
  lmat <- lvar[jtdex[windex,'j'],]*dtMat[windex,kdex==1]/dxMat[windex,kdex==1]

  lam1 <- as.vector(lmat) 
  
  tmp <- pmvnormAR1_approx(neach=10,xx=ggnow,psi,const=lam1,sig=1,lower=0,stages=c(2:4))
  ji  <- rep(jtiSamp[,'j'],nspec)
  si  <- rep(c(1:nspec),each=n1)
  
  qjs <- 1 - byFunctionRcpp(tmp,ji,si,matrix(0,nn,nspec),matrix(0,nn,nspec),MEAN=T)
  
  qjs
}


  
getProbZeroCount <- function(np,pgam,jdex,jtiSamp,mat,windex,ftimeG=firstTime){
  
  nr  <- nrow(pgam)*nbreak
  
  nnp     <- nrow(pgam)
  fjlong  <- rep(jdex,nspec*nbreak)
  fslong  <- rep(c(1:nspec),each=nnp*nbreak)
  
  pgam[pgam < 0] <- 0
  
  tmp <-  make2d(dxMat[windex,][-ftimeG,]*pgam*sampArea[windex,][-ftimeG,],nr=nr) 
  tmp[tmp < 0] <- 0
  
  ujs <-  byFunctionRcpp(as.vector(tmp),fjlong,fslong,
                         summat=mat*0,totmat=mat*0 ,MEAN=F)  # Pr(y = 0)
  
  #marginalize out gamma:
#  tmp <- dxMat[windex,1]*sampArea[windex,1]
#  ujs <- byFunctionRcpp(tmp,jtiSamp[,1],tmp*0+1,matrix(0,np),matrix(0,np),MEAN=T)
#  ujs <- 1/matrix(ujs,np,nspec)
  
  exp(-ujs)
}

probAbsentOld <- function(gam,beta,sigma){               #not used?
  
  #tree growth
    
  xb <- (xall[longSampFast,]%*%beta)[-ftime,]
  gg <- (make2d( groMat )[longSampFast,])[-ftime,]
    
  tmp <- getProbGroAbs(xb,gg,sigma)
  pjs <- tmp$prob                    # Pr(growth > 0|growth other spp)
  
  gg <- short2long(wgam,'STAGE')
  gg <- short2long(gg,'STAGE')
  
  tmp   <- omegaLambdaCDF(gg,lambda,psi,ovar,zeroIndex=absLong)
  qjs   <-  byFunctionRcpp(as.vector(tmp[-ftime,]),i=fastJ,j=fastS,
                           summat=matPlotSpec*0,totmat=matPlotSpec*0 ,MEAN=T)
    
#  tmp <- as.vector( make2d(dxMat*gam*sampArea)[longSampFast,] )[-flongSpec]
  
  tmp <-  make2d(dxMat*gam*sampArea)[longSampFast,] 
  gsj <-  byFunctionRcpp(as.vector(tmp[-ftime,]),i=fastJ,j=fastS,
                         summat=matPlotSpec*0,totmat=matPlotSpec*0 ,MEAN=T)
  u    <- exp(-gsj)                     # Pr(y = 0)
  
  logRegion <- log(pjs) - log( (pjs + (1 - pjs)*u) )  # Pr(tree growth > 0|growth other spp, y = 0)
  logLocal  <- log(qjs) - log( (qjs + (1 - qjs)*u) )  # Pr(pop growth > 0|growth other spp, y = 0)
  
  list(qjs = qjs, pjs = pjs, zeroPoisson = u, logLocal = logLocal, logRegion = logRegion)
}

omegaLambdaCDF <- function(gg,lambda,psi,ovar,zeroIndex=NULL){  # not log
  
  mm <- 1:5  # size class when not marginalizing over many
  
  mm <- 5
  nm <- length(mm)
  
  p1 <- matrix(NA,ntnp,nspec)
  d1 <- 1
  
  if(is.null(zeroIndex))zeroIndex <- rep(0,nrow(gg))
  
  for(s in 1:nspec){
    
  #  ss  <- which(longStageS[1:ntnpns] == s & zeroIndex == 0)  #species index
    
    ss  <- which(longStageS[1:ntnpns] == s)
    ds  <- dtLong[ss,1]
    ss  <- ss[ds > 0]
    ds  <- dtLong[ss,1]
    tt  <- c(1:ntnp)[ss - (s-1)*ntnp]
    
    ds[ds > 1] <- 1  #remove time increment from variance
    
    jj  <- longStageJ[1:(ntnpns)][ss]
    ls  <- lambda[jj,s]
    
    g1 <- gg[ss,]
    
  #  mm <- max(last[ss])
    
    v1 <- sqrt(rep(ds*ls*ovar,mm))
    
    if(length(ss) == 1)g1 <- matrix(g1,1)
    
    tmp <- matrix( pnorm(0,mean=as.vector(g1[,mm]),sd=v1,log=T,lower.tail=T) ,nrow(g1) ,nm)
    p1[tt,s] <- exp(rowSums(tmp))
    
  }
  
  #  p1[absMat[,kdex == 1] == 1]  <- NA   
  p1
}




plotSpec2fullFormat <- function(xx,jdex){
  
  # xplot is nplot X nspec
  # returns ntnp X nsnb
  
  rr <- xx[ jdex[,'j'], ]
  rr[,sdex]
  
}


getGamMatrixOld <- function( ffec=fecundMat,ggam=gam,ttheta=thetaMat,
                         gro=groMat,USEMU=F,index=c(1:ntnp), kstep = 2 ){
  
  dtdx <- dtMat[index,]/dxMat[index,]
  mmat <- gro*dtdx
  
  ntnpI <- length(index)
  
  if(USEMU){
    mu   <- xgroBeta
    mmat <- matrix(mu,ntnpG,nbreak*nspec)*dtMat[index,]/dxMat[index,]  #mean growth rate
  }
  
  mmat[absTreeMat[index,] == 1] <- 0
  
  #make kernel
  
  sdMat <- matrix( sqrt(diag(sigma)), ntnpI, nsnb,byrow=T )*dtdx
  
  #from same class
  kdist <- dxMat[index,]
  p0    <- pnorm(kdist,mmat,sdMat)  
  ptot  <- p0
  
  gg <- ggam*p0*(1 - ttheta)^dtMat[index,]
  
  for(k in 1:kstep){
    
    to   <- which(kdex > k)
    from <- which(kdex <= (nbreak-k))
    
    kdist[,to] <- kdist[,to] + dxMat[index,to]
    p0[,to]    <- pnorm(kdist[,to],mmat[,from],sdMat[,from]) - ptot[,to]  #fraction to + k
    gg[,to]    <- gg[,to] + ggam[,from]*p0[,to]*(1 - ttheta[,from])^dtMat[index,from]
    
    ptot[,to] <- ptot[,to] + p0[,to]
  }
  
  gg[tnextG,] <- gg[tnowG,]
  gg[tnextG,fstage] <- gg[tnextG,fstage] + ffec[tnowG,]*dtMat[index,][tnowG,fstage]
  gg
}
  

getGamMatrixRcpp <- function( ffec = fecundMat[wind,],ggam=gam[wind,],
                              ttheta=thetaMat[wind,],gro = groMat[wind,],
                              tnow=tnowG,tnext=tnextG,index=wind, kstep = 2){
  
  # index is in 1:ntnp
  # tnow, tnext is in 1:ntnpG
  
  groP <- gro
  groP[groP < 0] <- 0
  
  dtdx <- dtMat[index,]/dxMat[index,]
  mmat <- groP*dtdx                        #rho*dt/dx
  
  ntnpI <- length(index)
  
  mmat[absTreeMat[index,] == 1] <- 0
  
  #make kernel
  
  sdMat <- sqrt( matrix( diag(sigma), ntnpI, nsnb,byrow=T )*dtdx )
  
  #from same class
  kdist <- dxMat[index,]             
  p0    <- pnorm(kdist,mmat,sdMat)  
  ptot  <- p0
  
  all_powered <- (1 - ttheta)^dtMat[index,]
  gg <- ggam*p0*all_powered
  
  gg <- loop_in_getGamMatrix_cpp(kdex, 0:(length(kdex)-1), kstep, nbreak, kdist, (index-1), dxMat, 
                                 mmat, sdMat, p0, ptot, gg, ggam, all_powered)
  
  gg[tnext,] <- gg[tnow,]
  gg[tnext,fstage] <- gg[tnext,fstage] + ffec[tnow,]*dtMat[index,][tnow,fstage]
  gg
}


getGam <- function(ffec,ggam,ttheta,gro=groMat){
  
  gnow <- ggam*0
  
  gg <- groTrans(tnow,klast,gro)
  gf <- groTrans(tnow,fstage,gro)
  gg[gg > 1] <- 1
  gf[gf > 1] <- 1
  
  gnow[tnext,know] <- ggam[tnow,klast]*(gg - survTrans(tnow,klast,ttheta)) +
                      ggam[tnow,know]*(1 - gg - survTrans(tnow,know,ttheta) ) 
  gnow[tnext,fstage] <- ffec[tnow,]*dtMat[tnow,fstage] + ggam[tnow,fstage]*
                        (1 - gf - survTrans(tnow,fstage,ttheta) )
  
  gnow[ftime,] <-  gnow[ftime+1,]        #temporary: 1st yr set to 2nd yr
  
  gnow
}



  
omegaLambda1 <- function(gg1,gg2=gg1,lam1,lam2=lam1,psi1,psi2=psi1,
                         ovar1,ovar2=ovar1,zeroIndex=absLong,
                         lstageDt=longStageDt,lstageDx=longStageDx,longS=longSpec,longJ=longSpec){
  
  # density for MVN with AR(1) covariance
  # gg1 is ntnpns X nbreak
  # lam1 is ntnp X nspec
  # for subsample longJ <- lstageJ
  
  nr <- nrow(gg1)
  
  p1 <- p2 <- rep(0,nr)
  
  nn <- nr/nspec

  d1 <- 1
  
  if(is.null(zeroIndex))zeroIndex <- rep(0,nrow(gg1))
  
  oinv1 <- invertAR1(nbreak,psi1,1)
  ev1   <- eigen(oinv1, only.values = TRUE)$values
  oinv2 <- invertAR1(nbreak,psi2,1)
  ev2   <- eigen(oinv2, only.values = TRUE)$values
  
 # sindex <- rep(c(1:nspec),each=nn)
  
  l1 <- (as.vector(lam1)*lstageDt/lstageDx)[zeroIndex == 0]
  l2 <- (as.vector(lam2)*lstageDt/lstageDx)[zeroIndex == 0]
  
  g1 <- gg1[zeroIndex == 0,]
  g2 <- gg2[zeroIndex == 0,]
  
  tmp1 <- dmvnormAR1(xx=g1,psi=psi1,const=l1,sinv=oinv1,sig=ovar1,ev=ev1)  #AR1 with scaled variance
  tmp2 <- dmvnormAR1(xx=g2,psi=psi2,const=l2,sinv=oinv2,sig=ovar2,ev=ev2)  
  
  p1[zeroIndex == 0] <- tmp1
  p2[zeroIndex == 0] <- tmp2
  
  p1 <- matrix(p1,nn,nspec)
  p2 <- matrix(p2,nn,nspec)
  
  list(p1 = p1, p2 = p2)
}

###################################################
smoothAR1 <- function(ymat,rho){
   
  wt <- c(rho,1,rho)
  wt <- wt/sum(wt)
  
  k2 <- which(kdex < nbreak & kdex > 1)
  k1 <- k2 - 1
  k3 <- k2 + 1
  
  tmp <- ymat*0
  
  tmp[,k2] <- ymat[,k1]*wt[1] + ymat[,k2]*wt[2] + ymat[,k3]*wt[3]
  tmp[,kdex == 1] <- ymat[,kdex == 1]*(wt[1] + wt[2]) + ymat[,kdex == 2]*wt[3]
  
  tmp
}
    
gamLimits <- function(){
  
  #largest class by plot
  tmp <- y
  tmp[tmp > 0] <- 1
  tmp <- apply(tmp*wideKmat,1,max) + 1
  tmp <- by(tmp,jtdex[,'j'],max,na.rm=T)
  tmp <- tmp[jtdex[,'j']]
  tooBig <- matrix(tmp,ntnp,nsnb)
  tooBig[tooBig < wideKmat] <- 0
  tooBig[tooBig > 0] <- 1
  
  ymiss <- y0[,sdex]             #missing from plot == 1
  ymiss <- ymiss[jtdex[,'j'],]
  
  ga <- y/sampArea/dxMat         #mean intensity
  ga[is.na(ga)] <- 0
  
  ga <- smoothAR1(smoothAR1(smoothAR1(ga,.5),.5),.5)
  
  maxg <- apply(ga,2,max,na.rm=T)
 # mmax <- (.1 + matrix(maxg,ntnp,nsnb,byrow=T))*(1 - exp(-.1*sampArea*dxMat)) # max value ever observed
    
  vplot <- sqrt(1+y)/sampArea/dxMat         #approx variance
  vplot <- sqrt(ga+.01)
  logam <- ga - 2*sqrt(vplot)
  higam <- ga + 2*sqrt(vplot)
 # logam[y > 0] <- 0                         # obs size-species in plot
#  logam[ymiss == 0] <- 0                    # obs species in plot
  logam[ymiss == 1] <- -100                 # not obs species in plot
  
  maxBio <- 20*1000              #maximum biomass in a diameter class kg/ha
  tmp <- maxBio/biomass/dxMat
  higam[higam > tmp] <- tmp[higam > tmp]
#  higam <- higam*tooBig
  
  logam[logam > higam & higam > 0] <- .5*higam[logam > higam & higam > 0]
  logam[logam <= higam & higam <= 0] <- -1      #obs on plot but never for this size-species
  
#  wh <- which(jtdex[,'j'] %in% holdOutPlots)
  
#  if(nHoldOut > 0){
#    higam[wh,] <- matrix(maxg,length(wh),nsnb,byrow=T)
#  }
  higam[higam <= 1e-7] <- 1e-7
  
  # logam[,kdex == 1] <- 0     #all first stage could be recruitment
  
  colnames(logam) <- colnames(higam) <- colnames(y)
  
  ga[ga < logam] <- logam[ga < logam]
  ga[ga > higam] <- higam[ga > higam]
  
  list(gam = ga, logam = logam, higam = higam, tooBig = tooBig)
}

  
  

###################################### 
phiBootstrap <- function(pars){       
  
  bf  <- matrix(pars[-1],(length(pars)-1),1)
  sig <- pars[1]
  
  #  pnow <- pnew <- phi*0
  
  mu   <- fall[ij,]%*%bf                                    #new recruits unknown at time 1
  
  pzero <- pnorm(0,mu,sqrt(sig),lower.tail=T)
  
  mu[mu < 0] <- 0
  
  fm <- fecMat[ij,ss]
  da <- dtMat[ij,fstage[ss]]*sampArea[ij,fstage[ss]]
  
  gmu <- mu
  gmu[fm == 0] <- (pzero + (1 - pzero)*mu)[fm == 0]
  
  pr <- dpois( fm,(recruit[ij,ss]*mu)*da)
  pr[fm == 0] <- (pzero + (1 - pzero)*pr)[fm == 0]
  -sum(log(pr))
}

######################################

updatePhi <- function(){
  
  hi <- hiPhi[wind,]
  lo <- loPhi[wind,]
  
 # whi <- (2*fecMat[wind,] + 1)/(.1 + recruit[wind,])/dtmG[,fstage]/sampAreaG[,fstage]    #upper bound
 # whi[whi > 500] <- 500
 # hi[hi > whi] <- whi[hi > whi]
 # lo[lo > hi]  <- .8*hi[lo > hi]
  
 # phiG[phiG < lo] <- lo[phiG < lo]
 # phiG[phiG > hi] <- hi[phiG > hi]
  
  mu <- byFunctionRcpp(as.vector(xfecBeta[-lastTime,]),fastJ,fastS,
                                     summat=matJSG*0,totmat=matJSG*0 ,MEAN=T)
  xx <- byFunctionRcpp(as.vector(phiG[-firstTime,]),fastJ,fastS,
                       summat=matJSG*0,totmat=matJSG*0 ,MEAN=T)
  
  rjs <- exp( pmvnormCondRcpp(q=0,xx,mu,smat=fsigma,whichVar=c(1:nspec)) )  #not reproductive
  
  ww <- which(!is.finite(rjs),arr.ind=T)  #under/overflow
  if(length(ww) > 0){
    mw <- mu[ww]
    xw <- xx[ww]
    rjs[ ww[mw < xw,] ] <- 0   #hi phi means reproductive
    rjs[ ww[mw > xw,] ] <- 1   #lo phi means not
  }
  pjs <- 1 - (1 - babsentG)*(1 - rjs)       # either no trees or not reproductive
  
  
  
  da <- dtMat[wind,][-lastTime,fstage]*sampAreaG[-lastTime,fstage]
  ii <- (recruit[wind,]*xfecBeta + Imm[wind,])[-lastTime,]
  us <- exp( -ii*da)
  u  <- byFunctionRcpp(as.vector(us),fastJ,fastS,
                       summat=matJSG*0,totmat=matJSG*0 ,MEAN=T)  #Pr not obs given it can occur
  
  pfabsent <- log(pjs) - log( (pjs + (1 - pjs)*u) )
  
  if(nHolds > 0)pfabsent[holds,] <- log(pjs[holds,])
  pfabsent[is.na(pfabsent)] <- 0
  pfabsent <- exp(pfabsent)
  
  fl <- matrix( rbinom(msample*nspec,1,pfabsent), msample,nspec)   # np sample plots
  fl[f0[msamp,] == 0] <- 0
  fl[knownAbsent[msamp,] == 1] <- 1
  bnoRecG <- fl[jtind[,'j'],]                            # plot and yr 
  
  noPhi <- matrix( rbinom(msample*nspec,1,rjs), msample,nspec)
  noPhi[f0[msamp,] == 0] <- 0
#  noPhi[fl == 0] <- 0

  nojt <- noPhi[jtind[,'j'],]
  
  hi[nojt == 1] <- 0
  lo[nojt == 0] <- 0
  lo[nojt == 1] <- -500
  
  phiG[nojt == 0 & phiG < 0] <- lastPhiPos[wind,][nojt == 0 & phiG < 0]
  phiG[nojt == 1 & phiG > 0] <- lastPhiNeg[wind,][nojt == 1 & phiG > 0]
  
  
  psd <- abs(phiG)/10 + .01
  isd <- Imm[wind,]/20 + .01
  
  phiNew <- matrix( tnorm(ntnpnsG,as.vector(lo),as.vector(hi),
                          as.vector(phiG),psd), ntnpG,nspec)
  immNew <- matrix( tnorm(ntnpnsG,as.vector(loI[wind,]),as.vector(hiI[wind,]),
                          as.vector(Imm[wind,]),isd), ntnpG,nspec)
  
  pnow <- pnew <- phiG*0
  
  mu   <- xfecBeta[-lastTime,]                            #new recruits unknown at time 1
  pnow[-firstTime,] <- dmvnormZeroMean(phiG[-firstTime,] - mu ,smat=fsigma)
  pnew[-firstTime,] <- dmvnormZeroMean(phiNew[-firstTime,] - mu  ,smat=fsigma)
  
  prReg <- pnew - pnow ###
  
  posNow <- phiG
  posNew <- phiNew
  posNow[posNow < 0] <- 0
  posNew[posNew < 0] <- 0
  
  fm <- fecMat[wind,]
  da <- dtmG[,fstage]*sampAreaG[,fstage]   #[,fstage]
  
  frnew <- recruit[wind,]*posNew + immNew
  frnow <- recruit[wind,]*posNow + Imm[wind,]

  prData <- poissonRatio(fm,frnew*da,frnow*da)
  
  if(nholds > 0)prData[holdIndex,] <- 0
  
  prReg[-firstTime,] <- prReg[-firstTime,] + prData[-firstTime,]
  
  tmp <- getGamMatrixSample(f1=frnow,f2=frnew,
                            gam1=wgamG,gam2=wgamG,th1=thetaG,th2=thetaG,
                            gr1=groPos[wind,],gr2=groPos[wind,],index=wind, 
                            ntnpK=ntnpG,time1=tnowG,time2=tnextG,first=firstTime,
                            ktimes=4)
  predGamG <- tmp[[1]]
  predNew  <- tmp[[2]]
  
  prGam <- gaussianRatio(predNew[,fstage],predGamG[,fstage],wgamG[,fstage],gamVar[,fstage])
  
  a <- exp( prReg + prGam )
  
  z <- runif(ntnpnsG,0,1)
  ww <- which(z < a)
  phiG[ww] <- phiNew[ww]
  Imm[wind,][ww] <- immNew[ww]
  predGamG[,fstage][ww] <- predNew[,fstage][ww]
  
  list(phi = phiG, Imm = Imm, pfabsent = pfabsent, bnoRecPlot = fl, 
       noPhi = noPhi, bnoRecPlotYr = bnoRecG,predNow = predGamG,accept=length(ww)/ntnpnsG)
}


updateGammaOld <- function(){
  
  fecG <- fecundMat[wind,]
  
  tiny <- 1e-8
  
  wgamG[absMatG == 0 & wgamG < 0] <- lastGamPos[wind,][absMatG == 0 & wgamG < 0]
  wgamG[absMatG == 1 & wgamG > 0] <- lastGamNeg[wind,][absMatG == 1 & wgamG > 0]
  
  gamG[gamG > 0]  <- wgamG[gamG > 0]
  gamG[wgamG < 0] <- 0
  
  psd <- as.vector( abs(wgamG)*.01 ) + .001
  
  loG <- logam[wind,]
  hiG <- higam[wind,]
  hiG[absMatG == 1] <- 0
  loG[absMatG == 1] <- -50
  
  wgamG[wgamG > hiG] <- hiG[wgamG > hiG]
  wgamG[wgamG < loG] <- loG[wgamG < loG]
  
  gamG <- wgamG
  gamG[gamG < 0] <- 0
  
  lo <- as.vector(loG)
  hi <- as.vector(hiG)
  
  #proposal
  pam  <- matrix( tnorm(ntnpnsG*nbreak,lo,hi,as.vector(wgamG),psd),ntnpG,nsnb)  
  
  gamNew <- pam
  gamNew[gamNew < 0] <- 0
  
  gamG[absMatG == 1] <- gamNew[absMatG == 1] <- 0
  fecG[absMatG[,kdex == 1] == 1] <- 0
  
  #observations from SSD
  pgnow <-  dpois(y[wind,],as.matrix(sampByWidth[wind,]*gamG + 1e-10),log=T) #- log(tmp1)
  pgnew <-  dpois(y[wind,],as.matrix(sampByWidth[wind,]*gamNew + 1e-10),log=T) #- log(tmp2)
  
  pgnow[absMatG == 1] <- pgnew[absMatG == 1] <- 0
  
  if(nHolds > 0){
    pgnow[jtind[,'j'] %in% holds,] <- pgnew[jtind[,'j'] %in% holds,] <- 0
  }
  
  prSeed <- byFunctionRcpp(as.vector(pgnew - pgnow),wideJ,wideS*0+1,
                           matrix(0,msample,1),matrix(0,msample,1),MEAN=F)

  #predict data
  predGamG <- matrix( rpois(ntnpnsG*nbreak,sampByWidth[wind,]*gamG),ntnpG,nsnb)
  

  
  #total density
#  wmeanNew  <-byFunctionRcpp( as.vector(pam[,wstage]*dxMat[wind,wstage]),
#                              as.vector(wideJG[,wstage]),
#                              as.vector(wideSG[,wstage]),
#                              matJS*0,matJS*0,MEAN=T)
  
  # from SDM
  # prReg <- dmvnormZeroMean(wmeanNew - xgamBeta,gsigma) - 
  #          dmvnormZeroMean(wmean - xgamBeta,gsigma)
  
  prReg <- rep(0,msample)
  
  
  if(nHolds > 0)pgnow[holds,] <- pgnew[holds,] <- 0
  
  if(PDE){
    predNow <- getGam( fecG,gamG,thetaG,groPos[wind,] )
    predNew <- getGam( fecG,gamNew,thetaG,groPos[wind,] )
  }
  if(!PDE){
    predNow <- getGamMatrixRcpp( fecG,gamG,thetaG,groPos[wind,], index=wind)
    predNew <- getGamMatrixRcpp( fecG,gamNew,thetaG,groPos[wind,], index=wind )
  }
  
  predNow[tnextG,] <- predNow[tnextG,] - gamG[tnextG,]
  predNew[tnextG,] <- predNew[tnextG,] - gamNew[tnextG,]
  predNow[firstTime,] <- predNew[firstTime,] <-0
  
  ggnow <- short2long(predNow,'STAGE',n1=ntnpG,n2=ntnpnbG,n3=ntnpnsG)    #check n2, n3
  ggnew <- short2long(predNew,'STAGE',n1=ntnpG,n2=ntnpnbG,n3=ntnpnsG)

  tmp <- omegaLambda1(gg1=ggnow,gg2=ggnew,lam1=lvarG,psi1=psi,ovar1=1,
                      zeroIndex=absLongG,lstageDt=lstageGDt,lstageDx=lstageGDx,longJ=msamp[lstageJ])
  pnow <- tmp$p1
  pnew <- tmp$p2
  
  pnow[is.na(pnow)] <- 0
  pnew[is.na(pnew)] <- 0
  
  mm <- as.vector(wideJG[-firstTime,fstage])
  prGamma <- byFunctionRcpp(as.vector(pnew[-firstTime,] - pnow[-firstTime,]),mm,mm*0+1,
                            matrix(0,msample,1), matrix(0,msample,1),MEAN=F)
  
  aa <- exp(prSeed + prReg + prGamma)
  z  <- runif(msample,0,1)
  wa <- which(z < aa )
  
  wz <- which(wideJG %in% wa,arr.ind=T)
  
  gamG[wz]  <- gamNew[wz]
  wgamG[wz] <- pam[wz]
  predNow[wz] <- predNew[wz]
  
  list(gam = gamG, wgam = wgamG,mh = aa, predGam = predGamG, predModel = predNow,accept=length(wa)/msample)
}

getGamMatrixTwiceOld <- function(f1=fecG,f2=NULL,gam1=gamG,gam2=NULL,
                              th1=thetaG,th2=NULL,gr1=groPos[wind,],gr2=NULL, 
                              index=wind, 
                              ntnpK=ntnpG,ktimes=3){
  
  tmp2 <- numeric(0)
  
  xk   <- breakMat[index,]                          
  dtm  <- dtMat[index,]
  dxm  <- dxMat[index,]
  sigk <- sqrt(diag(sigma) )
  bk   <- matrix(sigk[sdex],ntnpK,nsnb,byrow=T)  # s.d. matrix
  
  gk1   <- gr1*dtMat[index,]                      # interval growth
  gg1   <- gam1*(1 - th1)^dtm
  mk1   <- xk + gk1                                 # kernel mean on x scale
  tmp1  <- gg1*dxm*dnorm(xk,mk1,bk)*dxm             # received from current stage
  
  if(!is.null(f2)){
    gk2   <- gr2*dtMat[index,]                      # interval growth
    gg2   <- gam2*(1 - th2)^dtm
    mk2   <- xk + gk2                               # kernel mean
    tmp2  <- gg2*dxm*dnorm(xk,mk2,bk)*dxm           # received from current stage
  }
  
  for(k in 1:(ktimes-1)){
    kfrom <- which(kdex <= (nbreak-k))  # contribution to k stages ahead
    kto   <- which(kdex > k)            # received k stages ahead
    
    kern       <- dnorm(xk[,kto],mk1[,kfrom],bk[,kfrom])*dxm[,kto]
    tmp1[,kto] <- tmp1[,kto] + gg1[,kfrom]*dxm[,kfrom]*kern
    if(!is.null(f2)){
      kern       <- dnorm(xk[,kto],mk2[,kfrom],bk[,kfrom])*dxm[,kto]
      tmp2[,kto] <- tmp2[,kto] + gg2[,kfrom]*dxm[,kfrom]*kern
    }
  }
  
  kfrom <- which(kdex <= (nbreak-ktimes))  # contribution to k stages ahead
  kto   <- which(kdex > ktimes)            # received k stages ahead
  
  kern       <- 1 - pnorm(xk[,kto],mk1[,kfrom],bk[,kfrom])
  tmp1[,kto] <- tmp1[,kto] + gg1[,kfrom]*dxm[,kfrom]*kern
  tmp1       <- tmp1/dxm     #from ha-1 to cm-1 ha-1
  tmp1[,fstage] <- tmp1[,fstage] + f1*dtm[,fstage]
  
  if(!is.null(f2)){
    kern       <- 1 - pnorm(xk[,kto],mk2[,kfrom],bk[,kfrom])
    tmp2[,kto] <- tmp2[,kto] + gg2[,kfrom]*dxm[,kfrom]*kern
    tmp2       <- tmp2/dxm     #from ha-1 to cm-1 ha-1
    tmp2[,fstage] <- tmp2[,fstage] + f2*dtm[,fstage]
  }
  
  list(predgam1 = tmp1, predgam2 = tmp2)
}
  


getGamMatrixOld <- function(feck=fecG,gamk=gamG,thek=thetaG,grok=groPos[wind,], 
                         index=wind, 
                         ntnpK=ntnpG,ktimes=3){
  
  tiny <- 1e-8
  xk   <- breakMat[index,]                          
  dtm  <- dtMat[index,]
  dxm  <- dxMat[index,]
  sigk <- sqrt(diag(sigma) )
  bk   <- matrix(sigk[sdex],ntnpK,nsnb,byrow=T)  # s.d. matrix
  
  gk   <- grok*dtMat[index,]                      # interval growth
  gg   <- gamk*(1 - thek)^dtm
  
  mk   <- xk - dxm/2 + gk                                 # kernel mean
  tmp  <- gg*dxm*dnorm(xk,mk,bk)*dxm               # received from current stage
  
  for(k in 1:(ktimes-1)){
    kfrom <- which(kdex <= (nbreak-k))  # contribution to k stages ahead
    kto   <- which(kdex > k)            # received k stages ahead
    kern  <- dnorm(xk[,kto],mk[,kfrom],bk[,kfrom])*dxm[,kto]
    tmp[,kto] <- tmp[,kto] + gg[,kfrom]*dxm[,kfrom]*kern
  }
  
  tmp <- tmp/dxm     #from ha-1 to cm-1 ha-1
  
  tmp[,fstage] <- tmp[,fstage] + feck*dtm[,fstage]
  tmp
}


getGamMatrixSampleOld <- function(f1=fecG,f2=NULL,gam1=gamG,gam2=NULL,
                              th1=thetaG,th2=NULL,gr1=groPos[wind,],gr2=NULL, 
                              index=wind, 
                              ntnpK=ntnpG,time1=tnowG,time2=tnextG,first=firstTime,
                              ktimes=4){

  xk    <- breakMat[index,]                          
  dtm   <- dtMat[index,]
  dxm   <- dxMat[index,]
  imat1 <-  xk*0
  tmp1  <-  xk*0
  
  mk1   <- xk - dxm/2 + gr1*dtMat[index,]   # kernel mean on x scale
  gx1   <- dxm*gam1*(1 - th1)^dtm
  
  ww    <- which(mk1 < xk)
  tmp1[ww]   <- gx1[ww]                     #stay in same class
  imat1[ww]  <- 1
  
  if(!is.null(f2)){
    imat2  <- xk*0
    tmp2   <- xk*0

    mk2   <- xk - dxm/2 + gr2*dtMat[index,]  
    gx2   <- dxm*gam2*(1 - th2)^dtm
    
    ww    <- which(mk2 < xk)
    tmp2[ww] <- gx2[ww]
    imat2[ww] <- 1
  }
  
  kto <- c(1:(ntnpG*nsnb))
  k1  <- c(1:ntnpG)
  
  for(k in 1:ktimes){
    
    kfrom <- c( 1:(ntnpG*(nsnb-k)) )
    kto   <- kto[-k1]
              
    ww    <- which(mk1[kfrom] < xk[kto] & imat1[kfrom] == 0)
    if(length(ww) > 0){
      tmp1[kto[ww]] <- tmp1[kto[ww]] + gx1[kfrom[ww]]
      imat1[kfrom[ww]]  <- 1
    }
    
    if(!is.null(f2)){
      ww    <- which(mk2[kfrom] < xk[kto] & imat2[kfrom] == 0)
      if(length(ww) > 0){
        tmp2[kto[ww]] <- tmp2[kto[ww]] + gx2[kfrom[ww]]
        imat2[kfrom[ww]]  <- 1
      }
    }
  }
  
  tmp1          <- tmp1/dxm     #from ha-1 to cm-1 ha-1
  tmp1[,fstage] <- tmp1[,fstage] + f1*dtm[,fstage]
  tmp1[time2,]  <- tmp1[time1,]
  tmp1[first,]  <- gam1[first,]
  
  if(!is.null(f2)){
    tmp2          <- tmp2/dxm     #from ha-1 to cm-1 ha-1
    tmp2[,fstage] <- tmp2[,fstage] + f2*dtm[,fstage]
    tmp2[time2,]  <- tmp2[time1,]
    tmp2[first,]  <- gam2[first,]
  }
  
  list(predgam1 = tmp1, predgam2 = tmp2)
}

getGamMatrixSample <- function(f1=fecG,f2=NULL,gam1=gamG,gam2=NULL,
                               th1=thetaG,th2=NULL,gr1=groPos[wind,],gr2=NULL, 
                               index=wind, 
                               ntnpK=ntnpG,time1=tnowG,time2=tnextG,first=firstTime,
                               ktimes=4){
  
  xk    <- breakMat[index[time1],]                          
  dtm   <- dtMat[index[time1],]
  dxm   <- dxMat[index[time1],]
  imat1 <-  xk*0
  tmp1  <-  xk*0
  nktot <- nrow(xk)
  
  mk1   <- xk - dxm/2 + gr1[time1,]*dtm   # kernel mean on x scale
  gx1   <- dxm*gam1[time1,]*(1 - th1[time1,])^dtm
  
  ww    <- which(mk1 < xk)
  tmp1[ww]   <- gx1[ww]                     #stay in same class
  imat1[ww]  <- 1
  
  if(!is.null(f2)){
    imat2  <- xk*0
    tmp2   <- xk*0
    
    mk2   <- xk - dxm/2 + gr2[time1,]*dtm  
    gx2   <- dxm*gam2[time1,]*(1 - th2[time1,])^dtm
    
    ww    <- which(mk2 < xk)
    tmp2[ww] <- gx2[ww]
    imat2[ww] <- 1
  }
  
  kto <- c(1:(nktot*nsnb))
  k1  <- c(1:nktot)
  
  for(k in 1:ktimes){
    
    kfrom <- c( 1:(nktot*(nsnb-k)) )
    kto   <- kto[-k1]
    
    ww    <- which(mk1[kfrom] < xk[kto] & imat1[kfrom] == 0)
    if(length(ww) > 0){
      tmp1[kto[ww]] <- tmp1[kto[ww]] + gx1[kfrom[ww]]
      imat1[kfrom[ww]]  <- 1
    }
    
    if(!is.null(f2)){
      ww    <- which(mk2[kfrom] < xk[kto] & imat2[kfrom] == 0)
      if(length(ww) > 0){
        tmp2[kto[ww]] <- tmp2[kto[ww]] + gx2[kfrom[ww]]
        imat2[kfrom[ww]]  <- 1
      }
    }
  }
  
  tmp1          <- tmp1/dxm     #from ha-1 to cm-1 ha-1
  tmp1[,fstage] <- tmp1[,fstage] + f1[time1,]*dtm[,fstage]*(1 - th1[time1,fstage])^dtm[,fstage]
  gam1[time2,]  <- tmp1
  
  if(!is.null(f2)){
    tmp2          <- tmp2/dxm     #from ha-1 to cm-1 ha-1
    tmp2[,fstage] <- tmp2[,fstage] + f2[time1,]*dtm[,fstage]*(1 - th1[time1,fstage])^dtm[,fstage]
    gam2[time2,]  <- tmp2
  }
  
  list(predgam1 = gam1, predgam2 = gam2)
}



getGamMatrixTwice <- function(f1=fecG,f2=NULL,gam1=gamG,gam2=NULL,
                              th1=thetaG,th2=NULL,gr1=groPos[wind,],gr2=NULL, 
                              index=wind, 
                              ntnpK=ntnpG,ktimes=3){
  tiny <- 1e-8
  tmp2 <- numeric(0)
  
  xk   <- breakMat[index,]                          
  dtm  <- dtMat[index,]
  dxm  <- dxMat[index,]
  sigk <- sqrt(diag(sigma) )
  bk   <- matrix(sigk[sdex],ntnpK,nsnb,byrow=T)  # s.d. matrix
  
  gk1   <- gr1*dtMat[index,]                      # interval growth
  gg1   <- gam1*(1 - th1)^dtm
  mk1   <- xk - dxm/2 + gk1                                 # kernel mean on x scale
 # tmp1  <- gg1*dxm*dnorm(xk,mk1,bk)*dxm             # received from current stage
  
  plast1 <- pcum1 <- pnorm(xk,mk1,bk)
  tmp1   <- gg1*dxm*plast1
  
  if(!is.null(f2)){
    gk2   <- gr2*dtMat[index,]                      # interval growth
    gg2   <- gam2*(1 - th2)^dtm
    mk2   <- xk - dxm/2 + gk2                        # kernel mean
 #   tmp2  <- gg2*dxm*pnorm(xk,mk2,bk)*dxm           # received from current stage
    plast2 <- pcum2 <- pnorm(xk,mk2,bk)
    tmp2   <- gg2*dxm*plast2
  }
  
  for(k in 1:ktimes){
    
    kfrom <- which(kdex <= (nbreak-k))  # contribution to k stages ahead
    kto   <- which(kdex > k)            # received k stages ahead
    
    ww1 <- which(plast1[,kfrom] < 1 - tiny,arr.ind=T)
    
    if(k < ktimes)pcum1[,kfrom][ww1]  <- pnorm(xk[,kto][ww1],mk1[,kfrom][ww1],bk[,kfrom][ww1])
    if(k == ktimes)pcum1 <- pcum1*0 + 1
    
    kern   <- pcum1[,kfrom] - plast1[,kfrom]
    plast1 <- pcum1
    tmp1[,kto] <- tmp1[,kto] + gg1[,kfrom]*dxm[,kfrom]*kern
    
    if(!is.null(f2)){
      
      ww2 <- which(plast2[,kfrom] < 1 - tiny,arr.ind=T)
      
      if(k < ktimes)pcum2[,kfrom][ww2]  <- pnorm(xk[,kto][ww2],mk2[,kfrom][ww2],bk[,kfrom][ww2])
      if(k == ktimes)pcum2 <- pcum2*0 + 1
      
      kern   <- pcum2[,kfrom] - plast2[,kfrom]
      plast2 <- pcum2
      tmp2[,kto] <- tmp2[,kto] + gg2[,kfrom]*dxm[,kfrom]*kern
    }
  }
  
  tmp1          <- tmp1/dxm     #from ha-1 to cm-1 ha-1
  tmp1[,fstage] <- tmp1[,fstage] + f1*dtm[,fstage]
  
  if(!is.null(f2)){
  
    tmp2          <- tmp2/dxm     #from ha-1 to cm-1 ha-1
    tmp2[,fstage] <- tmp2[,fstage] + f2*dtm[,fstage]
  }
  
  list(predgam1 = tmp1, predgam2 = tmp2)
}

updateGamma <- function(){
  
  tiny <- 1e-8
  
  loG <- logam[wind,]
  hiG <- higam[wind,]
  hiG[absMatG == 1] <- 0
  loG[absMatG == 1] <- -50
  loG[absMatG == 0 & loG < 0] <- 0
  
  lo <- as.vector(loG)
  hi <- as.vector(hiG)
  
  #proposal
  rg  <- rgamma(ntnpnbG*nspec,5,10)
 # psd <- as.vector( abs(wgamG)*.2 ) + rg
  
  psd <- as.vector(abs(wgamG)/(wideKmat[wind,] + rg))
  pam  <- matrix( tnorm(ntnpnsG*nbreak,lo,hi,as.vector(wgamG),psd),ntnpG,nsnb)*tooBig[wind,]  
  
#  wgamG[,kdex == nbreak] <- pam[,kdex == nbreak] <- 0
  
  gamNew <- pam
  gamNew[gamNew < 0] <- 0
  
 # gamG[absMatG == 1] <- gamNew[absMatG == 1] <- 0
  fecG[absMatG[,kdex == 1] == 1] <- 0
  
  #observations from SSD
  sm <- sampByWidth[wind,]
  
  prData <- poissonRatio(y[wind,],sm*gamNew + tiny,sm*gamG + tiny)
  if(nHolds > 0)prData[holdIndex,] <- 0         # no data on holdouts
  
  #predict data
  
  tmp <- getGamMatrixSample(f1=fecG,f2=fecG,gam1=wgamG,gam2=pam,
                            th1=thetaG,th2=thetaG,
                            gr1=groPos[wind,],gr2=groPos[wind,], 
                            index=wind, 
                            ntnpK=ntnpG,time1=tnowG,time2=tnextG,first=firstTime,
                            ktimes=4)
  predNow <- tmp[[1]]
  predNew <- tmp[[2]]
  
  ppred  <- predNow
  ppred[ppred < 0] <- 0
  
  #predict data one step ahead
  predYG <- matrix( rpois(ntnpnsG*nbreak,sm*ppred),ntnpG,nsnb)
  
  resNow <- predNow - wgamG
  resNew <- predNew - pam
  
  prGam <- gaussianRatio(resNew,resNow,0,gamVar)
  
  prGD <- byFunctionRcpp(as.vector(prGam + prData),as.vector(wideJG),as.vector(wideSG),
                        matJSG,matJSG,MEAN=F)
  
  #total plot basal area
  
  bNow <- gamG*biomass[wind,]*dxmG/1000
  bmat <- rowSums(bNow)
  bNow <- byFunctionRcpp(bmat,jtind[,'j'],rep(1,ntnpG),
                         matrix(0,msample,1),matrix(0,msample,1),MEAN=T)
  
  bNew <- gamNew*biomass[wind,]*dxmG/1000
  bmat <- rowSums(bNew)
  bNew <- byFunctionRcpp(bmat,jtind[,'j'],rep(1,ntnpG),
                         matrix(0,msample,1),matrix(0,msample,1),MEAN=T)
  
 # print(range(bNow))
  
  prB <- gaussianRatio(bNew,bNow,ecoBio[msamp,1],ecoBio[msamp,2]/ntot)
  
  mb <- 3*ecoBio[msamp,1]     
  lb <- .3*ecoBio[msamp,1] 
  
  ww <- which( bNew > mb | bNew < lb)
  if(length(ww) > 0)prB[ww] <- 1e-10
  
  
  tmp <- rowSums(prGD) + prB 
  
  predict <- predNow
 # predict[predict < 0] <- 0
  
  aa <- exp(tmp)
  z  <- runif(length(aa),0,1)
  aa[aa < z] <- 0
  aa[aa > 0] <- 1
  
 # aa <- aa[jtind[,'j'],]
 # aa <- aa[,sdex]
  
  wa <- which(aa == 1)
  wj <- which(jtind[,'j'] %in% wa)
  
  gamG[aa == 1,]    <- gamNew[aa == 1,]
  wgamG[aa == 1,]   <- pam[aa == 1,]
  predNow[aa == 1,] <- predNew[aa == 1,]
  
  predNow[predNow < 0] <- 0
  
  list(gam = gamG, wgam = wgamG,predGam = predict, predY = predYG,
       accept=sum(aa)/length(aa),bmass = bNow)
}


updateOmega <- function(){
  
  absTmp  <- babsentMat                         #absent and firstTime excluded
  absTmp[ftime,] <- 1
  wkind   <- which(absTmp == 0)
  
  kindex <- wideKmat[wkind]
  jindex <- wideJmat[wkind]
  sindex <- wideSmat[wkind]
  
  res <- wgam*0
  
  res <- (predGamModel - gam)^2
  
  SJL <- sum(1 - absTmp[,fstage])  #all stages present on same plots

  #sample tau
  
  kmat <- kappa[jtdex[,'j'],sdex]
  rr   <- res/kmat/dtMat
  
  resTime <- byFunctionRcpp(as.vector(rr[wkind]),rep(1,length(wkind)),
                            kindex,
                            matrix(0,1,nbreak),matrix(0,1,nbreak),MEAN=F)
  
  muk <- 100*exp(-.05*(1:nbreak))
  s1  <- SJL
  s2  <- muk*(s1 - 1)
  p1  <- s1 + SJL/2
  p2  <- s2 + 1/2/dxxdtt*resTime
  
  hi <- muk*3 + 1
  
  tau <- rtrunc(nn=nbreak,lo=.001,hi=hi,p1=rep(p1,nbreak),p2=p2,'igamma')
  tau[tau > hi] <- hi[tau > hi]
  tau[!is.finite(tau)] <- muk[!is.finite(tau)]
  
  tmp  <- (1 - absTmp)[wkind]
  KL   <- byFunctionRcpp(tmp,jindex,sindex,
                         matrix(0,nplot,nspec),
                         matrix(0,nplot,nspec),MEAN=F)
  
  kmat <- matrix(tau,ntnp,nsnb,byrow=T)
  
  rr <- res[wkind]/kmat[wkind]/dtMat[wkind]
  
  resTime <- byFunctionRcpp(rr,jindex,sindex,matrix(0,nplot,nspec),
                            matrix(0,nplot,nspec),MEAN=F)
  
  muk <- 1
  s1  <- 10 + KL
  s2  <- 1 + muk*(s1 - 1)
  p1  <- s1 + KL/2
  p2  <- s2 + 1/2/dxxdtt*resTime
  
  kappa <- matrix( 1/rgamma(nplot*nspec,p1,p2), nplot,nspec)
  kappa[,nspec] <- 1          #reference class
  kappa[kappa > 20] <- 20     #max value
  
  list(tau = tau, kappa = kappa)
}
  
  
  

updateOmegaQuick <- function(kexclude=NULL){
  
  #kexclude are size classes to omit from covariance
  
  mu  <- byFunctionRcpp(as.vector(wgam),as.vector(wideJTmat),
                        as.vector(wideSmat),matSpec*0,matSpec*0,MEAN=T) 
  
  mu <- rowMeans(wgam)
  tmp <- sum( (wgam[-ftime,] - mu[-ftime])*(wgam[-ltime,] - mu[-ftime]) )/sum( (wgam[-ftime,] - mu[-ftime])^2 )
  
 # predNow <- getGamMatrixRcpp(ffec=fecundMat,ggam=gam,
 #                                     ttheta=thetaMat,gro=groPos,index=c(1:ntnp))
  predNow <- predGamModel
  predNow[tnext,] <- predNow[tnext,] - wgam[tnext,]
  
  mu  <- rowMeans(predNow)
  rho <- sum( (predNow[-ftime,] - mu[-ftime])*(predNow[-ltime,] - mu[-ltime]) )/
    sum( (predNow[-ltime,] - mu[-ltime])^2 )
  
  if(rho > 0)rho <- 0
  
  muy <- predNow
  if(rho != 0)muy[-ltime,] <- rho*(wgam[-ltime,] - predNow[-ltime,]) + predNow[-ftime,]
  
  muy[babsentMat == 1] <- 0
  
  nn <- 1 - babsentMat
  
  if(!is.null(kexclude)){
    muy[wideKmat %in% kexclude] <- 0
    nn[wideKmat %in% kexclude] <- 0
  }
  
  my <- byFunctionRcpp(as.vector(muy^2),as.vector(wideJmat),
                       as.vector(wideSmat),matrix(0,nplot,nspec),matrix(0,nplot,nspec),MEAN=F)
  nn <- byFunctionRcpp(as.vector(nn),as.vector(wideJmat),
                       as.vector(wideSmat),matrix(0,nplot,nspec),matrix(0,nplot,nspec),MEAN=F)
  
  u1 <- 1 + nn/2
  u2 <- 1 + .5*my
  
  lvar <- matrix( rtrunc(nn=nplot*nspec,lo=.1,hi=5000,p1=as.vector(u1),p2=as.vector(u2),'igamma'), nplot,nspec)
  lvar <- lvar*dxMat[ftime,kdex==1]/dtMat[ftime,kdex==1]
  
  list(rho = rho, lvar = lvar)
}

  
updateOmegaFull <- function(){
#sample Omega
  
  if(PDE) predNow <- getGam(fecG,gamG,thetaG,gro=groPos[wind,])
  if(!PDE)predNow <- getGamMatrixRcpp(ffec=fecundMat,ggam=gam,
                                      ttheta=thetaMat,gro=groPos,index=c(1:ntnp))
  predNow[tnext,] <- predNow[tnext,] - wgam[tnext,]
  
  ggnow <- short2long(predNow,'STAGE',n1=ntnp,n2=ntnpnb,n3=ntnpns)
  
  lamProp <- lambda*0 + tnorm(nplot*nspec,0,500,lambda,.001)
  if(rbinom(1,0,.3) == 1)lamProp <- lambda
  lambda[,'other'] <- lamProp[,'other'] <- 1
  
  psiProp <- tnorm(1,loPsi,hiPsi,psi,.1)
  
  tmp   <- tnorm.mvtRcpp(c(ovar,psi), c(ovar,psi), oProp, c(0,loPsi), c(ohi,hiPsi))  
  ppsi  <- tmp[2]
  povar <- tmp[1]
  
  
  tmp <- omegaLambda1(gg1=ggnow,lam1=lambda[jtdex[,'j'],],lam2=lamProp[jtdex[,'j'],],
                      psi1=psi,psi2=ppsi,
                      ovar1=ovar,ovar2=povar,zeroIndex=absLong)
  pnow <- tmp$p1
  pnew <- tmp$p2
  
  if(nHolds > 0)pnow[jtdex[,'j'] %in% holds,] <- pnew[jtdex[,'j'] %in% holds,] <- 0
  
  pnow[absMat[,fstage] == 1] <- pnew[absMat[,fstage] == 1] <- 0
  
  a <- exp( sum(pnew,na.rm=T) - sum(pnow,na.rm=T) )
  z <- runif(1,0,1)
  if(z < a){
    ovar <- povar
    psi  <- ppsi
    lambda <- lamProp
  }

list(ovar = ovar, psi = psi,lambda = lambda)
}
  
short2long <- function(xx,dm='SPEC',n1=ntnp,n2=ntnpnb,n3=ntnpns){  
  #from (nt*nplot,nsnb) to (nt*nplot*nspec_or_nbreak,nspec_or_nbreak)
  
  stack3d( make3d(xx,n1=n1),dm,n1=n1,n2=n2,n3=n3 )
}

updateBetaMort <- function(xx,yy,bindex){
  
  mvar <- rep(0,nspec)
  mm   <- length(mnames)
  
  for(s in 1:nspec){
    
    wc <- c(1:mm)
    ws <- which(bindex[,s] == 1)
    XS <- xx[ws,]
    if( qr(XS)$rank < ncol(XS) ){    
      lose <- numeric(0)
      if( sum(XS[,'xeric']) == 0 )lose <- c(lose,grep('xeric',mnames))
      if( sum(XS[,'mesic']) == 0 )lose <- c(lose,grep('mesic',mnames))
      if( length(lose) == 0 )lose <- intPriorM   #if not full rank interactions are zero
      wc <- wc[-lose]
      XS <- XS[,wc]
      betaMort[lose,s] <- 0
    }  
      
    Y  <- yy[ws,s]
    
    V  <- msigma[j,j]*solve(crossprod(XS))
    v  <- 1/msigma[j,j]*crossprod(XS,Y)
    betaMort[wc,s] <- tnorm.mvtRcpp(betaMort[wc,s], V%*%v, V, loM[,s], hiM[,s])
    
    mvar[s] <- updateVariance(Y,XS%*%betaMort[wc,s],s1=mprior[1],s2=mprior[2],hi=10) 
  }
  
  list(betaMort = betaMort, mvar = mvar)
}

  
  
updateThetaPars <- function(){
  
  thetaMat[thetaMat < loth] <- loth[thetaMat < loth]         #change later?
  thetaMat[thetaMat > hith] <- hith[thetaMat > hith]
  
  tiny <- 1e-5
  thetaMat[thetaMat <= tiny] <- tiny
  
  logitNow <- logitNew <- plogit <- gall*0
  predNow  <- predNew <- gall*0
  
  logitTheta <- short2long( logit(thetaMat),n1=ntnpG,n2=ntnpnbG,n3=ntnpnsG )
  amat       <- short2long( absMat)
  
  pbeta    <- matrix( tnorm.mvtRcpp(as.vector(betaMort), as.vector(betaMort), 
                                    mProp, loM, hiM),pm,nspec)
  logitNow <- xall[,mnames]%*%betaMort
  logitNew <- xall[,mnames]%*%pbeta
  
  pnow   <- dnorm(logitTheta,logitNow,sdMat,log=T)
  pnew   <- dnorm(logitTheta,logitNew,sdMat,log=T)
  
  pnow[mortRegIndex,] <- pnew[mortRegIndex,] <- 0
  pnow[amat == 0] <- pnew[amat == 0] <- 0
  
  aa <- exp(colSums(pnew) - colSums(pnow))
  ww <- which( runif(nspec,0,1) < aa )  
  betaMort[,ww] <- pbeta[,ww]
  
  logitNow <- xall[,mnames]%*%betaMort
  res <- (logitTheta - logitNow)[mortRegIndex,]
  
  u1 <- mprior1 + nrow(res)/2
  u2 <- mprior2 + colSums( res^2 )/2
  varTheta <- 1/rgamma(nspec,u1,u2)
  
  list(beta = betaMort, varTheta = varTheta)
}

sampleTheta <- function(){
  
#  lo <- loth[wind,]
#  hi <- hith[wind,]
  
  lo <- .002
  hi <- hiTheta[wind,]
  
#  thetaG[thetaG < lo] <- lo[thetaG < lo]         #change later?
#  thetaG[thetaG > hi] <- hi[thetaG > hi]
  
#  tiny <- 1e-8
#  thetaG[thetaG <= tiny] <- tiny
  
  pu <- as.vector(thetaG)
  
  ss <- rep( rlnorm(ntnpG,log(.07),1), nsnb)
  ss <- ss*pu*(1 - pu)
  
  ss[whichGrow] <- ss[whichGrow]/(1 + as.vector( ydatG ) )
  
  thetaNew <- matrix( tnorm(ntnpG*nsnb,lo,hi,pu,ss), ntnpG, nsnb )
  
  list(thetaNow = thetaG, thetaNew = thetaNew)
}

sampleGrow <- function(){
  
  lg <- loGro[wind,]
  hg <- hiGro[wind,]
  hg[babsTreeMatG == 1] <- 0
  lg[babsTreeMatG == 1] <- -2
  lg[babsTreeMatG == 0] <- 0
  
#  FLAG <- F                   #change in current groG
#  ww <- which(groG < lg)
#  if(length(ww) > 0){
#    FLAG <- T
#    groG[ww] <- lg[ww]
#  }
  
 # groG[wwb & groG < 0]  <- lastGroPos[wind,][wwb & groG < 0]
 # groG[!wwb & groG > 0] <- lastGroNeg[wind,][!wwb & groG > 0]
  
 # sss <- abs(groG/100) + .001
  
  sss <- rep( rlnorm(ntnpG,log(.001),1), nsnb)
  
  sss[whichGrow] <- ygrSdG/7
  sss[sss > .05] <- .05
#  sss[sss < .001] <- .001
  
#  mu <- groG
#  mu[mu < lg] <- lg[mu < lg]
#  mu[mu > hg] <- hg[mu > hg]
  
  pgrow <- matrix( tnorm(ntnpnbG*nspec,lg,hg,groG,sss), ntnpG,nsnb)
  
  list(groNow = groG, groNew = pgrow)
}

updateDemogr <- function(){
  
  FLAG <- T
  
  tmp <- sampleTheta()
  thetaNow <- tmp$thetaNow
  thetaNew <- tmp$thetaNew
  logitNew <- logit(thetaNew)
  logitNow <- logit(thetaNow)
  
  tmp <- sampleGrow()
#  groNow <- tmp$groNow
  groNew <- tmp$groNew
#  FLAG   <- tmp$FLAG
#  if(FLAG)message('flag in updateDemogr')
  
  posNow <- groPosG
 # posNow <- groNow
#  posNow[posNow < 0] <- 0
  posNew <- groNew
  posNew[posNew < 0] <- 0
  
  # growth regression
  p1 <- p2 <- rep(0,ntnpnbG)
  xb <- xgroBeta[xFirstRegG,]
  
  p1[xLastRegG] <- dmvnormZeroMean( make2d(groG,nr=ntnpnbG)[xLastRegG,] - xb ,smat=sigma )
  p2[xLastRegG] <- dmvnormZeroMean( make2d(groNew,nr=ntnpnbG)[xLastRegG,] - xb ,smat=sigma )
  
  prGrowReg <- matrix(p2 - p1,ntnpG,nbreak) # by STAGE
  prGrowReg <- rowSums(prGrowReg)

  # theta regression
  vr  <- matrix(diag(msigma),ntnpG,nsnb,byrow=T)
  mu  <- matrix( xmortBeta,ntnpG,nsnb )
  tmp <- mu*0
  
  tmp[tnextG,] <- gaussianRatio(logitNew[tnextG,],logitNow[tnextG,],mu[tnowG,],vr[tnowG,])
  tmp <- byFunctionRcpp(as.vector(tmp),rep( c(1:ntnpG),nsnb),longSpecSG,
                         matJTS, matJTS,MEAN=F)
  prThetaReg <- rowSums(tmp)
  
  #SSD update
  if(!FLAG){
    predNew <- getGamMatrixSample( f1=fecG,gam1=wgamG,th1=thetaNew,gr1=posNew, 
                                   index=wind,
                                   ntnpK=ntnpG,time1=tnowG,time2=tnextG,first=firstTime,
                                   ktimes=4 )[[1]]
  }
  if(FLAG){
    tmp <- getGamMatrixSample( f1=fecG,f2=fecG,
                               gam1=wgamG,gam2=wgamG,th1=thetaNow,th2=thetaNew,
                               gr1=posNow,gr2=posNew, 
                               index=wind,
                               ntnpK=ntnpG,time1=tnowG,time2=tnextG,first=firstTime,
                               ktimes=4 ) 
    predGamG <- tmp[[1]]
    predNew  <- tmp[[2]]
  }
  
  tmp <- gaussianRatio(predNew,predGamG,wgamG,gamVar)
  tmp[,kdex == nbreak] <- 0
  prGamma <- rowSums( tmp )   
  
  # data ratio
 
  tmp <- gaussianRatio(groNew[whichGrow],groG[whichGrow],ygrMuG,(ygrSdG/growDataWt)^2)  #growth
  
  deathPr <- dtmG[whichGrow]
  gnow    <- (1 - thetaNow[whichGrow])^deathPr   #survival
  gnew    <- (1 - thetaNew[whichGrow])^deathPr
  tmp     <- tmp + binomialRatio(ydatG,ysurvG,gnew,gnow)*(1 + ydatG)*survDataWt
  tmp[!is.finite(tmp)] <- 0
  
  prData <- byFunctionRcpp( tmp,wideJTG[whichGrow],
                            tmp*0 + 1,matrix(0,ntnpG),matrix(0,ntnpG),MEAN=F) 
  prData[holdIndex] <- 0
 
  acrit <- prThetaReg + prGrowReg + prGamma + prData
  acrit[is.na(acrit)] <- 0
  
  aa <- exp(acrit)
  
  z <- runif(length(aa),0,1)
  ww <- which(z < aa)
  
  thetaNow[ww,] <- thetaNew[ww,]
  groG[ww,]   <- groNew[ww,]
  predGamG[ww,]  <- predNew[ww,]
  
  gall <- make2d(groG,nr=ntnpnbG,nc=nspec)
  
  list(thetaG = thetaNow, groG = groG, gall = gall,predModel = predGamG, accept = ww)
  
}


updateGrow <- function(){
  
  p1 <- p2 <- rep(0,ntnpnb)
  
  lg <- loGro
  hg <- hiGro
  lg[absTreeMat == 0] <- 0
  hg[absTreeMat == 1] <- 0
  lg[absTreeMat == 1] <- -2
  
  groMat[absTreeMat == 0 & groMat < 0] <- lastGroPos[absTreeMat == 0 & groMat < 0]
  groMat[absTreeMat == 1 & groMat > 0] <- lastGroNeg[absTreeMat == 1 & groMat > 0]
  
  groPos[groMat > 0] <- groMat[groMat > 0]
  
  sss <- abs(groMat/20) + .0000001
  
  sss[whichY] <- yGrowSd[whichY]
  
  #  psd <- rexp(ntot,sss)
  
  pgrow <- matrix( tnorm(length(groMat),lg,hg,groMat,sss), ntnp,nsnb)
  
  posNow <- groMat
  posNow[posNow < 0] <- 0
  posNew <- pgrow
  posNew[posNew < 0] <- 0
  
  #growth regression
  xb <- xall[-llong,]%*%beta
  
  pnow <- dmvnormZeroMean(make2d(groMat)[-flong,] - xb ,smat=sigma)
  pnew <- dmvnormZeroMean(make2d(pgrow)[-flong,] - xb ,smat=sigma)
  p1[-flong] <- pnow
  p2[-flong] <- pnew
  
  pnow <- matrix(p1,ntnp,nbreak)[-ftime,]
  pnew <- matrix(p2,ntnp,nbreak)[-ftime,]
  jmat <- matrix(jtdex[,'j'],ntnp,nbreak)
  
  ii <- list(plot = as.vector(jmat[-ftime,]))
  
  pnow1 <- byIndex(as.vector(pnow),ii,sum)
  pnew1 <- byIndex(as.vector(pnew),ii,sum)
  
  #gamma
  if(PDE){
    predNow <- getGam(fecundMat,gam,thetaMat,groPos)
    predNew <- getGam(fecundMat,gam,thetaMat,posNew)
  }
  if(!PDE){
    predNow <- getGamMatrixRcpp(fecundMat,gam,thetaMat,groPos)
    predNew <- getGamMatrixRcpp(fecundMat,gam,thetaMat,posNew)
  }
  
  predNow[tnext,] <- predNow[tnext,] - wgam[tnext,]
  predNew[tnext,] <- predNew[tnext,] - wgam[tnext,]
  predNow[ftime,] <- predNew[ftime,] <-0
  
  ggnow <- short2long(predNow,'STAGE')
  ggnew <- short2long(predNew,'STAGE')
  
  tmp <- omegaLambda1(ggnow,ggnew,lambda,lambda,psi,psi,ovar,ovar)
  pnow2 <- tmp$p1
  pnew2 <- tmp$p2
  
  pnow2[is.na(pnow2)] <- 0
  pnew2[is.na(pnew2)] <- 0
  
  tnow2 <- byFunctionRcpp( as.vector(pnow2[-ftime,]),as.vector(wideJmat[-ftime,fstage]),
                           as.vector(wideSmat[-ftime,fstage]),
                           matPlotSpec*0,matPlotSpec*0,MEAN=F) 
  tnew2 <- byFunctionRcpp( as.vector(pnew2[-ftime,]),as.vector(wideJmat[-ftime,fstage]),
                           as.vector(wideSmat[-ftime,fstage]),
                           matPlotSpec*0,matPlotSpec*0,MEAN=F) 
  pnow2 <- rowSums( tnow2 )
  pnew2 <- rowSums( tnew2 )  
  
  #individual growth
  pnow3 <- pnew3 <- groMat*0
  
  pnow3[whichY] <- dnorm(groMat[whichY],yGrowMu[whichY],yGrowSd[whichY],log=T)
  pnew3[whichY] <- dnorm(pgrow[whichY],yGrowMu[whichY],yGrowSd[whichY],log=T)
  
  pnow3[!is.finite(pnow3)] <- 0
  pnew3[!is.finite(pnew3)] <- 0
  
  pnow3 <- byFunctionRcpp( as.vector(pnow3[-ltime,]),as.vector(wideJmat[-ltime,]),
                           as.vector(wideSmat[-ltime,]),
                           matPlotSpec*0,matPlotSpec*0,MEAN=F) 
  pnew3 <- byFunctionRcpp( as.vector(pnew3[-ltime,]),as.vector(wideJmat[-ltime,]),
                           as.vector(wideSmat[-ltime,]),
                           matPlotSpec*0,matPlotSpec*0,MEAN=F) 
  
  pnow <- pnow1 + pnow2 + rowSums(pnow3)
  pnew <- pnew1 + pnew2 + rowSums(pnew3)
  
  a <- exp(pnew - pnow)
  w <- which( runif( nplot,0,1 ) < a)
  
  aindex <- which(wideJmat %in% w)
  
  groMat[aindex] <- pgrow[aindex]
  
  gall <- make2d(groMat)
  
  list(groMat = groMat, gall = gall, aindex = aindex)
  
}

  
ScInvWishFromSS <- function(SS,dff,delta=diag(1,nrow(SS)),
                            priorO=diag(1,nrow(SS)),priorOdf=(1+nrow(SS)),
                            priorDmu=rep(1,nrow(SS)),priorDvar=rep(10,nrow(SS)),
                            varLo=rep(1e-8,nrow(SS)),varHi=rep(10000,nrow(SS)) ){ 
  
  #varBound is permissible range of variances
  
  ww <- which(diag(delta) < 1e-5)
  if(length(ww) > 0)delta <- delta*10
  
  k <- nrow(SS)
  n <- dff - k - 1    #approximate n
  
  
  tmp <- getScWishMat(SS,dff,delta,priorO)
  IO  <- tmp$omegaInv
  O   <- tmp$omega
  sig <- tmp$sigma
  
  vars <- diag(sig)   # vars <- diag(SS)
  
  loDel <- sqrt( varLo/diag(O) )
  hiDel <- sqrt( varHi/diag(O) )
  
  ww <- which(loDel < 1e-6)
  if(length(ww) > 0)loDel[ww]  <- 1e-6
  ww <- which(hiDel < loDel)
  if(length(ww) > 0)hiDel[ww]  <- loDel[ww]*1.01
  
  psd <- diag(delta)/20
  
  propD <- diag( tnorm(k,loDel,hiDel,diag(delta),psd),k )
  
  ts <- IO*SS                
  t0 <- ts
  diag(t0) <- 0
  rs1 <- rowSums(t0/matrix(diag(delta),k,k,byrow=T))/diag(delta)
  rs2 <- rowSums(t0/matrix(diag(propD),k,k,byrow=T))/diag(propD)
  
  tsv <- diag(ts)/2/vars
  
  pnow <- -(n+1)*log(diag(delta)) - tsv - rs1 - ( (log(diag(delta)) - priorDmu)^2)/2/priorDvar
  pnew <- -(n+1)*log(diag(propD)) - tsv - rs2 - ( (log(diag(propD)) - priorDmu)^2)/2/priorDvar
  
  pnow <- sum(pnow)
  pnew <- sum(pnew)
  
  if( runif(1,0,1) < exp(pnew - pnow) ) diag(delta) <- diag(propD)
  
  sig <- delta%*%O%*%delta
  
  list(delta = delta, sigma = sig)
}

  
HB <- function(xx,yy,bb,sinv,minb,maxb,loPrior,hiPrior,avalue,bvalue){  #for stochastic approximation
  
  betaVar <- diag(100,length(bb))
  
  tp <- 0
  
  derivB <- t(xx)%*%(yy - xx%*%bb)%*%sinv
  bn     <- bb + avalue*derivB
  
  bn[bn < loPrior] <- loPrior[bn < loPrior]
  bn[bn > hiPrior] <- hiPrior[bn > hiPrior]
  
  p1 <- as.vector(bb)
  p2 <- as.vector(bn)
  
  dis <- dist2(p1,p2)
  IN <- F
  lob <- c(maxb - bn )
  hib <- c(bn  - minb)
  if( min(lob) >= 0 & max(hib) >= 0 )IN <- T
  
  if( dis < bvalue & IN ){
    bg <- bn
  } else{
    
  #  bg <- runif( length(bb),as.vector(minb),as.vector(maxb) )
    
    bg <- tnorm( length(bb), as.vector(minb),as.vector(maxb),0,4 )
    
  #  bg <- tnorm.mvtRcpp(bb*0,bb*0,betaVar,as.vector(minb),
  #                    as.vector(maxb),times=1)
    bg <- matrix(bg,nrow(bb),ncol(bb))
    dimnames(bg) <- dimnames(bb)
    tp <- 1
  }
  
  list(beta = bg, tp = tp)
}


dist2 <- function(x1,x2){
  
  sqrt(sum((x1 - x2)^2))
}

updateStocAppr <- function(xx,yy,bb,sigma,sinv=NULL,delta=diag(1,nspec),
                           lo,hi,loPrior,hiPrior,
                           tp=NULL,STOC=F,SSmat,VAR='MV',vprior=c(1,1),
                           varHi=rep(1000,nspec),avalue,bvalue){
  
  # VAR = 'MV'   - covariance matrix
  # VAR = 'UNI'  - covariance matrix is diagonal 
  # VAR = 'NONE' - no variance
  
  if(is.null(sinv))sinv <- chol2inv(chol(sigma))
  
  if(!STOC)bb <- bUpdateMVN_Rcpp(xx,yy,bb=bb,lo=lo,hi=hi,sigma)
  
  if(STOC){
    tmp <- HB( xx, yy ,bb,sinv,minb=lo,maxb=hi,loPrior,hiPrior,avalue,bvalue)
    bb  <- tmp$beta
    tp  <- tmp$tp
  }
  
  res <- yy - xx%*%bb
  
  if(VAR == 'MV'){
    
    SSx   <- crossprod(res)
    SSsum <- SSmat[nss+1,,] + SSx - SSmat[nss,,]
    SSmat[nss+1,,] <- SSmat[nss+1,,] + SSx - SSmat[nss,,]
    SSmat[2:nss,,] <- SSmat[1:(nss-1),,]
    SSmat[1,,]  <- SSx
    
    dff <- min(g*nrow(xx),ntnpnb) + nspec + 1
    tmp <- ScInvWishFromSS(SSsum,df=dff,priorDmu=rep(.1,nspec),priorDvar=rep(1,nspec),
                         varLo=rep(1e-8,nspec),varHi=varHi  )
    sigma <- tmp$sigma
    delta <- tmp$delta
    
  }
  if(VAR == 'UNI'){
    
    require(pscl)
      
    u1 <- vprior[1] + nrow(res)/2
    u2 <- vprior[2] + colSums( res^2 )/2
    varTheta <- 1/rgamma(nspec,u1,u2)
    
  #  varTheta <- rtrunc(nspec,varHi*0,varHi,u1,u2,'igamma')
    
    sigma    <- diag(varTheta)
    delta    <- numeric(0)
  }
    
  
  list(beta = bb, tp = tp, SSmat = SSmat, sigma = sigma, delta = delta)
}


hiLo <- function(lo,hi,loPrior,hiPrior,tp){
  
  lo <- lo - tmp$tp
  hi <- hi + tmp$tp
  lo[lo < loPrior] <- loPrior[lo < loPrior]
  hi[hi > hiPrior] <- hiPrior[hi > hiPrior]
  
  list(lo = lo, hi = hi)
}

treeVsPlotDens <- function(bins=NULL){
  
  if(is.null(bins))bins <- c(0,100,100000)
  tmp  <- y/sampArea
  plot <- byFunctionRcpp( as.vector(tmp),as.vector(wideJmat),as.vector(wideSmat),
                          matPlotSpec*0,matPlotSpec*0,MEAN=F)
  tree <- byFunctionRcpp( as.vector(y),as.vector(wideJmat),as.vector(wideSmat),
                          matPlotSpec*0,matPlotSpec*0,MEAN=F)
  tHist <- pHist <- numeric(0)
  for(s in 1:nspec){
    treeHist <- histWeight(plot[,s],plot[,s],bins=bins,removeZeros=F,
                         useMids=F)
    plotHist <- histWeight(plot[,s],tree[,s]*0+1,bins=bins,removeZeros=F,
                         useMids=F)
    tHist <- rbind(tHist,treeHist$dens)
    pHist <- rbind(pHist,plotHist$dens)
  }
  bins <- treeHist$bins
#  plot(treeHist$bins+1,plotHist$dens,type='s',col='black',log='x')
#  lines(treeHist$bins+1,treeHist$dens,type='s',col='red')
  
  list(bins = bins, tree = tHist, plot = pHist)
}
  
plotRate <- function(s,demName = 'grow',demog=growMeanDatAll,
                     temp = climAll[,'temp'],
                     prec = climAll[,'prec'],
                     maplon=range(plotLon),maplat=range(plotLat),LEGEND=F){
  
  require(stats)
  

  wc    <- which(sdex == s)
  nn    <- ytotDataAll[,wc]
  demog <- demog[,wc]
  
  www        <- which(nn == 0)
  nn[www]    <- NA
  demog[www] <- NA
  
  ytot  <- rowSums(nn,na.rm=T)
  dMean <- rowSums(demog,na.rm=T)/ytot
  
  ww <- which(is.finite(dMean))
  
  if(demName == 'surv')dMean <- 1 - dMean
  
  ww <- which(is.finite(dMean))
  
  xx <- temp[ww]
  yy <- dMean[ww]
  
  xseq   <- seq(min(xx) - 1,max(xx) + 1,length=500)
  
  ylab <- 'Growth rate (cm)'
  if(demName == 'surv')ylab <- 'Survival probability (1/yr)'
  
  plot(xx,yy,cex=.3,xlab='Temperature',ylab=ylab,xlim=c(-15,20))
  
  ma <- predict( smooth.spline(x=xx,y=yy,nknots=20) ,xseq)$y
  lines(xseq,ma,lwd=5,col='grey')
  lines(xseq,ma,lwd=3,lty=2)
  
  
  ma <- predict( smooth.spline(x=xx,y=yy,nknots=20) ,xseq)$y
  lines(xseq,ma,lwd=5,col='grey')
  lines(xseq,ma,lwd=3,lty=2)
  
  if(demName == 'grow'){
    wtAve <- sum(xx*yy)/sum(yy)
    wtSd  <- sqrt( sum(xx^2*yy)/sum(yy) - wtAve^2 )
    ym    <- dnorm(xseq,wtAve,wtSd)
    lines(xseq,ym/max(ym)/2,lwd=4,col='white')
    lines(xseq,ym/max(ym)/2,lwd=2,lty=3)
    lines(xseq,(ym/max(ym)/2),lwd=2)
  }
  
  
  xx <- prec[ww]
  xseq   <- seq(min(xx) - 1,max(xx) + 1,length=500)
  
  plot(xx,yy,cex=.3,xlab='Precipitation',ylab=ylab,xlim=c(500,2200))
  
  ma <- predict( smooth.spline(x=xx,y=yy,nknots=20) ,xseq)$y
  lines(xseq,ma,lwd=5,col='grey')
  lines(xseq,ma,lwd=3,lty=2)
  
  ma <- predict( smooth.spline(x=xx,y=yy,nknots=20) ,xseq)$y
  lines(xseq,ma,lwd=5,col='grey')
  lines(xseq,ma,lwd=3,lty=2)
  
  if(demName == 'grow'){
    wtAve <- sum(xx*yy)/sum(yy)
    wtSd  <- sqrt( sum(xx^2*yy)/sum(yy) - wtAve^2 )
    ym    <- dnorm(xseq,wtAve,wtSd)
    lines(xseq,ym/max(ym)/2,lwd=4,col='grey')
    lines(xseq,ym/max(ym)/2,lwd=2,lty=3)
    lines(xseq,(ym/max(ym)/2),lwd=2)
  }
  
  if(LEGEND)legend('topright',c('model','spline'),lty=c(1,2),lwd=c(2,2),bty='n')
  
}



getClimEnvelop <- function(climRef,climMat,lonLat,
                           climVars = c('temp','prec'),
                           delta = c('+ 5','* .7'),lty=1,
                           POINTS=F,CONTOURS=F,col='black',lwd=2,FILL=T){
  
  # delta - character, not formula format
  
  op <- substr(delta,1,1)
  dl <- as.numeric( substr(delta,3,10) )
  
  nref <- nrow(climMat)
  
  wpoints <- numeric(0)
  
  if(length(lwd) == 1)lwd <- rep(lwd,nrow(climRef))
  if(length(col) == 1)col <- c(1:nrow(climRef))
  
  for(j in 1:nrow(climRef)){
    
    wj <- c(1:nref)
    
    for(k in 1:length(climVars)){
      
      print(k)
      x  <- climRef[j,climVars[k]]
      y  <- climMat[,climVars[k]]
      
      hi <- lo <- x
      
      if(op[k] == '+'){
        hi <- x + dl[k]
        lo <- x - .1
      }
      if(op[k] == '-'){
        lo <- x - dl[k]
        hi <- x + .1
      }
      if(op[k] == '*'){
        if(dl[k] > 1){
          hi <- x*dl[k]
          lo <- x*.99
        }
        if(dl[k] < 1){
          lo <- x*dl[k]
          hi <- x*1.01
        }
      }
      
      wk <- which(y < lo | y > hi)
      wj <- wj[!wj %in% wk]
    }
    
    wpoints <- append(wpoints,list(wj))
    
    zz <- rep(0,nref)
    zz[wj] <- 1
    
    wsea <- which(topo$z < 0,arr.ind=T)
    xsea <- topo$x[wsea[,1]]
    ysea <- topo$y[wsea[,2]]
    zsea <- rep(min(zz)-1,length(xsea))
    
    lon <- c(lonLat[,'lon'],xsea)
    lat <- c(lonLat[,'lat'],ysea)
    zzz <- c(zz,zsea)
    
    
    if(POINTS)points(lonLat[wj,'lon'],lonLat[wj,'lat'],col=col)
    if(CONTOURS){
      values2contour(xx=lon,yy=lat,
                           zzz,nx=80,ny=80,col=col[j],lty=lty,
                           zlevs=c(.5,100),lwd=lwd[j],add=T,fill=FILL)
    }
    
  }
  
  invisible(wpoints)
}



sizeSpeciesDistribution <- function(specNames,plotIndex,bins2Plot,yscale=20){
  
  nss <- length(specNames)
  nbb <- length(bins2plot)
  
  sIndex <- match(specNames,specs)
  yy <- y[plotIndex,]
  area <- sampArea[plotIndex,]
  dxx  <- dxMat[plotIndex,]
  bmat <- breakMat[plotIndex,]
  
  pnames <- sort( unique(plotnames[jtdex[plotIndex,'j']]) )
  np     <- length(pnames)
  
  perHa <- matrix(0,nss,nbb)
  rownames(perHa) <- specNames
  
  for(s in 1:nss){
    
    ws     <- which(sdex == sIndex[s])
    if(length(ws) == 0)next
    bs     <- bmat[,ws]
    mids   <- (bs[,-nbreak] + bs[,-1])/2
    bounds <- cbind(mids[,1]*0,mids)
    mids[,nbreak-1] <- bs[,nbreak-1] + bs[,nbreak-1] - bs[,nbreak-2]
    mids <- cbind(mids,mids[,nbreak-1] + 500)
    
    wt <- yy[,ws]/area[,ws]/length(plotIndex)
    
    w0 <- which(wt > 0,arr.ind=T)
    
    tmp <- histWeight(mids[w0],wt[w0],bins=bins2plot,removeZeros=T,useMids=F)
    
    perHa[s,] <- tmp$wtHist
  }
  
  
  ncol <- 100
#  colseq <- mapColors(ncol)
  
  greyColors   <- colorRampPalette(c('white','orange','brown','black'))
  
  colseq <- greyColors(ncol)
  
  scale <- 10^seq(log10(.1),log10(100),length.out=ncol)
  
  ww   <- as.matrix(expand.grid(c(1:nss),c(1:nbb)))
  
  icol <- findInterval(perHa,scale,all.inside=T)
  coli <- colseq[icol]
  
  dlast <- diff(bins2plot)[nbb-1]
  blast <- c(bins2plot[-1],bins2plot[nbb]+dlast)
  xleft <- bins2plot[ww[,2]]
  xright <- blast[ww[,2]]
  ybottom <- ww[,1] 
  ytop    <- ww[,1] + 1
  
  tmp <- mapSetup(range(bins2plot),c(1,nss),yscale)
   
  par(bty='n',las=1,pin=c(tmp[1],tmp[2]*7))
  plotSetup(xtic=bins2plot,ytic=.5+seq(1,nss,by=1),xvals = bins2plot, 
            yvals = specNames, xlabel='Diameter (cm)',ylabel=' ',fc='white',
            endFactor=c(.01,.05))
  
  rect(xleft,ybottom,xright,ytop,col=coli,border=coli)
  
  invisible(perHa)
}


getMainInt <- function(isIntPrior = isIntPriorX, varNames = xnames){
  
  mainNames <- character(0)
  mainCols <-  numeric(0)
  
  ww <- which(varNames %in% colnames(isIntPrior))
  if(length(ww) > 0){
    
    intNames  <- varNames[ww]
    wcol      <- which(colnames(isIntPrior) == intNames)
    wmain     <- rownames( which(isIntPrior[,wcol] == 1,arr.ind=T) )
    mainNames <- matrix(wmain,length(intNames),2,byrow=F)

    mainCols  <- matrix( match(mainNames,varNames),length(intNames),2)
    colnames(mainNames) <- colnames(mainCols) <- intNames
  }
  list( mainNames = mainNames, mainCols = mainCols )
}
    


predictX_plotBA <- function(predictYnames = c('gall','theta','phi'),priorMux = .5,
                            priorVx = 1){

#  njj <- ntnpnb
#  if(vshort %in% c('temp','prec') | vk %in% c('therm','deficit'))njj <- ntnp
  
#  priorxj <- priorMuX[1:njj,]
  v  <- Vi <- rep(0,ntnpnb)

  
  for(jj in 1:length(predictYnames)){
    
    tmp <- getRegParts(predictYnames[jj])
    b   <- tmp$b
    sig <- tmp$sig
    sigIv <- tmp$sigIv
    xx  <- tmp$x
    yy  <- tmp$y
    ix  <- tmp$ix
    iy  <- tmp$iy
    
    if(!'plotBA' %in% colnames(xx))next
    
    ws <- match('plotBA',colnames(xx))         #also works for phi: first ix rows is 1st break

    yx     <- t(yy[iy,] - xx[ix,-ws]%*%b[-ws,])
    v[ix] <- v[ix] + t(b[ws,])%*%sigIv%*%yx
    
    Vi[ix] <- Vi[ix] + b[ws,]%*%sigIv%*%matrix(b[ws,],nspec,1)
    
  }
  
  V  <- 1/(Vi + 1/priorVx)
  mu <- V*(v + priorMux/priorVx)
  
  #   plot(xPrior[c(1:ntnp)[-ftime],vk],mu[-ftime],cex=.2)
  #   abline(0,1)
  #   title(predictYnames[jj])
  
  
  tnorm(length(v),0,1,mu,sqrt(V))
  
  
} 

predictXLinear <- function(predictXnames=c('plotBA','temp','therm'),
                           predictYnames = c('wmeantot','gall','theta','phi'),
                           priorMuX = xPrior,
                           priorVX = 1/ntnp ){
  
  # wgam  - abundance
  # gall  - growth
  # phi   - fecundity
  # theta - mortality
  # priorIX - prior inverse covariance matrix
  
  # 'temp', 'therm', 'prec' - 1 value for each ntnp: xall[1:ntnp,**] gets all
  # 'plotBA' - does not repeat: there are ntnpnb unique values
  #          - fall contains maximum plotBA (at ground surface)
  
  # vname: plotBA has ntnpnb values
  #        climate has ntnp   values for phi
  #                    nplot  values for wmeantot
  #                    ntnpnb values for theta, gall
  
  nxx <- length(predictXnames)
  
  predX <- matrix(0,ntnpnb,nxx)
  colnames(predX) <- predictXnames
  
  
  for(kk in 1:nxx){
    
    vk  <- predictXnames[kk]
    vshort <- 'nothing'
    
    ww <- grep('prec',vk)
    if(length(ww) > 0)vshort <- 'prec'
    ww <- grep('temp',vk)
    if(length(ww) > 0)vshort <- 'temp'
    
    njj <- ntnpnb
    if(vshort %in% c('temp','prec') | vk %in% c('therm','deficit'))njj <- ntnp
    
    priorxj <- priorMuX[1:njj,]
    v  <- Vi <- matrix(0,njj,1)

    for(jj in 1:length(predictYnames)){

      tmp <- getRegParts(predictYnames[jj])
      b   <- tmp$b
      sig <- tmp$sig
      sigIv <- tmp$sigIv
      xx  <- tmp$x
      yy  <- tmp$y
      ix  <- tmp$ix
      iy  <- tmp$iy
      
      if(!vk %in% colnames(xx))next
      
      tmp <- getVv(vname=vk,rname=predictYnames[jj],xx,yy,b,sigIv,ix,iy )
      
  #    iz <- c(1:length(tmp[[1]]))
      if(length(tmp$v) == nplot)iz <- jtdex[,'j']
      if(length(tmp$v) == ntnp)iz <- c(1:ntnp)[-ftime]
      
      v[ iz ]      <- v[ iz ] + tmp$v[ iz ]
      Vi[iz] <- Vi[iz] + tmp$Vi[1]
      
    }
      
    V  <- 1/(Vi + 1/priorVX)
    mu <- V*(v + priorxj[,vk]/priorVX)
    
 #   plot(xPrior[c(1:ntnp)[-ftime],vk],mu[-ftime],cex=.2)
 #   abline(0,1)
 #   title(predictYnames[jj])
                                        
    
    tmp <- tnorm(length(v),0,1,mu,sqrt(V))
    
    if(vshort %in% c('temp','prec') | vk %in% c('therm','deficit'))predX[ rep(c(1:ntnp),nbreak),vk] <- tmp
    if(vk == 'plotBA')predX[,vk] <- tmp
    
  } 
  predX
}


predictXFactor <- function(predictXnames=c('xeric','mesic'),
                           predictYnames = c('wmeantot','gall','theta','phi'),
                           priorMuX = xPrior,
                           priorWt = 10 ){
  
  nxx <- length(predictXnames)
  
  namesX <- predictXnames
  namesI <- character(0)
  
  XINT <- F
  tmp <- getMainInt(isIntPrior=isIntPriorX, varNames = xnames)
  mainNames <- tmp$mainNames
  mainCols  <- tmp$mainCols
  
  if(length(mainCols) > 0)XINT <- T
  
  nnn <- nagg <- nplot
  ixx <- ftime
  if(XINT){
    nnn  <- ntnpnb   #interaction with a climate variable
    ixx  <- 1:nnn
    nagg <- nplot     #aggregate over ntnp plot-intervals
  }
  
  nowX   <- xall[ixx,predictXnames]
  priorX <- xPrior[ixx,predictXnames]
  nowX[nowX[,'xeric'] == -1 & nowX[,'mesic'] == 1,] <- priorX[nowX[,'xeric'] == -1 & nowX[,'mesic'] == 1,]
  propX  <- nowX*0
  
  pX <- nowX[1:nplot,]*0
  
  tmp  <- t( rmultinom(nplot,1,c(1/3,1/3,1/3)) )
  pX[tmp[,3] == 1,2] <- 1
  pX[tmp[,1] == 1,1] <- -1
  
  wxer <- which(priorX[ftime,'xeric'] == -1)
  wmes <- which(priorX[ftime,'mesic'] == 1)
  
  xer <- rbinom(length(wxer),1,.5)
  mes <- rbinom(length(wmes),1,.5)
  
  pX[wxer,'xeric'] <- -xer
  pX[wxer,'mesic'] <- 0
  pX[wmes,'mesic'] <- mes
  pX[wmes,'xeric'] <- 0
  
  if(nnn == ntnpnb)propX <- pX[rep(jtdex[,'j'],nbreak),]  
  
  colnames(propX) <- predictXnames
  
  predX <- xall[ftime,namesX]
  xm    <- matrix(0,nagg,1)
  pdiff <- rep(0,nagg)
  
  for(jj in 1:length(predictYnames)){
    
    rname <- predictYnames[jj]
    tmp <- getRegParts(rname)
    b   <- tmp$b
    sig <- tmp$sig
    xx  <- tmp$x
    yy  <- tmp$y
    ix  <- tmp$ix
    iy  <- tmp$iy
    nxx <- nrow(xx)
    pnow <- pnew <- rep(0,nxx)
    
    #   if(!vk %in% colnames(xx))next
    
    nxx   <- nrow(xx)
    
    if(nxx == ntnp)  iz <- jtdex[,'j']
    if(nxx == ntnpnb)iz <- rep(jtdex[,'j'],nbreak)
    if(nxx == nplot) iz <- ftime
    
    xprop <- xx
    xprop[,predictXnames] <- propX[iz,]
    
    for(jjj in 1:ncol(mainNames)){
      xprop[,colnames(mainNames)[jjj] ] <- xprop[,mainNames[1,jjj]]*xprop[,mainNames[2,jjj]]
    }
    
    xb <- xx[ix,]%*%b  
    xp <- xprop[ix,]%*%b
    
    pnow[iy] <- dmvnormZeroMean( yy[iy,] - xb ,smat=sig ) 
    pnew[iy] <- dmvnormZeroMean( yy[iy,] - xp ,smat=sig )
    
    if(nxx == ntnp & nagg == nplot) iz <- jtdex[,'j']
    #     if(nxx == ntnp & nagg == ntnp)  iz <- c(1:ntnp)
    #     if(nxx == nplot & nagg == ntnp) iz <- jtdex[ftime,'j']
    #     if(nxx == ntnp & nagg == ntnpnb)iz <- stuff
    if(nxx == ntnpnb)iz <- rep( jtdex[,'j'],nbreak )
    
    if(nxx == nplot) iz <- 1:nplot
    pdiff    <- pdiff + byFunctionRcpp(pnew - pnow,iz,iz*0+1,xm,xm,MEAN=F)
  }
  
  rnow <- rowSums(nowX)
  rnew <- rowSums(propX)
  rpri <- rowSums(priorX)
  
  pnow <- pnew <- rep(1,nagg)/(1 + priorWt)
  pnow[rnow == rpri] <- priorWt/(1 + priorWt)
  pnew[rnew == rpri] <- priorWt/(1 + priorWt)
  
  iz <- rep( jtdex[,'j'],nbreak )
  xm <- matrix(0,nplot,1)
  pdiff <- pdiff + byFunctionRcpp( log(pnew) - log(pnow), iz, iz*0+1,xm,xm,MEAN=F)
  
  a <- exp(pdiff)
  z <- runif(nplot,0,1)
  
  ww <- which(z < a)
  
  newX <- nowX[ftime,]
  newX[ww,predictXnames] <- pX[ww,]
  
  # predX[ww,predictXnames] <- propX[ww,]
  
  xall[,namesX] <- pX[rep( jtdex[,'j'],nbreak ),] 
  
  if( XINT ){      
    
    for(jjj in 1:ncol(mainNames)){
      xall[,colnames(mainNames)[jjj] ] <- xall[,mainNames[1,jjj]]*xall[,mainNames[2,jjj]]
    }
    
    #   xall[,namesX[3]] <- xall[,namesX[1]]*xall[,namesI[1]]
    #   xall[,namesX[4]] <- xall[,namesX[2]]*xall[,namesI[2]]
  }
  
  list(predX = newX, xall = xall)
}

  
  
#################################
predictX <- function(predictXnames,
                     predictYnames = c('wmeantot','gall','theta','phi'),
                     priorMuX = xPrior,
                     priorWtX = 10,xUpdate,x2update = c('xeric','mesic') ){
  
  nxx   <- length(predictXnames)
  longj <- rep( jtdex[,'j'],nbreak )
  longb <- rep(1:ntnp,nbreak)
  
  namesX <- predictXnames
  namesI <- character(0)
  
  XINT    <- F
  XFACTOR <- F
  if(length(ixnames) > 0)XINT <- T
  if('xeric' %in% xnames)XFACTOR <- T     #only xeric, mesic are factors
  
  np    <- length(predictXnames)
  nxx   <- ntnp*length(predictXnames)
  sdd   <- rexp(nxx,1/.05)
  tmp   <- tnorm(nxx,0,1,xUpdate[1:ntnp,predictXnames] ,sdd)
  propX <- matrix(tmp,ntnp,np)
  propX <- propX[longb,]
  colnames(propX) <- predictXnames
  
  if(XFACTOR){          
    
    tmp <- myrmultinom(1,xPropHydro) 
    pX   <- xUpdate[ftime,c('xeric','mesic')]*0
    pX[tmp[,1] == 1,'xeric'] <- -1
    pX[tmp[,3] == 1,'mesic'] <- 1
    
    propX[,c('xeric','mesic')] <- pX[longj,c('xeric','mesic')] 
  }
  
  xm    <- matrix(0,nplot,1)      #accept at plot (not plot-year)
  pdiff <- rep(0,nplot)
  
  for(jj in 1:length(predictYnames)){
    
    rname <- predictYnames[jj]
    
    tmp <- getRegParts(rname)
    b   <- tmp$b
    sig <- tmp$sig
    xx  <- tmp$x
    yy  <- tmp$y
    ix  <- tmp$ix
    iy  <- tmp$iy
    isIntPrior <- tmp$isIntPrior
    noIntPrior <- tmp$noIntPrior
    intPrior   <- tmp$intPrior
    nxx <- nrow(xx)
    pnow <- pnew <- rep(0,nxx)
    
    pnames <- predictXnames[predictXnames %in% colnames(xx)]
    
    nxx   <- nrow(xx)
    
    if(nxx == nplot) iz <- 1:nplot
    if(nxx == ntnp)  iz <- 1:ntnp
    if(nxx == ntnpnb)iz <- 1:ntnpnb

    xprop <- xx
    xprop[,pnames] <- propX[iz,pnames]
    
    if(XINT){
      for(jjj in 1:length(intPrior)){
        wjj <- rownames(isIntPrior)[isIntPrior[,jjj] == 1]
        xprop[,intPrior[jjj]] <- xprop[,wjj[1]]*xprop[,wjj[2]]
      }
    }
    
    xb <- xx[ix,]%*%b  
    xp <- xprop[ix,]%*%b
    
    pnow[iy] <- dmvnormZeroMean( yy[iy,] - xb ,smat=sig ) 
    pnew[iy] <- dmvnormZeroMean( yy[iy,] - xp ,smat=sig )
    
    if(nxx == nplot) iq <- c(1:nplot)
    if(nxx == ntnp)  iq <- jtdex[,'j']
    if(nxx == ntnpnb)iq <- longj
    
    pdiff <- pdiff + byFunctionRcpp(pnew - pnow,iq,iq*0+1,xm,xm,MEAN=F)
  }
  
  if(XFACTOR){
    rnow <- rowSums(xUpdate[ftime,c('xeric','mesic')])
    rnew <- rowSums(pX)
    rpri <- rowSums(priorMuX[ftime,c('xeric','mesic')])
    
    pwt <- priorWtX['xeric']
    
    pnow <- pnew <- rep(1,nplot)
    pnow[rnow == rpri] <- pwt
    pnew[rnew == rpri] <- pwt
  }

  pr <- pdiff + log(pnew) - log(pnow)
  a  <- exp(pr)
  z  <- runif(nplot,0,1)
  ww <- which(z < a)
  
  wlong  <- which(longj %in% ww)
  wshort <- which(jtdex[,'j'] %in% ww)
  
  xUpdate[wlong,predictXnames]  <- propX[wlong,]
  xall[,x2update] <- xUpdate[,x2update]  #UPDATES MAIN BUT NOT INTERACTIONS
  
  if(XINT){
    for(jjj in 1:length(intPriorX)){
      wjj <- rownames(isIntPriorX)[isIntPriorX[,jjj] == 1]
      xUpdate[,intPriorX[jjj]] <- xUpdate[,wjj[1]]*xUpdate[,wjj[2]]
  #    xall[,intPriorX[jjj]]    <- xUpdate[,intPriorX[jjj]]
      xall[,intPriorX[jjj]]    <- xall[,wjj[1]]*xall[,wjj[2]]
    }
  }
  
  pnames <- predictXnames[predictXnames %in% fnames]
  fall   <- xall[1:ntnp,fnames]
  
  pnames <- predictXnames[predictXnames %in% gnames]
  xgall  <- xall[ftime,gnames]
  
  list(xall = xall, fall = fall, xgall = xgall, xUpdate = xUpdate)
}

    
getRegParts <- function(rname){
  
  if(rname == 'phi'){
    ix  <- fecRegIndexF
    iy  <- fecRegIndexL
    xx  <- fall
    yy  <- phi
    b   <- betaFec
    sig <- fsigma
    sigIv <- finv
    isIntPrior <- isIntPriorF 
    noIntPrior <-  noIntPriorF
    intPrior   <- intPriorF 
  }
  if(rname == 'theta' | rname == 'gall'){
    ix <- groRegIndexL
    iy <- groRegIndexF 

    if(rname == 'theta'){
      yy  <- logitTheta
      b   <- betaMort
      xx <- xall[,mnames]
      sig <- msigma
      sigIv <- minv
      isIntPrior <- isIntPriorM 
      noIntPrior <-  noIntPriorM
      intPrior   <- intPriorM 
    }
    if(rname == 'gall'){
      yy  <- gall
      b   <- beta
      xx  <- xall
      sig <- sigma
      sigIv <- sinv
      isIntPrior <- isIntPriorX 
      noIntPrior <-  noIntPriorX
      intPrior   <- intPriorX 
    }
  }
  if(rname == 'wmeantot'){
    ix  <- iy <- c(1:nplot)
    xx  <- xgall
    yy  <- wmeantot
    b   <- betaGam
    sig <- gsigma
    sigIv <- ginv
    isIntPrior <- isIntPriorG 
    noIntPrior <-  noIntPriorG
    intPrior   <- intPriorG 
  }
  list(b = b, sig = sig, sigIv = sigIv, x = xx, y = yy, ix = ix, iy = iy,
       isIntPrior = isIntPrior, noIntPrior = noIntPrior, intPrior = intPrior )
}

       
  
getVv <- function(vname = 'temp',rname = 'theta',
                  xx,yy,b,sigIv,
                  ix=1:nrow(xx),iy=1:nrow(yy)){
  
  # rname: phi only contributes to first stage, 1:ntnp
  #        wmeantot only contributes to 1:nplot
  #        theta, gall to 1:ntnpnb
  # vname: plotBA has ntnpnb values
  #        climate has ntnp   values for phi
  #                    nplot  values for wmeantot
  #                    ntnpnb values for theta, gall
    
    ws <- match(vname,colnames(xx))
    
    yx <- t(yy[iy,] - xx[ix,-ws]%*%b[-ws,])
    v0 <- t(b[ws,])%*%sigIv%*%yx
    
    v  <- rep(0,nrow(xx))
    
    ww <- grep('temp',vname)
    if(length(ww) > 0)vname <- 'temp'

    ww <- grep('prec',vname)
    if(length(ww) > 0)vname <- 'prec'
 
    v <- matrix(0,ntnp,1)
    
 #   if(vname %in% c('therm','deficit'))v <- matrix(0,ntnp,1)  
    if(vname == 'plot BA')             v <- matrix(0,ntnpnb,1)  
    
    if(vname %in% c('temp','prec','therm')){ 
      
      if(rname %in% c('theta','gall')){ # aggregate over stage; do nothing for wmeantot
      
        intnp   <- 1:ntnp
        intnpnb <- rep(intnp,nbreak)
      
        xm <- matrix(0,ntnp,1)
    
        tmp <- byFunctionRcpp( as.vector(v0),intnpnb[iy],ix*0+1,xm, xm,MEAN=F)
        nn  <- byFunctionRcpp( as.vector(v0)*0 + 1,intnpnb[iy],ix*0+1,xm, xm,MEAN=F)
        v   <- tmp/nn
        v[!is.finite(v)] <- 0
      }
      
      if(rname == 'phi')v[ix] <- v0
         
    }
    
    if(vname == 'plotBA')  v[ix,] <- v0
    if(rname == 'wmeantot')v <- v0
    
    Vi <- b[ws,]%*%sigIv%*%matrix(b[ws,],nspec,1)
    
    list(v = v, Vi = Vi)
}
      
    
getRegIndex <- function(){
  
  tmp <- y
  tmp[tmp > 0] <- 1
  tmp <- apply(tmp*wideKmat,1,max) + 1
  tmp <- by(tmp,jtdex[,'j'],max,na.rm=T)
  tmp <- tmp[jtdex[,'j']]
  
  tmp <- rep(tmp,nbreak)                #last class
  cmp <- rep(c(1:nbreak),each=ntnp)     #class in xall
  
  WLAST <- tmp >= cmp
  
  jj <- rep(jtdex[,'j'],nbreak )
  
  groRegIndexL <-  which(WLAST & !jj %in% holdOutPlots & !c(1:ntnpnb) %in% llong)
  groRegIndexF <-  which(WLAST & !jj %in% holdOutPlots & !c(1:ntnpnb) %in% flong)
  
  regFirstIndex <- regLastIndex <- rep(0,ntnpnb)
  regFirstIndex[groRegIndexL] <- 1
  regLastIndex[groRegIndexF] <- 1
  
  list(groRegIndexL = groRegIndexL, groRegIndexF = groRegIndexF, 
       regFirstIndex = regFirstIndex, regLastIndex = regLastIndex)
}

climateMap <- function(maplon,maplat,lon,lat,topo,climVec,zlevs,colorRamp,nx,ny,
                       xbox,xscale,labside='left',legside='bottomright',
                       endLabels=c('100','4000'),addContour,contourLabel=F,
                       mapscale=5){
  
  lonLat <- cbind(lon,lat)
  xr     <- apply(lonLat,2,range)
  
  wx <- which(topo$x >= xr[1,1] & topo$x <= xr[2,1])
  wy <- which(topo$y >= xr[1,2] & topo$y <= xr[2,2])
  
  xx <- topo$x[wx]
  yy <- topo$y[wy]
  zz <- topo$z[wx,wy]
  
  xs <- seq(1,length(xx),by=1)
  ys <- seq(1,length(yy),by=1)
  
  xx <- xx[xs]
  yy <- yy[xs]
  zz <- zz[xs,ys]
  
  minz <- min(zz,na.rm=T)
  
  tmp <- oceanValues(minz=minz,minc=0,xx,yy,zz,lonLat,climVec)
  cVec <- tmp$clim
  cLoc <- tmp$lonLat
  cMin <- min(climVec)
  
#  wo <- which(topo$z < -10,arr.ind=T)
#  ocean <- cbind(topo$x[wo[,1]],topo$y[wo[,2]])
  
  colF   <- colorRampPalette(colorRamp)
  col <- colF(length(zlevs))
  
  regMap(x=topo$x,y=topo$y,z=topo$z,IMAGE=F,xl=maplon,yl=maplat,
         mapscale=mapscale,lineCol='grey',axes=F,lwd=2) 
  
  values2contour(xx=cLoc[,'lon'],yy=cLoc[,'lat'],
                 z=cVec,nx=nx,ny=ny,col=col,lty=1,labcex=1,
                 zlevs=zlevs,lwd=2,fill=T,add=T,drawlabels=F)
  
  regMap(x=topo$x,y=topo$y,z=topo$z,IMAGE=F,xl=maplon,yl=maplat,
         mapscale=mapscale,lineCol='white',axes=F,ADD=T,lwd=1) 
  
  
  if(!missing(addContour)){
    values2contour(xx=cLoc[,'lon'],yy=cLoc[,'lat'],
                 z=cVec,nx=nx,ny=nx,
                 col='black',lty=1,labcex=1,
                 zlevs=addContour,lwd=2,fill=F,add=T,drawlabels=contourLabel)
  }
  
#  dx <- diff(range(topo$x))/70
#  symbols(ocean[,1],ocean[,2],squares=rep(dx,nrow(ocean)),inches=F,fg='white',
#          bg='white',add=T)
  
  mapMask(cLoc[,'lon'],cLoc[,'lat'],dx=.1,dy=.1,whiteOutDist=.3,col='white')
  
  map('state',interior=F,col='grey',add=T,lwd=4) 
  map('state',interior=F,col='white',add=T,lwd=2)
  
 # rect(xbox[1,1],xbox[1,2],xbox[2,1],xbox[2,2],col='grey',border='white')
  colorLegend(xscale[,1],xscale[,2],scale=c(1:length(col)),cols=col,
              labside=labside,endLabels=endLabels)
  
}


getGibbsSummaries <- function(chainList,ng){
  
  
  for(k in 1:length(chains)){
    
    tmp <- chainList[[k]]
    if(is.matrix(tmp))  chainList[[k]] <-tmp[1:ng,]
    if(!is.matrix(tmp)) chainList[[k]] <-tmp[1:ng]
  }
  
  nmatrix <- matrix(nsample,nplot,nspec)
  nmatrix[nmatrix == 0] <- 1
  
  npmatrix <- matrix(npsample,ntnp,nspec)
  npmatrix[npmatrix == 0] <- 1
  
  prAbsentMu <- prAbsent/nmatrix
  prAbsRecMu <- prAbsRec/nmatrix
  prAbsFecMu <- prAbsFec/nmatrix
  
  prAbsTreeMu   <- prAbsTree/nmatrix
  prAbsTreeNCMu <- prAbsTreeNC/nmatrix
  
  prAbsentMu[!is.finite(prAbsentMu)]       <- 0  # no pop growth
  prAbsTreeMu[!is.finite(prAbsTreeMu)]     <- 0  # no tree growth
  prAbsTreeNCMu[!is.finite(prAbsTreeNCMu)] <- 0  # no tree growth w/o comp
  prAbsRecMu[!is.finite(prAbsRecMu)]       <- 0  # no recruitment
  prAbsFecMu[!is.finite(prAbsFecMu)]       <- 0  # no recruit potential
  
  outfile <- outFile(outfolder=outfolder,'mcmc.txt')
  
  tmp <- processMCMC(chains = chains, sums = sums, burnin = burnin,outfile=outfile,
                     outfolder=savefolder,
                     chainList = chainList, chainSum = chainSum,PLOTS=F)
  postSummary <- tmp$postSummary
  postMuSe    <- tmp$postMuSe 
  parMeans    <- tmp$parMeans
  parSds      <- tmp$parSds
  kappaMu     <- postMuSe$kappa[,1:nspec]
  kappaSe     <- postMuSe$kappa[,-c(1:nspec)]
  tauPost     <- postSummary[grep('tau',rownames(postSummary)),]
  tauMu       <- tauPost[,'estimate']
  
  predGamMu <- predGamSum/predGamNo
  predGamMu[predGamMu < 0] <- 0
  fecMu     <- fecundMatSum/npmatrix
  phiMu     <- phiSum/npmatrix
  ImmMu     <- ImmSum/npmatrix
  gamMu     <- gamSum/matrix(npsample,ntnp,nsnb)
  gamSd     <- sqrt(gamSum2/matrix(npsample,ntnp,nsnb) - gamMu^2)
  groMu     <- groMatSum/matrix(npsample,ntnp,nsnb)
  thetaMu   <- thetaMatSum/matrix(npsample,ntnp,nsnb)
  
  ws <- grep('sigma',rownames(postSummary))
  smean <- matrix( postSummary[ws,1],nspec,nspec )
  colnames(smean) <- rownames(smean) <- specs
  cmeanGro <- cov2cor(smean)
  
  ws <- grep('fsigma',rownames(postSummary))
  kmean <- matrix( postSummary[ws,1],nspec,nspec )
  colnames(kmean) <- rownames(kmean) <- specs
  cmeanFec <- cov2cor(kmean)
  
  ws <- grep('gsigma',rownames(postSummary))
  gmean <- matrix( postSummary[ws,1],nspec,nspec )
  colnames(gmean) <- rownames(gmean) <- specs
  cmeanGam <- cov2cor(gmean)
  
  list( prAbsentMu  = prAbsentMu,prAbsRecMu = prAbsRecMu, prAbsRec = prAbsRec,
         prAbsFecMu = prAbsFecMu, prAbsFec = prAbsFec,prAbsTreeMu = prAbsTreeMu,
         prAbsTree = prAbsTree,prAbsTreeNCMu = prAbsTreeNCMu,postSummary = postSummary,
         postMuSe = postMuSe,parMeans = parMeans,parSds = parSds,kappaMu = kappaMu,
         kappaSe = kappaSe,tauPost = tauPost,tauMu = tauMu,predGamMuJ = predGamMu,
         fecMu = fecMu,phiMu = phiMu,ImmMu = ImmMu,gamMuJ = gamMu,gamSdJ = gamSd,
         groMu = groMu,thetaMu = thetaMu,cmeanGro = cmeanGro,cmeanFec = cmeanFec,
         cmeanGam = cmeanGam, chainList = chainList )
}


getPosteriorScores <- function(trim=-100,minVar=.001){
  
  xpredMu <- xPred/npred
  xpredVr <- xPred2/npred - xpredMu^2
  
  xpredVr[xpredVr < minVar] <- minVar
  
  pc  <- predictXnames[!predictXnames %in% c('xeric','mesic')]
  pc  <- c(pc,'plotBA')
  nc  <- length(pc)
  tmp <- getScoreNorm(xPrior[1:ntnp,pc],xpredMu[1:ntnp,pc],xpredVr[1:ntnp,pc])
  
  tmp[tmp < trim] <- trim
  
  score <- byFunctionRcpp(as.vector(tmp),
                          rep(jtdex[,'j'],nc),
                          rep(1:nc,each=ntnp),
                          matrix(0,nplot,nc),matrix(0,nplot,length(pc)),MEAN=T)
  colnames(score) <- pc
  
  tmp <- getScoreNorm(xPrior[,'plotBA'],xpredMu[,'plotBA'],xpredVr[,'plotBA'])
  tmp[tmp < -300] <- -300
  tmp <- byFunctionRcpp(tmp,rep(jtdex[,'j'],nbreak),
                        rep(c(1:nbreak),each=(ntnp)),
                        matrix(0,nplot,nbreak),matrix(0,nplot,nbreak),MEAN=T)
  tmp[tmp < -300] <- -300
  scoreBAbySize <- tmp
  scoreBAbySize[,1] <- 0
  
  size4BA <- 20     # diameter (cm) for BA prediction
  tmp <- wideKmat
  tmp[breakMat > size4BA] <- 0
  ki <- apply(tmp,1,max)[ftime]
  
  score[,'plotBA'] <- scoreBAbySize[cbind(c(1:nplot),ki)]
  
  list(score = score, xpredMu = xpredMu, xpredVr = xpredVr)
}


getX2climate <- function(vname='temp'){
  
  # from xall to climate variable
  # uncenter, unstandardize
  
  cmat   <- matrix(0,nplot,1)
  wc     <- grep(vname,xnames)[1]
  tmp    <- byFunctionRcpp(xall[1:ntnp,wc],jtdex[,'j'],rep(1,ntnp),cmat,cmat,MEAN=T)
  if( vname %in% c('deficit') )tmp    <- 1 - tmp  # in reverseSign
  wr     <- xrange[,xnames[wc]]
  tmp*(diff(wr)) + wr[1]
  
}

packPlotMat <- function(mat,rname){
  rownames(mat) <- rname
  mat
}

combinePlots <- function(mat,ROUND=3,sortCols=F,sortRows=T){
  
  #retain columns in mat, aggregate rows by plotname REGION-PLOT
  #rows combined by idex
  
  nn <- nrow(mat)
  nc <- ncol(mat)
  
  tmp <- matrix( unlist(strsplit(rownames(mat),'-')),nrow(mat),2,byrow=T) 
  
  allnames <- unique(tmp[,2])
  nall     <- length(allnames)
  
  idex <- match(tmp[,2],allnames)
  
  ii <- rep(idex,nc)
  jj <- rep(c(1:nc),each=nn)
  
  mm <- matrix(0,nall,nc)
  
  tmp <- byFunctionRcpp( as.vector(mat),ii,jj,mm,mm,MEAN=T)
  colnames(tmp) <- colnames(mat)
  rownames(tmp) <- allnames
  
  tmp <- round(tmp,ROUND)
  
  if(sortRows)tmp <- tmp[order(rownames(tmp)),]
  if(sortCols)tmp <- tmp[,order(colnames(tmp))]
  tmp
}


getRegionChains <- function(includeSections,chainList){
  
  kchains <- numeric(0)
  for(r in 1:length(includeSections)){
    ww <- which(names(chainList) == includeSections[r])
    kchains <- append(kchains,list(chainList[[ww]]))
    #  if(names(chainList)[r] %in% includeSections)kchains <- append(kchains,list(chainList[[r]]))
  }
  names(kchains) <- includeSections
  kchains
}

UJSDM_fia_gibbs <- function(x,y,breaks=NULL,
                        thetaPriorSample=NULL,
                        thetaPriorSpecies=NULL,
                        betaPrior=NULL,ng=100,burnin=1,
                        nugget = diag(diag(cov(y)))*.000001,
                        TYPE=NULL,ZEROINFL=F,maxBreaks = 50,
                        compCols=c(1:ncol(y))){
  
  # x          - n by p design, first column intercept (ones)
  # y          - n by S response, continuous and discrete depending on TYPE
  
  # betaPrior  - beta prior truncated at zero
  # ng         - no. iterations
  # maxBreaks  - if too many classes, aggregated to this value
  # breaks     - NULL means ordinal data, breaks fitted
  # wmiss      - cbind(row,col) of missing values
  # TYPE       - 'ordinal','discrete','continuous','composition','presenceAbsence'
  # nugget     - may be needed if large S, small n
  # compCols   - columns that sum to one for composition data
  
  # thetaSamplePrior, thetaSpeciesPrior - used if ZEROINFL
  # thetaSamplePrior  - n by 2 beta prior values (q1, q2)
  # thetaSpeciesPrior - S by 2 beta prior values (q1, q2)
  
  #      continuous: c(-Inf,0,U,Inf)     censored at U
  #     composition: c(-Inf,0,1)         fractions
  #        discrete: c(-Inf,0,...,U,Inf) U is largest count minus 1 or censored value
  #         ordinal: c(-Inf,0,...,Inf)   no. of classes + 1
  # presenceAbsence: c(-Inf,0,Inf)   no. zeros and ones
  
  types <- c('continuous','discrete','composition','ordinal','presenceAbsence')
  tvec  <- paste0(types,collapse=", ")
  
  if(is.null(TYPE))    stop( paste('specify TYPE: ', tvec,sep="") )   
  if(!TYPE %in% types) stop( paste('TYPE must be one of: ',tvec,sep="") )
  
  if(burnin >= ng)stop( 'burnin must be greater than no. MCMC steps, ng' )
  
  # missing values in x
  wmiss <- which(!is.finite(x),arr.ind=T)
  nmiss <- length(wmiss)
  
  if(length(wmiss) > 0){         #initialize missing values with means
    xmean    <- colMeans(x,na.rm=T)
    x[wmiss] <- xmean[wmiss[,2]]
    xprior   <- x[wmiss]
    warning( paste(nmiss,' missing values in x will be imputed',sep='') )
  }
  
  #check design
  checkInt <- range(x[,1])
  if(checkInt[1] != 1 | checkInt[2] != 1)stop( paste('1st X must be intercept (ones)') )
  colnames(x)[1] <- 'intercept'
  
  tmp <- checkDesign(x)
  
  maxy <- max(y)
  
  print(tmp)
  if(tmp$rank < tmp$p)stop( 'x not full rank' )
  
  S <- ncol(y)
  q <- ncol(x)
  n <- nrow(y)
  
  snames <- colnames(y)
  xnames <- colnames(x)
  
  if(is.null(snames))snames <- paste('S',1:S,sep='-')
  if(is.null(xnames))xnames <- paste('x',1:q,sep='-')
  
  snames <- sub('_','-',snames)
  xnames <- sub('_','-',xnames)
  
  colnames(y) <- snames
  colnames(x) <- xnames
  
  updateBeta <- UJSDM_updateBetaNoPrior
  loBeta <- hiBeta <- NULL
  
  if(!is.null(betaPrior)){
    
    loBeta <- matrix(-Inf,q,S)
    xx     <- matrix( unlist(strsplit(betaPrior$negBeta,'_')),ncol=2,byrow=T )
    loBeta[ cbind( match(xx[,2],xnames),match(xx[,1],snames) ) ] <- 0
    loBeta <- as.vector(loBeta)
    
    hiBeta <- matrix(Inf,q,S)
    xx     <- matrix( unlist(strsplit(betaPrior$posBeta,'_')),ncol=2,byrow=T )
    hiBeta[ cbind( match(xx[,2],xnames),match(xx[,1],snames) ) ] <- 0
    hiBeta <- as.vector(hiBeta)
    
    updateBeta <- UJSDM_updateBetaPrior
  }                 
  
  if(ZEROINFL & min(y) > 0)ZEROINFL <- F        #must be zeros for ZEROINFL
  fixTheta <- NULL
  loSpec   <- hiSpec <- numeric(0)
  
  w <- z <- y
  
  if(TYPE == 'presenceAbsence'){
    
    if(min(z) == 0)z <- z + 1    #zeros are first class
    breaks <- c(-Inf,0,Inf)
    ncut   <- length(breaks)
    cuts   <- matrix(breaks,S,ncut,byrow=T)
    
    UJSDM_updateW <- UJSDM_updateWdiscrete
    w                <- UJSDM_Y2W(y,breaks,cuts,lohi=1)  
    
    tmp <- speciesCountSummary(z,aggregate=F,TYPE=TYPE)
    classBySpec <- tmp$classBySpec
  }  
  
  if(TYPE == 'ordinal'){
    if(min(y) == 0)y <- y + 1    #zeros are first class
    breaks <- c(-Inf,c(0:(max(y)-2)),Inf)
    z <- y
    ZEROINFL <- F
  }
  
  if(TYPE == 'discrete'){
    if(is.null(breaks))breaks <- c(-Inf,seq(0,(max(y)-1)),Inf)
    if(length(breaks) > 50){
      warning('breaks created')
      breaks <- round( seq(min(y),(max(y) - 1),length=maxBreaks), 0)
      breaks <- c(-Inf,breaks,Inf)
      z      <- matrix( findInterval(y,breaks),n,S)
    }
  }
  
  if(TYPE == 'ordinal' | TYPE == 'discrete'){
    
    ncut <- length(breaks)
    cuts <- matrix(breaks,S,ncut,byrow=T)
    
    UJSDM_updateW <- UJSDM_updateWdiscrete
    w <- UJSDM_Y2W(y,breaks,cuts,lohi=1)  
    
    tmp <- speciesCountSummary(z,aggregate=F,TYPE=TYPE)
    if(min(z) == 0)z <- z + 1
    classBySpec  <- tmp$classBySpec
    classAll     <- tmp$classAll
    lowestClass  <- tmp$lowestClass
    highestClass <- tmp$highestClass
    
    if(TYPE == 'ordinal'){
      
      fixTheta <- infTheta <- numeric(0)
      
      if(max(lowestClass) > 1){                #shift lowest class
        wl <- which(lowestClass > 1)
        for(kk in wl){
          ks <- 3:(lowestClass[kk]+1)
          fixTheta <- rbind(fixTheta, cbind(rep(kk,length(ks)),ks ) )
        }
      }
      if(min(highestClass) < (ncut-1)){       #find highest class
        wl <- which(highestClass < (ncut-1))
        for(kk in wl){
          ks <- (highestClass[kk]+1):(ncut-1)
          infTheta <- rbind(infTheta, cbind(rep(kk,length(ks)),ks ) )
        }
      }
    }
  }
  
  if(TYPE == 'composition'){
    
    UJSDM_updateW <- UJSDM_updateWcomp
    
    if(is.null(breaks)){
      breaks <- c(-Inf,0,1)
      breaks <- matrix(breaks,S,3,byrow=T)
    }
    cuts <- breaks
    
    z   <- y*0 + 1
    z[y > 0] <- 2
    z[y < 0] <- 2     #neg values must be predictors
    
    w  <- y
    ww <- which(z == 1)
    w[ww] <- runif(length(ww),-.1,0)
    
    tmp <- speciesCountSummary(z,aggregate=F,TYPE)
    classBySpec <- tmp$classBySpec
    ncut        <- ncol(classBySpec) + 1
  }
  
  if(TYPE == 'continuous'){
    
    UJSDM_updateW <- UJSDM_updateWcont
    
    #censored values set to 2nd to last cut
    if(!is.null(breaks) & !is.matrix(breaks)){
      breaks <- matrix(breaks,n,length(breaks),byrow=T)
    }
    if(is.null(breaks)){
      breaks <- matrix(c(-Inf,0,Inf),n,length(breaks),byrow=T)
    }
    
    cuts <- breaks
    ncut <- length(breaks)
    
    wy0 <- which(breaks == 0,arr.ind=T)[,1]      #specs with zero class
    wy1 <- as.matrix( which(y == 0,arr.ind=T) )  #lower threshold (usually 0)
    wy1 <- wy1[wy1[,2] %in% wy0,]
    
    w[wy1] <- runif(nrow(wy1),-maxy/2,0)
    
    if(length(wy1) > 0)loSpec <- sort(unique(wy1[,2]))
    
    i2  <- as.matrix( which(y == breaks[3],arr.ind=T) )     #censored above
    wy3 <- numeric(0)
    
    if(length(i2) > 0){
      w[i2]  <- runif(nrow(i2),breaks[3],breaks[3]+1)
      z      <- w
      wy3    <- i2
      hiSpec <- sort(unique(wy3[,2]))
    }
    
    #classes defined for discrete/continuous
    z <- y
    z[y == 0] <- 1
    z[y > 0]  <- 2
    z[y < 0]  <- 2                       #negative values must be predictors
    if(ncol(breaks) > 3){                #censored class
      z3 <- matrix(breaks[,3],n,S,byrow=T)
      z[y > z3] <- 3
    }
    
    tmp <- speciesCountSummary(z,aggregate=F)
    classBySpec <- tmp$classBySpec
  }
  
  tmp <- UJSDM_specBreakI(z)
  i1  <- tmp$i1  
  i2  <- tmp$i2
  
  #observed zeros
  y0 <- y
  y0[y > 0] <- 1
  y0 <- 1 - y0
  
  b <- y*0
  thetag  <- 1                                     # initial Pr fail to observe when present
  wabs   <- which(w <= 0,arr.ind=T)                # if not ZEROINF, absent where not seen
  wabs   <- wabs[wabs[,2] %in% compCols,]
  wpres  <- numeric(0)                             # present but not observed
  
  whereZero <- which(y == 0,arr.ind=T)
  noZeros   <- which(!c(1:n) %in% sort(unique(whereZero[,1])))
  
  b[wabs] <- 1
  
  if(ZEROINFL){
    sumb <- b*0
    b[whereZero] <- 1
    w[whereZero] <- tnorm(nrow(whereZero),-1,0,0,.4)
    thetag <- .1
    wabs   <- which(b == 1,arr.ind=T)              # absent -- required for sampleW
    wpres  <- which(b == 0 & y == 0,arr.ind=T) 
    
    if(is.null(thetaPriorSample) & is.null(thetaPriorSpecies)){
      stop('must have thetaPrior if ZEROINFL')
    }
    if(!is.null(thetaPriorSample)){
      updateZero <- UJSDM_updateZeroSample
      thetaPrior <- thetaPriorSample
    }
    if(!is.null(thetaPriorSpecies)){
      updateZero <- UJSDM_updateZeroSampleSpecies
      thetaPrior <- thetaPriorSpecies
    }
  }
  wfill <- rbind(wabs,wpres)
  
  clim <- Inf
  
  cnames <- paste('C',1:ncut,sep='-')
  
  XX  <- crossprod(x)
  IXX <- solve(XX)
  WX  <- crossprod(x,w)
  WIX <- IXX%*%WX
  
  tg       <- cutg <- cuts
  sg       <- nugget*100 + diag(.1,S)
  sigmaDf  <- n - q + S - 1
  alpha    <- matrix(0,q,S)
  
  alpha <- matrix( updateBeta(WIX,IXX,sg,alpha,loBeta,hiBeta),q,S)
  bg    <- UJSDM_A2B(alpha,sg)
  xb    <- x%*%alpha
  
  colnames(y) <- rownames(classBySpec) <- rownames(sg) <- colnames(sg) <- snames
  colnames(x) <- xnames
  
  colnames(classBySpec) <- cnames[1:ncol(classBySpec)]
  
  print(classBySpec)
  
  rgibbs <- matrix(0,ng,S*S)
  bgibbs <- matrix(0,ng,S*q)
  
  colnames(rgibbs) <- as.vector( outer(snames,snames,paste,sep='_') )
  sgibbs <- lgibbs <- rgibbs
  colnames(bgibbs) <- as.vector( t(outer(snames,xnames,paste,sep='_')) )
  
  
  if(TYPE == 'ordinal'){
    cgibbs <- matrix(0,ng,(ncut-3)*S)
    rownames(cutg) <- snames
    colnames(cgibbs) <- as.vector( outer(snames,cnames[-c(1,2,ncut)],paste,sep='_') )
  }
  
  ypred  <- ypred2 <- sumb <- y*0
  sumDev <- 0   #for DIC
  sMean  <- sg*0
  ntot   <- 0
  
  pbar <- txtProgressBar(min=1,max=ng,style=1)
  
  for(g in 1:ng){
    
    WX  <- crossprod(x,w)                                #covariance
    WIX <- IXX%*%WX
    SS  <- crossprod(w) - t(WX)%*%WIX + nugget
    SI  <- solve(SS)
    sg  <- solve( rwish(sigmaDf,SI) )
    
    alpha <- bg <- matrix( updateBeta(WIX,IXX,sg,alpha,loBeta,hiBeta),q,S )
    muw   <- x%*%bg
    
    
    if(TYPE == 'ordinal' | TYPE == 'presenceAbsence') bg    <- UJSDM_A2B(alpha,sg)  
    
    if(TYPE == 'ordinal'){
      tg    <- UJSDM_updateTheta(w,z,tg,fixTheta,infTheta,lohi=100)
      cutg  <- UJSDM_theta2cuts(tg,sg)                 # correlation scale
      cgibbs[g,] <- cutg[,-c(1,2,ncut)]
    }
    
    w <- UJSDM_updateW(w,muw,sigma = sg,y = y,b = b,theta = tg,loHi = 100,
                       i1 = i1, i2 = i2,
                       wabs = wabs, wfill = wfill, wy1 = wy1, wy3 = wy3, 
                       noZeros = noZeros,
                       loSpec = loSpec, hiSpec = hiSpec,compCols=compCols)
    
    if(ZEROINFL){
      thetag <- updateZero(thetaPrior,y0,b)
      b      <- UJSDM_sampleB(muw,w,y,sg,thetag)          # latent indicator of absence
      wabs   <- which(b == 1,arr.ind=T)                   # absent -- required for sampleW
      wpres  <- which(b == 0 & y == 0,arr.ind=T)          # present but not observed
      wfill  <- rbind(wabs,wpres)
      sumb   <- sumb + b
    }
    
    if(nmiss > 0){
      x[wmiss] <- imputX_MVN(x,w,alpha,wmiss,sinv,xprior)[wmiss]
      XX  <- crossprod(x)
      IXX <- solve(XX)
    }
    
    cg <- cov2cor(sg)
    
    rgibbs[g,] <- cg
    bgibbs[g,] <- bg
    lgibbs[g,] <- cor(w)
    sgibbs[g,] <- sg
    
    if(g > burnin){
      ntot <- ntot + 1
      if(!ZEROINFL)yp <- myrmvnorm(n,xb,sg)
      
      if(ZEROINFL){
        lo <- y*0
        hi <- y*0 + maxy
        lo[wabs] <- -maxy
        hi[wabs] <- 0
        yp <- tnorm.mvtMatrixRcpp(avec=w, muvec=muw, smat=sg,lo=lo,hi=hi)
      }
      yp[yp < 0] <- 0
      
      if(TYPE == 'ordinal' | TYPE == 'presenceAbsence')yp <- UJSDM_w2y(yp,sg)
      
      ypred  <- ypred + yp
      ypred2 <- ypred2 + yp^2
      sumDev <- sumDev + sum(deviance(w,x,alpha,s=sg,LIKE='mvnorm'))
      sMean  <- sMean + sg
    }
    
    setTxtProgressBar(pbar,g)
  }
  
  sMean <- sMean/ntot
  
  chains <- list( rgibbs = rgibbs, lgibbs = lgibbs, 
                  sgibbs = sgibbs, bgibbs = bgibbs ) 
  
  tmp <- processPars(bgibbs)$summary
  bMu <- matrix(tmp[,'estimate'],q,S)
  bSe <- matrix(tmp[,'se'],q,S)
  aMu <- bMu
  
  yMu <- ypred/ntot
  ySd <- sqrt(ypred2/ntot - yMu^2)
  cMu <- cuts
  cSe <- numeric(0)
  
  # note: on latent w scale
  meanDev <- sumDev/ntot
  pd  <- meanDev - sum(deviance(yMu,x,aMu,sMean,LIKE='mvnorm'))
  DIC <- 2*pd + meanDev
  
  score <- sum( getScoreNorm(y,yMu,ySd^2) )  # gaussian w
  
  if(TYPE == 'ordinal'){
    tmp <- processPars(cgibbs)$summary
    cMu <- matrix(tmp[,'estimate'],S,ncut-3)
    cSe <- matrix(tmp[,'se'],S,ncut-3)
    colnames(cMu) <- colnames(cSe) <- cnames[-c(1,2,ncut)]
    aMu <- UJSDM_B2A(bMu,sMean)
    rownames(cMu) <- rownames(cSe) <- colnames(aMu) <- snames
    chains <- c(chains,list(cgibbs = cgibbs))
    
    zMu <- yMu*0
    for(s in 1:S)zMu[,s] <- findInterval(yMu[,s],c(-Inf,0,cMu[s,],Inf))
    yMu <- zMu
    ySe <- numeric(0)
  }
  
  if(TYPE == 'presenceAbsence'){
    zMu <- yMu*0
    for(s in 1:S)zMu[,s] <- findInterval(yMu[,s],c(-Inf,0,Inf))
    yMu <- zMu - 1
    ySe <- numeric(0)
  }
  
  tmp <- processPars(rgibbs)$summary
  rMu <- matrix(tmp[,'estimate'],S,S)
  rSe <- matrix(tmp[,'se'],S,S)
  
  tmp <- processPars(lgibbs)$summary
  lMu <- matrix(tmp[,'estimate'],S,S)
  lSe <- matrix(tmp[,'se'],S,S)
  
  tmp <- processPars(sgibbs)$summary
  sMu <- matrix(tmp[,'estimate'],S,S)
  sSe <- matrix(tmp[,'se'],S,S)
  
  
  betaCov <- t(bMu)%*%cov(x)%*%bMu
  
  prAbs <- numeric(0)
  if(ZEROINFL)prAbs <- sumb/ntot
  
  colnames(bMu) <- colnames(bSe) <-
    rownames(rMu) <- colnames(rMu) <- rownames(rSe) <- 
    colnames(rSe) <- snames
  rownames(bMu)   <- rownames(bSe) <- xnames
  
  absProb <- sumb/ng
  
  list(classBySpec = classBySpec, TYPE = TYPE, absProb = absProb,
       chains = chains, breaks = breaks, x = x, y = y,
       cutMu = cMu, cutSe = cSe, betaMu = bMu, betaSe = bSe, 
       corMu = rMu, corSe = rSe, sigMu = sMu, sigSe = sSe,
       wcorMu = lMu, wcorSe = lSe, betaCov = betaCov,
       yMu = yMu, ySd = ySd, prAbs = prAbs, compCols = compCols,
       DIC = DIC, score = score)
}


mapDemogr <- function(nn,demName='y',specI,sizeRange=c(0,500),
                      mapscale=5,presProb=.6,
                      zrange,mapx=maplon,mapy=maplat,PRES=F,insetPos=NULL){
  
  # specI - species
  # sizeRange - in cm
  # presProb - draw line where species is predicted to have positive rates
  # mapx, mapy - bounds of map
  # PRES - draw probability line
  
  
  if(demName == 'growDat')demog  <- growMeanDatAll
  if(demName == 'growEst')demog  <- growMeanEstAll
  if(demName == 'reproEst')demog <- reproEstAll
  if(demName == 'reproDat')demog <- reproDatAll
  if(demName == 'fecnEst')demog  <- phiEstAll
  if(demName == 'survDat')demog  <- survDatAll
  if(demName == 'survEst')demog  <- survEstAll
  if(demName == 'y')demog  <- ytotAllData
  
  colF   <- colorRampPalette(c('white','wheat',
                               'tan','brown','black'))  #order by fit

  
  ws <-  which(specs %in% specI)
    
  if(demName %in% c('reproDat','reproEst','fecnEst')){
    wc    <- ws
    dMean <- demog[,wc]
    pres  <- 'rep'
    
  } else {
    wc    <- which(specs[sdex] %in% specI)
    nn <- nn[,wc]
    demog <- demog[,wc]
    bk    <- breaksTot[,kdex][,wc]
    
    www        <- which(bk < sizeRange[1] | bk > sizeRange[2])
    nn[www]    <- NA
    demog[www] <- NA
    www        <- which(nn == 0)
    nn[www]    <- NA
    demog[www] <- NA
    
    ytot  <- rowSums(nn,na.rm=T)
    dMean <- rowSums(demog,na.rm=T)/ytot
    pres <- 'tree'
  }
  
  if(length(wc) == 0)stop('species not present')
  
  if(demName == 'survDat')dMean <- 1 - dMean
  
  ww <- which(is.finite(dMean))
  
  dMean[dMean > zrange[2]] <- zrange[2]
  
  zlevs <- c(seq(zrange[1],zrange[2],length=100),zrange[2]*1.2)
  regMap(topo$x,topo$y,topo$z,IMAGE=F,
         xl=mapx, yl = mapy, mapscale=mapscale,lineCol='grey')
  
  values2contour(lonLatAll[ww,1],lonLatAll[ww,2],
                 dMean[ww],nx=50,ny=50,col=colF(length(zlevs)),
                 zlevs=zlevs,lwd=2,add=T,drawlabels=F,fill=T)
  
  if(PRES){
    if(pres == 'tree'){
      z0 <- presTreeNCAll[,ws]
      values2contour(lonLatAll[,1],lonLatAll[,2],
                     z0,nx=60,ny=60,col=c('white','wheat'),
                     zlevs=c(-1,presProb),lwd=c(1,2),add=T,drawlabels=F,fill=T)
      
      z1 <- presPopAll[,ws]
      values2contour(lonLatAll[,1],lonLatAll[,2],
                     z1,nx=60,ny=60,col=c('white','brown'),
                     zlevs=c(presProb,presProb),lwd=c(5,3),add=T,drawlabels=F)
    } else {
      if(demName == 'reproEst')z0 <- presReprAll[,ws]
      if(demName == 'reproDat')z0 <- presReprAll[,ws]
      if(demName == 'fecnEst') z0 <- presFecAll[,ws]
      
      values2contour(lonLatAll[,1],lonLatAll[,2],
                     z0,nx=60,ny=60,col=c('white','wheat'),
                     zlevs=c(-1,presProb),lwd=c(1,2),add=T,drawlabels=F,fill=T)
    }
  }
  
  mapMask(lonLatAll[,1],lonLatAll[,2],dx=.15,dy=.15,whiteOutDist=.35,col='white')
  
  map('state',add=T,col='white',lwd=5,interior=F)
  map('state',add=T,col='grey',lwd=1,interior=T)
  
  if(!is.null(insetPos))add.scatter(hist2Add(dMean[ww],xlim=zrange),posi=insetPos,
                                     ratio=.11,inset=c(.03,.13),bg.col='grey')
  
  
  ksn <- paste(sizeRange[1],sizeRange[2],sep='-')
  
  
}

specsInLonLat <- function(mlon,mlat,lon,lat,yy,minPlot=20){
  
  ww <- which(lon > mlon[1] &
                lon < mlon[2] &
                lat > mlat[1] &
                lat < mlat[2])
  
  yy <- yy[ww,]
  yy[yy > 0] <- 1
  
  tmp <- colSums(yy)
  
  tmp[tmp >= minPlot]
}


