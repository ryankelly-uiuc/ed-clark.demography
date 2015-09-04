
library(mvtnorm)
library(Matrix)
library(SparseM)


getSeedProb2 <- function(u,fd=0,fmat,wtree=c(1:ntree),wyr=c(1:nyr)){  

  # fmat - ntree by nyr matrix of fecundities
  # wtree - which trees to include
  # wyr can be 1 yr or c(1:nyr)
  # fd is coefficient, multiplied by BA (bap), for background seed

  tiny <- 1e-6

  ww <- which(seedmat[,'t'] %in% wyr)
  yy  <- matrix(NA,nrow(seedmat),2)
  colnames(yy) <- c('lambda','lnL')

  if(!is.matrix(fmat))fmat <- as.matrix(fmat)

  fec <- numeric(0)

  bap <- fd*baPlotYr[cbind(seedmat[,'j'],seedmat[,'t'])]
  bap[is.na(bap)] <- tiny


  kern <- getKern(u,distall[,wtree],tinyKern=1e-10)
  kern[kern < tiny] <- 0

  fecmt  <- fmat[wtree,]
  fecmt[is.na(fecmt)] <- 0
  w0     <- which(rowSums(fecmt,na.rm=T) > 1)
  if(length(w0) > 0){
  
    fecmt  <- fecmt[w0,]
    nc     <- length(w0)
    tmp <- getSeedProb_cpp(ntrap,nc,nyr,fecmt,kern)
  }

 # dpois(seedmat[,'count'],ss*tmp,log=T)
}



seedPred <- function(fmat,kern,fdg){

  kern%*%fmat + fdg*BASeedMat
}


getSeedProb3 <- function(fmat,kern,seed,fdg=0,spred=NULL){ 

   if(is.null(spred))spred <- as.matrix( seedPred(fmat,kern,fdg) )
   dpois(seed,spred*trAreaMat + 1e-10,log=T)
}



##########################################

XtimesB2 <- function(xx,b,tindex=c(1:n),YR=F){

   b1 <- 0
   if(YR)b1 <- b[ cbind(yrIndex[tindex],treemat[tindex,'spec']) ]

   bb <- t( b[pCol,treemat[tindex,'spec']])
   rowSums(xx[tindex,pCol]*bb) + b1
}

XtimesBMV <- function(xx,z,blist,tindex=c(1:n),yrList=yrList){
  
  wi <- which(treeRows(matr=0) | treeRows(sex=0))
  
  for(k in 1:length(blist)){
    YR     <- F
    b      <- get(blist[k])
    if(blist[k] == 'bd' & 'YRD' %in% yrList)YR <- T
    if(blist[k] == 'bh' & 'YRH' %in% yrList)YR <- T
    if(blist[k] == 'bf' & 'YRF' %in% yrList)YR <- T

    z[,k] <- XtimesB2(xx,b,tindex,YR=YR)

    if(blist[k] == 'bd')z[wi,ynames == 'dinc'] <- XtimesB2(xx,bdjuv,tindex=wi,YR=YR)
    if(blist[k] == 'bh')z[wi,ynames == 'hinc'] <- XtimesB2(xx,bhjuv,tindex=wi,YR=YR)
    if(blist[k] == 'bf')z[wi,ynames == 'fecn'] <- 0
  }
  z
}
  


XtimesB <- function(x,b){

   b <- t(b[,treemat[,'spec']])
   rowSums(x*b)
}

getXY <- function(xx,yy){

  solve( crossprod(xx) ) %*% crossprod(xx,yy) 

}

getKmat <- function(light,kappa){

#  for(k in 1:ny)Kmat[,k] <- light*kappa[k]
  for(k in 1:ny)Kmat[,k] <- light/(light + kappa[k])

  Kmat
}

y2z <- function(y,Kmat){  # yK(KK')^(-1), both as 2 columns
  
  cbind(y/Kmat)

}

setupBeta <- function(xx,yy,posPriors,negPriors){

  q    <- length(znames)
  beta <- matrix(0,q,1)

  yy <- yy/.2

  bp          <- getXY(xx[,pCol],yy)
  beta[pCol,] <- bp

  if(length(yrCol) > 0){
      wc <- which(colSums(xx[,yrCol]) > 0)
      bt <- getXY(xx[,yrCol[wc]],yy - xx[,pCol]%*%bp)
      beta[yrCol[wc],] <- bt
  }

  ww  <- which(znames %in% negPriors & beta > 0)
  if(length(ww) > 0)beta[ww,] <- -.001
  ww  <- which(znames %in% posPriors & beta < 0)
  if(length(ww) > 0)beta[ww,] <- .001

  res <- yy - xx%*%beta

  sg <- crossprod(res)/length(yy)

  beta <- matrix(beta,q,nttype)
  rownames(beta) <- znames
  colnames(beta) <- treeNames
  list(beta = beta, sg = sg)
}
  

posNegClimPars <- function(pnlist,vnames){

  zn <- xnames
  wi <- grep('X',xnames)
  if(length(wi) > 0)zn <- xnames[-wi]
  
  nv <- length(vnames)

  for(kk in 1:nv){
    wt <- grep(vnames[kk],pnlist)
    wx <- grep(vnames[kk],zn)
    if(length(wx) == 0)next
    if(length(wt) > 0 & length(wx) > 0) pnlist[wt] <- vnames[wx] 
  }
  pnlist
}



initY <- function(znames,posX=NULL,negX=NULL,zeroX=NULL,xx,yy,windex,loB,hiB){
  
  if( length(posX)  > 0) {
  #  pos <- posNegClimPars( pnlist=posX,vnames=posX )
    loB[znames %in% posX,] <- 0
  }
  if( length(negX)  > 0){
  #  neg <- posNegClimPars(negX,c('temp','prec','pdsi'))
    hiB[znames %in% negX,] <- 0  
  }
  if( length(zeroX) > 0){
    loB[znames %in% zeroX,] <- -.01  
    hiB[znames %in% zeroX,] <- .01  
  }
  
  tmp  <- setupBeta(xx=xx[windex,],yy=yy[windex],c(posX,zeroX),c(negX,zeroX))
  beta <- tmp$beta
  sg   <- tmp$sg
  
  res <- yy - xx%*%beta
  res <- vec2Mat(tin[windex,1],tin[windex,2],res[windex],nr=ntree,nc=nyr)
  res <- rowMeans(res,na.rm=T)
  avar <- advar <- var(res,na.rm=T)
  
  list(loB = loB, hiB = hiB, beta = beta, sig = sg, avar = avar)
}


truncatedPriors <- function( posList = NULL, negList = NULL, zeroList = NULL ){

  phi <- loPhi <- hiPhi <- numeric(0)
  
  bg  <- matrix(0,ncol(X),nttype)
  rownames(bg) <- colnames(X)
  colnames(bg) <- treeNames

  loB <- bg*0 - 50
  hiB <- bg*0 + 50
  
  beta <- loBeta <- hiBeta <- numeric(0)
  
  ww  <- which(treeRows())

  if('dinc' %in% ynames){
    
    yy  <- y[,'dinc']
    
    tmp  <- initY(znames,posX=posList$dinc,negX=negList$dinc,
                         zeroX=zeroList$dinc,xx=X,yy,windex=ww,loB,hiB)
    bd   <- bg <- tmp$beta
    if(!'dinc' %in% yrEffect)bd[-pCol,] <- bg[-pCol,] <- 0
    loBd <- tmp$loB
    hiBd <- tmp$hiB
    sdg  <- sg <- tmp$sig/10
    advar <- advar <- tmp$avar
    
    if('years' %in% xnames & 'dinc' %in% yrEffect){
      loBd[yrCol,] <- -.01
      hiBd[yrCol,] <- .01
    }
    
    beta    <- append(beta,list(bd = bg))
    loBeta  <- append(loBeta,list(lo_bd = loBd))
    hiBeta  <- append(hiBeta,list(hi_bd = hiBd))
  }
    
  if(HT){
    
    yy  <- y[,'hinc']
    
    tmp  <- initY(znames,posList$hinc,negList$hinc,
                         zeroList$hinc,X,yy,windex=ww,loB,hiB)
    bh   <- tmp$beta
    if(!'hinc' %in% yrEffect)bh[-pCol,] <- 0
    loBh <- tmp$loB
    hiBh <- tmp$hiB
    shg  <- tmp$sig/10
    ahvar <- tmp$avar
    
    if('years' %in% xnames & 'hinc' %in% yrEffect == 'fecn'){
      loBh[yrCol,] <- -.01
      hiBh[yrCol,] <- .01
    }
    
    beta    <- append(beta,list(bh = bh))
    loBeta  <- append(loBeta,list(lo_bh = loBh))
    hiBeta  <- append(hiBeta,list(hi_bh = hiBh))
  }
  
  if('fecn' %in% ynames){
    
    yy  <- y[,'fecn']
    
    tmp  <- initY(znames,posX=posList$fecn,negX=negList$fecn,
                         zeroX=zeroList$fecn,X,yy,windex=ww,loB,hiB)
    bf   <- tmp$beta
    if(!'fecn' %in% yrEffect)bf[-pCol,] <- 0
    lo_bf <- tmp$loB
    hi_bf <- tmp$hiB
    sfg  <- tmp$sig/10
    afvar <- tmp$avar
    if('diam' %in% xnames)lo_bf['diam',] <- 0
    
    if('years' %in% xnames & 'fecn' %in% yrEffect == 'fecn'){
      lo_bf[yrCol,] <- -.01
      hi_bf[yrCol,] <- .01
    }
    
    beta    <- append(beta,list(bf = bf))
    loBeta  <- append(loBeta,list(lo_bf = lo_bf))
    hiBeta  <- append(hiBeta,list(hi_bf = hi_bf))

    phi   <- matrix(0,ncol(M),1)  #maturation parameters
    nm    <- ncol(M)
    rownames(phi) <- mnames
    loPhi <- rep(-20,nm)
    hiPhi <- rep(20,nm)
    loPhi[mnames %in% posList$fecn] <- 0
    hiPhi[mnames %in% negList$fecn] <- 0
    hiPhi[mnames == 'intercept'] <- -7
    loPhi[mnames == 'diam'] <- .5

    phi <- matrix(glm(Q ~ M[,-1], family=binomial())$coefficients,ncol=1)
    phi[phi < loPhi] <- loPhi[phi < loPhi]
    phi[phi > hiPhi] <- hiPhi[phi > hiPhi]
    wq  <- which(mnames %in% negList$fecn & phi > 0)
    if(length(wq) > 0)phi[wq] <- -.1
    wq  <- which(mnames %in% posList$fecn & phi < 0)
    if(length(wq) > 0)phi[wq] <- .1
    rownames(phi) <- mnames
  }
  
  ssg <- aag <- rep(1,length(ynames))
  sg <- diag(ssg)
  ssg[ynames == 'dinc'] <- sdg
  ssg[ynames == 'fecn'] <- sfg
  aag[ynames == 'dinc'] <- advar
  aag[ynames == 'fecn'] <- afvar
  if(HT){
    ssg[ynames == 'hinc'] <- shg
    aag[ynames == 'hinc'] <- ahvar
  }
  ag <- diag(aag)
  colnames(ag) <- colnames(sg) <- rownames(ag) <- rownames(sg) <- ynames

  list(beta = beta, loB = loBeta, hiB = hiBeta, phi = phi, loPhi = loPhi, 
       hiPhi = hiPhi, sg = sg, ag = ag)
}



getKern <- function(u,dij,tinyKern=0){
  kk <- u/pi/(u + dij^2)^2
  kk[is.na(kk) | kk < tinyKern] <- 0
  Matrix(kk)
}

updateU <- function(minU=4,maxU=500){

  propu <- tnorm(1,minU,maxU,ug,rexp(1,1/4))
  propd <- tnorm(1,0,5,fdg,rexp(1,100))       #LDD

  kernProp <- getKern(propu,distall,tinyKern=1e-8)

  pnow <- sum( getSeedProb3(sparseFecn,sparseKern,seedMat,fdg),na.rm=T) + 
          dnorm(ug,priorU,priorVsd,log=T)
  pnew <- sum( getSeedProb3(sparseFecn,kernProp,seedMat,propd),na.rm=T) +  
          dnorm(propu,priorU,priorVsd,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a){
    ug  <- propu
    fdg <- propd
  }
  list(ug = ug, fdg = fdg)
}

  
getSeedProb <- function(u,fmat,fd=0,wtree=c(1:ntree),wyr=c(1:nyr),GETF=F){  

  # fmat - ntree by nyr matrix of fecundities
  # wtree - which trees to include
  # wyr can be 1 yr or c(1:nyr)
  # fd is coefficient, multiplied by BA (bap), for background seed

  tiny <- 1e-5

  ww <- which(seedmat[,'t'] %in% wyr)
  yy  <- matrix(NA,nrow(seedmat),2)
  colnames(yy) <- c('lambda','lnL')

  if(!is.matrix(fmat))fmat <- as.matrix(fmat)

  fec <- numeric(0)
  if(GETF){
    fec     <- matrix(0,ntree,nyr)
    qmat    <- fec*0
    qmat[tin] <- Q
    diamMax <- apply(diamMat,1,max,na.rm=T)
  }

  bap <- baPlotYr
  bap[is.na(bap)] <- tiny

  for(t in wyr){

     ws <- which(seedmat[,'t'] == t)

     si <- seedmat[ws,'sindex']
     sa <- trapArea[seedmat[ws,'j']]

     w0 <- which(fmat[wtree,t] > 0)
     wt <- wtree[ w0 ]

     la <- tiny + fd*bap[cbind(seedmat[ws,'j'],rep(t,length(ws)))]

     if(length(wt) > 0){
       kern <- getKern(u,distall[si,wt])
       la   <- la + kern%*%col2Mat(fmat[wt,t]) 
     }
     
     yy[ws,'lambda']  <- as.vector(la)
     yy[ws,'lnL'] <- dpois(seedmat[ws,'count'],sa*yy[ws,'lambda'],log=T)

     if(GETF){

        for(j in 1:nplot){

            wtt <- which(treeData[wtree,'j'] == j & diamMax > minMatrDiam)
            if(length(wtt) == 0)next
            wss <- which(seedmat[,'j'] == j & seedmat[,'t'] == t)
            if(length(wss) == 0)next

            stgrid <- as.matrix( expand.grid(seedmat[wss,'sindex'],wtt) )
            dst    <- distall[stgrid]
            osd    <- stgrid[order(dst,decreasing=F),2]
            osd    <- osd[!duplicated(osd)]

            if(length(osd) > 20)osd <- osd[1:15]

            kern <- getKern(u,distall[seedmat[wss,'sindex'],])
            kern <- kern[,osd]

            ft <- invMat( crossprod(kern)  + diag(tiny,length(osd)) ) %*% 
                       crossprod(kern,seedmat[wss,'count']/sa[seedmat[wss,'j']])

            ft[ft < minF] <- minF
            ft[ft > maxF] <- maxF
            fec[osd,t] <- as.vector(ft)
        }
        fec[fec == 0 & qmat == 1] <- minF
     }
  }

  yy <- cbind(seedmat,yy)[ww,]
  list(seedProb = yy, fec = fec)
}

treeRows <- function(surv=1,matr=NULL,spec=NULL,reg=NULL,plt=NULL,years=NULL,sex=NULL){

  ws <- wg <- wm <- wy <- wr <- wp <- rep(T,n)

  if(!is.null(spec)) ws <- treemat[,'spec'] == spec
  if(!is.null(sex))  wg <- female[treemat[,'tindex']] == sex
  if(!is.null(matr)) wm <- Q == matr
  if(!is.null(years))   wy <- treemat[,'t'] %in% years
  if(!is.null(reg))  wr <- treemat[,'reg'] == reg
  if(!is.null(plt))  wp <- treemat[,'j'] == plt
  (is.finite(survVec) & survVec == surv) & ws & wg & wm & wy & wr & wp

}


initQ <- function(){

  q  <- matrix(NA,ntree,nyr)
  q[tin] <- rbinom(n,1,getRho(M,phi))
  q[is.na(q)] <- 0
  qsum <- t( apply(q,1,cumsum) )
  qsum[qsum > 1] <- 1
  qsum[diamMat < minMatrDiam] <- 0
  qsum[diamMat > maxMatrDiam] <- 1
  qsum[tin]

}

initFem <- function(){

  fem    <- vec2Mat(tin[,1],tin[,2],treemat[,'monoec'],nc=nyr)[,nyr]
  fem[sample(ntree,round(ntree/3,0))] <- 1
  fem[isMale] <- 0
  fem[isFem]  <- 1
  fem[is.na(fem)] <- 1
  fem
}

initF <- function(kern,lambda){

  invMat(crossprod(kern))%*%crossprod(kern,count)

}
  

 
getRho <- function(m,phi){ invlogit(m%*%phi) }  #probit maturation

getS <- function(k,s) crossprod(t(k))*s


getDelta <- function(theta1,theta2,dinc1,dinc2,phi,obsVec,femStatus){

  notSeen <- rep(1,ntree)
  w1 <- which(is.finite(obsVec) & obsVec == 0)     #obs to be immature
  if(length(w1) > 0)notSeen[w1] <- 1 - verror*femStatus[w1]
  d1theta <- phi[mnames == 'diam']*dinc1*theta1*(1 - theta1)
  d2theta <- phi[mnames == 'diam']*dinc2*theta2*(1 - theta2)

   d1theta*notSeen/(d1theta*notSeen + d2theta)
}

propFecnNow <- function(fmat){  #propose from current fecn

  ff <- fmat
  ww <- which(is.finite(ff))
  nw <- length(ww)
  
  pv <- log(ff[ww])/50
  pv[!is.finite(pv) | pv < .1] <- .1

  ff[ff == 0] <- minF

  ff[ww] <- exp( tnorm(nw,log(minF),log(maxF),log(ff[ww]),pv) )
  Matrix(ff)
}

propFecnReg <- function(){ #propose from regression

  ll <- Kmat[,2]
  xmat      <- matrix(NA,ntree,nyr)
  xmat[tin] <- (xbdfMat[,2] + ag[tin[,1]])*ll

  ww        <- which(is.finite(xmat))
  nw        <- length(ww)

  fmat[ww]  <- exp( tnorm(nw,log(minF),log(maxF),xmat[ww],sqrt(sg)) )
  fmat
}

propFecnSeed <- function(u,t,ffix,wfix,wpred){  #predict from seed and wfix trees

  tiny  <- 1e-10

  jpred <- treeData[wpred,'j']
  ws <- which(seedmat[,'t'] == t & seedmat[,'j'] %in% jpred)
  sa <- trapArea[seedmat[ws,'j']]

  kern2 <- getKern(u,distall[seedmat[ws,'sindex'],wfix])
  kern1 <- as.matrix( getKern(u,distall[seedmat[ws,'sindex'],wpred]) )

  p2 <- seedmat[ws,'count']/sa - kern2%*%col2Mat(ffix)
  p2[p2 < 0] <- 0

  fpred <- as.numeric( invMat( crossprod(kern1) + diag(tiny,length(wpred)) )%*%t(kern1)%*%p2 )
  fpred[fpred > maxF] <- maxF
  fpred[fpred < minF] <- minF
  fpred
}

getDiamError <- function(){

  u1 <- dd1 + nrow(diamObsIndex)/2
  u2 <- dd2 + .5*sum( (diamMat[diamObsIndex] - diamObsMat[diamObsIndex])^2, na.rm=T )
  1/rgamma(1,u1,u2)

}

XtimesBMVold <- function(X,bd,bf,index=c(1:n)){

  x1 <- rowSums( X[index,]*t( bd[,treemat[index,'spec']] ))
  x2 <- rowSums( X[index,]*t( bf[,treemat[index,'spec']] ),1,sum)
  cbind(x1,x2)
}


getMonod <- function( light,kg,xx=xbdfMat,svar=NULL,index=c(1:n) ){

  Kmat <- getKmat(light,kg)
  mu   <- (xx[index,] + ag[tin[index,1],])*Kmat[index,]

  jj <- 1:ny
  for(j in 1:ny){
    sigmaCols[,jj] <- Kmat[,j]*Kmat*matrix(svar[j,],n,ny,byrow=T)
    jj <- jj + ny
  }

  list(mu = mu, sigmaCols = sigmaCols[index,])
}



updateD <- function(){  #diameter growth
  
  
  wsCol  <- which(ynames == 'dinc')
  wsCol  <- wsCol*ny - (ny - wsCol)

  pnow <- pnew <- rep(0,n)

  accept <- matrix(0,ntree,nyr)

  smat <- vec2Mat(tin[,1],tin[,2],survVec) 
  sindex <- which(is.na(smat),arr.ind=T)  #dead tree years

  D1now <- diamMat[ firstIndex ]
  D1new <- tnorm(ntree,diamFirstTime[,1],diamFirstTime[,2],D1now,.1)  #first yr
  dinew <- tnorm(n,minDincMat[tin],maxDincMat[tin],dincMat[tin],rexp(n,30))

  dimat  <- matrix(NA,ntree,nyr)
  dimat[tin] <- dinew
  DnewMat <- dinc2diam(D1new,firstTime,dimat)
  dimat   <- diam2dinc(diam=DnewMat,minInc=minDinc,firstTime=firstTime)

  dlagNew <- lagGrowthRcpp(dimat,firstTime,lastTime,deathPeriod,deathDiscount)

  xmat <- matrix(0,ntree,nyr)
  
 # if(!DANDF)xmat[tin] <- xbgMat + ag[tin[,1],1]

  sdmat <- matrix(sqrt(sg),ntree,nyr)
  
  ynow <- ynew <- y
  ynow[,'dinc'] <- dincMat[tin]
  ynew[,'dinc'] <- dimat[tin]
  
  wi  <- which( treeRows(matr=0) | treeRows(matr=1,sex=0) )  #immature diameter
  wm  <- which( treeRows(matr=1,sex = 1) )                   #mature
  
  tmp      <- getMonod( light,kg,xbdfMat,svar=svar )  
  pnow[wi] <- dnorm(ynow[wi,'dinc'],tmp$mu[wi,'dinc'],sqrt(tmp$sigmaCols[wi,wsCol]),log=T)
  pnew[wi] <- dnorm(ynew[wi,'dinc'],tmp$mu[wi,'dinc'],sqrt(tmp$sigmaCols[wi,wsCol]),log=T)
  
  p1 <- conditionalNorm(ynow[wm,],tmp$mu[wm,],tmp$sigmaCols[wm,],wsCol)
  mu <- p1$mu
  vr <- p1$vr
  pnow[wm] <- dnorm(ynow[wm,'dinc'],mu,sqrt(vr),log=T)
    
  p1 <- conditionalNorm(ynew[wm,],tmp$mu[wm,],tmp$sigmaCols[wm,],wsCol)
  mu <- p1$mu
  vr <- p1$vr
  pnow[wm] <- dnorm(ynew[wm,'dinc'],mu,sqrt(vr),log=T)
  
  

  pnowMat <- pnewMat <- diamMat*0

  pnowMat[diamObsIndex]   <- dnorm(diamMat[diamObsIndex],diamObsMat[diamObsIndex],
                                   sqrt(sgDiam),log=T) 
  if(DINCDATA)pnowMat[dincObsIndex] <- pnowMat[dincObsIndex] + 
                                       dnorm(dincMat[dincObsIndex],dincObsMat[dincObsIndex],
                                       sqrt(sgDinc),log=T) 
  if(DENDDATA)pnow[dendObsIndex]   <- pnowMat[dendObsIndex] + 
                                      dnorm(dincMat[dendObsIndex],dendObsMat[dendObsIndex],
                                      sqrt(sgDend),log=T) 
  pnowMat[diamPriorIndex] <- pnowMat[diamPriorIndex] + 
                                    dnorm(diamMat[diamPriorIndex],diamLastMu,diamLastSd,log=T) 
  
  pnowMat[tin] <- pnowMat[tin] + pnow

  tmp <- survDincDiam(dlagMat,diamMat)
  liveG <- tmp$lg
  deadG <- tmp$dg
  liveD <- tmp$ld
  deadD <- tmp$dd
  wl    <- tmp$liveIndex
  wd    <- tmp$deadIndex

  tmp <- survLikelihood(dincSurv,diamSurv,liveG,deadG,liveD,deadD)
  pnowMat[wl] <- pnowMat[wl] + tmp$live
  pnowMat[wd] <- pnowMat[wd] + tmp$dead

  
  pnewMat[diamObsIndex]  <- dnorm(DnewMat[diamObsIndex],diamObsMat[diamObsIndex],
                                  sqrt(sgDiam),log=T) 
  if(DINCDATA)pnewMat[dincObsIndex]  <- pnewMat[dincObsIndex] + 
                                        dnorm(dimat[dincObsIndex],dincObsMat[dincObsIndex],
                                        sqrt(sgDinc),log=T) 
  if(DENDDATA)pnewMat[dendObsIndex]  <- pnewMat[dendObsIndex] + 
                                        dnorm(dimat[dendObsIndex],dendObsMat[dendObsIndex],
                                        sqrt(sgDend),log=T) 
  pnewMat[diamPriorIndex] <- pnewMat[diamPriorIndex] + 
                                dnorm(DnewMat[diamPriorIndex],diamLastMu,
                                diamLastSd,log=T) 
  pnewMat[tin] <- pnewMat[tin] + pnew
#  pnew[tin] <- pnew[tin] + dnorm(dimat[tin],xmat[tin],sdmat[tin],log=T)

  tmp <- survDincDiam(dlagNew,diamMat)
  liveG <- tmp$lg
  deadG <- tmp$dg
  liveD <- tmp$ld
  deadD <- tmp$dd
  wl    <- tmp$liveIndex
  wd    <- tmp$deadIndex

  tmp <- survLikelihood(dincSurv,diamSurv,liveG,deadG,liveD,deadD)
  pnewMat[wl] <- pnewMat[wl] + tmp$live
  pnewMat[wd] <- pnewMat[wd] + tmp$dead

  pnowMat[sindex] <- pnewMat[sindex] <- 0

  pnow <- rowSums(pnowMat,na.rm=T)
  pnew <- rowSums(pnewMat,na.rm=T)

  a <- exp(pnew - pnow)
  z <- runif(ntree,0,1)
  diamMat[z < a,] <- DnewMat[z < a,]
  dincMat[z < a,] <- dimat[z < a,]
  accept[z < a,]  <- 1

  list(diamMat = diamMat, dincMat = dincMat, accept = accept)
}

size2inc <- function(sizeMat,minInc=0){

  sm <- sizeMat

  if(!is.matrix(sm)){
    sizeMat <- matrix(NA,ntree,nyr)
    sizeMat[tin] <- sm
    colnames(sizeMat) <- yrvec
  }

  slast <- lastTimeMat(sizeMat,firstTime,INC=T,minInc)

  xl <- sizeMat - slast
  xx <- cbind(xl[,-1],xl[,ncol(xl)])
  xx[xx < minInc] <- minInc
  colnames(xx) <- colnames(sizeMat)
  xx
}

inc2size <- function(firstSize,firstTime=rep(1,length(firstDiam)),incMat){    
   #firstDiam is vector of diams at firstTime

  nr   <- length(firstSize)
  nyr  <- ncol(incMat)

  dtmp <- incMat
  dtmp[is.na(dtmp)] <- 0
  dtmp <- cbind(rep(0,nr),dtmp)
  dtmp[ cbind(c(1:nr),firstTime) ] <- firstSize

  dtmp <- t( apply(dtmp,1,cumsum) )
  dtmp <- dtmp[,-(nyr+1)]
  colnames(dtmp) <- yrvec
  dtmp
}




dinc2diam <- function(firstDiam,firstTime=rep(1,length(firstDiam)),dincMat){    
   #firstDiam is vector of diams at firstTime

  nr   <- length(firstDiam)
  nyr  <- ncol(dincMat)

  dtmp <- dincMat
  dtmp[is.na(dtmp)] <- 0
  dtmp <- cbind(rep(0,nr),dtmp)
  dtmp[ cbind(c(1:nr),firstTime) ] <- firstDiam

  dtmp <- t( apply(dtmp,1,cumsum) )
  dtmp <- dtmp[,-(nyr+1)]
  dtmp[dtmp < 0] <- 0
  colnames(dtmp) <- yrvec
  dtmp
}


propFecnSeed2 <- function(wfix,wprop){  #predict from seed and wfix trees

  tiny  <- 1e-10

  kern2 <- sparseKern[,wfix]
  kern1 <- sparseKern[,wprop]
  p2    <- seedMat/trAreaMat - kern2%*%sparseFecn[wfix,] - fdg*BASeedMat
  p2[p2 < 0] <- 0
  p2[is.na(p2)] <- 0

  fpred <-  invMat( crossprod(kern1) + diag(tiny,length(wprop))) %*%t(kern1)%*%p2 
  fpred[fpred > maxF] <- maxF
  fpred[fpred < minF] <- minF

  tnorm(length(fpred),minF,maxF,as.vector(fpred),10)
}


getDelta2 <- function(theta1,theta2,dinc1,dinc2,phi,femStatus){

  tiny <- 10e-10
  huge <- 1 - tiny

  d1theta <- phi[mnames == 'diam']*dinc1*theta1*(1 - theta1)
  d2theta <- phi[mnames == 'diam']*dinc2*theta2*(1 - theta2)

  if(length(obsImm) > 0)d1theta[obsImm] <- d1theta[obsImm]*(1 - verror*femStatus[obsImm[,1]])

  delta <- d1theta/(d1theta + d2theta)
  delta[delta < tiny] <- tiny

 # delta[delta < huge & (d1theta + d2theta) < tiny] <- tiny

  delta[delta > huge] <- huge
  delta
}

updateFemale <- function(){

  q1 <- rowSums(qmat,na.rm=T)
  p1 <- sum(female[q1 > 0])
  p2 <- sum(1 - female[q1 > 0])

  rbeta(1,prFemA + p1,prFemB + p2)

}


updateF2 <- function(){

  tiny <- 1e-10
  gFraction <- .1
  
  wfe <- which(ynames == 'fecn')

  accept <- pnow <- pnew <- mmat <- censorMat*0

  wmin <- which( (is.na(fmat) & qmat == 1) | fmat == 0 & qmat == 1 )

  if(length(wmin) > 0){
    fmat[wmin] <- minF
    sparseFecn[wmin] <- minF
  }

  # propose gender
  PROPG <- 0
  if(DIOEC)PROPG <- rbinom(1,1,.3)  # PROPG == 0 means retain current gender

  femNew <- female
  if(PROPG == 1){
    ss <- sample(1:ntree,round(ntree*gFraction,0))
    femNew[ss] <- rbinom(length(ss),1,.5)
  }

  femNew[isMale] <- 0
  femNew[isFem]  <- 1
  femNew[treeData[,'monoec'] == 1] <- female[treeData[,'monoec'] == 1] <- 1
  femNewMat <- matrix(femNew,ntree,nyr)

  #propose maturation time
  qnew <- qmat <- qmat*censorMat

  firstOne <- qmat - lastTimeMat(qmat,firstTime)  #first imputed mature
  firstOne <- which(firstOne == 1,arr.ind=T)
  lastZero <- firstOne
  lastZero[,2] <- lastZero[,2] - 1

  lastZero[lastZero[,2] < 1,2] <- 1

  ncs <- nrow(firstOne)
  change <- sample(ncs,ncs/2)

  qnew[firstOne[change,]] <- 0              #half immature
  qnew[lastZero[c(1:ncs)[-change],]] <- 1   #half mature

  #starts mature
  mstart  <- which(qmat[ cbind(1:ntree,firstTime) ] == 1)
  qnew[ cbind(mstart,firstTime[mstart])] <- sample(c(0,1),length(mstart),replace=T)

  #ends immature
  iend  <- which(qmat[ cbind(1:ntree,lastTime) ] == 0)
  qnew[ cbind(iend,lastTime[iend]) ] <- sample(c(0,1),length(iend),replace=T)

  qmat[indexMature]   <- qnew[indexMature] <- 1
  qmat[indexImmature] <- qnew[indexImmature] <- 0

  qnowVec <- qmat[tin]
  qnewVec <- qnew[tin]

  wmNow <- which( treeRows() & female[treemat[,'tindex']] == 1 & qnowVec == 1 )
  wmNew <- which( treeRows() & femNew[treemat[,'tindex']] == 1 & qnewVec == 1 )

  tiNow <- tin[wmNow,]
  tiNew <- tin[wmNew,]

  #propose fecundity

  closeStat   <- rowSums(qmat[closeTreesAll,],na.rm=T)*femNew[closeTreesAll]
  sampTrees   <- closeTreesAll[closeStat > 0]
  sampleClose <- sample(sampTrees,length(sampTrees)/2)

  fnew     <- propFecnNow(fmat*qnew*femNewMat)
  propSeed <- propFecnSeed2(wfix=c(1:ntree)[-sampleClose],wprop=sampleClose)

  fnew[sampleClose,] <- propSeed
  fnew <- fnew*qnew*femNewMat
  fnew[is.na(fnew)] <- 0
  fnew[fnew < minF & qnew == 1] <- minF

  # contribution from regression
  tmp <- getMonod( light,kg,xbdfMat,svar=svar )
  tmp <- conditionalNorm(y,tmp$mu,tmp$sigmaCols,wfe)
  mu  <- tmp$mu
  vr  <- tmp$vr

  pnow[tiNow] <- dnorm(log(fmat[tiNow]),mu[wmNow],sqrt(vr[wmNow]),log=T)
  pnew[tiNew] <- dnorm(log(fnew[tiNew]),mu[wmNew],sqrt(vr[wmNew]),log=T)

  mmat[tin] <- getRho(M,pg)
  mmat2     <- nextTimeMat(mmat,lastTime)
  dincNext  <- nextTimeMat(dincMat,lastTime)
  deltaNow  <- getDelta2(mmat,mmat2,dincMat,dincNext,pg,female)
  deltaNew  <- getDelta2(mmat,mmat2,dincMat,dincNext,pg,femNew)

  matnow <- qmat*log(deltaNow) + (1 - qmat)*log(1 - deltaNow)
  matnew <- qnew*log(deltaNew) + (1 - qnew)*log(1 - deltaNew)

  matnow[is.na(matnow)] <- 0
  matnew[is.na(matnew)] <- 0

  pnow <- pnow + matnow
  pnew <- pnew + matnew

  ptnow <- sumBy2Index(pnow[ tin ],treemat[,'j'],treemat[,'t'],nplot,nyr)
  ptnew <- sumBy2Index(pnew[ tin ],treemat[,'j'],treemat[,'t'],nplot,nyr)

  snow <- getSeedProb3(sparseFecn,sparseKern,seedMat,fdg)
  snow <- sumBy2Index(snow[ seedmat[,c('sindex','t') ] ],seedmat[,'j'],seedmat[,'t'],nplot,nyr)

  snew <- getSeedProb3(Matrix(fnew),sparseKern,seedMat,fdg)
  snew <- sumBy2Index(snew[ seedmat[,c('sindex','t') ] ],seedmat[,'j'],seedmat[,'t'],nplot,nyr)

  pnowSum <- ptnow + seedwt*snow
  pnewSum <- ptnew + seedwt*snew

  if(PROPG == 0){  
     a <- exp(pnewSum - pnowSum)
     z <- matrix(runif(nyr*nplot,0,1),nplot,nyr)
     ww <- which(z < a,arr.ind=T)
     wt <- which(treemat[,'j'] %in% ww[,1] & treemat[,'t'] %in% ww[,2])
  }
  if(PROPG == 1){

     q1 <- rowSums(qmat,na.rm=T)
     q2 <- rowSums(qnew,na.rm=T)

     anow <- rowSums(pnowSum,na.rm=T) + dbinom(sum(female[q1 > 0]),ntree,prFemg,log=T)
     anew <- rowSums(pnewSum,na.rm=T) + dbinom(sum(femNew[q2 > 0]),ntree,prFemg,log=T)
     a    <- exp(anew - anow)
     z    <- runif(nplot,0,1)
     ww <- which(z < a,arr.ind=T)
     wt <- which(treemat[,'j'] %in% ww)
  }
  if(length(wt) > 0){
    fmat[tin[wt,]] <- fnew[tin[wt,]]
    qmat[tin[wt,]] <- qnew[tin[wt,]]
    female[treeData[,'j'] %in% ww] <- femNew[treeData[,'j'] %in% ww]
    accept[tin[wt,]] <- 1
  }

  list(fmat = Matrix(fmat),fecn = fmat[tin], qmat = qmat, 
       Q = qmat[tin], female = female, accept = accept)
}


updateF <- function(){  # fecundity and maturation

  tiny <- 1e-10
  gFraction <- .1  #fraction to propose new gender

  PROPG <- 0
  if(DIOEC)PROPG <- rbinom(1,1,.8)  # PROPG == 0 means retain current gender
  PROPREG <- rbinom(1,1,.2)

  femNew <-rbinom(ntree,1,.5)
  femNew[isMale] <- 0
  femNew[isFem]  <- 1

  femNew[treeData[,'monoec'] == 1] <- female[treeData[,'monoec'] == 1] <- 1

  propNow <- propFecnNow()
  propReg <- propFecnReg()

  if(PROPG == 0)femNew <- female
  if(PROPG == 1){
    ss <- sample(1:ntree,round(ntree*(1 - gFraction),0))
    femNew[ss] <- female[ss]
  }

  fecn[fecn < minF] <- minF

  qmat <- fmat <- mmat <- xmat <- smat <- matrix(NA,ntree,nyr)
  qmat[tin] <- Q
  fmat[tin] <- fecn
  xmat[tin] <- xbdfMat[,2] + ag[tin[,1],1]
  sdMat <- matrix(sqrt(sg),ntree,nyr)

  wm <- which(treeRows(matr=1))

  if(FONLY){

    xmatf <- xmat
    tmp <- getMonod( light,kg,bf,bf,svar=sg,index=wm )
    mu  <- tmp$mu
    vr  <- tmp$sigma
    xmat[tin[wm,]] <- mu
    sdMat[tin[wm,]] <- sqrt(vr)
  }

  if(DANDF){
    xmatd <- xmatf <- xmat

    xx <- cbind(dincMat[tin], log(fmat[tin]) )
    xx[xx[,2] == -Inf,2] <- log(1)

    tmp <- getMonod( light,kg,bd,bf,svar=svar,index=wm )
    tmp <- conditionalBiVarNorm(xx[wm,],tmp$mu,tmp$sigma,2)
    mu  <- tmp$mu
    vr  <- tmp$vr
    xmat[tin[wm,]] <- mu
    sdMat[tin[wm,]] <- sqrt(vr)
    
  }

  mmat[tin] <- getRho(M,pg)
  smat[tin] <- survVec

  

  fmat[diamMat < minMatrDiam] <- 0
  qmat[diamMat < minMatrDiam] <- 0
  qmat[diamMat > maxMatrDiam] <- 1

  fnewAll <- fmat
  qnewAll <- qmat

  pnowTot <- pnewTot <- rep(0,nplot)

  accept <- matrix(0,ntree,nyr)

  seedPlotYr <- matrix(NA,nplot,nyr)
  
  for(t in 1:nyr){

     stayImm <- which(lastImm >= t & smat[,t] == 1)
     stayMat <- which(firstMat <= t & smat[,t] == 1)

     fmat[firstTime > t,t] <- NA
     qmat[firstTime > t,t] <- NA

     qnow <- qnew <- qmat[,t]
     fnow <- fnew <- fmat[,t]

     fnow[is.na(fnow) & qnow == 1] <- minF

     wi   <- which(is.finite(fnow) & firstTime <= t & smat[,t] == 1)
     if(length(wi) == 0)next
     wobs <- which(sexObs[,t] == 0)    #obs not mature

     m1 <- mmat[,t]
     d1 <- dincMat[,t]

     if(t == 1){
         wq <- which(qmat[,t+1] == 1)
         m2 <- mmat[,t+1]
         d2 <- dincMat[,t+1]
     }
     if(t == nyr){
         wq <- which(qmat[,t-1] == 0 & diamMat[,t] > minMatrDiam)
         m2 <- mmat[,t]
         d2 <- dincMat[,t]
     }
     if(t > 1 & t < nyr){
         wq <- which(qmat[,t-1] == 0 & qmat[,t+1] == 1 & diamMat[,t] > minMatrDiam)
         m2 <- mmat[,t+1]
         d2 <- dincMat[,t+1]
     }

     wq <- wq[wq %in% wi]     #make transition (t-1, t+1)

     obsVec <- sexObs[,t]

     if(length(wq) > 0)qnew[wq] <- rbinom(length(wq),1,.5)     #propose status transition trees

     if(length(stayImm) > 0)qnow[stayImm] <- qnew[stayImm] <- 0
     if(length(stayMat) > 0)qnow[stayMat] <- qnew[stayMat] <- 1

     wm <- which(qnow == 1)
     fnow[is.finite(qnow) & qnow == 1 & fnow < minF] <- minF
##############################
     fnew[wi] <- propNow[wi,t]
     
     wprop <- which(qnew == 1 & femNew == 1 & closeTrap > 30)  #far trees
     if(length(wprop) > 0 & PROPREG == 1)fnew[wprop] <- propReg[wprop,t]

     wprop <- which(qnew == 1 & femNew == 1 & closeTrap < 15)  #close trees
     if(length(wprop) > 0){
       wfix  <- which(qnew == 1 & femNew == 1 & closeTrap > 15)
       ffix  <- matrix(propReg[wfix,t],ncol=1)
       fnew[wprop] <- propFecnSeed(ug,t,ffix,wfix,wprop)
     }

     fnewAll[,t] <- fnew
     qnewAll[,t] <- qnew
##############################
     wnow <- which(qnow == 1 & female == 1)
     wnew <- which(qnew == 1 & femNew == 1)

     pnow <- pnew <- rep(0,ntree)

     pnow[wnow] <- dnorm(log(fnow[wnow]),xmat[wnow,t],sdMat[wnow,t],log=T) * .1*closeTrap[wnow]
     pnew[wnew] <- dnorm(log(fnew[wnew]),xmat[wnew,t],sdMat[wnew,t],log=T) * .1*closeTrap[wnew]

     deltaNow <- getDelta(m1,m2,d1,d2,pg,obsVec,female)
     deltaNew <- getDelta(m1,m2,d1,d2,pg,obsVec,femNew)

     matnow <- qnow*log(deltaNow) + (1 - qnow)*log(1 - deltaNow)
     matnew <- qnew*log(deltaNew) + (1 - qnew)*log(1 - deltaNew)

     matnow[is.na(matnow)] <- 0
     matnew[is.na(matnew)] <- 0

     pnow <- pnow + matnow
     pnew <- pnew + matnew

     ptnow <- sumByIndex(pnow,treeData[,'j'],nplot)
     ptnew <- sumByIndex(pnew,treeData[,'j'],nplot)

     fmat[,t] <- fnow

     psnow <- getSeedProb(ug,fmat,fdg,wtree=wnow,wyr=t)$seedProb
     psnow <- sumByIndex(psnow[,'lnL'],psnow[,'j'],nplot)
     seedPlotYr[,t] <- psnow
     posum <- ptnow + seedwt*psnow

     fmnew <- fmat
     fmnew[,t] <- fnew

     psnew <- getSeedProb(ug,fmnew,fdg,wtree=wnew,wyr=t)$seedProb
     psnew <- sumByIndex(psnew[,'lnL'],psnew[,'j'],nplot)
     pesum <- ptnew + seedwt*psnew

     if(PROPG == 0){
       a <- exp(pesum - posum)
       z <- runif(nplot,0,1)
       wz <- which(z < a)
       wz <- which(treeData[,'j'] %in% wz)

       fmat[wz,t] <- fnew[wz]
       qmat[wz,t] <- qnew[wz]
       accept[wz,t] <- 1
     }
     if(PROPG == 1){
       pnowTot <- pnowTot + posum
       pnewTot <- pnewTot + pesum
     }
   }

  if(PROPG == 1){
       a <- exp(pnewTot - pnowTot)
       z <- runif(nplot,0,1)
       wz <- which(z < a)
       wz <- which(treeData[,'j'] %in% wz)

       fmat[wz,] <- fnewAll[wz,]
       qmat[wz,] <- qnewAll[wz,]
       female[wz] <- femNew[wz]
       accept[wz,] <- 1
  }

  list(fecn = fmat[tin], Q = qmat[tin], female = female, accept = accept)
}

updateKappa <- function(){

     kprop <- tnorm.mvtRcpp(kg,kg,propK,lo=c(0,0),hi=c(2,2),times=1)

     wm   <- treeRows(matr=1)
     know <- getKmat(light,kg)
     knew <- getKmat(light,kprop)

     pnow <- sum( dnorm(y[wm,1],z[wm,1]*know[wm,1],sqrt(rg[1]),log=T) ) + 
             sum( dnorm(y[wm,2],z[wm,2]*know[wm,2],sqrt(rg[4]),log=T) ) 
     pnew <- sum( dnorm(y[wm,1],z[wm,1]*knew[wm,1],sqrt(rg[1]),log=T) ) + 
             sum( dnorm(y[wm,2],z[wm,2]*knew[wm,2],sqrt(rg[4]),log=T) ) 

     wi  <- which( treeRows(matr=0) | treeRows(matr=1,sex=0) )

     pnow <- pnow + sum(dnorm( y[wi,1],z[wi,1]*know[wi,1],sqrt(rg[1]),log=T))
     pnew <- pnew + sum(dnorm( y[wi,1],z[wi,1]*knew[wi,1],sqrt(rg[1]),log=T))

     a <- exp(pnew - pnow)
     z <- runif(1,0,1)
     if(z < a)kg <- kprop
     kg
}
  

updateLight <- function(){

  wm <- which( treeRows(matr=1,sex=1,years=c(1:(nyr-1))) )    #mature and female
  wi <- which( treeRows(matr=0,years=c(1:(nyr-1))) |  
                  treeRows(sex=0,years=c(1:(nyr-1))) )

  propL <- tnorm(n,loL,hiL,light,light*(1 - light))

  pnow <- pnew <- propL*0

   tmp <- getMonod(light,kg,xbdfMat,svar)
   mu  <- tmp$mu
   vr  <- tmp$sigmaCols
   if(ny == 2)pnow[wm] <- dbivarNormFromCols(y[wm,],mu[wm,],vr[wm,])
   if(ny == 3)pnow[wm] <- dtrivarNormFromCols(y[wm,],mu[wm,],vr[wm,])
   if(nyi == 1)pnow[wi] <- dnorm(y[wi,wyi],mu[wi,wyi],sqrt(vr[wi,wyi]),log=T )
   if(nyi == 2)pnow[wi] <- dbivarNormFromCols(y[wi,wyi],mu[wi,wyi],vr[wi,wyi]) 

  tmp  <- getMonod(propL,kg,xbdfMat,svar)
   mu  <- tmp$mu
   vr  <- tmp$sigmaCols
   if(ny == 2)pnew[wm] <- dbivarNormFromCols(y[wm,],mu[wm,],vr[wm,])
   if(ny == 3)pnew[wm] <- dtrivarNormFromCols(y[wm,],mu[wm,],vr[wm,])
   if(nyi == 1)pnew[wi] <- dnorm(y[wi,wyi],mu[wi,wyi],sqrt(vr[wi,wyi]),log=T )
   if(nyi == 2)pnew[wi] <- dbivarNormFromCols(y[wi,wyi],mu[wi,wyi],vr[wi,wyi]) 

  a <- exp(pnew - pnow)
  z <- runif(n,0,1)

  light[z < a] <- propL[z < a]
  light
}


getKappa <- function(){  #before z

  wm <- which( treeRows(matr=1,sex=1,years=c(1:(nyr-1))) )    #mature and female
  wi <- which( treeRows(matr=0,years=c(1:(nyr-1))) |  
                  treeRows(sex=0,years=c(1:(nyr-1))) )

   pk <- diag( (kg - loK)/(hiK - loK) )*propK

   kprop <- tnorm.mvtRcpp(kg,kg,pk,loK,hiK,times=1)

   pnow <- pnew <- light*0
   
   tmp <- getMonod(light,kg,xbdfMat,svar)
   mu  <- tmp$mu
   vr  <- tmp$sigma
   if(ny == 2)pnow[wm]  <- dbivarNormFromCols(y[wm,],mu[wm,],vr[wm,]) 
   if(ny == 3)pnow[wm]  <- dtrivarNormFromCols(y[wm,],mu[wm,],vr[wm,]) 
   if(nyi == 1)pnow[wi] <- dnorm(y[wi,wyi],mu[wi,wyi],sqrt(vr[wi,wyi]),log=T )
   if(nyi == 2)pnow[wi] <- dbivarNormFromCols(y[wi,wyi],mu[wi,wyi],vr[wi,wyi]) 
   pnow <- sum(pnow) + sum(dnorm(kg,c(kmu,kmu),sqrt(kvar),log=T) )
   
   tmp <- getMonod(light,kprop,xbdfMat,svar)
   mu  <- tmp$mu
   vr  <- tmp$sigma
   if(ny == 2)pnew[wm]  <- dbivarNormFromCols(y[wm,],mu[wm,],vr[wm,]) 
   if(ny == 3)pnew[wm]  <- dtrivarNormFromCols(y[wm,],mu[wm,],vr[wm,]) 
   if(nyi == 1)pnew[wi] <- dnorm(y[wi,wyi],mu[wi,wyi],sqrt(vr[wi,wyi]),log=T )
   if(nyi == 2)pnew[wi] <- dbivarNormFromCols(y[wi,wyi],mu[wi,wyi],vr[wi,wyi]) 
   pnew <- sum(pnew) + sum(dnorm(kprop,c(kmu,kmu),sqrt(kvar),log=T))
   
  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)kg <- kprop
  kg
}


bgMatrix <- function(x){

  cnames <- c('bd','bf','bjuv')

  meanList <- sdList <- numeric(0)

  kj <- 0

  for(k in 1:length(cnames)){

     bmeans <- bsds <- numeric(0)

     wk <- grep(paste(cnames[k],'_',sep=''),rownames(x))
     if(length(wk) == 0)next

     kj <- kj + 1

     xk <- x[wk,]

     if(nttype > 1){

       for(kk in 1:nttype){
         wkk <- grep(treeNames[kk],rownames(xk))
         xkk <- xk[wkk,]
         bmeans <- cbind(bmeans,xkk[,1])
         bsds   <- cbind(bsds,xkk[,2])
       }
       colnames(bmeans) <- colnames(bsds) <- treeNames
       rownames(bmeans) <- znames

     }

     mn <- paste(cnames[k],'Means',sep='_')
     ms <- paste(cnames[k],'Sds',sep='_')
     assign( mn,bmeans)
     assign( ms,bsds)
     meanList <- append(meanList,list(get(mn)))
     sdList <- append(sdList,list(get(ms)))
     names(meanList)[kj] <- names(sdList)[kj] <- cnames[k]
  }

  list(means = meanList, sds = sdList)
}

updateBetaUV <- function(){          #growth of juveniles
  
  loB <- loBd
  hiB <- hiBd
  
  ws <- which(ynames == 'dinc')
  savar <- svar[ws,ws] + avar[ws,ws]
   
  for(k in 1:nttype){

    if(treeBySpec[k] < minTrees | length(ww) <= nrow(bjuv))next
    
    wi <- which( treeRows(spec=k,matr=0,years=c(1:(nyr-1))) | 
                 treeRows(spec=k,matr=1,sex=0,years=c(1:(nyr-1))))  
#    wm <- which( treeRows(spec=k,years=c(1:(nyr-1))) )   
    
    b1 <- bjuv[pCol,k]                        #fixed effects
    y1 <- y[wi,'dinc']
    y1 <- y1/Kmat[wi,'dinc']
    
    
    if(length(yrCol) > 0){
      if(!'dinc' %in% yrEffect | length(wi) < nyr*2)bjuv[yrCol,] <- 0
    }
    
    
    if(length(yrCol) > 0 & 'dinc' %in% yrEffect){
      
      y2 <- y1 - X[wi,pCol]%*%b1                       #remove fixed effects
      wc <- which(colSums(X[wi,yrCol]) > 0)
      yc <- yrCol[wc]
      b2 <- bjuv[yc,k]
      
      b2 <- bUpdateNorm(X[wi,yc],y2,b2,
                    priorB[yc,k],priorIVB[yc,yc],
                    lo=loB[yc,k],hi=hiB[yc,k],savar)
      
      b2 <- b2 - mean(b2)                              #sum to zero
      bjuv[yc,k] <- b2
      y1 <- y1 - X[wi,yc]%*%b2                       #remove yr effects
    }         
    
    b1     <- bUpdateNorm(X[wi,pCol],y1,b1,
                      priorB[pCol,k],priorIVB[pCol,pCol],
                      lo=loB[pCol,k],hi=hiB[pCol,k],savar)
    bjuv[pCol,k] <- b1
    
  }
  bjuv
}
    
    


updateMissX <- function(){          #does not yet include interactions

  missColX <- unique(missingX[,2])

  if(DANDF) sinv <- invMat(svar)
  if(!DANDF)sinv <- 1/sg

  wm <- which( treeRows(matr=1,sex=1) )    #mature and female
  wi <- which( treeRows(matr=0) )   

  z  <- y

  z <- y2z( y, Kmat )

  for(k in missColX){

    missK <- missingX[missingX[,2] == k,]

    Vi <- rep(1/priorXV[k],n)
    v  <- priorX[k]*Vi

    bk <- cbind( bd[k,treemat[,'spec']],bf[k,treemat[,'spec']])            #mature

    qk <- z[wm,] - XtimesBMV(X[,-k],bd[-k,],bf[-k,],index=wm) - ag[tin[wm,1],]

    bsx <- cbind( bk[wm,1]*sinv[1] + bk[wm,2]*sinv[2], bk[wm,1]*sinv[2] + bk[wm,2]*sinv[4] )

    Vi[wm] <- Vi[wm] + bsx[,1]*bk[wm,1] + bsx[,2]*bk[wm,2]
    v[wm]  <- v[wm]  + bsx[,1]*qk[,1] + bsx[,2]*qk[,2]

   #immature
    qk     <- z[wi,1] - XtimesB(X[,-k],bjuv[-k,])[wi] - ag[tin[wi,1],1]
    Vi[wi] <- Vi[wi] + bjuv[k,treemat[wi,'spec']]^2/sg
    v[wi]  <- v[wi]  + bjuv[k,treemat[wi,'spec']]*qk/sg

    mu <- (v/Vi)[missK[,1]]
    vr <- 1/Vi[missK[,1]]

    X[missK] <- tnorm(nrow(missK),priorXrange[1,k],priorXrange[2,k],mu,sqrt(vr))

  }

  X
}



updateBetaMVOld <- function(ylist=ynames,MATR=T){
  
  if(!MATR)ylist <- ylist[ylist != 'fecn']
  wj    <- match(ylist,ynames)
  jlist <- blist[wj]
  nj    <- length(wj)
  if(!MATR)jlist <- paste(jlist,'juv',sep='')
  outlist <- vector('list',length=nj)
  names(outlist) <- jlist

  bdd <- bff <- bd*0
  
  bk <- matrix(0,nrow(bd),nj)
  colnames(bk) <- ynames[wj]
  lo <- hi <- bk
  
  wye <- match(yrEffect,ynames)
  wye <- wye[ynames[wye] %in% ylist]
  nye <- length(wye)

  for(k in 1:nttype){

    if(MATR) wm <- which( treeRows(spec=k,matr=1,sex=1,years=c(1:(nyr-1))) )    #mature and female
    if(!MATR)wm <- which( treeRows(spec=k,matr=0,years=c(1:(nyr-1))) | 
                          treeRows(spec=k,sex=0,years=c(1:(nyr-1))) )
      
    if(length(wm) <= nrow(bf))next

    y1 <- y2z( y[wm,], Kmat[wm,] )

    
    for(kk in wj){
      bk[,ylist[kk]] <- get(jlist[kk])[,k]
      if(nye == 0)bk[-pCol,] <- 0
      if(nye > 0){
        if(wye != kk)bk[-pCol,] <- 0
      }
      lo[,ylist[kk]] <- get(paste('lo',blist[kk],sep='_'))[,k]
      hi[,ylist[kk]] <- get(paste('hi',blist[kk],sep='_'))[,k]
    }
    
    b1 <- matrix(bk[pCol,],ncol=nj)

    if(nye > 0){

       wc <- which(colSums(X[wm,yrCol]) > 0)
       yc <- yrCol[wc]

       y2 <- y1[,wj] - X[wm,pCol]%*%b1                #remove fixed effects
       b2 <- bk[yc,] 
       
       if( nye > 1 ){
         b2[,wye] <- bUpdateMVN_Rcpp(X[wm,yc],y2[,wye],b2[,wye],
                                     lo[yc,wye],hi[yc,wye],(svar+avar)[wye,wye])
       }
       if(nye == 1){
         b2[,wye] <- bUpdateNorm(X[wm,yc],matrix(y2[,wye]),b2[,wye],
                              lo=lo[yc,wye],hi=hi[yc,wye],sigma=(svar+avar)[wye,wye])
       }
         
       b2 <- b2 - matrix( colMeans(b2),length(yc),nj,byrow=T )
       bk[yc,] <- b2
       
       
       y1[,wj] <- y1[,wj] - X[wm,yc]%*%b2  #remove yr effects
    }    
    
    if(nj == 1)b1 <- bUpdateNorm(X[wm,pCol],matrix(y1[,wj]),b1,
                             lo=lo[pCol,],hi=hi[pCol,],sigma=(svar+avar)[wj,wj])
    
    if(nj > 1) b1 <- bUpdateMVN_Rcpp(X[wm,pCol],y1[,wj],b1,
                                     lo[pCol,],hi[pCol,],(svar+avar)[wj,wj])
    bk[pCol,] <- b1
    
    if('bd' %in% jlist | 'bdjuv' %in% jlist)bd[,k] <- bdjuv[,k] <- bk[,'dinc']
    if('bh' %in% jlist | 'bhjuv' %in% jlist)bh[,k] <- bhjuv[,k] <- bk[,'hinc']
    if('bf' %in% jlist)bf[,k] <- bk[,'fecn']
  }
  
  for(k in 1:length(jlist))outlist[[k]] <- get(jlist[k]) 
  
  outlist
}


updateBetaMV <- function(MATR=T){
  
  ylist <- ynames
  if(!MATR)ylist <- ylist[ylist != 'fecn']
  wj    <- match(ylist,ynames)
  jlist <- blist[wj]
  nj    <- length(wj)         #NO. RESPONSES
  if(!MATR)jlist <- paste(jlist,'juv',sep='')
  outlist <- vector('list',length=nj)
  names(outlist) <- jlist

  bdd <- bff <- bd*0
  

  bk <- matrix(bd[,1],length(xnames),ncol=nj)
  rownames(bk) <- rownames(bd)

  colnames(bk) <- ynames[wj]
  lo <- hi <- bk

  
  wye <- match(yrEffect,ynames)       #responses that have yr effects
  wye <- wye[ynames[wye] %in% ylist]
  nye <- length(wye)                  #no. responses having yr effects

  for(k in 1:nttype){

    if(MATR) wm <- which( treeRows(spec=k,matr=1,sex=1,years=c(1:(nyr-1))) )    #mature and female
    if(!MATR)wm <- which( treeRows(spec=k,matr=0,years=c(1:(nyr-1))) | 
                          treeRows(spec=k,sex=0,years=c(1:(nyr-1))) )

    if(length(wm) <= nrow(bf))next

    z <- y2z( y[wm,], Kmat[wm,] )

    for(kk in wj){

      bk[,ylist[kk]] <- get(jlist[kk])[,k]
      if(nye == 0)bk[-pCol,] <- 0
      if(!kk %in% wye)bk[-pCol,] <- 0
      lo[,ylist[kk]] <- get(paste('lo',blist[kk],sep='_'))[,k]
      hi[,ylist[kk]] <- get(paste('hi',blist[kk],sep='_'))[,k]
    }
    
    b1 <- bk[pCol,]
    if(!is.matrix(b1)){
       b1 <- matrix(bk[pCol,],ncol=nj)
       rownames(b1) <- rownames(bk)[pCol]
    }

    if(nye > 0){

       wc <- which( colSums(X[wm,yrCol]) > 0 )
       yc <- yrCol[wc]

       y2 <- z[,wj] - X[wm,pCol]%*%b1                #remove fixed effects
       b2 <- bk[yc,] 
       if(!is.matrix(b2))b2 <- matrix(b2,ncol=1)
       
       if( nye > 1 ){
         b2[,wye] <- bUpdateMVN_Rcpp(X[wm,yc],y2[,wye],b2[,wye],
                                     lo[yc,wye],hi[yc,wye],(svar+avar)[wye,wye])
       }
       if(nye == 1){
         b2[,wye] <- bUpdateNorm(X[wm,yc],matrix(y2[,wye]),b2[,wye],
                              lo=lo[yc,wye],hi=hi[yc,wye],sigma=(svar+avar)[wye,wye])
       }
         
       b2 <- b2 - matrix( colMeans(b2),length(yc),nj,byrow=T )
       bk[yc,] <- b2
       
       z[,wj] <- z[,wj] - X[wm,yc]%*%b2  #remove yr effects
    }    
    
  #  if(nj == 1)b1 <- bUpdateMVN_Rcpp(xx=X[wm,pCol],yy=z[,wj],bb=b1,
  #                                   lo=lo[pCol,],hi=hi[pCol,],sigma=(svar+avar)[wj,wj])

    
    b1 <- updateBetaMVN(xx=X[wm,pCol],yy=z[,wj],bb=b1,lo=lo[pCol,],hi=hi[pCol,],
                                     sigma=(svar+avar)[wj,wj],
                                      INT=INT,isIntPrior=isIntPriorX,intPrior=intPriorX,
                                      noIntPrior=noIntPriorX,SIGMA=F)$beta
    bk[pCol,] <- b1
    
    if('bd' %in% jlist | 'bdjuv' %in% jlist)bd[,k] <- bdjuv[,k] <- bk[,'dinc']
    if('bh' %in% jlist | 'bhjuv' %in% jlist)bh[,k] <- bhjuv[,k] <- bk[,'hinc']
    if('bf' %in% jlist)bf[,k] <- bk[,'fecn']
  }
  
  for(k in 1:length(jlist))outlist[[k]] <- get(jlist[k]) 
  
  outlist
}



updateBeta <- function(){

  bg <- bd
  lo <- loBd
  hi <- hiBd

  if(FONLY){
    bg <- bf
    lo <- loBf
    hi <- hiBf
  }
  if(DANDF)bg <- bjuv
   
  for(k in 1:nttype){

    if(treeBySpec[k] < minTrees)next

    wa <- which( treeRows(spec=k) )          #all
    wi <- which( treeRows(spec=k,matr=0) )          #immature or male
    wm <- which( treeRows(spec=k,matr=1,sex=1) )    #mature and female

    if(DONLY)ww <- wa                        #univariate dinc is ww
    if(FONLY)ww <- wm
    if(DANDF)ww <- wi                        #multivariate is ww (immature) and wm
      
    if(length(ww) <= length(bg))next

    b1  <- bg[pCol,k]                        #univariate and immature for multivariatebd
    y1  <- y[ww,1] - ag[tin[ww,1],1]

    if(DANDF){                               #mature only
        y1m <- y[wm,] - ag[tin[wm,1],]
        b1m <- cbind(bd[pCol,k],bf[pCol,k])
    }

    if(length(yrCol) > 0){
      wc <- which(colSums(X[ww,yrCol]) > 0)
      yc <- yrCol[wc]
      b2 <- bg[yc,k]
      y1 <- y[ww,1] - ag[tin[ww,1],1] - X[ww,yc]%*%b2
      if(DANDF){
         wcm <- which(colSums(X[wm,yrCol]) > 0)
         ycm <- yrCol[wcm]
         b2m <- cbind(bd[ycm,k],bf[ycm,k])
         y1m <- y[wm,] - ag[tin[wm,1],] - X[wm,ycm]%*%b2m
      }
    }

    b1     <- bUpdateNorm(X[ww,pCol],y1,b1,
                      priorB[pCol,k],priorIVB[pCol,pCol],
                      lo=lo[pCol,k],hi=hi[pCol,k],sg)
    bg[pCol,k] <- b1

    if(DANDF){
       b1m <- bUpdateMVN_Rcpp(X[wm,pCol],y1m,b1m,cbind(loB[pCol,k],loF[pCol,k]),
                                            cbind(hiB[pCol,k],hiF[pCol,k]),svar+avar)
       bd[pCol,k] <- b1m[,1]
       bf[pCol,k] <- b1m[,2]
    }

    if(length(yrCol) > 0){
        y2 <- y[ww,1] - ag[tin[ww,1],1] - X[ww,pCol]%*%b1
        b2 <- bUpdateNorm(X[ww,yc],y2,b2,
                     priorB[yc,k],priorIVB[yc,yc],
                     lo=lo[yc,1],hi=hi[yc,1],sg)
        bg[yc,k] <- b2

        if(DANDF){
            y2m <- y[wm,] - ag[tin[wm,1],] - X[wm,pCol]%*%b1m
            b2m <- bUpdateMVN_Rcpp(X[wm,ycm],y2m,b2m,cbind(loB[ycm,k],loF[ycm,k]),
                                               cbind(hiB[ycm,k],hiF[ycm,k]),svar+avar)
           bd[ycm,k] <- b2m[,1]
           bf[ycm,k] <- b2m[,2]
        }
     }         
  }

 if(DONLY | FONLY)return(bg)
 
 list(bd = bd, bf = bf, bjuv = bg)
}

nobsByIndiv <- function(wm=1:ntree){
  
  aa <- diamMat*0
  aa[tin[wm,]] <- 1
  list(ni = rowSums(aa), imat = aa)
}

updateA <- function(){   #random effect and variance, univariate

  yy  <- y2z(y,Kmat)
  res <- yy - xbdfMat  #does not have random effects

  mtree <- which( treeRows(matr=1,sex=1,years=c(1:(nyr-1))) )    #mature and female
  itree <- which( treeRows(matr=0,years=c(1:(nyr-1))) |  
                  treeRows(sex=0,years=c(1:(nyr-1))) )
  
  if(nyi == 1){
    siI <- matrix(1/svar[wyi,wyi],nyi,nyi)
    aiI <- matrix(1/avar[wyi,wyi],nyi,nyi)
  }
  if(nyi > 1){
    siI <- invMat(svar[wyi,wyi])
    aiI <- invMat(avar[wyi,wyi])
  }

  siM <- invMat(svar)
  aiM <- invMat(avar)
  
  ni <- nobsByIndiv(itree)$ni
  nm <- nobsByIndiv(mtree)$ni

  ni[is.na(ni)] <- 0
  nm[is.na(nm)] <- 0

  simatM <- matrix(siM,ntree,ny^2,byrow=T)
  simatI <- matrix(siI,ntree,nyi^2,byrow=T)

  # bigV
  if(nyi == 1)V1 <- 1/(simatI*matrix(ni,ntree,1) + matrix(aiI[wyi,wyi],ntree,1) )
  if(nyi == 2)V1 <- invertcol2(simatI*matrix(ni,ntree,4) + matrix(aiI[wyi,wyi],ntree,4,byrow=T) )
  
  if(ny == 2)V2  <- invertcol2(simatM*matrix(nm,ntree,ny^2) + matrix(aiM,ntree,ny^2,byrow=T))
  if(ny == 3)V2  <- invertcol3(simatM*matrix(nm,ntree,ny^2) + matrix(aiM,ntree,ny^2,byrow=T))$I
    
  # small v
  jjI  <- c(1:nyi)
  jjM  <- c(1:ny)
  vI  <- VvI <- matrix(0,ntree,nyi)
  vM  <- VvM <- matrix(0,ntree,ny)
  
  for(k in 1:ny){      
    rr      <- diamMat*0
    rr[tin] <- res[,k]
    if(k %in% wyi)vI  <- vI + simatI[,jjI]*rowSums(rr*(1 - qmat),na.rm=T)    #immature
    vM  <- vM + simatM[,jjM]*rowSums(rr*qmat*femMat,na.rm=T)   #mature & female
    jjI <- jjI + nyi
    jjM <- jjM + ny
  }

  # Vv  (mean)
  jjI  <- c(1:nyi)
  jjM  <- c(1:ny)
  for(k in 1:ny){
    if(k %in% wyi)VvI  <- VvI + V1[,jjI]*vI[,k]
    VvM  <- VvM + V2[,jjM]*vM[,k]
    jjI <- jjI + nyi
    jjM <- jjM + ny
  } 
  
  wm <- which(rowSums( qmat*femMat,na.rm=T) > 1)
  wi <- c(1:ntree)[-wm]
  
 # mu <- rep
    
  for(k in 1:ny){

     if(k %in% wyi & nyi == 1)tmp1 <- list( mu = VvI, vr = V1 )
     if(k %in% wyi & nyi == 2)tmp1 <- conditionalBiVarNorm(x=ag,mu=VvI,sigMat=V1,k)
     
     if(ny == 2)tmp2 <- conditionalBiVarNorm(x=ag,mu=VvM,sigMat=V2,k)
     if(ny == 3)tmp2 <- conditionalTriVarNorm(x=ag,mu=VvM,sigMat=V2,k)
          
     if(k %in% wyi)ag[wi,k] <- tnorm(length(wi),-tVals[k],tVals[k],tmp1$mu[wi],sqrt(tmp1$vr[wi]))
     ag[wm,k] <- tnorm(length(wm),-tVals[k],tVals[k],tmp2$mu[wm],sqrt(tmp2$vr[wm]))
  }
       
  if(ny == 1){
    u1 <- ai1 + length(wtree)/2
    u2 <- ai2 + .5*sum(ag[wtree]^2)
    avar <- 1/rgamma(1,u1,u2)
  }
  if(ny > 1){

   py <- xbdfMat[mtree,] + ag[tin[mtree,1],]
   ns <- length(mtree)

   sdelta <- adelta <- diag(1,ny)
   
   tmp <- updateScInvWish(yz=yy[mtree,],predy=py,delta=sdelta,priorDmu=rep(1,ny),
                          priorDvar=rep(.001,ny),varLo=loSigma,varHi=hiSigma)
   sdelta <- tmp$delta
   svar   <- tmp$sigma
 
  # svar <- updateWishart(yy[mtree,],py,priorS,ny+1)
   
   tmp <- updateScInvWish(ag[wm,],ag[wm,]*0,delta=adelta,priorDmu=rep(1,ny),
                          priorDvar=rep(.1,ny),varLo=loSigma,varHi=hiSigma)

 # avar <- updateWishart(ag[wm,],ag[wm,]*0,priorS,ny+1)

   adelta <- tmp$delta
   avar   <- tmp$sigma
  }

  colnames(avar) <- rownames(avar) <- colnames(svar) <- rownames(svar) <- ynames
       
  list(ag = ag, avar = avar, svar = svar, adelta = adelta, sdelta = sdelta)
}




demDeviance <- function(wi,wm,predy,sigmaCols){
  
  if(length(imCols) == 1)m1 <- dnorm(y[wi,imCols],py[wi,imCols],scols[wi,isCols],log=T)
  if(length(imCols) == 2)m1 <- dbivarNormFromCols(y[wi,imCols],py[wi,imCols],scols[wi,isCols]) 
  if(ny == 2)m2 <- dbivarNormFromCols(y[wm,],py[wm,],scols[wm,])
  if(ny == 3)m2 <- dtrivarNormFromCols(y[wm,],py[wm,],scols[wm,]) 
  regLik <- sum(m1) + sum(m2)

  -2*regLik
}



updatePhi <- function(){   # maturation parameters

  tiny <- 1e-15

  prop <- matrix(tnorm.mvtRcpp(pg,pg,proppg,lo=loP,hi=hiP,times=1),ncol=1)

  wm  <- which( treeRows() )
  qm  <- Q[wm]

  tnow <- getRho(M[wm,],pg)
  tnew <- getRho(M[wm,],prop)

  tnow[tnow < tiny] <- tiny
  tnow[tnow > (1 - tiny)] <- 1 - tiny
  tnew[tnew < tiny] <- tiny
  tnew[tnew > (1 - tiny)] <- 1 - tiny  

  pnow <- qm*log(tnow) + (1 - qm)*log(1 - tnow) + dmvnorm(t(pg),priorPhi,priorVPhi,log=T)
  pnew <- qm*log(tnew) + (1 - qm)*log(1 - tnew) + dmvnorm(t(prop),priorPhi,priorVPhi,log=T)

  a <- exp(sum(pnew,na.rm=T) - sum(pnow,na.rm=T))
  z <- runif(1,0,1)
  if(z < a)pg <- prop
  pg
}

updateVerror <- function(){  # observation error on status

  svec <- sexObs[tin]
  oimm <- which(svec == 0 & Q == 1)  # not recognized as repro
  omat <- which(svec == 2 & Q == 1)  # recognized as reproductive  CHECK, ZEROS HERE

  p1 <- v1 + length(omat)
  p2 <- v2 + length(oimm)

  rbeta(1,p1,p2)
}


predData <- function(bjuv,bd,bf,svar,avar,ag,sg,kg=NULL,SAMPLE=F){

    wi <- which(treeRows(matr=0) | treeRows(matr=1,sex=0))
    wm <- which(treeRows(matr=1,sex=1))

    py <- matrix(NA,n,2)
    samp <- NULL
    if(SAMPLE)samp <- py

      py[wi,1] <- (xbgMat[wi] + ag[tin[wi,1],1])*Kmat[wi,1]
      tmp  <- getMonod(light,kg,bd,bf,svar,wm)
      py[wm,] <- tmp$mu
      pvar    <- tmp$sigma

      if(SAMPLE){
       samp[wi,1] <- tnorm(length(wi),minDinc,maxDinc,py[wi,1],sqrt(sg)*Kmat[,1])
       samp[wm,]  <- rbvnormFromCols(py[wm,],pvar,lo=0,hi=100)
      }

    list(mu = py, z = samp)
}

predModel <- function(ename,pname,quant=NULL,SPEC=1,add=F,col=1,FUN,ylim=NULL,
                      xlab=pname,ylab=ename,lite=mean(light,na.rm=T)){

  FUN <- match.fun(FUN)

  nsim <- 2000
  ww <- which(chains == ename)
  vv <- chainList[[ww]]

     kgvals <- chainList[[which(chains == 'kg')]]
     kgcol  <- 1
     if(ename == 'bf')kgcol=2

  xx <- X
  if(ename == 'pg')xx <- M

  ecols <- colnames(xx)

  mmeans <- colMeans(xx)

  ns <- 100
  xpred <- matrix(mmeans,ns,ncol(xx),byrow=T)
  colnames(xpred) <- ecols
  xpred[,pname] <- seq(min(xx[,pname]),max(xx[,pname]),length=ns)

  if(!is.null(quant)){

    for(k in 1:length(quant)){
      wk <- which(ecols == names(quant)[k])
      qk <- quant[[k]]
      xpred[,wk] <- quantile(xx[,wk],qk)
    }
  }

  gi <- sample(burnin:ng,nsim,replace=T)

  mp <- mrange <- numeric(0)
  kk <- 1

  for(k in SPEC){
    mpred <- matrix(NA,nsim,ns)
    kcol  <- c(1:ncol(xx))
    if(ename == 'bf')kcol  <- grep(treeNames[k],colnames(vv))
    
    for(i in 1:nsim){
      mu <- xpred%*%vv[gi[i],kcol]
      if(ename == 'bf')mu <- mu*getKmat(lite,kgvals[gi[i],])[kgcol]

      mpred[i,] <- FUN(xpred%*%vv[gi[i],kcol])
    }
    mp <- append(mp,list(apply(mpred,2,quantile,c(.5,.025,.975))) )
    if(is.null(ylim))ylim <- range(c(mrange,mp[[kk]]))

    if(!add & kk == 1)plot(xpred[,pname],mp[[kk]][1,],type='l',lwd=3,ylim=ylim,xlab=xlab,ylab=ylab)
    lines(xpred[,pname],mp[[kk]][1,],lwd=3,col=col)
    lines(xpred[,pname],mp[[kk]][2,],lwd=2,lty=2,col=col)
    lines(xpred[,pname],mp[[kk]][3,],lwd=2,lty=2,col=col)
    title(treeNames[kk])
    kk <- kk + 1
  }

  invisible(list(x = xpred, y = mp) )
}

mapPlot <- function(jplot='CW_118',tplot=NULL,zname='diam',add=F,col=1,zscale=1,
                    fill=F,RELATIVE=F,specLegend=F){

  wp <- which(treeData[,'j'] %in% jplot)
  if(length(wp) < 2)return(numeric(0))


  specColors <- mapColors(nttype)

  seedRange <- specCodes <- 0

  utj <- numeric(0)

  for(j in jplot)utj <- rbind(utj,UTMhull[[j]])
  
  ran <- apply(utj,2,range)
    
  x <- UTMtree[wp,'x']
  y <- UTMtree[wp,'y']
  colVec <- specColors[ treeData[wp,'spec'] ]
  sym <- 'circles'

  if(is.null(tplot))tplot <- c(1:nyr)

  if(zname == 'diam'){
     z <- diamMat[wp,tplot]
     if(length(tplot) > 1)z <- rowMeans(z,na.rm=T)
  }
  if(zname == 'dinc'){
     z <- dincMat[wp,tplot]
     if(length(tplot) > 1)z <- rowMeans(z,na.rm=T)
  }
  if(zname == 'fecn'){
     z <- sqrt(fmat[wp,tplot])
    if(length(tplot) > 1)z <- rowMeans(z,na.rm=T)
  }
  if(zname == 'seed'){
     ww <- which(seedmat[,'j'] %in% jplot & seedmat[,'t'] %in% tplot)
     ss <- seedmat[ww,]
     x <- UTMseed[ss[,'sindex'],'x']
     y <- UTMseed[ss[,'sindex'],'y']
     z <- ss[,'count']

     seedRange <- range(z,na.rm=T)

     if(length(tplot) > 1){

         jkmat <- xmat <- ymat <- matrix(0,max(seedmat[,'k']),nplot)

         for(t in tplot){
           wt <- which(ss[,'t'] == t)
           jk <- cbind(ss[wt,'k'],ss[wt,'j'])
           jkmat[jk] <-  jkmat[jk] + ss[wt,'count']
           xmat[jk]  <- UTMseed[ss[wt,'sindex'],'x']
           ymat[jk]  <- UTMseed[ss[wt,'sindex'],'y']
         }

         z <- (jkmat/length(tplot))[jk]
         x <- xmat[jk]
         y <- ymat[jk]
    }

     colVec <- rep(col,length(z))
     sym <- 'squares'
  }

  if(max(z,na.rm=T) == 0)return()

  if(RELATIVE)z <- z/max(z,na.rm=T)

  z <- z*zscale*max(apply(ran,2,diff))/100

  mapSpecies(x, y, z,
           mapx=ran[,1], mapy=ran[,2], add=add, sym=sym, colVec = colVec,fill=fill)

  if(specLegend){
     wj <- sort( unique(treeData[treeData[,'j'] %in% jplots,'spec']) )
     un <- grep('UNKN',treeNames)
     if(length(un) > 0)wj <- wj[wj != un]

     legend('topleft',treeNames[wj],text.col=specColors[wj],cex=.6,bg='white',bty='white')
  }


 # if(!add)title(paste(zname,plotnames[jplot]))

  invisible( list(xyz = cbind(x,y,z), seedRange = seedRange, specCodes = colVec) )
}


mapSeedPred <- function(jplot,tplot, u, fecn, wide = 10, add=F,colVec=1, spec=NULL){  #  wide - grid in meters

  wp <- which(treeData[,'j'] %in% jplot)
  if(!is.null(spec))wp <- which(treeData[,'j'] %in% jplot & treeData[,'spec'] == spec)
  ntt <- length(wp)
  if(ntt == 0)return()

  tmp <- plotUTMrange(jplot)
  xr  <- x2 <- tmp$xrange
  yr  <- yr <- tmp$yrange
  ran <- tmp$urange
  xyz <- tmp$utmxyz

  xseq <- seq(ran[1,1],ran[2,1],by=wide)
  yseq <- seq(ran[1,2],ran[2,2],by=wide)
  grid <- as.matrix( expand.grid(xseq,yseq) )

  z <- rep(0,nrow(grid))

  si <- c(1:n)
  if(!is.null(spec))si <- which(treemat[,'spec'] == spec)

   for(i in 1:ntt){

     wi <- which(treemat[,'tindex'] == wp[i] & treemat[,'t'] %in% tplot)
     if(length(wi) == 0)next

     fi <- fecn[wi]
     ni <- length(fi)
     if(ni == 0)next

     di <- distmat(xyz[i,'x'],xyz[i,'y'],grid[,1],grid[,2])
     z <- z + rowMeans( getKern(u,di)%*%matrix(fi,1) )
   }

  zmat <- matrix(z,length(xseq),length(yseq))

  kk <- length(colVec)

  for(j in 1:length(colVec)){
    contour(xseq,yseq,zmat,add=add,col=colVec[j],lwd=2*(kk - j + 1),levels=c(1,10,100,1000))
  }

  invisible(list(x = xseq, y = yseq, z = zmat))
}




treePlot <- function(xmat,ymat,add=F,xlab=' ',ylab=' ',cols='black',
            xlim=range(xmat,na.rm=T),ylim=range(ymat,na.rm=T),nPerBin=NULL,logy=F){  #i by t matrices of x and y

  if(length(cols) == 1)cols <- rep(cols,nrow(xmat))

  if(!add){
    if(!logy)plot(xmat[1,],ymat[1,],type='l',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,col=cols)
    if(logy) plot(xmat[1,],ymat[1,],type='l',xlim=xlim,ylim=ylim,log='y',xlab=xlab,ylab=ylab,col=cols)
  }
  for(i in 1:nrow(xmat)){
    lines(xmat[i,],ymat[i,],col=cols[i])
  }

  if(!is.null(nPerBin)){

      nbb <- nPerBin/nrow(xmat)
      nbb <- seq(0,1,by=nbb)
      if(max(nbb) < 1)nbb <- c(nbb,1)
      bins <- quantile(xmat,nbb,na.rm=T)
      nbin <- length(bins)
      mids <- (bins[-nbin] + bins[-1])/2
      mids <- mids[mids > 0]
      nbin <- length(mids)

      qplot <- c(.5,.025,.1,.9,.975)
      ci <- matrix(NA,length(qplot),nbin)

     for(k in 1:(nbin-1)){
        ci[,k]  <- quantile(ymat[xmat >= mids[k] & xmat <= mids[k+1]],qplot,na.rm=T)
     }
     lines(mids,ci[1,],col='white',lwd=6)
     lines(mids,ci[1,],lwd=2)
     for(k in 2:length(qplot)){
        lines(mids,ci[k,],col='white',lwd=4,lty=2)
        lines(mids,ci[k,],lwd=1,lty=2)
     }
  }
}

u2d <- function(u){  #dispersal parameter to mean distance

  sqrt(u)*pi/2
}
d2u <- function(d){  #mean distance to dispersal parameter

 (2*d/pi)^2
}



lagGrowth <- function(dimat=dincMat,timelag,discount){   

  #weighted ave grow of previous timelag years, discount < 1

  dd   <- matrix(0,ntree,nyr)
  dtm  <- dt0 <- dimat
  cc   <- dincMat*0
  sw    <- 0

  fMat <- matrix(firstTime,ntree,nyr)
  timeMat  <- matrix(1:nyr,ntree,nyr,byrow=T)

  for(t in 1:timelag){

    wi <- which(timeMat - fMat >= (t-1),arr.ind=T)
    wt  <- t^(-discount)
    dd[wi]  <- dd[wi] + wt*dtm[wi]
    cc[wi]  <- cc[wi] + wt

    if(t < timelag)dtm <- lastTimeMat(dtm,firstTime,INC=F,minDinc)

  }
  dd <- dd/cc
  dd[dd < minDinc] <- minDinc
  dd[dd > maxDinc] <- maxDinc

  dd
}

survDincDiam <- function(dd=dlagMat,dmat=diamMat){

  wd <- which(survMat*censorMat == 0,arr.ind=T)
  wl <- which(survMat*censorMat == 1,arr.ind=T)

  liveGrow <- findInterval(dd[wl],dincSeq,all.inside=T)
  deadGrow <- findInterval(dd[wd],dincSeq,all.inside=T)

  liveDiam <- findInterval(dmat[wl],diamSeq,all.inside=T)
  deadDiam <- findInterval(dmat[wd],diamSeq,all.inside=T)

  list(lg = liveGrow, dg = deadGrow, ld = liveDiam, dd = deadDiam, liveIndex = wl, deadIndex = wd)
}

survLikelihoodOld <- function(parG,parD,liveG,deadG,liveD,deadD){  # log likelihood, inputs are binned

  tiny   <- 1e-10

  tlive <- parG[liveG] + parD[liveD] - 2
  tdead <- 1 - exp( parG[deadG] + parD[deadD] - 2)

  tdead[tdead < tiny] <- tiny
  
  list(live = tlive, dead = log(tdead) )

}

survLikelihoodOld2 <- function(dincSurv,diamSurv,liveG,deadG,liveD,deadD){  # log likelihood, inputs binned

  tiny   <- 1e-4

  tlive <- 1 - exp(-dincSurv[liveG] - diamSurv[liveD])
  tdead <- exp(-dincSurv[deadG] - diamSurv[deadD])

  tdead[tdead < tiny] <- tiny
  tlive[tlive < tiny] <- tiny
  
  list(live = log(tlive), dead = log(tdead) )

}

survLikelihood <- function(dincSurv,diamSurv,liveG,deadG,liveD,deadD){  # log likelihood

  tiny   <- -10

  tlive <- diamSurv[liveD] * log(dincSurv[liveG])
  tdead <- log(1 - dincSurv[deadG]^diamSurv[deadD])

  tdead[tdead < tiny] <- tiny
  tlive[tlive < tiny] <- tiny
  
  list(live = tlive, dead = tdead )

}


updateDeathPars <- function(){  #growth and diameter as competing risks

 # maxdvar <- 6
 # maxdiam <- .1

  maxdvar <- 1
  maxdiam <- 3

  accept <- 0

  d1   <- dincSurv[-1]
  d2   <- dincSurv[-ndinc]
  mids <- (d1 + d2)/2
  lo   <- c(0,mids)
  hi   <- c(mids,maxdvar)
  hi[ndinc-1] <- lo[ndinc] <- maxdvar
  sp   <- (hi - lo)^runif(ndinc,.5,3)
  dincNew <- tnorm(ndinc,lo,hi,dincSurv,sp)
  dincNew[1] <- 0
  dincSurv[ndinc] <- dincNew[ndinc] <- maxdvar

  d1   <- diamSurv[-1]
  d2   <- diamSurv[-ndiam]
  mids <- (d1 + d2)/2

  lo   <- c(1,mids)
  hi   <- c(mids,maxdiam)

  sp   <- hi - lo
  diamNew <- tnorm(ndiam,lo,hi,diamSurv,sp)
  diamNew[1] <- diamSurv[1] <- 1

  tmp   <- survDincDiam(dlagMat,diamMat) #bin growth and diameter
  liveG <- tmp$lg
  deadG <- tmp$dg
  liveD <- tmp$ld
  deadD <- tmp$dd

  tmp1  <- survLikelihood(dincSurv,diamSurv,liveG,deadG,liveD,deadD) 
  pnow  <- sum(tmp1$live) + sum(tmp1$dead)

  tmp2  <- survLikelihood(dincNew,diamNew,liveG,deadG,liveD,deadD) 
  pnew  <- sum(tmp2$live) + sum(tmp2$dead)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a){
    dincSurv <- dincNew
    diamSurv <- diamNew
    accept   <- 1
  }

  list(dincSurv = dincSurv, diamSurv = diamSurv, accept = accept)
}

  

sampleDeathYr <- function(){

  dead <- rep(NA,length(deadTrees))
  k <- 0

  dinc <- dincMat

  for(i in deadTrees){

    k <- k + 1

    if(is.finite(treeDeath[i,'censoryr']))next
    ti <- treeDeath[i,'t1']:treeDeath[i,'t2']
    ni <- length(ti)

    di <- dincSurv[findInterval(dlagMat[i,ti],dincSeq,all.inside=T)]  #time-specific death rate
    dd <- diamSurv[findInterval(diamMat[i,ti],diamSeq,all.inside=T)]
    st <- di^dd                                                   #pr survive interval
    st[is.na(st)] <- 1

    f1 <- matrix(st,ni,ni,byrow=T)
    diag(f1) <- 1 - st
    f1[lower.tri(f1,diag=F)] <- 1
    ff <- apply(f1,2,cumprod)[ni,]   
    deadyr <- ti[ rmultinom(1,1,ff) == 1]

    survMat[i,deadyr] <- 0
    survMat[i,years > deadyr] <- NA
    survMat[i,firstTime[i]:(deadyr-1)] <- 1

    dead[k] <- deadyr

    dinc[i,deadyr:nyr] <- minDinc/2
  }

  deadIndex = cbind(deadTrees,dead)

  ww <- which(is.na(dead))                       #exclude censored individuals
  if(length(ww) > 0)deadIndex <- deadIndex[-ww,]

  survMat <- survMat

  list(survMat = survMat, survVec = survMat[tin], deadIndex = deadIndex, dincMat = dinc)
}


getShade <- function(bspec,ball){

  tmp <- dateColumn(bspec,'diam',yrvec)
  bdiam <- tmp$x
  byr   <- tmp$yr

  tmp <- dateColumn(ball,'diam',yrvec)
  adiam <- tmp$x
  ayr   <- tmp$yr

  hall <- adiam
  hb   <- bdiam

  if(length(htAllomFile) > 0){                 #ht differences

    hall <- allomConvert(adiam,ball[,'species'],allomFile=paste(dataPath,htAllomFile,sep=''),
                         codeColumn='species',
                        'htInt','htSlope',defaultSpec='other')
    hb   <- allomConvert(bdiam,bspec[,'species'],allomFile=paste(dataPath,htAllomFile,sep=''),
                         codeColumn='species',
                        'htInt','htSlope',defaultSpec='other')
  }
  
  rangeh <- t( apply(hall,1,range,na.rm=T) )
    wt  <- which(is.finite(rangeh[,1]))
    tmp <- interpRows(hall[wt,],INCREASING=T,minVal=rangeh[,1],maxVal=rangeh[,2])
    hall[wt,] <- tmp

  rangeh <- t( apply(hb,1,range,na.rm=T) )
    wt  <- which(is.finite(rangeh[,1]))
    tmp <- interpRows( row2Mat(hb[wt,]),INCREASING=T,minVal=rangeh[,1],maxVal=rangeh[,2],tinySlope=.1)
    hb[wt,] <- tmp
  

  first <- matrix(ball[,'censinyr'],nrow(hall),ncol(hall))
  last  <- matrix(apply(ball[,c('censoryr','deathyr')],1,min,na.rm=T),nrow(hall),ncol(hall))
  ymat  <- matrix(yrvec,nrow(hall),ncol(hall),byrow=T)
  hall[first > ymat | last < ymat] <- 0

  shade <- matrix(NA,nrow(bspec),nyr)

  maxh <- apply(hall,1,max,na.rm=T)
  
  for(i in 1:nrow(bspec)){

     hi <- hb[i,]
     ci <- which(maxh > min(hi,na.rm=T))
     
     dij  <- distmat(bspec[i,'UTMx'],bspec[i,'UTMy'],ball[ci,'UTMx'],ball[ci,'UTMy'])
     wij  <- which(dij < 15 & dij > 0)

     if(length(wij) == 0){  #no taller neighbors
       shade[i,] <- 0
       next
     }
     hij  <- hall[ci[wij],]                                    #neigbors

     hmat <- matrix(hi,length(wij),nyr,byrow=T)  #focal individual
     hij[hij < hmat] <- NA

     wt <- matrix(dij[wij]^2,nrow(hmat),nyr) + 1
     cc <- colSums( (hij - hmat)/wt,na.rm=T)
     shade[i,] <- 1 - exp(-cc)
     shade[i,is.na(hb[i,])] <- NA
  }

  shade <- round(shade,5)

  list(shadeIndex = shade, htmat = round(hb,2))
}



 ###########################################
inDEM <- function(pname=NULL,rname=NULL,path){   #old           

#pname has REGION_PLOT, e.g. DF_BW, path = datafiles/DEM or datafiles/TCI
#rname is REGION, e.g. 'DF', 'CW'
 
 if(!is.null(pname))region <- strsplit(pname,'_')[[1]][1]
 if(!is.null(rname))region <- rname
 fname  <- paste(path,region,'.txt',sep='')
 header <- read.table(fname, nrow = 6)

 for(i in 1:6)assign(as.character(header$V1[i]), header$V2[i])

   dat <- read.table(fname, skip = 6, na.strings = '-9999')
   colnames(dat) <- seq(from = xllcorner, by = cellsize, length = ncols)
   rownames(dat) <- seq(to = yllcorner, by = - cellsize, length = nrows)

 dat
}

inTCI <- function(pname,path='datafiles/TCI/'){

 region <- strsplit(pname,'_')[[1]][1]
 fname  <- paste(path,region,'_TCI.txt',sep='')
 header <- read.table(fname, nrow = 6)

 tmp
}


getMinMaxX <- function(xnames){

  minMaxDiam  <- c(0,100)     # cm
  minMaxLight <- c(0,1)       # fraction full exposure
  minMaxTemp  <- c(0,16)    # degrees C
  minMaxP     <- c(500,3500)  # mm
  minMaxPDSI  <- c(-5,5)
  if( length(grep('wetness',hydroPath)) > 0 ) minMaxHydro <- c(-5,5)
  if( length(grep('flow',hydroPath)) > 0 )    minMaxHydro <- c(1,4)
  if( length(grep('gradient',hydroPath)) > 0 )minMaxHydro <- c(1,4)
  if( length(grep('tci',hydroPath)) > 0 )     minMaxHydro <- c(0,1e+7)

  minMaxX <- matrix(0,length(xnames),2)
  rownames(minMaxX) <- xnames
  colnames(minMaxX) <- c('min','max')

  for(k in 1:length(xnames)){

    if(xnames[k] == 'diam')               minMaxX[k,]  <- minMaxDiam
    if(xnames[k] == 'light')              minMaxX[k,] <- minMaxLight
    if(length(grep('prec',xnames[k])) > 0)minMaxX[k,] <- minMaxP
    if(length(grep('temp',xnames[k])) > 0)minMaxX[k,] <- minMaxTemp
    if(length(grep('pdsi',xnames[k])) > 0)minMaxX[k,] <- minMaxPDSI
    if(xnames[k] == 'hydro')              minMaxX[k,] <- minMaxHydro
  }
  minMaxX
}

nearestHydro <- function(utmx,utmy,hydroFile){

   nn <- length(utmx)

   if( !file.exists(hydroFile) ){
     warning( paste(hydroFile,'missing') )
     return(rep(NA,nn))
   }

   tmp <- read.asciigrid.r(hydroFile,as.image=T, plot.image=F)
   hydrox <- tmp$x
   hydroy <- tmp$y
   hydroz <- tmp$z

   TCI <- rep(NA,nn)

   for(i in 1:nn){

     minX    <- which.min( (utmx[i] - hydrox)^2 )
     minY    <- which.min( (utmy[i] - hydroy)^2 )
     TCI[i]  <- hydroz[minX,minY]

      if(is.na(TCI[i]) | TCI[i] == 0){
        searchX <- c((minX - 10):(minX + 10))
        searchY <- c((minY - 10):(minY + 10))
        wx <- which(searchX < 0 | searchX > nrow(hydroz))
        wy <- which(searchY < 0 | searchY > ncol(hydroz))
        if(length(wx) > 0)searchX <- searchX[-wx]
        if(length(wy) > 0)searchY <- searchY[-wy]
        
        search  <- expand.grid(searchX,searchY)
        wf <- which(!is.na(hydroz[cbind(search[,1],search[,2])]) & 
                           hydroz[cbind(search[,1],search[,2])] > 0)
        search <- search[wf,]

        dista <- t(distmat(utmx[i],utmy[i],hydrox[search[,1]],hydroy[search[,2]]))
        closest <- apply(dista,1,which.min)
        TCI[i] <- hydroz[search[closest,1],search[closest,2]]
      }
    }

    wf <- which(!is.na(TCI) & TCI > 0)
    wn <- which(is.na(TCI) | TCI == 0)

    if(length(wf) > 0){
      dista   <- t(distmat(utmx[wn],utmy[wn],utmx[wf],utmy[wf]))
      closest <- apply(dista,1,which.min)
      TCI[wn] <- TCI[wf[closest]]
    }

  TCI
}
     

hydroByTree <- function(treeData,UTMdata){

  nn    <- nrow(treeData)
  hydro <- rep(NA,nn)

  tpath <- paste(dataPath,hydroPath,sep='')

  for(j in 1:nplot){

    wj      <- which(treeData[,'plot'] == plotnames[j])
    if(length(wj) == 0)next

    region <- strsplit(plotnames[j],'_')[[1]][1]
    fname  <- paste(tpath,region,'.asc',sep='')

    UTM    <- row2Mat(UTMdata[wj,])

    hydro[wj] <- nearestHydro(UTMdata[wj,1],UTMdata[wj,2],fname)
  }
  hydro
}


plotPlot <- function(x,jindex,xname=' '){

  plot(jitter(jindex),x,cex=.2)
  
  text(c(1:nplot),5,plotnames,srt=90)
  title(xname)

}


#diam2mass <- function(d,alloB,alloL){  #allo has c(int,slope) on log10 scale

#  b <- 10^alloB[,1] *d^alloB[,2]
#  l <- 10^alloL[,1] *d^alloL[,2]
#  l[is.finite(l) & l > b] <- b[is.finite(l) & l > b]*.9
  
#  list(biom = b, leaf = l)
#}

stemMass2diam <- function(stem,leaf,alloB){

  10^( (log10(stem + leaf) - alloB[,1])/alloB[,2] )
} 

mass2diam <- function(mass,allo){

  10^((log10(mass) - allo[,1])/allo[,2])
}




getNames <- function(genus=NULL,omit=NULL){

  seeds <- trees <- mono <- character(0)

  traits <- vector('list',2)
  names(traits) <- c('gm1000Seed','kgM3Wood')

  treeNames <- read.table(paste(dataPath,treeCodeFile,sep=''),header=T)
  trees     <- as.character( treeNames[treeNames[,'genus'] %in% genus,'code'] )

  if(!is.null(omit)){
      wo <- grep(omit,trees)
      if(length(wo) > 0)trees <- trees[-wo]
  }


  if('gm1000Seed' %in% colnames(treeNames))traits$gm1000Seed <- treeNames[,'gm1000Seed']
  if('kgM3Wood' %in% colnames(treeNames))traits$kgM3Wood <- treeNames[,'kgM3Wood']

  if('fecn' %in% ynames){
    seedNames <- read.table(paste(dataPath,seedCodeFile,sep=''),header=T)
    seeds     <- seedNames[seedNames[,'genus'] %in% genus,'code']
    mono      <- treeNames[treeNames[,'genus'] %in% genus,'monoecious']

    for(k in 1:length(excludeSeedCode)){
       wk <- grep(excludeSeedCode[k],seeds)
       if(length(wk) > 0)seeds <- seeds[-wk]
    }
  }

  seeds <- as.character(seeds)
  trees <- as.character(trees)

  if(length(trees) == 0)stop('genus not in treeCodeFile')
  if(length(seeds) == 0)stop('genus not in seedCodeFile')

  print(as.character(trees))
  print(as.character(seeds))
  
  #default for allometric equations
  df <- 'other'
  if(genus == 'Acer')    df <- 'acerRubr'
  if(genus == 'Betula')  df <- 'betuAlle'
  if(genus == 'Carya')   df <- 'caryGlab'
  if(genus == 'Fraxinus')df <- 'fraxAmer'
  if(genus == 'Magnolia')df <- 'magnFras'
  if(genus == 'Nyssa')   df <- 'nyssSylv'
  if(genus == 'Pinus')   df <- 'pinuEchi'
  if(genus == 'Quercus') df <- 'querAlba'
  if(genus == 'Ulmus')   df <- 'ulmuAmer'

  list(treeNames = trees, seedNames = seeds, mono = mono, traits = traits, defaultAllomSpec = df)
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


inCovariates <- function(vnames,minMax=NULL){

  qnames <- vnames[vnames != 'years']

  cyears <- c( (yrvec[1] - climLag):max(yrvec) )

  n <- nrow(treemat)
  p <- length(vnames)
  x <- matrix(1,n,length(qnames))

  yrmat <- yrIndex <- numeric(0)

  if(!is.null(minMax) & !is.matrix(minMax)){
    minMax <- row2Mat(minMax)
    rownames(minMax) <- vnames
  }

  climSummary <- numeric(0)
     qj <- numeric(0)
  
  graphics.off()

  par(mfrow=c(2,2),bty='n')

  for(k in 1:p){

     qj <- qc <- numeric(0)

     if(vnames[k] == 'intercept') next

     if(vnames[k] %in% censusVariables){
         wc    <- match(vnames[k],censusVariables)
         x[,k] <- get(censusVariables[wc])
     }

     CK <- numeric(0)
     for(kk in 1:length(climateNames)){
       if( length(grep(climateNames[kk],vnames[k])) > 0)CK <- kk
     }
     
     kp <- 0

     if( vnames[k] %in% climateNames | length(CK) > 0 ){

         for(j in 1:nplot){

           cj <- regVec[j]

           wj <- which(treemat[,'j'] == j)
           if(length(wj) == 0)next

           months <- NULL
           if(length(grep('-',vnames[k])) > 0)months   <- unlist(strsplit(vnames[k],'-'))[2]

           nss <- length(subplotNames[[j]])
           if(nss == 0)nss <- 1

           for(kk in 1:nss){
             
             slowdown <- apply(matrix(100,1000,1000),1,sum)

             tmp  <- getClimateData(x=climateNames[CK],months,
                                  whichplot=plotnames[j],
                                  whichsub=subplotNames[[j]][kk],years=cyears,
                                  climatePath=climatePath)

             if(is.null(tmp))next

             clim <- tmp$climvar
             mvec <- tmp$monthVec

             clim[clim == -99] <- NA

             if(nss == 1){             #no subplots
               wss <- which(treeData[,'j'] == j)
               if(length(wj) == 0)next
             }

             if(nss > 1){             #there are subplots
                wss <- which(treeData[,'j'] == j & 
                             treeData[,'subplot'] == as.character(subplotNames[[j]][kk]))
                wj <- which(treemat[,'tindex'] %in% wss)
                if(length(wj) == 0)next
             }

             tyr <- yrvec[treemat[wj,'t']]
             if(max(mvec) > 6)tyr <- tyr - climLag   #previous year climLag = 1

             ti   <- match(tyr,clim[,'yr'])
             ww   <- which(!is.finite(ti))
             if(length(ww) > 0)message( paste('missing ',vnames[k],' in ',k,sep='') )
             x[wj,k]   <- clim[ti,2]

             if(kp == 0){
               mm <- range(clim[,2],na.rm=T)
               if(!is.null(minMax))mm <- minMax[k,]
               plot(0,0,xlim=range(yrvec),ylim=mm,xlab='yr',ylab=vnames[k])
               title(vnames[k])
               abline(h=0,lty=2)
             }

             lines(yrvec,clim[clim[,'yr'] %in% yrvec,2],col=cj)
             kp <- kp + 1

           }

             qc <- signif(c(mean(clim[,2],na.rm=T),range(clim[,2],na.rm=T)),3)
             names(qc) <- paste(vnames[k],c('mean','min','max'))
             qj <- rowBind(qj,qc,plotnames[j])
        }

         climSummary <- cbind(climSummary,qj)


        legend('topleft',names(plotgroups),text.col=c(1:length(plotgroups)),bty='n')
    #    dev.print(device=postscript,file=paste(outMain,'climateVars.ps',sep=''),width=6,horizontal=F)
         dev.copy2pdf(file=paste(outMain,'climateVars.pdf',sep=''))
     }

  #   print(climSummary)
     
     if(!vnames[k] %in% censusVariables & length(CK) == 0 & vnames[k] != 'years'){
        x[,k] <- get(vnames[k])
     }

     if(vnames[k] == 'years'){
        
        rmat <- matrix(0,n,nreg)
        rmat[cbind(c(1:n),treemat[,'reg'])] <- 1
        tmat <- matrix(0,n,nyr) 
        tmat[cbind(c(1:n),treemat[,'t'])] <- 1

        yrmat <- numeric(0)
        for(r in 1:nreg)yrmat <- cbind(yrmat, matrix(rmat[,r],n,nyr)*tmat)
        ynames <- outer(yrvec,c(1:nreg),paste,sep='-')
        colnames(yrmat) <- ynames
        
        ysum  <- colSums(yrmat)
        yrmat <- yrmat[,ysum > 0]

      }
  }############end k loop

  colnames(x) <- qnames
  p     <- ncol(x)
  yrCol <- numeric(0)

  missing <- which(is.na(x),arr.ind=T)
  if(length(missing) > 0){
      xmean <- colMeans(x,na.rm=T)
      x[missing] <- xmean[missing[,2]]
  }

  if(length(yrmat) > 0){              #yr effects
     yrIndex <- which(t(yrmat) == 1,arr.ind=T)[,1] + ncol(x)
     x <- cbind(x,yrmat)
     yrCol <- c((p+1):ncol(x))
  }

  meanX <- colMeans(x)
  sdX   <- apply(x,2,sd)
 # xstand <- (x - matrix(meanX,nrow(x),ncol(x),byrow=T))/matrix(sdX,nrow(x),ncol(x),byrow=T)
 # if(STANDARDIZE)x[,2:p] <- xstand[,2:p]

  if(!is.null(minMax)){           #standardize to minMax scale
    for(k in 1:p){
      if(rownames(minMax)[k] %in% c('intercept','years'))next
      x[,k] <- (x[,k] - minMax[k,1])/(minMax[k,2] - minMax[k,1])
    }
  }

  if('hydro' %in% vnames & reverseHydro)x[,'hydro'] <- -x[,'hydro']   #wet sites have large hydro

  if(length(negPrior) > 0){
    x[,negPrior] <- 1 - x[,negPrior]
  }
      
   list(x = x, missing = missing, p = p, yrCol = yrCol, yrIndex = yrIndex, meanX = meanX, sdX = sdX)
}



inTreePlotData <- function(){
  
  lf <- list.files(dataPath)
  if(!areaFile %in% lf)stop( paste(areaFile,' not in ',path) )
  if(!datFile %in% lf) stop( paste(datFile,' not in ',path) )

  plotArea <- read.table(paste(dataPath,areaFile,sep=''),header=T,row.names=1)
  colnames(plotArea) <- matrix( unlist( strsplit( colnames(plotArea) ,'X') ),ncol=2,byrow=T)[,2]

  if( as.numeric(max(colnames(plotArea))) < max(yrvec) ){
    misscol <- nyr - ncol(plotArea)
    mmat <- matrix( plotArea[,ncol(plotArea)],nrow(plotArea),ncol=misscol)
    colnames(mmat) <- yrvec[(nyr - misscol): nyr][-1]
    plotArea <- cbind(plotArea,mmat)
  }

  plotArea <- plotArea[,colnames(plotArea) %in% yrvec]

  tmp      <- read.table(paste(dataPath,datFile,sep=''),header=T,row.names=1)
  plotData <- tmp[match(plotnames,rownames(tmp)),]

  if(!'seedArea' %in% colnames(plotData))
      warning('seedArea missing from plotData, contains seed trap area')

  list(plotData = plotData, plotArea = plotArea)
}


filename2plot <- function(files){   
  # plot names for data files: removes version no. and active.txt
  # names is REGION_PLOT_VERSION_active.txt

  f1 <- matrix(unlist(strsplit(files,'_active')),ncol=2,byrow=T)[,1]  #remove 'active.txt'
  f2 <- matrix(unlist(strsplit(f1,'_')),ncol=3,byrow=T)[,1:2]
  if(!is.matrix(f2))f2 <- matrix(f2,1)
  paste(f2[,1],f2[,2],sep='_')

}


inSeedTraps <- function(omitTrap='trapnum',firstMonth = 1){

  allplots <- unlist(plotgroups)
  nplot    <- length(allplots)

  path   <- paste(dataPath,seedPath,sep='') 
  afiles <- list.files(path,pattern='active.txt')

  fplots <- filename2plot(afiles)

  seedMat <- perM2 <- numeric(0)
  allCodes <- character(0)

  trapPlotYr <- matrix(NA,nplot,nyr)
  rownames(trapPlotYr) <- plotnames
  colnames(trapPlotYr) <- yrvec

  for(j in 1:nplot){

    wf     <- which(fplots == plotnames[j])
    if(length(wf) == 0)next

    ff     <- afiles[wf]
    backup <- grep("~",ff)
    if(length(backup) > 0)ff <- ff[-backup]

    fname  <- paste(path,ff,sep='')
    b      <- read.table(fname,header=T,sep="\t",fill=T)

    ww <- numeric(0)

    if(!is.null(omitTrap)){
      checkb <- b[,omitTrap]
      if(is.matrix(checkb)) ww <- which(is.na(rowSums(checkb)))
      if(!is.matrix(checkb))ww <- which(is.na(checkb))
      if(length(ww) > 0)b <- b[-ww,]
      if(nrow(b) == 0)next
    }

    allCodes <- unique(c(allCodes,colnames(b)))

    uname <- paste(path,allplots[j],'_UTM.txt',sep='')
    uf   <- read.table(uname,header=T)

    if(allplots[j] == 'CW_218')b <- b[b[,'trapnum'] <= 20,]

    snames <- seedNames[seedNames %in% colnames(b)]
    if(length(snames) == 0 & !plotnames[j] %in% colnames(mtree))next

    if(length(snames) > 0){
      count <- b[,colnames(b) %in% snames]
      if(length(snames) > 1)count <- rowSums(count,na.rm=T)
    }
    if(length(snames) == 0)count <- rep(0,nrow(b))

    trap  <- b[,'trapnum']
    ntrap <- max(trap,na.rm=T)
    cyr   <- b[,'year']
    cmo   <- b[,'month']
    
    yri   <- match(cyr,yrvec)
    yri[cmo < firstMonth] <- yri[cmo < firstMonth] - 1   #collected before Sept assigned to previous yr
    yri[yri < 1] <- 1

    trapByYr <- matrix(NA,ntrap,nyr)
 
    seeds <- tapply(count,list(trap=trap,yri=yri),sum,na.rm=T)
    trapByYr[,as.numeric(colnames(seeds))] <- seeds

    ks <- as.numeric(rownames(seeds))
    ut <- uf[match(uf[,'trapnum'],ks),c('UTMx','UTMy')]
    ik <- which(ks %in% rownames(ut))
    if(length(ik) < length(ks)){
      message( paste('missing UTMs: ',ks[!ks %in% ik],sep='') )
      trapByYr <- trapByYr[ik,]
    }

    jk <- cbind( rep(j,nrow(trapByYr)),  ks[ik], ut)
    colnames(jk) <- c('j','k','x','y')
    trapByYr <- cbind(jk,trapByYr)
    seedMat  <- rbind(seedMat,trapByYr)

    meanseed <- colSums(trapByYr,na.rm=T)/colSums(trapByYr*0+1,na.rm=T)/trapArea[j]
    perM2    <- rowBind(perM2, round(meanseed,2)[-c(1:4)], allplots[j])
    rss      <- colSums(trapByYr*0+1,na.rm=T)[-(1:4)]

    trapPlotYr[j,] <- rss
  }

  colnames(seedMat)[-c(1:4)] <- yrvec
  colnames(perM2) <- yrvec

  list(seedMat = seedMat[,-c(1:4)], perM2 = perM2, seedData = as.matrix(seedMat[,1:4]),
       UTMseed = seedMat[,c('x','y')],trapPlotYr = trapPlotYr)
}


checkSeedUTM <- function(j,uf){
   
    utree <- UTMtree[treeData[,'j'] == j,]
    xtree <- treeData[treeData[,'j'] == j,c('x','y')]

    tsamp <- sample(1:nrow(utree),300,replace=T)

    dx <- matrix(uf[,'x'],nrow(uf),length(tsamp)) - 
          matrix(xtree[tsamp,'x'],nrow(uf),length(tsamp),byrow=T) 
    ux <- matrix(utree[tsamp,'x'],nrow(uf),length(tsamp),byrow=T) + dx
    ux <- rowMeans(ux)

    dy <- matrix(uf[,'y'],nrow(uf),length(tsamp)) - 
          matrix(xtree[tsamp,'y'],nrow(uf),length(tsamp),byrow=T) 
    uy <- matrix(utree[tsamp,'y'],nrow(uf),length(tsamp),byrow=T) + dy
    uy <- rowMeans(uy)

   newf <- uf
   newf[,c('UTMx','UTMy')] <- round(cbind(ux,uy),1)
   write.table(newf,paste(seedPath,plotnames[j],'_UTMa.txt',sep=''),quote=F,row.names=F)

   plot(utree[,'x'],utree[,'y'],cex=.3)
   points(uf[,'UTMx'],uf[,'UTMy'],col=2)
   title(plotnames[j])

   dev.print(device=postscript,file=paste('seed',plotnames[j],'map.ps',sep=''),width=6,horizontal=F)

   invisible(newf)
}



getTreeMat <- function(){

  nv   <- length(censusVariables)
  minv <- 0
  maxv <- 10000
  INCR <- T
  dflt <- NULL

  vnames <- paste(censusVariables,'Vec',sep='')
  unames <- paste(censusVariables,'Mat',sep='')

  for(k in 1:nv){

     kb <- get( unames[k] )

     if(!censusVariables[k] %in% noInterpolation){

        minVal <- minv
        maxVal <- maxv
        INCVAL <- INCR
        dff     <- dflt
        tiny    <- .01
        if(censusVariables[k] == 'canopy'){
          minVal <- 0
          INCVAL <- F
          dff <- 0
        }
        if(censusVariables[k] == 'diam'){
           minVal <- 0
           maxVal <- 300
           tiny <- minDinc
        }
        if(censusVariables[k] == 'ht'){
           minVal <- 0
           maxVal <- 70
           tiny <- -10
           hall <- allomConvert(diamMat,treeNames[treeData[,'spec']],
                                allomFile=paste(dataPath,htAllomFile,sep=''),
                                codeColumn='species',
                                'htInt','htSlope',defaultSpec=defaultAllomSpec)
           kb[is.na(kb) & !is.na(hall)] <- hall[is.na(kb) & !is.na(hall)]
        }
        kb <- interpRows(kb,startIndex=firstTime,INCREASING=INCVAL,
                  minVal=minVal,maxVal=maxVal,defaultValue=dff,tinySlope=tiny)
        
        if(censusVariables[k] == 'ht'){  #if seedlings, some will have ht, but not diam
        
          tmp <- allomConvert(kb,treeNames[treeData[,'spec']],
                            allomFile=paste(dataPath,htAllomFile,sep=''),
                            codeColumn='species',
                            'htInt','htSlope',defaultSpec=defaultAllomSpec, invert=T)
          
          diamMat[is.na(diamMat) & !is.na(tmp)] <- tmp[is.na(diamMat) & !is.na(tmp)]
          tmp <- mat2Vec( diamMat ,GETINDEX=T,lastTime=lastTime)
          tin <- tmp$index
          
          assign( 'diamVec' , tmp$x )
        }
        
        
        if(censusVariables[k] == 'canopy') kb <- round(kb,0)
     }

     tmp <- mat2Vec( kb ,GETINDEX=T,lastTime=lastTime)
     tin <- tmp$index

     assign( vnames[k] , tmp$x )
  }



  ii <- c(1:nrow(tin))

  # cull measurements before first time, after lasttime
  ww <- numeric(0)
  for(i in 1:ntree){
     ww <- c(ww, which(tin[,1] == i & tin[,2] < firstTime[i]) )
     ww <- c(ww, which(tin[,1] == i & tin[,2] > treeDeath[i,'t2']) )
  }

  if(length(ww) > 0)tin <- as.matrix(tin[-ww,])
  tin <- as.matrix(tin)

  ii <- ii[!ii %in% ww]
     
  jplot  <- match(treeData[,'plot'],plotnames)
  j      <- matrix(jplot,ntree,nyr)[tin]
  reg    <- matrix(treeData[,'reg'],ntree,nyr)[tin]
  spec   <- matrix(treeData[,'spec'],ntree,nyr)[tin]
  monoec <- matrix(treeData[,'monoec'],ntree,nyr)[tin]

  out <- cbind(reg,spec,monoec,tin,j)

  for(k in 1:nv)out <- cbind(out,get(vnames[k])[ii])
  colnames(out) <- c('reg','spec','monoec','tindex','t','j',censusVariables)

  diamFirst <- diamMat[ cbind(c(1:ntree),firstTime) ]

  smat <- matrix(treeDeath[,'t2'],ntree,nyr)
  tmat <- matrix(1:nyr,ntree,nyr,byrow=T)
  smat[is.na(smat) | smat > tmat] <- 1
  smat[smat != 1] <- 0
  survVec <- smat[tin]

  list(treemat = out, surv = survVec)
}


diam2dinc <- function(diam,minInc=0,firstTime=firstTime,ntree=nrow(diam),nyr=ncol(diam)){
  
  #if diam is a vector requires index matrix tin in main program

  diamMat <- diam

  if(!is.matrix(diam)){
    diamMat <- matrix(NA,ntree,nyr)
    diamMat[tin] <- diam
  }

  dlast <- lastTimeMat(diamMat,firstTime,INC=T,minInc)

  xl <- diamMat - dlast
  xx <- cbind(xl[,-1],xl[,ncol(xl)])
  colnames(xx) <- colnames(diam)
  xx
}

cleanFirstTime <- function(x,firstTime){

  f <- matrix(firstTime,ntree,nyr)
  m <- matrix(c(1:nyr),ntree,nyr,byrow=T)
  x[m < f] <- NA
  x
}
  



getSeedMat <- function(){

  mm <- matrix(1:nyr,nplot,nyr,byrow=T)
  mm[match(rownames(seedM2),plotnames),] <- mm[match(rownames(seedM2),plotnames),]* (seedM2*0 + 1)
  fy <- apply(mm,1,min,na.rm=T)
  fy[!plotnames %in% rownames(seedM2)] <- nyr + 1

  tmp       <- mat2Vec( seedMat ,GETINDEX=T )
  count     <- tmp$x
  seedIndex <- as.matrix(tmp$index)
  colnames(seedIndex) <- c('sindex','t')

  j <- matrix( seedData[,'j'],ntrap,nyr )[seedIndex]
  k <- matrix( seedData[,'k'],ntrap,nyr )[seedIndex]
  
  ss <- cbind(seedIndex,j,k,count)
  ww <- which(fy[ss[,'j']] <= ss[,'t'])
  ss[ww,]
}


inDendroBand <- function( path = paste(dataPath,dendPath,sep='') ){

  afiles <- list.files(path,pattern='_active.txt')
  backup <- grep("~",afiles)
  if(length(backup) > 0)afiles <- afiles[-backup]
  plot  <- unlist(strsplit(afiles,'_active.txt') ) 

  nfile <- length(plot)

  dendMat <- dimat <- nMat <- matrix(NA,ntree,nyr)

  for(k in 1:nfile){

    kplot <- match(plot[k],plotnames)
    wtree <- which(treeData[,'j'] == kplot)
    treek <- treeData[wtree,]
    if(nrow(treek) == 0)next

    b <- read.table(paste(path,afiles[k],sep=''),header=T)

    ws <- which(b[,'site'] == plotnames[kplot] & b[,'species'] %in% treeNames)
    wb <- which(b[,'ID'] %in% treek[,'ID'])

    if(length(ws) > length(wb)){
      print('unknown increment cores')
      print(afiles[k])
      print(b[ws[!ws %in% wb],])
    }

    if(length(wb) == 0)next

    b <- b[wb,]
    wk <- match(b[,'ID'],treek[,'ID'])       

    tmp        <- dateColumn(b,'diam',yrvec)
    diamB      <- tmp$x      
    dimat[wtree[wk],] <- diamB

    tmp        <- dateColumn(b,'dinc',yrvec)
    dincB      <- tmp$x      
    dendMat[wtree[wk],] <- dincB
}

 dendMat[dendMat < 0] <- 0

 list(dendMat = dendMat, diamMat = dimat)
}


inTreeIncr <- function(){

  path   <- paste(dataPath,ringPath,sep='')
  afiles <- list.files(path,pattern='_active.txt')
  backup <- grep("~",afiles)
  if(length(backup) > 0)afiles <- afiles[-backup]
  files  <- unlist(strsplit(afiles,'_active.txt') ) 

  parsf  <- matrix( unlist(strsplit(files,'_')) ,ncol=3,byrow=T)
  plot   <- paste(parsf[,1],parsf[,2],sep='_')
  years  <- parsf[,3]

  nfile <- length(years)

  dincMat <- nMat <- matrix(0,ntree,nyr)

  for(k in 1:nfile){

    kplot <- match(plot[k],plotnames)
    treek <- treeData[treeData[,'j'] == kplot,]
    if(nrow(treek) == 0)next

    print(afiles[k])

    b <- read.table(paste(path,afiles[k],sep=''),header=T)

    wb <- which(b[,'ID'] %in% treek[,'ID'])
    if(length(wb) == 0)next

    b <- b[wb,]

    bid <- sort(unique(b[,'ID']))
    nb  <- length(bid)
    for(i in 1:nb){

      bi <- b[b[,'ID'] == bid[i] & b[,'yr'] %in% yrvec,]

      wr <- which(treeData[,'j'] == kplot & treeData[,'ID'] == bid[i])
      wc <- match(bi[,'yr'],yrvec)

      dincMat[wr,wc] <- dincMat[wr,wc] + bi[,'cm']
      nMat[wr,wc]    <- nMat[wr,wc] + 1

      wlo <- which(2*bi[,'cm'] < minDinc)
      if(length(wlo) > 0){
         message("growth rate < minDinc")
         print(afiles[k])
         print(bi[wlo,'cm'])
      }

   }
  }

  dincMat <- (dincMat/nMat) * 2         #radial measurements to diameter growth for increment cores 
  dincObs <- which(nMat > 0,arr.ind=T)

  list(dincMat = dincMat, dincObs = dincObs)
}


distanceObjects <- function(){


  distall <- distmat(UTMtree[,1],UTMtree[,2],UTMseed[,1],UTMseed[,2])

  closeTrap <- round(apply(distall,2,min,na.rm=T),1)  #distance to closest trap

  closeTrees <- matrix(NA,nrow(seedData),10)
  dmax       <- apply(diamMat,1,max,na.rm=T)

  for(k in 1:nrow(seedData)){

    dk <- distall[k,]
    dk[dmax < minMatrDiam] <- 1e+10
    ord <- order(dk,decreasing=F)[1:10]
    closeTrees[k,] <- ord
  }

  for(j in 1:nplot){

    ws <- which(seedData[,'j'] == j)
    wt <- which(treeData[,'j'] != j)
    if(length(ws) > 0 & length(wt) > 0)distall[ws,wt] <- NA
  }

  list(closeTrap = closeTrap, distall = distall, closeTrees = closeTrees)
}


sexChange <- function(char=NULL,num=NULL){

  if(is.null(char)){
    if(!is.matrix(num))num <- as.matrix(num)
    x <- matrix(character(0),nrow(num),ncol(num))
    x[num == 1] <- 'un'
    x[num == 2] <- 'n'
    x[num == 3] <- 'r'
    x[num == 4] <- 'f'
    x[num == 5] <- 'm'
    colnames(x) <- colnames(num)
  }
  if(is.null(num)){
    if(is.matrix(char)) x <- matrix(NA,nrow(char),ncol(char))
    if(!is.matrix(char))x <- rep(NA,length(char))
    x[char == 'un'] <- 1
    x[char == 'n'] <- 2
    x[char == 'r'] <- 3
    x[char == 'f'] <- 4
    x[char == 'm'] <- 5
    colnames(x) <- colnames(char)
  }
  x
}

matureStatus <- function(){  #requires sexMat

  #last observed immature

   immCode <- sexChange(char='n')
   matCode <- sexChange(char=c('r','f','m'))
   unCode  <- sexChange(char='un')

   tmat <- matrix(c(1:nyr),ntree,nyr,byrow=T)

   imMat <- matrix(0,ntree,nyr)                  #not repro
   imMat[sexMat == immCode] <- 1
   imMat <- -(t( apply(imMat,1,cumsum) ) - 1)*tmat
   lastImm <- apply(imMat,1,max)
   lastImm[lastImm == nyr] <- NA

   maMat  <- imMat*0                             #repro
   maMat[sexMat %in% matCode] <- 1
   maMat <- t(apply(maMat,1,cumsum) )
   maMat[maMat == 0] <- 9999
   firstMat <- apply(maMat*tmat,1,min)
   firstMat[firstMat > nyr] <- NA
  
   sexObs <- sexMat*0
   sexObs[sexMat == unCode] <- -1
   sexObs[sexMat == immCode] <- 0
   sexObs[sexMat %in% matCode] <- 1
 #  sexObs[sexMat == sexChange(char='m')] <- 3
 #  sexObs[sexMat == sexChange(char='f')] <- 4

   list(lastImm = lastImm, firstMat = firstMat, sexObs = sexObs)
}

plotUTMrange <- function(jplot){  # jplot - index for plots to include

    utmxyz <- numeric(0)

    xr <- yr <- numeric(0)
    for(j in jplot){

     utmxyz <- rbind(utmxyz,UTMtree[treeData[,'j'] == j,])
     xr <- range(xr,UTMhull[[j]][,1],na.rm=T)
     yr <- range(yr,UTMhull[[j]][,2],na.rm=T)
    }

    list(urange = cbind(xr,yr), utmxyz = utmxyz)
}


#ynamess <- function(){

#  ynames <- character(0)
 # if('dinc' %in% ynames)ynames <- 'D'
 # if('hinc' %in% ynames)  ynames <- c(ynames,'H')
 # if('fecn' %in% ynames)ynames <- c(ynames,'F')
 # ynames
#}

  

inTreeCensus <- function(GETSHADE=F){

  require(grDevices)

  requireTree <- c('ID','tag','species','censinyr','growinyr','deathyr','censoryr','UTMx','UTMy')
  requireLing <- c('species','deathyr','censoryr')
  utmNames    <- c('UTMx','UTMy')
  deathNames  <- c('birthyr','censinyr','growinyr','deathyr','censoryr')

  kname <- c('open','sex','canopy','diam')
  if(HT)kname <- c(kname,'ht')
  kyr   <- c('openYr','sexYr','canYr','sampleYr','sampleYr')

  htObserve <- numeric(0)

  UTM <- treeData <- treeDeath <- treeVar <- shade <- htMat <- numeric(0)
  treeSpec <- character(0)
  UTMhull <- vector( 'list',nplot )
  names(UTMhull) <- plotnames
  subplotNames <- UTMhull

  nv <- length(censusVariables)

  for(k in 1:nv)assign(censusVariables[k],numeric(0))

  tpath <- paste(dataPath,treePath,sep='')
  lpath <- paste(dataPath,lingPath,sep='')

  tfiles <- list.files(tpath,pattern='_active.txt')   #in tree directory
  backup <- grep("~",tfiles)
  if(length(backup) > 0)tfiles <- tfiles[-backup]
  tplots <- filename2plot(tfiles)

  lplots <- character(0)     
  if(LING){
    lfiles <- list.files(lpath,pattern='_active.txt')   #in seedling directory
    backup <- grep("~",lfiles)
    if(length(backup) > 0)lfiles <- lfiles[-backup]
    lplots <- filename2plot(lfiles)
  }

  baTable <- numeric(0)

  for(j in 1:nplot){

    btree <- bling <- ball <- hull <- numeric(0)

    wt     <- which(tplots == plotnames[j])
    wl     <- which(lplots == plotnames[j])

    if(length(wt) == 0 & length(wl) == 0)next

    print(plotnames[j])

    if(length(wt) > 0){          #########tree files

      ff     <- tfiles[wt]
      fname  <- paste(tpath,ff,sep='')
      ball   <- read.table(fname,header=T,sep="\t",fill=T)
      btree  <- ball[ball[,'species'] %in% treeNames,]

      #tree plots have UTM for each tree    
      missingCol(btree,requireTree,'stop')     
      btree    <- missingCol(btree,requireLing,action='add')
      btree    <- missingCol(btree,'subplot','add',value='Z')
      btree    <- missingCol(btree,'type','add',value='tree')
      btree    <- missingCol(btree,'birthyr','add',value=NA)

      dcols <- grep('diam',colnames(ball))
      hull <- ball[,c('UTMx','UTMy')]

      if(length(dcols) < 2)btree <- ball <- numeric(0)  #at least 2 censuses

      if(length(dcols) >= 2){

        hull <- ball[,c('UTMx','UTMy')]

        dmat   <- pi*(ball[,dcols]/2)^2 /10000/plotArea[plotnames[j],nyr]
        bmax   <- apply(dmat,1,max,na.rm=T)
        bmax[!is.finite(bmax)] <- 0

        baj <- t( round(as.matrix( unlist( by(bmax,ball[,'species'],sum,na.rm=T) ) ),5) )

        baTable <- appendMatrix(baTable,baj)
        rownames(baTable)[nrow(baTable)] <- plotnames[j]

        print(fname)

        #rm trees without UTMs
        ww <- which( is.na(rowSums(btree[,utmNames])) | btree[,'deathyr'] <= yrvec[1] )
        if(length(ww) > 0)btree <- btree[-ww,]
       }
    }

    if(length(wl) > 0){          ###############seedling files

      ff     <- lfiles[wl]
      fname  <- paste(lpath,ff,sep='')
      bmat   <- read.table(fname,header=T,sep="\t",fill=T)
      bling  <- bmat[bmat[,'species'] %in% treeNames,]

      print(fname)

      missingCol(bling,requireLing,'stop')     

      nb <- NA
      if(nrow(bling) > 0){
        nb <- 1:nrow(bling)
        ww <- grep('diam',colnames(bling))
        bling[,ww] <- bling[,ww]*.1       #seedling diams from mm to cm
        if(HT){
          ww <- grep('ht',colnames(bling))
          bling[,ww] <- bling[,ww]*.1       #seedling ht from cm to m
        }
      }

      bling <- missingCol(bling,'ID','add',value=nb)
      bling <- missingCol(bling,'tag','add',value=nb)
      bling <- missingCol(bling,requireTree,'add')
      bling <- missingCol(bling,'subplot','add',value='Z')
      bling <- missingCol(bling,'type','add',value='ling')
      if(!'birthyr' %in% colnames(bling)){
        if('plantyr' %in% colnames(bling))colnames(bling)[colnames(bling) == 'plantyr'] <- 'birthyr'
        if(!'plantyr' %in% colnames(bling))missingCol(bling,'birthyr','add',value='NA')
      }    

      #seedling plots have separate UTM file

      subs <- sort(unique(bling[,'subplot']))
      nsub <- length(subs)
      si   <- match(bling[,'subplot'],subs)
      subplotNames[[j]] <- as.character(subs)

      ufile <- paste(lpath,plotnames[j],'_UTM.txt',sep='')
      umat  <- read.table(ufile,header=T,sep="\t",fill=T)

      t2u   <- match(subs[si],umat[,'subplot'])
      uj    <- umat[t2u,c('UTMx','UTMy')]
      bling[,utmNames] <- uj
      hull <- rbind(hull,uj)
    }###########################################33

    hull <- hull[is.finite(hull[,1]) & is.finite(hull[,2]),]
    UTMhull[[j]] <-  hull[chull(hull),] 

    notree <- noling <- F
    if( is.null(nrow(btree)) )notree <- T
    if( is.null(nrow(bling)) )noling <- T
    if(!notree)if(nrow(btree) == 0)notree <- T
    if(!noling)if(nrow(bling) == 0)noling <- T

    if(notree & noling)next

    #merge tree and seedling diams and hts

    sampleYr <- numeric(0)

    ltyp <- c('btree','bling')
    kvar <- 'diam'
    knam <- 'diamB'
    if(HT){
      kvar <- c(kvar,'ht')
      knam <- c(knam,'htB')
    }

    diamB <- htB <- numeric(0)
    sexYr <- canYr <- openYr <- numeric(0)

    for(l in 1:2){                

      lb  <- get(ltyp[l])
      if(length(lb) == 0)next
      sexYr <- c(sexYr,dateColumn(lb,'sex',yrvec)$yr)
      canYr <- c(canYr,dateColumn(lb,'canopy',yrvec)$yr)
      if(ltyp[l] == 'bling')openYr <- c(openYr,dateColumn(lb,'open',yrvec)$yr)

      for(k in 1:length(kvar)){
        tmp        <- dateColumn(lb,kvar[k],yrvec)
        assign( knam[k],appendMatrix(get(knam[k]),tmp$x) )
        sampleYr   <- c(sampleYr,tmp$yr)
      }
    }
    sexYr   <- sort(unique(sexYr))
    canYr   <- sort(unique(canYr))
    openYr   <- sort(unique(openYr))
    sampleYr   <- sort(unique(sampleYr))
    sampleTime <- match(sampleYr,yrvec)

    kcols <- character(0)

    for(k in 1:length(kname)){

      if(length(btree) > 0){
          btree <- missingCol(btree,paste(kname[k],get(kyr[k]),sep=''),'add',NA)
          kcols <- c(kcols,colnames(btree)[grep(kname[k],colnames(btree))])
      }
      if(length(bling) > 0){
          bling <- missingCol(bling,paste(kname[k],get(kyr[k]),sep=''),'add',NA)
          kcols <- c(kcols,colnames(bling)[grep(kname[k],colnames(bling))])
      }
    }
    kcols <- sort(unique(kcols))
    kcols <- c(union(requireTree,requireLing),'birthyr','type','subplot',kcols)

    b <- numeric(0)
    if(length(btree) > 0)b <- btree[,kcols]
    if(length(bling) > 0)b <- rbind(b,bling[,kcols])
   
    xx <- diamB
    if(HT)xx <- cbind(diamB,htB)
    ww <- which( is.na(rowMeans(xx,na.rm=T)) )              # trees missing diam and ht
    
    if(length(ww) > 0){
        b <- b[-ww,]
        if(nrow(b) == 0)next
        diamB <- diamB[-ww,]
        if(HT)htB   <- htB[-ww,]
    }

    t2 <- t1 <- sampleTime[ match(b[,'deathyr'],sampleYr) ]              #determine death interval
    t2[is.finite(b[,'censoryr'])] <- NA                                  #censored
    
    wi <- which(is.finite(t2))
    if(length(wi) > 0){
      t1     <- t2
      td     <- match(t2[wi],sampleTime) - 1
      td[td < 1] <- 1
      t1[wi] <- sampleTime[ td ]
    }
    wi <- which((t2 - t1) < 1)
    if(length(wi) > 0)t2[wi] <- t1[wi] + 1
    
    deathInt <- cbind(t1,t2)

   # trees not dead or censored (lost)
    lastSamp <- matrix(c(1:nyr),nrow(diamB),nyr,byrow=T)*(diamB*0 + 1)

    if(HT){
      lastSampH <- matrix(c(1:nyr),nrow(htB),nyr,byrow=T)*(htB*0 + 1)
      lastSamp  <- apply(cbind(lastSamp,lastSampH),1,max,na.rm=T)
    }
    lastSamp  <- apply(lastSamp,1,max,na.rm=T)

    wm <- which( is.na(deathInt[,2]) & is.na(b[,'censoryr']) )

    lostSince     <- rep(NA,nrow(b))
    lostSince[wm] <- lastSamp[wm]
    lostSince[lostSince == max(sampleTime)] <- NA

    if(is.finite(lostTime)){
      wl <- which( (nyr - lostSince) > lostTime )
      if(length(wl) > 0){
        yl <- lostSince[wl]
        deathInt[wl,] <- cbind(yl,sampleTime[ findInterval(yl,sampleTime) + 1] )
        b[wl,'deathyr'] <- yrvec[deathInt[wl,2]]
      }
    }

    firstYr   <- apply(b[,c('birthyr','censinyr','growinyr')],1,min,na.rm=T)
    dt <- b[,'deathyr']  - firstYr
    dc <- b[,'censoryr'] - firstYr
    dn <- max(yrvec) - firstYr
    firstTime <- match(firstYr,yrvec)
    ww <- which(dt < minYears | dc < minYears | dn < minYears)
    if(length(ww) > 0){
        b         <- b[-ww,] 
        diamB     <- row2Mat( diamB[-ww,] )
        if(HT)htB <- row2Mat( htB[-ww,] )
        deathInt  <- row2Mat( deathInt[-ww,] )
        if(nrow(b) == 0)next
    }

    if(HT)htObs <- htB

    if(GETSHADE){

      shj <- matrix(NA,nrow(diamB),ncol(diamB))

      wling <- which(b[,'type'] == 'ling')
      if(length(wling) > 0){         
         if(HT){
           h0  <- row2Mat( htB[wling,] )
           d2h <- .15 + 8.1*diamB[wling,]
           h0[!is.finite(h0) & is.finite(d2h)] <- d2h[!is.finite(h0) & is.finite(d2h)]
           htB[wling,] <- interpRows(h0,firstTime[wling],minVal=0)
         }
         
         wlite <- grep('open',colnames(b))         #gap treatment as light index

         if(length(wlite) > 0){
            sh <- dateColumn(b,'open',yrvec)$x
            sh[sh == 0] <- .03
            sh[sh == 1] <- .95
            shj <- sh
         }
      }

     wtree <- which(b[,'type'] == 'tree')
      if(length(wtree) > 0){         
        tmp <- getShade(b[wtree,],ball)
        if(HT)htB[wtree,] <- tmp$htmat
        shj[wtree,] <- tmp$shadeIndex
      }
      shade <- rbind(shade,shj)
    }

    if(HT)htMat <- rbind(htMat,htB)

    spec   <- match(b[,'species'],treeNames)
    monoec <- female <- rep(0,nrow(b))

    if('fecn' %in% ynames){
      monoec[mono[spec]] <- 1
      female <- monoec
      female[monoec == 0] <- 0
    }

    wreg <- match( unlist(strsplit(plotnames[j],'_'))[1] , names(plotgroups) ) #region for plot
    reg  <- rep(wreg,nrow(b))

    for(k in 1:nv){
      kb <- dateColumn(b,censusVariables[k],yrvec)$x
      if(censusVariables[k] == 'sex')kb <- sexChange(char=kb,num=NULL)
      if(censusVariables[k] == 'ht')kb <- htB
      assign(censusVariables[k], appendMatrix( get(censusVariables[k]), kb ) )
    }   

    wc <- which( utmNames %in% colnames(b))
    if(length(wc) < 2)warning('UTMs missing from file')

    ut     <- as.matrix( b[,utmNames[wc]] )
    UTM    <- rbind( UTM,   ut  )
    
    treeD    <- b[,treeID]
    plot     <- rep(plotnames[j],nrow(b)) 
    jp       <- rep(j,nrow(b))
    subplot  <- b[,'subplot']
    type     <- b[,'type']
    treeD    <- cbind(reg,jp,plot,subplot,spec,monoec,female,type,treeD)
    treeData <- rbind( treeData, treeD )
    treeSpec <- c(treeSpec,as.character(b[,'species']) )

  #  birthyr <- rep(NA,length(jp))
    if(length(jp) == 1)jdead <- row2Mat( c(unlist(b[,deathNames]),deathInt) )
    if(length(jp) >  1)jdead <- cbind(as.matrix(b[,deathNames]),deathInt )
    colnames(jdead) <- c(deathNames,'t1','t2')

    treeDeath <- rbind(treeDeath,jdead)

    if(HT)htObserve <- rbind(htObserve,htObs)

 }################end plot loop


  tab <- tabulate(treeData[,'spec'],nbin=nttype)
  ww  <- which(tab < minTrees & tab > 0)

  if(length(ww) > 0){
    ws <- which(treeData[,'spec'] %in% ww)
    treeData  <- treeData[-ws,]
    treeDeath <- treeDeath[-ws,]
    treeSpec  <- treeSpec[-ws]
    UTM       <- UTM[-ws,]
    shade     <- shade[-ws,]

    for(k in 1:length(censusVariables)){
      kname <- paste(censusVariables[k],'Mat',sep='')
      assign( censusVariables[k], get(censusVariables[k])[-ws,] )
    }
  }

  treeDeath <- as.matrix(unlist(treeDeath))

  newNames  <- unique(treeNames[ treeData[,'spec'] ] )

  wn   <- which(treeNames %in% newNames)
  wnew <- match(treeNames[treeData[,'spec']],treeNames[wn])

  treeNames <- treeNames[wn]
  treeData[,'spec'] <- wnew

  colnames(treeData)[colnames(treeData) == 'jp'] <- 'j'

  stems <- table(treeSpec,treeData[,'plot'])
  wa    <- which(treeData[,'type'] == 'tree')
  wj    <- which(treeData[,'type'] == 'ling')
  stemsTree  <- table(treeSpec[wa],treeData[wa,'plot'])
  stemsLing    <- table(treeSpec[wj],treeData[wj,'plot'])

  vlist <- vector('list',nv)
  for(k in 1:nv)vlist[[k]] <- get(censusVariables[k])
  names(vlist) <- censusVariables

  if(!'sex' %in% censusVariables){
   sex <- rep(1,length(diam))
   vlist[[nv+1]] <- sex
  }

  colnames(UTM)[colnames(UTM) == 'UTMx'] <- 'x'
  colnames(UTM)[colnames(UTM) == 'UTMy'] <- 'y'
 # colnames(UTM)[colnames(UTM) == 'UTMz'] <- 'z'

  list( varList = vlist, UTM = UTM, UTMhull = UTMhull, treeData = treeData, htMat = htMat,
        treeDeath = as.matrix(treeDeath), stems = stems, shade = shade, treeNames = treeNames,
        subplotNames = subplotNames, htObs = htObserve, baTable = baTable,
        stemsTree = stemsTree, stemsLing = stemsLing)
}


dateColumn <- function(b,vname,yrvec=NULL,nyr=length(yrvec)){    

  #returns matrix where columns are years for vname, rows from b
  
  nr <- nrow(b)
  x  <- matrix(NA,nr,nyr); colnames(x) <- yrvec
  
  wc <- grep(vname,colnames(b))
  if(length(wc) == 0)return( list(x = x,y = 0) )
  
  cc <- matrix(unlist(strsplit(colnames(b)[wc],vname)),ncol=2,byrow=T)[,2] 
  yr <- as.numeric(cc) 
  if(is.na(yr[1]))return( list(x = x,y = 0) )

  if(is.null(yrvec))yrvec <- yr
  nyr <- length(yrvec)
  
#  yr <- yr[yr %in% yrvec]

  wy <- match(yr,yrvec)
  wc <- wc[is.finite(wy)]
  wy <- wy[is.finite(wy)]

  wc <- wc[order(yrvec[wy])]
  wy <- wy[order(yrvec[wy])]
  
  if(length(wc) == 0){
    warning('no obs within year range',immediate.=T)
    return( list(x = x, yr = yr) )
  }

  btmp <- matrix(unlist(b[,wc]),nr,length(wc))

  x[,wy]  <- btmp
  
  yr <- yrvec[wy]
  
  tmp <- colSums(x,na.rm=T)
  tt  <- as.numeric(names(tmp)[tmp > 0])
  yr  <- yr[yr %in% tt]

 # if(length(wc) == 1){ x <- matrix(x,ncol=1); colnames(x) <- yr }
 # yrange  <- c(yr[1],max(yr))
  list(x = x, yr = yr)
}
  


seedLikeF <- function(propsd){ #sample fecn

  atab   <- matrix(0,mplot,nt)
  tabnow <- tabnew <- matrix(0,mplot,max(ntrap))

  tmp   <- propFecn(propsd)
  fnow  <- tmp[,1]
  fnew  <- tmp[,2]

  for(t in 1:length(yrvec)){ 

    tnow <- which(treeRows(matr=1,yr=t,sex=1))

    if(length(tnow) == 0)next

    treenow <- row2Mat(treemat[tnow,])
 
    jplot   <- unique(treenow[,"j"]) 
    seednow <- seedmat[seedmat[,"t"] == t & seedmat[,"j"] %in% jplot,]

    pnow  <- getprob1(fnow[tnow],upar,treenow[,"tindex"],seednow[,"sindex"],seednow)$like 
    pnew  <- getprob1(fnew[tnow],upar,treenow[,"tindex"],seednow[,"sindex"],seednow)$like 

    tabnow[ cbind(seednow[,'j'],seednow[,'k']) ] <- pnow
    tabnew[ cbind(seednow[,'j'],seednow[,'k']) ] <- pnew

    pbynow <- rowSums( row2Mat(tabnow[jplot,]) )
    pbynew <- rowSums( row2Mat(tabnew[jplot,]) )

    a  <- exp(pbynew - pbynow)
    z  <- runif(length(a),0,1)
    wa <- which(z < a)
    if(length(wa) > 0){
      for(j in jplot[wa]){
        wr <- treeRows(plt=j,yr=t)
        fnow[wr] <- fnew[wr]
        atab[j,t] <- atab[j,t] + 1
      }
    }
  }
  list(fnow = fnow, accept = atab)
}

mapHydro <- function(xutm,yutm,dataPath='datafiles/',
                     mapPath='datafiles/hydro/wetness/',reg='DF',TOPO=F,
                     buffer=0,conInt=NULL,reverseColor=F,BOUND=F,
                     xfactor=1,nlevs=5){
  
  #reg - region (DF, CW, MH)
  
  BOUND <- F
  FLAG  <- T
  
  xr <- range(xutm)
  yr <- range(yutm)
  xr <- c(xr[1] - buffer,xr[2] + buffer)
  yr <- c(yr[1] - buffer,yr[2] + buffer)
  
  
  ncol <- 100
  colseq <- terrain.colors(ncol)
  if(reverseColor) colseq  <- rev(colseq)

  sc <- max( c(diff(xr),diff(yr)) )/5     #m per inch
  mapwide <- diff(xr)


  mapSetup(xr,yr,scale=sc)
  plot(xr,yr,xlab='long',ylab='lat',cex=.1)

  mfile <- paste(mapPath,reg,'.asc',sep='')
  
  print(mfile)
  
  if(!file.exists(mfile)){
    warning(paste('missing map file:',mfile,'map omitted'))
    return( list(FLAG = F) )
  }

 # BOUND <- F
 # bfile <- paste(dataPath,boundaryFile,sep='')
 # if(file.exists(bfile)){
 #   boundaries <- read.table(bfile,header=T)
 #   wj         <- which(boundaries[,'plot'] %in% jplots)
 #   if(length(wj) > 0)BOUND <- T
 # }

  tmp <- read.asciigrid.r(mfile,as.image=T, plot.image=F)
  x <- tmp$x
  y <- tmp$y
  z <- tmp$z*xfactor
  

  z  <- z[x >= xr[1] & x <= xr[2],y >= yr[1] & y <= yr[2]]
  xl <- x[x >= xr[1] & x <= xr[2]]
  yl <- y[y >= yr[1] & y <= yr[2]]
  
  if(length(z) == 0){
    warning('no data in map coordinates')
    return( list(FLAG = F))
  }

  zrange <- range(z,na.rm=T)
  zlevs  <- seq(zrange[1],zrange[2],length=(nlevs+1))

  z[z < zrange[1]] <- zrange[1]
  z[z > zrange[2]] <- zrange[2]

  zr <- round( log10(diff(zrange) ), 1)
  
  print(zrange)

  image(xl,yl,z,add=T,col=colseq)
  
  
  # return values
  xygrid <- expand.grid(xl,yl)
  zgrid  <- as.vector(z)
  dxy    <- distmat(xygrid[,1],xygrid[,2],xutm,yutm)
  mind   <- apply(dxy,1,which.min)
                 

  if(TOPO){
  
    mfile <- paste(paste('datafiles/DEM/',reg,'.asc',sep=''))
  
    if(!file.exists(mfile)){
      warning(paste('missing map file:',mfile,'map omitted'))
      return( list(FLAG = F) )
    }
  
    tmp <- read.asciigrid.r(mfile,as.image=T, plot.image=F)
    x <- tmp$x
    y <- tmp$y
    z <- tmp$z*xfactor
  
    z  <- z[x >= xr[1] & x <= xr[2],y >= yr[1] & y <= yr[2]]
    xl <- x[x >= xr[1] & x <= xr[2]]
    yl <- y[y >= yr[1] & y <= yr[2]]
  
    zrange <- range(z,na.rm=T)
    zlevs  <- seq(zrange[1],zrange[2],length=(nlevs+1))
  
    z[z < zrange[1]] <- zrange[1]
    z[z > zrange[2]] <- zrange[2]
  
    zr <- round( log10(diff(zrange) ), 1)
    
    
    if(is.null(conInt)){
      conInt <- c(diff(zrange)/5,diff(zrange)/2)
    }
  
    for(k in 1:length(conInt)){
      lev <- seq(round(zrange[1],zr),round(zrange[2],zr),by=conInt[k])
      contour(xl,yl,z,add=T,levels=lev,col='brown',lwd=k/2,drawlabels=F) 
   }
  }

    if(BOUND){
      wj <- which(boundaries[,'plot'] == jplots[j])
      ss <- unique(boundaries[wj,'section'])
      for(jj in ss){
        wjk <- which(boundaries[,'plot'] == jplots[j] & boundaries[,'section'] == jj)
        #     polygon(boundaries[wjk,'UTMx'],boundaries[wjk,'UTMy'],border=1,lwd=3)
      }
    }
  
  xl1 <- xr[1] + diff(xr)/2
  lines( c(xl1-100,xl1+100),mapwide/20+c(yl[1],yl[1]),lwd=5,lend=1)
  text(xl1,yl[1],'200 m')
  
  
  values <- list(xy = xygrid[mind,],z = zgrid[mind], FLAG=T)
  
  invisible(values)
}



mapRegions <- function(jplots,zname = 'meters',mapPath,conInt = NULL, 
                       out=F, STANDARDIZE=F, TOPO=T, buffer=200, PLOTHULL=F,
                       reverseColor=F,BOUND=F){

  # plots - e.g., c('DF_EW','DF_EE'), i.e., in the same region
  
  graphics.off()

  reg <- matrix( unlist(strsplit(jplots[1],'_')), ncol=2,byrow=T)[1,1]
  np <- length(jplots)

  nlevs <- 100

    ncol <- 100
    colseq <- terrain.colors(ncol)
    if(reverseColor) colseq  <- rev(colseq)
 
  xfactor <- 1

  if(zname == 'ft'){
      xfactor <- .3048
      zname   <- 'elev'
  }
    

    utmPlot <- utmPlotHull <- utm <- numeric(0)
    for(j in 1:np){
      utm  <- rbind(utm,UTMhull[[jplots[j]]])
    }
    if(length(utm) == 0){
        warning(paste('no UTMs for region'))
        return()
    }
  
  
   xr <- range(utm[,1])
   yr <- range(utm[,2])
   xr <- c(xr[1] - buffer,xr[2] + buffer)
   yr <- c(yr[1] - buffer,yr[2] + buffer)
  
   sc <- max( c(diff(xr),diff(yr)) )/5     #m per inch
   mapwide <- diff(xr)
  
   mapHydro(utm[,1],utm[,2],dataPath='datafiles/',mapPath,reg,buffer=10,
            conInt=conInt,BOUND=T,TOPO=TOPO)


    if(PLOTHULL | BOUND){

     for(j in 1:np){

       if(length(UTMhull[[jplots[j]]]) == 0)next

       xy <- matrix( as.matrix(UTMhull[[jplots[j]]]), ncol=2)

       if(nrow(xy) == 1){
          text(xy[1],xy[2],jplots[j])
          next
       }

       if(PLOTHULL){
      #   polygon(xy[,1],xy[,2],border='white',col='white',lwd=2)
         polygon(xy[,1],xy[,2],border='white',lwd=2)
       }
       
       wt <- apply(xy,2,max) 
       text(wt[1],wt[2],jplots[j],pos=3)
       
      }
    }

    xl1 <- xr[1] + diff(xr)/2

    lines( c(xl1-100,xl1+100),mapwide/20+c(yr[1],yr[1]),lwd=5,lend=1)
    text(xl1,yr[1],'200 m')

    xll <- c( xr[2] - sc/10, xr[2] + sc/10)
    yll <- c( yr[1] , yr[1] + sc/2 )
 #   colorLegend(xll,yll,ytic=round(zrange,1),
 #        scale=seq(zrange[1],zrange[2],length=ncol),cols=colseq,labside='left')

    pname <- zname
    if(zname == hydroPath)pname <- gsub('/','_',pname) 

    if(out) dev.copy2pdf(file=paste(outMain,jplots[1],'_',pname,'Map.pdf',sep=''))
  
}

read.asciigrid.r <- function (fname, as.image = FALSE, plot.image = FALSE) {

    t = file(fname, "r")
    l5 = readLines(t, n = 6)
    l5s = strsplit(l5, "\\s+", perl = T)
    xllcenter = yllcenter = xllcorner = yllcorner = as.numeric(NA)
    for (i in 1:6) {
        fieldname = casefold(l5s[[i]][1])
        if (length(grep("ncols", fieldname))) 
            ncols = as.numeric(l5s[[i]][2])
        if (length(grep("nrows", fieldname))) 
            nrows = as.numeric(l5s[[i]][2])
        if (length(grep("xllcorner", fieldname))) 
            xllcorner = as.numeric(l5s[[i]][2])
        if (length(grep("yllcorner", fieldname))) 
            yllcorner = as.numeric(l5s[[i]][2])
        if (length(grep("xllcenter", fieldname))) 
            xllcenter = as.numeric(l5s[[i]][2])
        if (length(grep("yllcenter", fieldname))) 
            yllcenter = as.numeric(l5s[[i]][2])
        if (length(grep("cellsize", fieldname))) 
            cellsize = as.numeric(l5s[[i]][2])
        if (length(grep("nodata_value", fieldname))) 
            nodata.value = as.numeric(l5s[[i]][2])
    }
    if (is.na(xllcorner) && !is.na(xllcenter)) 
        xllcorner = xllcenter - 0.5 * cellsize
    else xllcenter = xllcorner + 0.5 * cellsize
    if (is.na(yllcorner) && !is.na(yllcenter)) 
        yllcorner = yllcenter - 0.5 * cellsize
    else yllcenter = yllcorner + 0.5 * cellsize
    map = scan(t, as.numeric(0), quiet = TRUE)
    close(t)
    if (length(as.vector(map)) != nrows * ncols) 
        stop("dimensions of map do not match that of header")
    map[map == nodata.value] = NA
    if (as.image) {
        img = matrix(map, ncols, nrows)[, nrows:1]
        img = list(z = img, x = xllcorner + cellsize * ((1:ncols) - 
            0.5), y = yllcorner + cellsize * ((1:nrows) - 0.5))
        if (plot.image) {
            image(img, asp = 1)
            return(invisible(img))
        }
        else return(img)
    }
    df = data.frame(map)
    names(df) = colname

   df
}

postPlotBySpec <- function(bname,unames,ylim=c(-20,20),log=F){

  ww <- which(chains == bname)

  if(unames[1] == 'years'){
    
    rvals <- as.numeric(matrix( unlist(strsplit(rownames(bd)[-pCol],'-')),ncol=2,byrow=T)[,2])
    rvals <- unique(rvals)

    par(mfcol=c(length(rvals),2),bty='n')

    for(r in 1:nreg){

      yc <- paste('-',r,'_',sep='')

      wc <- grep(yc,colnames(chainList[[ww]]))
      if(length(wc) == 0)next

      for(k in 1:nttype){

        wk <- intersect(wc,grep(treeNames[k],colnames(chainList[[ww]])))
        ci <- apply(chainList[[ww]][,wk],2,quantile,c(.5,.025,.975))
        yk <- as.numeric( matrix( unlist(strsplit(colnames(ci),yc)),ncol=2,byrow=T)[ ,1])
        if(k == 1)plot(yk,ci[1,],xlim=range(yrvec),ylim=ylim,lwd=2,type='l',xlab='Year',ylab='Year effect')
        lines(yk,ci[1,],lwd=2,col=k)
        for(kk in 2:3)lines(yk,ci[kk,],lty=2,col=k)
     }
     title(regions[r])
     abline(h=0,lty=2)
   }
   return()
  }

  unames <- unames[unames != 'years']
  cc <- ceiling(length(unames)/2)

  par(mfrow=c(cc,cc),bty='n')

  for(j in 1:length(unames)){

    xvals <- yvals <- matrix(NA,512,nttype)

    wc <- grep(unames[j],colnames(chainList[[ww]]))

    for(k in 1:nttype){

      wk <- intersect(wc,grep(treeNames[k],colnames(chainList[[ww]])))
      tmp <- density( chainList[[ww]][,wk] )
      xvals[,k] <- tmp$x
      yvals[,k] <- tmp$y
      if(log)yvals[,k] <- log10(tmp$y)
    }

    ylim <- c(0,max(yvals,na.rm=T))
    if(log)ylim[1] <- round(ylim[2] - 2,0)

    plot(xvals[,1],yvals[,1],xlim=range(xvals),ylim=ylim,type='l',lwd=2,xlab=unames[j])
    for(k in 1:nttype)lines(xvals[,k],yvals[,k],col=k,lwd=2)
    title(bname)
    if(j == 1)legend('topright',treeNames,text.col=c(1:nttype))
  }
}



getPlotNames <- function(dataPath = 'datafiles/',datFile='plotDataDuke.txt',
                         INCLUDE=F,treeOnly=F,regions=c('CW','DF','MH')){
  
  # INCLUDE - include only plots having include==T in plotDataDuke.txt

  tmp      <- read.table(paste(dataPath,datFile,sep=''),header=T)
  
  r2i <- paste(regions,'_',sep='')
  wr  <- numeric(0)
  for(j in 1:length(regions)){
    wr  <- c(wr,grep(r2i[j],tmp[,'plotname']))
  }
  wr <- sort(wr)
  
  if(INCLUDE){
    ww <- tmp[,'include']
    if(treeOnly)ww <- ww & tmp[,'type'] == 'tree'
    wr <- wr[ww[wr]]
  }

   plotnames <- as.character( tmp[wr,'plotname'] )

   ri   <- matrix( unlist(strsplit(plotnames,"_")),ncol=2,byrow=T)[,1]
   regs <- unique(ri)

   regs <- regs[regs %in% regions]
   plotnames <- plotnames[ri %in% regions]
   ri   <- ri[ri %in% regions]
 
   nr   <- length(regs)

   plotgroups <- vector('list',nr)
   names(plotgroups) <- regs

   for(k in 1:nr){
      plotgroups[[k]] <- plotnames[ri == regs[k]]
   }

   list(plotgroups = plotgroups, plotnames = plotnames, data = tmp[wr,])
}




allometric <- function(int,slope,d){

  10^( int + slope*log10(d) )

}

diam2MaxCan <- function( dataPath='datafiles/' ){

  allomCoeff <- read.table(paste(dataPath,crownAllomFile,sep=''),header=T)

  canMax <- diamMat*0

  for(k in 1:nttype){

    wk <- which( treeData[,'spec'] == k )
    coeff <- allomCoeff['all',]
    ws <- which(rownames(allomCoeff) == treeNames[k])
    if(length(ws) > 0)coeff <- allomCoeff[ws,]
    canMax[wk,] <- allometric(coeff[1,1],coeff[1,2],diamMat[wk,])
  }

  canMax
}

 
baPerHa <- function(diamCm,parea){  #diam in cm, area of plot in ha, returns m2/ha

  pi*(diamCm/2)^2 /parea/10000

}
  

getDens <- function(pname,tname,chain){

   nd <- 512
   nc <- ncol(chain)

   xt <- yt <- matrix(NA,nc,nd)

   for(kk in 1:nc){

      tmp <- density(chain[,kk],n = nd)
      xt[kk,]  <- tmp$x
      yt[kk,]  <- tmp$y
   }
   rownames(xt) <- rownames(yt) <- paste(tname,pname,colnames(chain),sep='-')
   list(x = xt, y = yt)
}



