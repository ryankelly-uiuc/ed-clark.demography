

removeRows <- function(ww,...){  # ... are names of matrices
	xx <- list(...)
	nx <- length(xx)
	out <- numeric(0)
	
	for(k in 1:nx)out <- append(out, list(get(xx[[k]])[-ww,]))
	names(out) <- xx
	out
}

outFile <- function(outfolder=character(0),file){
  paste(outfolder,file,sep='')
}

conditionalMVNVecRcpp <- function(xx, mu, sigma, cdex, p=ncol(mu)){
  
  gdex <- (1:p)[-cdex]-1
  conditionalMVNVec_cpp((cdex-1), gdex, xx, mu, sigma)
  
}


gaussianRatio <- function(q1,q2,mu,variance){ 
  #returns log ratio of 2 gaussians (e.g., metropolis)
  
  -.5*( (q1^2 - q2^2 -2*mu*(-q1 + q2))/variance)
}

binomialRatio <- function(nn,yy,t1,t2){
  
  yy*(log(t1) - log(t2)) - (nn - yy)*(log(1 - t1) - log(1 - t2))
}

poissonRatio <- function(y,lambda1,lambda2){
  
  y*(log(lambda1) - log(lambda2)) - lambda1 + lambda2
  
}
  

comp_outindex <- function(index_in, total) return((1:total)[-index_in]-1)

pmvnormCondRcpp <- function(q=0,xx,mu,smat,whichVar=c(1:nrow(smat))){  
  
  if( !identical(dim(xx),dim(mu)) )stop( "dimensions disagree in pmvnormCondRcpp" )
  
  pj <- mu*0  
  gindexes      <- sapply(whichVar, comp_outindex, total = ncol(xx))
  pj[,whichVar] <- pmvnormCond_cpp(q, whichVar-1 , gindexes, xx, mu, smat)
  pj  
}



bUpdateMVN_Rcpp <- function(xx,yy,bb = NULL,lo = NULL, hi = NULL,sigma,times=1){
  
  require(mvtnorm)
  
  nc <- 1
  if(is.matrix(yy))nc <- ncol(yy)
  if(!is.null(bb))dm <- dimnames(bb)
  if(is.null(bb)) dm <- list(colnames(xx),colnames(yy))
  
  testv <- try(chol(crossprod(xx)),T)
  if(inherits(testv,'try-error')){
    message('X not full rank')
    return( bb ) 
  }
  
  cx   <- chol2inv(testv)
  mu   <- matrix(cx %*% crossprod(xx,yy),1)
  vv   <- kronecker(sigma,cx)
  
  if(is.null(lo)){
    b <- matrix( rmvnorm(1,mu,vv),ncol=nc)
    dimnames(b) <- dm
    return(b)
  }
  
  lo <- as.vector(lo)
  hi <- as.vector(hi)

  testvv <- try(chol(vv),T)
  if(inherits(testvv,'try-error')){
    message('kronecker singular')
    return( bb ) 
  }

  tmp <- tnorm.mvtRcpp(as.vector(bb),mu,vv,lo,hi,times=times)
  bb <- matrix(tmp[times,],ncol=nc)
  dimnames(bb) <- dm
  bb
}


makeXY <- function(xdata,ydata,xnames,ynames,factors=NULL,standard=T,complete=F,
                   minLevels=10,combineFactors=NULL){
  
  #xnames columns in xdata for design matrix; should include 'intercept', unless
  #       all factor levels included; if the latter, each is a full intercept
  #factors are names in x taken to be factors
  #interactions have this form: xname1Xxname2--identified by 'X'
  #ynames columns in ydata for response (univarite will have one value)
  # if 'complete' only complete obs, otherwise NAs are in data;
  #factors cannot be NA
  #minLevels - observations with too few factor levels removed
  
  
  y <- as.matrix( ydata[,ynames] )
  ww <- rowSums(y)
  
  y <- y[ww > 0,]
  xdata <- xdata[ww > 0,]
  
  n <- nrow(xdata)
  x <- numeric(0)
  
  if(!is.null(factors)){
    
    for(j in 1:length(factors)){
      
      xj <- as.character( xdata[,factors[j]] )
      xj[is.na(xj)] <- 'bogus'

      jtab <- table(xj)
      wmin  <- which(jtab < minLevels | names(jtab) == 'bogus')
      wkeep <- which(jtab >= minLevels & names(jtab) != 'bogus')
      
      fnames <- names(jtab)[wkeep]
      other  <- names(jtab)[wmin]
      
      xj[xj %in% other] <- 'other'
      fnames <- c(fnames,'other')
      
      nf     <- length(fnames)
      f      <- matrix(0,n,nf)
      xx <- match(xj,c(fnames))
     
      
      f[ cbind(c(1:n),xx) ] <- 1
      cj <- paste(factors[j],fnames,sep='_')
      colnames(f) <- cj
   #   if('intercept' %in% xnames & j == 1){
        f <- as.matrix(f[,-1],n)
   #     colnames(f) <- cj[-1]
   #   }
      
      x <- cbind(x,f)
    }
    
    otherCols <- grep('other',colnames(x))
    if(!is.null(combineFactors)){
      otherCols <- otherCols[ !otherCols %in% grep(combineFactors,colnames(x)) ]
    }
    otherRows <- unique( which(x[,otherCols] == 1,arr.ind=T)[,1] )
    
    x <- x[-otherRows,]
    x <- x[,-otherCols]
    y <- y[-otherRows,]
  }

  
  vnames <- xnames[xnames %in% colnames(xdata) & !xnames %in% factors]
  p      <- length(vnames)
  
  if(p == 1){
    xmean <-  mean(xdata[,vnames],na.rm=T)
    xsd   <-  sd(xdata[,vnames],na.rm=T)
    names(xmean) <- names(xsd) <- vnames
  }
  if(p > 1){
    xmean <- colMeans(xdata[,vnames],na.rm=T)
    xsd   <- apply(xdata[,vnames],2,sd,na.rm=T)
  }
  
  if(p > 0){
    xcenter <- matrix( xdata[,vnames] - matrix(xmean[vnames],n,p,byrow=T), n,p)
    xstand  <- matrix( xcenter/matrix(xsd[vnames],n,p,byrow=T), n, p)
    colnames(xcenter) <- colnames(xstand) <- vnames
    
    if(!standard)x <- cbind(x,xdata[,vnames])
    if(standard) x <- cbind(x,xstand)
    
    ww <- grep('X',xnames)
    
    if(length(ww) > 0){
      
      for(j in ww){
        inames <- unlist( strsplit(xnames[j],'X') )
        if(!inames[1] %in% factors){
          xi1 <- x[,inames[1]]
          name1 <- inames[1]
        }
        if(inames[1] %in% factors){
          xi1 <- x[,grep(inames[1],colnames(x))]
          name1 <- colnames(x)[grep(inames[1],colnames(x))]
        }
        if(!inames[2] %in% factors){
          xi2 <- x[,inames[2]]
          name2 <- inames[2]
        }
        if(inames[2] %in% factors){
          xi2 <- x[,grep(inames[2],colnames(x))]
          name2 <- colnames(x)[grep(inames[2],colnames(x))]
        }
        xi <- matrix(xi1*xi2,n)
        colnames(xi) <- paste( name1,name2,sep='X')
        x <- cbind(x,xi)
      }
    }
  }
  
  n <- nrow(x)
  
  if('intercept' %in% xnames){
    intercept <- matrix(1,n,1)
    x <- cbind(intercept,x)
    colnames(x)[1] <- 'intercept'
  }
  
  x <- as.matrix(x)
  wmiss <- which(!is.finite(x),arr.ind=T)
  
  xx <- x
  xx[is.na(xx)] <- 0
  rank <- qr(xx)$rank
  
  if(rank < ncol(x))warning('x not full rank')
  
  
  list(x = x, y = y, wmiss = wmiss)
}



makeAR1 <- function(nd,rho,sigma){   # construct AR(1) correlation and covariance matrices of dimension nd
  
  x    <- diag(nd)
  cmat <- rho^abs(row(x) - col(x))
  vmat <- sigma*cmat
  
  list(cormat = cmat, varmat = vmat)
}

invertAR1old <- function(nn,phi,sig){  #dimension, correlation, variance
  
  x1 <- 1/(1 - phi^2)/sig
  x  <- diag( c(x1, rep( (1 + phi^2)/(1 - phi^2)/sig, nn-2), x1) )
  x[ row(x) == col(x) - 1 ] <- x[ row(x) == col(x) + 1 ] <-  -phi/(1 - phi^2)/sig 
  x
}

invertAR1 <- function(nd,phi,sig){  #dimension, correlation, variance
  
  x  <- diag( c(1, rep( (1 + phi^2), nd-2), 1) )
  x[ row(x) == col(x) - 1 ] <- x[ row(x) == col(x) + 1 ] <-  -phi 
  x/sig
}




xy2area <- function(xy){  # convex hull and polygon area
  
  require(splancs)
  
  xy <-  xy[is.finite(xy[,1]) & is.finite(xy[,2]),] 
  hull  <- as.matrix(xy[chull(xy),])
  area  <- areapl(hull)
  
  list(hull = hull, area = area)
}



pasteVector <- function(vec,sep=''){     #paste a vector to create one variable

  x <- vec[1]
  if(length(vec) > 1)for(k in 2:length(vec))x <- paste(x,vec[k],sep=sep)
  x
}



aggregateSequence <- function(index,dat,action='mean',minObs=NULL){    
  # index - column index for aggregation
  # dat is n by t
  
  nc <- length(index)
  if(!is.matrix(dat)){
    dat <- matrix(dat,ncol=nc)
  }
  nn <- nrow(dat)
	
	allDates <- sort(unique(index))
  nm <- length(allDates)
  
	newDat   <- matrix(NA,nn,nm)
	rownames(newDat) <- rownames(dat)
  colnames(newDat) <- allDates
	
	jk <- 0
	for(j in allDates){
		jk <- jk + 1
		wj <- which(index == j)
		mj <- dat[,wj]
    
    if(!is.null(minObs))if(length( which(is.finite(mj)) ) < minObs)next

    if(action == 'mean'){
		  if(is.matrix(mj)) newDat[,jk] <- rowMeans(mj,na.rm=T)
		  if(!is.matrix(mj))newDat[,jk] <- mean(mj,na.rm=T)
    }
    if(action == 'sum'){
		  if(is.matrix(mj)) newDat[,jk] <- rowSums(mj,na.rm=T)
		  if(!is.matrix(mj))newDat[,jk] <- sum(mj,na.rm=T)
    }
	}
	list(dates = allDates, data = newDat)
}

gapFillSequence <- function(y,baseline,fitwidth,x=NULL,tooLo=-Inf,tooHi=Inf){
	
	#y        - sequence with gaps
	#baseline - reference sequence
	#fitwidth - number of preceding values for calibration
	#x        - additional predictors
	
	nf <- length(y)
	if(length(x) > 0 & !is.matrix(x))x <- matrix(x,ncol=1)
	
	xmat <- matrix(NA,nrow(x),ncol=fitwidth)
		
   ki <- c(fitwidth:nf)
   kj <- ki
		
   for(k in 1:fitwidth){
       xmat[ki,k] <- baseline[kj]
		 kj <- kj - 1
	}
	if(length(x) > 0)xmat <- cbind(xmat,x)
	
	wr   <- which(is.finite(rowSums(xmat)) & is.finite(y))
	wr   <- wr[y[wr] > tooLo]              #bad values
	wr   <- wr[y[wr] < tooHi]
	xx   <- xmat[wr,]
	yy   <- y[wr]

	xcol <- c(1:ncol(xx))
	if(sum(xx[,1]) == sum(xx[,2])){
		xcol <- c(1:ncol(xx))[-1]
	}
	b  <- solve(crossprod(xx[,xcol]))%*%crossprod(xx[,xcol],yy)
		
	predy <- xmat[,xcol]%*%b
		
	wf    <- which(!is.finite(y))
	y[wf] <- predy[wf]
	y
}

daysSinceDate <- function(month0,day0,yr0,mo,da,yr,nineteen=F){
	
	#vectors for mo, da, yr since initial date 0

	dd <- clark_mdy.date(mo,da,yr,nineteen)
	dd - clark_mdy.date(month0, day0, yr0,nineteen) + 1
}


clark_mdy.date <- function (month, day, year, 
                            nineteen = TRUE, fillday = FALSE, 
                            fillmonth = FALSE){
                            	
    #days since 1/1/1960

    temp <- any( (month != trunc(month)) | (day != trunc(day)) | 
        (year != trunc(year)))
    if (!is.na(temp) && temp) {
        warning("Non integer input values were truncated in mdy.date")
        month <- trunc(month)
        day   <- trunc(day)
        year  <- trunc(year)
    }
    if (nineteen) 
      year <- ifelse(year < 100, year + 1900, year)
     temp  <- numeric(length(month + day + year))
     month <- month + temp
     day   <- day + temp
     year  <- year + temp
    if (fillmonth) {
        temp <- is.na(month)
        month[temp] <- 7
        day[temp]   <- 1
    }
    if (fillday) 
        day[is.na(day)] <- 15
        month[month < 1 | month > 12] <- NA
        day[day < 1] <- NA
        year[year == 0] <- NA
    year   <- ifelse(year < 0, year + 1, year)
    tyear  <- ifelse(month > 2, year, year - 1)
    tmon   <- ifelse(month > 2, month + 1, month + 13)
    julian <- trunc(365.25 * tyear) + trunc(30.6001 * tmon) + 
        day - 715940
    temp <- trunc(0.01 * tyear)
    save <- ifelse(julian >= -137774, 
                   julian + 2 + trunc(0.25*temp) - temp, 
                   julian)
    year   <- ifelse(month == 12, year + 1, year)
    month  <- ifelse(month == 12, 1, month + 1)
    day    <- 1
    tyear  <- ifelse(month > 2, year, year - 1)
    tmon   <- ifelse(month > 2, month + 1, month + 13)
    julian <- trunc(365.25 * tyear) + trunc(30.6001 * tmon) + day - 715940
    temp   <- trunc(0.01 * tyear)
    save2  <- ifelse(julian >= -137774, julian + 2 + trunc(0.25*temp) - temp, julian)
    temp   <- as.integer(ifelse(save2 > save, save, NA))
    temp
}

yearMonthVecCens <- function(yrvec){   
	
	#returns days to end of each month for years in yrvec
	
	allYears <- sort(unique(yrvec))
	nyr <- length(allYears)
	
	mNames   <- c('jan','feb','mar','apr','may','jun','jul','aug',
	              'sep','oct','nov','dec')
	mDays <- rep(31,12)
	mDays[c(4,6,9,11)] <- 30
	mDays[2] <- 28
	names(mDays) <- mNames
	endMonth <- cumsum(rep(mDays,times=nyr))
	endYr    <- cumsum(rep(sum(mDays),nyr))
	list(endMonth = endMonth,endYr = endYr)
}

bigObjects <- function(nn=10){  #find big objects

  #find large objects
  zz <- sapply(ls(pos = 1), function(x)
                 object.size(get(x, envir = globalenv())))

  as.matrix(rev(sort(zz))[1:nn])
}


logit2Prob <- function(y){   #multivar logit to fractions
	
  if(!is.matrix(y))y <- matrix(y,nrow=1)
  
  zs   <- rowSums(exp(y))
  z1   <- 1/(1 + zs)
  zm   <- exp(y)/ (1 + zs)
  cbind(zm,z1)
  
}

prob2Logit <- function(y){     #fractions to multivar logit      
  
  log(y[,-r]/(1 - rowSums(y[,-r])))
  
}
################################
stateSpaceMVNlogit <- function(){
	#z is the observed proportions of plot/time by species
	
#	krons <- kronecker(diag(1,2),sg)

   omat <- matrix(obsError,n,ns,byrow=T)
   accept <- 0
	
	for(t in 1:nt){
				
		ti <- t + tindex
		propy <- matrix(rnorm(n*ns,yg[ti,],.3),n,ns)
		
		pnow <- rowSums(dnorm(yg[ti,],yobs[ti,],sqrt(obsError),log=T))
		pnew <- rowSums(dnorm(propy,yobs[ti,],sqrt(obsError),log=T))
		
		if(t > 1){
			mu <- xtime[ti-1,]%*%bg + yg[ti-1,]%*%ag
			pnow <- pnow + diag(-(yg[ti,] - mu)%*%sinv%*%t(yg[ti,] - mu)/2)
			pnew <- pnew + diag(-(propy - mu)%*%sinv%*%t(propy - mu)/2)
		}	
		if(t < nt){
			mu1 <- xtime[ti,]%*%bg + yg[ti,]%*%ag
			mu2 <- xtime[ti,]%*%bg + propy%*%ag
			pnow <- pnow + diag(-(yg[ti+1,] - mu1)%*%sinv%*%t(yg[ti+1,] - mu1)/2)
			pnew <- pnew + diag(-(yg[ti+1,] - mu2)%*%sinv%*%t(yg[ti+1,] - mu2)/2)
		}	
		a <- exp(pnew - pnow)
		za <- runif(n,0,1)
		wp <- which(za < a)
		if(length(wp) > 0){
			yg[ti[wp],] <- propy[wp,]
			accept <- accept + length(wp)
		}
   }
   list(yg = yg, accept = accept/nrow(yg))
}

obsErrorMVNlogit <- function(){
	
	u1 <- s1 + .5*n*ns*nt
	u2 <- s2 + .5*sum( (yg - yobs)^2 )
	1/rgamma(1,u1,u2)
}		


vec2Mat <- function(ir,ic,x,nr=max(ir),nc=max(ic)){   #vector x to matrix x[ir,ic]

  xmat <- matrix(NA,nr,nc)
  xmat[cbind(ir,ic)] <- x
  xmat

}

mat2Vec <- function(x,GETINDEX=F,lastTime=NULL){  #matrix x to vector, after lastTime removed

  xx  <- t(x)
  nr <- nrow(xx)
  nc <- ncol(xx)

  rowCol <- index <- numeric(0)

  if(is.null(colnames(xx)))colnames(xx) <- 1:nc
  if(is.null(rownames(xx)))rownames(xx) <- 1:nr

  
  ir <- as.vector( matrix(c(1:nr),nr,nc) )
  ic <- as.vector( matrix(c(1:nc),nr,nc,byrow=T) )
  y <- xx[cbind(ir,ic)]

  if(GETINDEX){
    index  <- expand.grid(col = c(1:nr),row = c(1:nc) )[,c(2,1)]

    if(!is.null(lastTime)){

      wl <- which(is.finite(lastTime))
      if(length(wl) > 0){
         wk <- numeric(0)
         for(k in wl)wk <- c(wk, which(index[,1] == k & index[,2] > lastTime[k]) )

         if(length(wk) > 0){
            y <- y[-wk]
            index <- index[-wk,]
         }
       }
    }
  }

  list(x = y, index = index)

}


interp <- function(y,INCREASING=F,minVal=-Inf,maxVal=Inf,defaultValue=NULL,
                   tinySlope=NULL){  #interpolate vector x

  if(is.null(defaultValue))defaultValue <- NA

  tiny <- .00001
  if(!is.null(tinySlope))tiny <- tinySlope

  y[y < minVal] <- minVal
  y[y > maxVal] <- maxVal

  n  <- length(y)
  wi <- which(is.finite(y))

  if(length(wi) == 0)return(rep(defaultValue,n))
  if(length(wi) == 1)ss <- tiny

  xx  <- c(1:n)
  z  <- y

  if(wi[1] != 1) wi <- c(1,wi)
  if(max(wi) < n)wi <- c(wi,n)

  ss <- diff(z[wi])/diff(xx[wi])

  ss[is.na(ss)] <- 0

  if(length(ss) > 1){
    if(length(ss) > 2)ss[1] <- ss[2]
    ss[length(ss)] <- ss[length(ss)-1]
  }
  if(INCREASING)ss[ss < tiny] <- tiny

  if(is.na(y[1]))  z[1] <- z[wi[2]] - xx[wi[2]]*ss[1]
  if(z[1] < minVal)z[1] <- minVal
  if(z[1] > maxVal)z[1] <- maxVal

  for(k in 2:length(wi)){

     ki <- c(wi[k-1]:wi[k])
     yk <- z[wi[k-1]] + (xx[ki] - xx[wi[k-1]])*ss[k-1]
     yk[yk < minVal] <- minVal
     yk[yk > maxVal] <- maxVal
     z[ki] <- yk
  }
  z
}
  

interpRows <- function(x,startIndex=rep(1,nrow(x)),endIndex=rep(ncol(x),nrow(x)),
                       INCREASING=F,minVal=-Inf,maxVal=Inf,
                       defaultValue=NULL,tinySlope=.001){  
  #interpolate rows of x subject to increasing

  nn  <- nrow(x)
  p  <- ncol(x)
  xx <- c(1:p)

  if(length(minVal) == 1)minVal <- rep(minVal,nn)
  if(length(maxVal) == 1)maxVal <- rep(maxVal,nn)

  ni   <- rep(NA,nn)
  flag <- numeric(0)

  z <- x

  for(i in 1:nn){
    if(startIndex[i] == endIndex[i]){
      z[i,-startIndex[i]] <- NA
      next
    }
    z[i,startIndex[i]:endIndex[i]] <- interp(x[i,startIndex[i]:endIndex[i]],
                                             INCREASING,minVal[i],maxVal[i],
                                             defaultValue,tinySlope)
  }
  
  z
}



byFunction <- function(x,i,j,mat,FUN=sum){  
  
  # for 2-D, ... can be i,j
  # mat is max(i) by max(j)
  
  FUN <- match.fun(FUN)
  
  g <- interaction(i,j)
  split(x, g) <- lapply( split(x,g), FUN, na.rm=T )
  mat[cbind(i,j)] <- x
  mat
  
}


byIndex <- function(xx,INDICES,FUN,coerce=F,...){  
  
#INDICES is list, each same length as  x
  
#  fun <- match.fun(FUN)
  
  nl <- length(INDICES)
  
  tmp  <-  unlist(by( as.vector(xx),INDICES,FUN,...) ) 
  nd   <- dim(tmp)
  tmp  <- array(tmp,dim=nd, dimnames=dimnames(tmp))
  
  tmp[is.na(tmp)] <- 0
  
  if(!coerce)return(tmp)
  
  dname <- dimnames(tmp)
  mk    <- rep(0,length(nd))
  
  for(k in 1:length(nd))mk[k] <- max(as.numeric(dimnames(tmp)[[k]]))
  
  wk <- which(mk > nd)
  if(length(wk) > 0){
    tnew  <- array(0,dim=mk)
    if(length(dim(tnew)) == 1)tnew <- matrix(tnew,dim(tnew),1)
    for(k in wk){
      newk <- c(1:mk[k])
      mat  <- match(dimnames(tmp)[[k]],newk)
      if(k == 1){
        tnew[mat,] <- tmp
        rownames(tnew) <- 1:nrow(tnew)
      }
      if(k == 2){
        tnew[,mat] <- tmp
        colnames(tnew) <- c(1:ncol(tnew))
      }
      tmp <- tnew
    }
  }
  tmp
}
 





sumByIndex <- function(x,pvec,imax){  

  #x - vector to sum, vector of indices,max index
  #    if x is a matrix, then one value of pvec per row of x
  #returns sum by index, inclusive from 1:imax
  
  if(is.matrix(x))x <- rowSums(x,na.rm=T)

     svec <- as.matrix( by(x,pvec,sum,na.rm=T) )

     psum <- rep(0,imax)
     psum[ match(rownames(svec),c(1:imax)) ] <- svec
     names(psum) <- c(1:imax)
     psum
}

sumBy2Index <- function(x,pvec1,pvec2,imax1,imax2){   #

  smat <- matrix(NA,imax1,imax2)

  svec <- by( x, list(i = pvec1, j =  pvec2), sum, na.rm=T) 
  svec <- as( unlist(svec),'matrix' )

  ir <- match(rownames(svec),c(1:imax1))
  ic <- match(colnames(svec),c(1:imax2))
  ii <- as.matrix( expand.grid(ir,ic) )
  smat[ii] <- as.vector(svec)
  smat
}



inData <- function(filename, xnames = NULL, ynames = NULL, tname = NULL, 
                   iname = NULL, na.rm = F, INTERCEPT = F){  
                   	
  #read in data file, return design matrix x, response y
  #xnames, ynames, tname, iname are column headings in filename
  #time indicator t
  #individual indicator i

  data <- read.table(filename,header=T)
  
  if(is.atomic(xnames)){
    x    <- data[,xnames]
    if(!is.matrix(x))x <- as.matrix(x)
    if(INTERCEPT){
      intercept <- rep(1,nrow(data))
  	   x <- cbind(intercept,x)
  	   colnames(x) <- c('intercept',xnames)
    }
  }
  if(is.atomic(ynames)){
    y <- data[,ynames]
    if(!is.matrix(y))y <- as.matrix(y) 
    y  <- matrix(y,nrow(data),length(ynames))
    colnames(y) <- ynames
  }
  
  wf <- c(1:nrow(data))
  
  if(na.rm){
  	 wf <- which(is.finite(rowSums(x)) & is.finite(rowSums(y)))
  	 x  <- x[wf,]
  	 y  <- y[wf,]
  }
  
  z  <- list(x = x, y = y)
  
  if(is.atomic(tname))z$t <- data[wf,tname]
  if(is.atomic(iname))z$i <- data[wf,iname]
  z
}

treebytime <- function(iindex,tindex,x,nr=max(iindex),nc=max(tindex)){ #matrix for x with individuals by time
	
  y <- matrix(NA,nr,nc)
  y[cbind(iindex,tindex)] <- x
  y
}

pRate <- function(par){

  g <- par[1]
  r <- par[2]
  c <- par[3]
  s <- par[4]
  
  f <- g*(1 - exp(-r*x))
  -sum(dnorm(y,c + f,s*f,log=T))
}

byHour <- function(q){  #for capacitance data (Ward et al. 2013)
	
	hrSeq <- sort(unique(b[,'hour']))
	dySeq <- sort(unique(b[,'JD']))
	nh    <- length(hrSeq)
	nd    <- length(dySeq)
	
	qmat <- matrix(NA,nd,nh)
	colnames(qmat) <- hrSeq
	rownames(qmat) <- dySeq
	
	di <- match(b[,'JD'],dySeq)
	hi <- match(b[,'hour'],hrSeq)
	qmat[cbind(di,hi)] <- q
	list(q = qmat, days = dySeq, hours = hrSeq)
	
}


############################### map species

mapSpecies <- function(x,y,z,mapx=range(x),mapy=range(y),
                       scale=0,add=F,sym='circles',colVec=rep(1,length(x)),fill=F){
	
   fillCol <- NA
   if(fill)fillCol <- colVec
   if(scale > 0)mapSetup(mapx,mapy,scale)
   if(sym == 'circles')symbols(x,y,circles=z/10,inches=F,xlim=mapx,ylim=mapy,fg=colVec,bg=fillCol,lwd=2,add=add)
   if(sym == 'squares')symbols(x,y,squares=z/10,inches=F,xlim=mapx,ylim=mapy,fg=colVec,bg=fillCol,lwd=2,add=add)
}

myECDF <- function(x){
  n  <- length(x)
  xx <- sort(unique(x))                 
  fx <- cumsum(tabulate(match(x,xx))/n) 
  list(x = xx, fx = fx)                 
}

samplePlots <- function(mapx,mapy,wide,nplot,mapscale=1,PLOTIT = T){

  yt      <- seq(mapy[1],mapy[2],by=wide)      #y grid locations
  yl      <- length(yt)
  xt      <- seq(mapx[1],mapx[2],by=wide)      #x grid locations
  xl      <- length(xt)
  mapgrid <- cbind(rep(xt,each=yl),rep(yt,xl))       #x and y locations
  loc     <- mapgrid[sample(yl*xl,nplot,replace=F),] #samples from grid

  sl <- .5*wide                             #plot edges
  xbound <- cbind((loc[,1]-sl),(loc[,1]+sl))
  ybound <- cbind((loc[,2]-sl),(loc[,2]+sl))
  xindex <- c(1,2,2,1,1)
  yindex <- c(1,1,2,2,1)
  
  specname  <- sort(unique(treedata[,'species']))
  nspec     <- length(specname)

 # plot.data  <- numeric(0)                   #list of observations
  tableSpec <- matrix(0,nspec,nplot)
  rownames(tableSpec) <- specname

  if(PLOTIT)plot(-1000,0,xlim=mapx,ylim=mapy)

  for(i in 1:nplot){

  # extract trees on sample plot i and store them in table.spec
    xt <- !is.na(cut(treedata[,'x'],breaks=xbound[i,],exclude=NA))
    yt <- !is.na(cut(treedata[,'y'],breaks=ybound[i,],exclude=NA))
    
    tmp <- treedata[xt & yt,]
    if(nrow(tmp) > 0){
      tableSpec[,i] <- table(tmp[,'species'])
      if(PLOTIT)symbols(tmp[,'x'],tmp[,'y'],circles=tmp[,'dbh']/20,inches=F,add=T)
    }

  # draw a box around each plot
    if(PLOTIT){
      xvec <- xbound[i,xindex]
      yvec <- ybound[i,yindex]
      lines(xvec,yvec)
    }
  }

  tableSpec
}


appendData <- function(oldfile,newfile,oldDates,newDates){  #append new plot data
    	
    rold <- rownames(oldfile)
    rnew <- rownames(newfile)

    	if(!is.matrix(newfile)){
    		newfile <- matrix(newfile,length(newfile),1)
    		colnames(newfile) <- newDates
    	}

    	if(length(oldfile) == 0)return(newfile)
    	
    	if(!is.matrix(oldfile)){
    		oldfile <- matrix(oldfile,length(oldfile),1)
    		colnames(oldfile) <- oldDates
    	}
    	
      allDates <- unique(c(oldDates,newDates))
    	
      newMat  <- matrix(NA,nrow(newfile),length(allDates))
      colnames(newMat) <- allDates

      wwc     <- match(colnames(newfile), allDates)
      newMat[,wwc] <- newfile
      newMat  <- matrix(newMat,nrow(newfile),length(allDates))
      colnames(newMat) <- allDates
      rownames(newMat) <- rnew

      oldMat <- matrix(NA,nrow(oldfile),length(allDates))
      colnames(oldMat) <- allDates
      wwc    <- match(colnames(oldfile), allDates)
      oldMat[,wwc] <- oldfile
      oldMat <- matrix(oldMat,nrow(oldfile),length(allDates))
      colnames(oldMat) <- allDates
      rownames(oldMat) <- rold
      
      rbind(oldMat,newMat)
}


missingCol <- function(x,colName,action='warn',value=NA){   

# check for missing columns in vector colName
# action = 'warn','stop','add'

  nx <- nrow(x)
  
  newCols <- character(0)
  
  checkx <- which( duplicated(colnames(x)) )
  checkc <- which( duplicated(colName ) )
  
  if(length(checkc) > 0){
    warning('duplicated columns in species list')
  }

  for(j in 1:length(colName)){

    if(colName[j] %in% colnames(x))next

    if(action == 'warn')warning( paste('missing',colName[j]) )
    if(action == 'stop')   stop( paste('missing',colName[j]) )

    if(action == 'add'){

      if(length(value) == 1) new <- rep(value,nx)
      if(length(value) == nx)new <- value
      if(!length(value) %in% c(1,nx))stop('value must be 1 or nrow(x)')
     
      x   <- cbind(x,new)
      colnames(x)[ncol(x)] <- colName[j]
      newCols <- c(newCols,colName[j])
    }
  }
#  x <- x[,colName]
  invisible(x)
}

hist2Add <- function(xx,xlim=NULL,title=''){
  
  # add histogram to existing plot using add.scatter, library ade4
  
  if(is.null(xlim))xlim <- quantile(xx,c(.001,.999),na.rm=T)
  
  xx[xx > xlim[2]] <- xlim[2]
  opar <- par('mar','yaxt','plt')
  on.exit(par(opar))
  par(mar = rep(.1,4),yaxt='n',plt=par('plt'))
  hist(xx,xlab='',ylab='',main=title,col='brown',proba=T,nclass=60,xlim=xlim)
}

KLdivergence <- function(probs,q,NORM=T){  # discrete K-L divergence
  
  # probs, q are normalized
  # if matrices, rows are distributions
  # if !NORM, then must be normalized
  
  if(!is.matrix(probs)){
    probs <- matrix(probs,1)
    q     <- matrix(q,1)
  }
  if(!NORM){
    p1 <- rowSums(probs)
    q1 <- rowSums(q)
    w0 <- which(p1 == 0)
    probs <- probs/matrix(p1,nrow(probs),ncol(probs))
    q     <- q/matrix(q1,nrow(probs),ncol(probs))
  }
  
  rat <- probs/q
  prr <- log(rat)
  prr[rat == 0 | !is.finite(prr)] <- 0
  
  list(KL = rowSums(prr*probs), zero = w0)
}


#######################################
appendMatrix <- function(m1,m2,fill=NA,SORT=F){  # matches matrices by column names

   if(length(m1) == 0){
     if(!is.matrix(m2))m2 <- matrix(m2,nrow=1)
     return(m2)
   }
   if(length(m2) == 0){
     if(!is.matrix(m1))m1 <- matrix(m1,nrow=1)
     return(m1)
   }

   c1 <- colnames(m1)
   c2 <- colnames(m2)
   r1 <- rownames(m1)
   r2 <- rownames(m2)
   n1 <- nrow(m1)
   n2 <- nrow(m2)

   allc <-  unique( c(c1,c2) ) 
   if(SORT)allc <- sort(allc)

   nr <- n1 + n2
   nc <- length(allc)

   if(is.null(r1))r1 <- paste('r',c(1:n1),sep='-')
   if(is.null(r2))r2 <- paste('r',c((n1+1):nr),sep='-')
   new <- c(r1,r2)

   mat1 <- match(c1,allc)
   mat2 <- match(c2,allc)

   out <- matrix(fill,nr,nc)
   colnames(out) <- allc
   rownames(out) <- new

   out[1:n1,mat1] <- m1
   out[(n1+1):nr,mat2] <- m2
   out
}
   


############################################
row2Mat <- function(vec){

  if(is.matrix(vec))return(vec)
  vn  <- names(vec)
  vec <- matrix(vec,1)
  colnames(vec) <- vn
  vec
}
############################################
col2Mat <- function(vec,namecol=NULL){

  if(is.matrix(vec))return(vec)
  vn  <- names(vec)
  vec <- matrix(vec,ncol=1)
  rownames(vec) <- vn
  colnames(vec) <- namecol
  vec
}
##############################################
rowBind <- function(matnow,row2add,rowName){

  if(length(matnow) == 0){
    matnow <- row2add
    if(!is.matrix(row2add)){
        matnow <- matrix(matnow,nrow=1)
        if(!is.null(names(row2add)))colnames(matnow) <- names(row2add)
    }
    rownames(matnow) <- rowName
    return(matnow)
  }

  if(!is.matrix(row2add))row2add <- t(as.matrix(row2add))

  matnow <- rbind(matnow,row2add)
  lastrows <- c( (nrow(matnow) - nrow(row2add) + 1):nrow(matnow) )
  rownames(matnow)[lastrows] <- rowName
  matnow
}


variancePrior <- function(mu,wt,maxFactor=100){

  s1 <- wt
  s2 <- mu*(s1 - 1)
  lo <- 1e-10
  hi <- maxFactor*mu

  list(mu = mu, s1 = s1, s2 = s2, lo = lo, hi = hi)
}

#############################################

values2contour <- function(xx,yy,z,nx=100,ny=100,lty=1,labcex=.7,
                           col='black',lwd=1,zlevs=NULL,add=F,fill=F,
                           drawlabels=F){    

  # contours where x,y is not a uniform grid, requires 'spatial' library

  require(mapdata)
 # require(MBA)

  xyzmat <- cbind(xx,yy,z)
  
  colnames(xyzmat) <- c('x','y','z')
  
 # surf  <- mba.surf(xyz=xyzmat,no.X=nx,no.Y=ny,h=7,sp=F,extend=F)$xyz.est
  
  surf  <- myMBA.surf(xyz=xyzmat,no.X=nx,no.Y=ny,h=7,sp=F,extend=F)$xyz.est
  
  if(is.null(zlevs)){
    zlevs <- signif(seq(min(z), max(z), length=3),1)
  }
  contour(surf, levels=zlevs,lwd=lwd,lty=lty,col=col,add=T,labcex=labcex,
          drawlabels=drawlabels)
  
  if(fill){
    zl <- zlevs
    if(length(zl) == 1)stop('fill.contour() requires at least 2 contour lines')
    .filled.contour(surf$x,surf$y,surf$z,levels=zl,col=col)
  }
  invisible(surf)
}


nearestNeighbors <- function(xy1,xy2,k=1){ #find k nearest neighbots
  
  require(RANN)
  
  nn2(xy1,xy2,k)$nn.idx
}

#####################################################

crosscorByRow <- function(xmat,ymat,lag=ncol(xmat),BOOTSTRAP=F,PLOT=F){  
  #cross correlation for each row of xmat[i,] vs ymat[i,]

  xmat <- row2Mat(xmat)
  ymat <- row2Mat(ymat)

  nn <- nrow(xmat)
  xx <- c(-lag:lag)
  nc <- length(xx)
  yy <- matrix(NA,nn,nc)
  ii <- numeric(0)

  ciMean <- numeric(0)

  for(i in 1:nn){
    di <- xmat[i,]
    fi <- ymat[i,]
    wi <- which(is.finite(fi) & is.finite(di) & fi > 0)
    if(length(wi) < lag/2)next
    if(var(fi[wi]) == 0)next

  #  xxx <- xx[wi]
    xxx <- wi - mean(wi)
    ddd <- di[wi]
    fff <- fi[wi]
    di  <- lm(ddd ~ xxx)$residuals
    fi  <- lm(fff ~ xxx)$residuals

    cross <- ccf(di,fi,type='correlation',plot=F)
    cx <- cross$lag
    cy <- cross$acf

    cy <- cy[cx %in% xx]
    cx <- cx[cx %in% xx]
    yy[i,match(cx,xx)] <- cy
    ii <- c(ii,i)
  }

  nk <- length(ii)   #sample size for good series
  ci <- matrix(NA,3,ncol(yy))
  rownames(ci) <- c('50%','2.5%','97.5%')

  if(nk == 1)ci[1,] <- yy[ii,]
  if(nk > 1 & nk < 10)ci[1,] <- apply(yy,2,mean,na.rm=T)
  if(nk > 10){
    ci <- apply(yy,2,quantile,c(.5,.025,.975),na.rm=T)
 #   ci[1,] <- apply(yy,2,mean,na.rm=T)
    colnames(ci) <- xx
  }
  ciMean <- ci

  if(BOOTSTRAP & nk > 10){

    nboot <- 2000
    mu <- matrix(NA,nboot,nc)
    for(g in 1:nboot){
      isamp <- sample(ii,nk,replace=T)
      mu[g,] <- apply(yy[isamp,],2,mean,na.rm=T)
    }

    ciMean <- apply(mu,2,quantile,c(.5,.025,.975),na.rm=T)
    colnames(ciMean) <- xx
  }


  if(PLOT){
    par(bty='n')
    plot(xx,ci[1,],type='l',lwd=2,ylim=c(-.6,.6),ylab='Correlation',xlab='Lag',col=2)
    abline(h=0,lwd=2,col='grey')
    abline(v=0,lwd=2,col='grey')

    for(j in 1:3)lines(xx,ci[j,],lty=2)
    for(j in 1:2)lines(xx,ciMean[j,],lty=2,col=2,lwd=2)

    text(xx[1],.5,paste('n = ',nk),pos=4)
  }

  list(lag = xx, ci = ci,  ciMean = ciMean, n = nk)
}

PDForPS <- function(PDF,fname){

  if(PDF)dev.copy2pdf(file=paste(fname,'.pdf',sep=''))
  if(!PDF)dev.print(device=postscript,file=paste(fname,'.ps',sep=''),width=6,horizontal=F)
}

####################################################

points2grid <- function(xx,yy,grid){
  
  # grid is 2 cols (xgrid,ygrid)
  
  require(RANN)
  
  tmp <- nn2(grid,cbind(xx,yy),k=1)[[1]]
  out <- grid[tmp,]
  tt  <- unique(tmp)
  fraction <- length(tt)/nrow(grid)
  
  list(outValues = out, gridIndex = tmp, gridFraction = fraction)
}
  

points2contour <- function(x,y,q=NULL,xlabel=NULL,ylabel=NULL,main=NULL,
                           xp=NULL,yp=NULL,levs=NULL,add=F,maxz=Inf,col='black',fill=F){  

  #creates contours for (x,y) at density q
  
  if(is.null(q))q <- 10
  if(is.null(xp))xp <- range(x,na.rm=T)
  if(is.null(yp))yp <- range(y,na.rm=T)

  xr    <- range(x,na.rm=T)
  yr    <- range(y,na.rm=T)
  xd <- (xr[2] - xr[1])/q
  yd <- (yr[2] - yr[1])/q

  xgrid <- seq(xr[1],xr[2],by=xd)
  ygrid <- seq(yr[1],yr[2],by=yd)

  xf <- cut(x,xgrid)
  yf <- cut(y,ygrid)

  z <- table(xf,yf)
  z[z > maxz] <- maxz

  xmids <- (xgrid - xd/2)[-1]
  ymids <- (ygrid - yd/2)[-1]
  
  d     <- max(c(diff(xr),diff(yr)),na.rm=T)
  if(is.null(q))q <- d/20
  
  
  if(is.null(levs))levs <- signif(seq(min(z),max(z),length.out=4),1)
  
  lwdd <- 2
  if(length(levs) > 1)lwdd <- seq(2,length(levs),by=1)
  
  

  if(!add)image(xmids,ymids,z,xlab=xlabel,ylab=ylabel,xlim=c(xp[1],xp[2]),ylim=c(yp[1],yp[2]))
  contour(xmids,ymids,z,add=T,levels=levs,lwd=lwdd,col=col)
  if(fill).filled.contour(xmids,ymids,z,levels=levs,col=col)
  if(!is.null(main))title(main)
  
  invisible( list(x = xmids, y = ymids, z = z) )
}

####################################################

histf <- function(vec,minv,maxv)hist(vec,breaks=seq(minv,maxv,by=.02))


####################################################

myrmultinom <- function(size,p){  

  # n multinomial r.v. for a n by ncol(p) matrix of probs
  # each row of p is a probability vector
  # size is one integer or a length-n vector of integers

  p <- row2Mat(p)

  n     <- nrow(p)
  J     <- ncol(p)

  if(length(size) == 1)size <- rep(size,n)

  jord  <- sample(J,J)    #randomize order

  p <- row2Mat(p[,jord])

  y <- yy  <- matrix(0,n,J)
  sizej <- size
  sumj  <- rep(0,n)
  dpj   <- rep(1,n)
  pj    <- p
  wj    <- c(1:n)

  for(j in 1:(J-1)){
    a     <- round(pj[wj,1],10)
    y[wj,j] <- rbinom(length(wj),sizej[wj],a)
    sumj  <- sumj + y[,j]
    sizej <- size - sumj
    dpj   <- dpj - p[,j]
    pj    <- matrix(p[,c((j+1):J)]/dpj,nrow(p))
    wj    <- which(sumj < size,arr.ind=T) 
  }

  if(n == 1)y[,J] <- size - sum(y)
  if(n > 1) y[,J] <- size - rowSums(y)

  yy[,jord] <- y
  yy

}

####################################################

truncpars <- function(x,lo,hi){        

  #JS Clark
  #fit truncated multivariate normal to posterior x, known lo and hi

  if(is.vector(x))x <- matrix(x,length(x),1)  
  px <- ncol(x)          #dimension of original x

  ww <- which(!is.finite(x),arr.ind=T)
  if(length(ww) > 0){
    ww <- unique(ww[,2])
    wk <- which(!c(1:ncol(x)) %in% ww,arr.ind=T)
    x <- x[,wk]
  }
  muvec <- apply(x,2,mean,na.rm=T)
  sig   <- cov(x,use="complete.obs")
  pk <- ncol(x)
  nn <- nrow(x)

  ngg   <- 1000
  nkeep <- ngg - 300
  nk    <- 0
  mug   <- rep(0,pk)
  cvg   <- rep(0,pk^2)

  for(g in 1:ngg){
    tmp   <- truncmvtnorm(x,muvec,sig,lo,hi)
    muvec <- tmp$mu
    sig   <- tmp$sig

    if(g > nkeep){
      print(muvec)
      nk       <- nk + 1
      mug      <- mug + muvec
      cvg      <- cvg + as.vector(sig)
    }
  }
  mvec <- mug/nk
  covmat <- matrix(cvg/nk,pk,pk)

  if(length(ww) > 0){
    mnew     <- rep(NA,px)
    mnew[wk] <- mvec
    covnew   <- matrix(NA,px,px)
    covnew[wk,wk] <- covmat
    mvec   <- mnew
    covmat <- covnew
  }

  list(mu = mvec, cm = covmat)
}

####################################################

trunclogis <- function(n,lo,hi,bpars,xvars,wp){ 

  #truncated logistic
  # bpars - parameter vector for logit
  # xvars - variables corresponding to bpars, e.g., c(1,x1,x2)
  # wp    - which variable and parameter in the logit model 
  #         (not 1, because bpars[1] is intercept)
  #         (xvars[wp] is not used)

  xlo     <- xvars
  xlo[wp] <- lo
  xhi     <- xvars
  xhi[wp] <- hi
  sl      <- sum(xlo*bpars)
  sh      <- sum(xhi*bpars)
  lflo    <- exp(sl)/(1 + exp(sl))
  lfhi    <- exp(sh)/(1 + exp(sh))

  z <- runif(n,lflo,lfhi)
  (log(z/(1 - z)) - sum(xvars[-wp]*bpars[-wp]))/bpars[wp]
}

####################################################

truncmvtnorm <- function(x,muvec,sig,lo,hi){
	
	require(mvtnorm)

  #sample from a truncated normal

  if(is.vector(x))x <- matrix(x,length(x),1)
  if(length(sig) == 1)sig <- matrix(sig,1,1)
  y  <- x*0
  pk <- ncol(x)
  n  <- nrow(x)

  for(j in 1:pk){

    if(j == 1){
      t <- muvec[1]
      w <- sqrt(sig[1,1])
    }
    if(j > 1){

       svec <- c(1:(j-1))
       vmat <- sig[j,svec] %*% solve(sig[svec,svec])
       ymu  <- y[,svec] - matrix(rep(muvec[svec],n),n,(j-1),byrow=T)
       t    <- t(muvec[j] +  vmat %*% t(ymu))
       w    <- as.numeric(sqrt( sig[j,j] - vmat %*% c(sig[svec,j]) ))
    }

    up <- pnorm(x[,j],t,w) - pnorm(lo[j],t,w)
    do <- pnorm(hi[j],t,w) - pnorm(lo[j],t,w)

    add <- w*qnorm(up/do)
    add[!is.finite(add)] <- 0
    y[,j] <- t + add

  }

  muy   <- colMeans(y)
  muvec <- myrmvnorm(1,muy,sig/n)

 #the covariance matrix

  mumat <- matrix(muvec,n,pk,byrow=T)
  sy    <- crossprod(y - mumat) 
   ss   <- solve(sy)
   df   <- n 
  alpha <- myrmvnorm(df,rep(0,pk),ss)
  th    <- solve(crossprod(alpha))

  list(mu = muvec, sig = th)

} 

####################################################

logit <- function(x){log(x/(1-x))}  #logit

####################################################

invlogit <- function(x, log = FALSE){  #inverse logit

  if(log)return(-log(1 + exp(-x)))
  1/(1 + exp(-x))

}

##############################
invMat <- function(SS,NEARPD=F){  #matrix inversion, if NEARPD find closest PD matrix

    require(Matrix)
    
    testv <- try(chol(SS),T)

    if( inherits(testv,'try-error') ){
       message('near pos definite used in invMat')
       if(NEARPD){
         require(Matrix)
         SS     <- as.matrix( nearPD(SS)$mat )
         testv <- try(chol(SS),T)
       }
    }

    chol2inv(testv)
}

mydmvnorm <- function(xx,mu,sigma,sinv=NULL,log=FALSE){
  
  xx <- xx - mu
  if(!is.matrix(xx))xx <- matrix(xx,1)
  
  if(!is.null(sinv)){
    distval <- diag( xx%*%sinv%*%t(xx) )
    ev      <- eigen(sinv, only.values = TRUE)$values
    logdet  <- -sum(log(ev))
  }
  
  if(is.null(sinv)){
    testv <- try(chol(sigma),T)
    if(inherits(testv,'try-error')){
      tiny <- min(abs(xx))/100 + 1e-5
      sigma <- sigma + diag(diag(sigma + tiny))
      testv <- try(chol(sigma),T)
    }
    covm <- chol2inv(testv)
    distval <- rowSums((xx %*% covm) * xx)
    ev <- eigen(sigma, only.values = TRUE)$values 
    if(min(ev) < 0)ev <- nearPD(sigma)$eigenvalues
    logdet   <- sum(log( ev ))
  }
  
  logret <- -(ncol(xx) * log(2 * pi) + logdet + distval)/2
  if(log)return(logret)
  exp(logret)
}
  
################################################
mydmvnormOld <- function(xx,mu,sigma,log=FALSE){

  #mv normal density

    if (is.vector(xx))xx <- matrix(xx, ncol = length(xx))
    if (is.vector(mu))mu <- matrix(mu, ncol = length(xx))

 #   zz   <- sweep(xx, 2, mu)
    zz <- xx - mu
    
    ss <- colSums(zz^2)%*%sigma

    testv <- try(chol(sigma),T)
    if(inherits(testv,'try-error')){
       message('error in mydmvnorm')
       return(mu)
    }

    cov <- chol2inv(testv)
    distval <- rowSums((zz %*% cov) * zz)
    names(distval) <- rownames(zz)

  #  distval <- mahalanobis(zz, mu, sigma)
    logdet   <- sum(log( eigen(sigma, symmetric = TRUE, only.values = TRUE)$values ))
    if(is.na(logdet))return(logdet)
    logretval <- -(ncol(zz) * log(2 * pi) + logdet + distval)/2
    if(log)return(logretval)
    exp(logretval)
}

  


bUpdateMVN <- function(xx,yy,bb = NULL,lo = NULL, hi = NULL,sigma){

  #  yy ~ MVN(xx%*%bb,sigma)
  #   xx - design matrix (n X p)
  #   yy - response matrix (n X q)
  #   bb - parameter matrix to update (p X q)
  #   lo - lower truncation (p X q)
  #   hi - upper truncation (p X q)
  #   if truncation (lo and/or hi provided) bb must also be provided

  require(mvtnorm)

  testv <- try(chol(crossprod(xx)),T)
  if(inherits(testv,'try-error')){
    message('X not full rank')
    return( bb ) 
  }
 
  cx   <- chol2inv(testv)
  mu   <- matrix(cx %*% crossprod(xx,yy),1)
  vv   <- kronecker(sigma,cx)
  
  if(is.null(lo)){
    b <- matrix( rmvnorm(1,mu,vv),ncol=2)
    dimnames(b) <- dimnames(bb)
    return(b)
  }

  lo <- as.vector(lo)
  hi <- as.vector(hi)

  b <- matrix( tnorm.mvt(as.vector(bb),mu,vv,lo,hi,times=3), ncol=2)
  dimnames(b) <- dimnames(bb)
  b
}


##############################################
rbvnormFromCols <- function(muMat,sigMat,lo=NULL,hi=NULL){  
	
	#one random vector per row for bivariate normal, based on cholesky
	#sigMat has 4 cols (s11,s12,s12,s22)
	
	n   <- nrow(muMat)
	rho <- sigMat[,2]/sqrt(sigMat[,1]*sigMat[,4])
	
        if(is.null(lo)){
	  z1 <- rnorm(n)
	  z2 <- rnorm(n)
        }
        if(!is.null(lo)){
          l1 <- (lo - muMat[,1])/sqrt(sigMat[,1])
          h1 <- (hi - muMat[,1])/sqrt(sigMat[,1])
          z1 <- tnorm(n,l1,h1,0,1)

          l2 <- ((lo - muMat[,2])/sqrt(sigMat[,4]) - rho*z1)/sqrt(1 - rho^2)
          h2 <- ((hi - muMat[,2])/sqrt(sigMat[,4]) - rho*z1)/sqrt(1 - rho^2)
          z2 <- tnorm(n,l2,h2,0,1)
        }
	
	cbind(muMat[,1] + sqrt(sigMat[,1])*z1,
	      muMat[,2] + sqrt(sigMat[,4])*(rho*z1 + sqrt(1 - rho^2)*z2) )
}

dbvnormFromCols <- function(y,muMat,sigMat){	
  
	sinv <- invertcol2(sigMat)
	z    <- y - muMat
	
	q <- y*0
	q[,1] <- z[,1]*sinv[,1] + z[,2]*sinv[,2]
	q[,2] <- z[,1]*sinv[,3] + z[,2]*sinv[,4]
	q     <- q[,1]*z[,1] + q[,2]*z[,2]
	
	logdet <- log(sigMat[,1]*sigMat[,4] - 2*sigMat[,2])
	
	-(log(2 * pi) + logdet + q)/2
	
}
	


pmvnormApprox <- function(q,mu=rep(0,nrow(sigma)),sigma ,times=100){
  
  #'joint' - pr that each vector < q and
  #'marginal' - pr each each element marginally, i.e., pnorm(q,mu,sqrt(diag(sigma)) )

if(!is.matrix(mu))mu <- matrix(mu,1)
if(length(q) == 1)q <- rep(q,nrow(sigma))

nr <- nrow(mu)
nc <- nrow(sigma)

avec <-  mu
z    <- mu*0

joint <- rep(0,nr)

for(i in 1:times){
  
  zi <- mu*0
  
  for(k in 1:nc){
    
    #   r <- tnorm.mvtRcpp(mu, mu, sigma, lo= -Inf, hi=Inf, whichSample=k, times=1) 
    
    tmp <- conditionalMVNVecRcpp(avec,mu,sigma,k)
    muk <- tmp$mu
    sgk <- tmp$vr
    
    if(length(muk) == 0)next
    
    avec[,k] <- rnorm(nr,muk,sqrt(sgk))
    wj  <- which(avec[,k] < q[k])
    zi[wj,k] <- zi[wj,k] + 1
    z[wj,k]  <- z[wj,k] + 1
  }
  joint[rowSums(z) == nc] <- joint[rowSums(z) == nc] + 1
}

 list(joint = joint/times, marginal = z/times)
}


####################################################

myrmvnorm <- function (nn, mu, sigma){

    sigsvd <- svd(sigma)
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    
    retval <- matrix(rnorm(nn * ncol(sigma)), nn) %*% retval
    retval + mu
}


####################################################

mypmvnorm <- function (lower, upper, mu,sigma){

   corr <- cov2cor(sigma)
   lower <- (lower - mu)/sqrt(diag(sigma))
   upper <- (upper - mu)/sqrt(diag(sigma))
   mean <- rep(0, length(lower))
   RET  <- mvt(lower = lower, upper = upper, df = 0,
                corr = corr, delta = mu, maxpts = 25000, abseps = 0.001,
                releps = 0)
   return(RET$value)
}
####################################################

myqmvnorm <- function (p, interval = c(-10, 10), tail = c("lower.tail", "upper.tail",
    "both.tails"), mean = 0, corr = NULL, sigma = NULL, maxpts = 25000,
    abseps = 0.001, releps = 0, ...)
{
    if (length(p) != 1 || (p <= 0 || p >= 1))
        stop(sQuote("p"), " is not a double between zero and one")
    tail <- match.arg(tail)
    dim <- length(mean)
    if (is.matrix(corr))
        dim <- nrow(corr)
    if (is.matrix(sigma))
        dim <- nrow(sigma)
    lower <- rep(0, dim)
    upper <- rep(0, dim)
    args <- checkmvArgs(lower, upper, mean, corr, sigma)
    dim <- length(args$mean)
    pfct <- function(q) {
        switch(tail, both.tails = {
            low <- rep(-abs(q), dim)
            upp <- rep(abs(q), dim)
        }, upper.tail = {
            low <- rep(q, dim)
            upp <- rep(Inf, dim)
        }, lower.tail = {
            low <- rep(-Inf, dim)
            upp <- rep(q, dim)
        }, )
        pmvnorm(lower = low, upper = upp, mean = args$mean, corr = args$corr,
            sigma = args$sigma, abseps = abseps, maxpts = maxpts,
            releps = releps) - p
    }
    if (tail == "both.tails") {
        interval[1] <- 0
        interval <- abs(interval)
    }
    qroot <- uniroot(pfct, interval = interval, ...)
    names(qroot)[1:2] <- c("quantile", "f.quantile")
    qroot
}


####################################################

fit.tnorm <- function(parvec){  

  #fit parameters for truncated normal
  #vector parvec includes mean, diagonal, & offdiagonals for cov matrix

  mu <- parvec[mi]
  sigmat <- matrix(0,length(mu),length(mu))
  sigmat <- diag(parvec[p2])
  sigmat[si] <- parvec[p3]
  sigmat[cbind(si[,2],si[,1])] <- parvec[p3]

  c1     <- mydmvnorm(x,mu,sigmat,log=T)
  c2     <- log(mypmvnorm(lo,hi,mu,sigmat))
  -sum(c1 - c2)

}

patch2patch <- function(from,to,q,dx,sigma){

   n2   <- length(to)
   n1   <- length(from)
   r    <- matrix(q[to],n2,n1)*exp(-(dx/sigma)^2)

   if(n1 > 1) p <- r/matrix(colSums(r),n2,n1,byrow=T)
   if(n1 == 1)p <- r/sum(r)
   p
}

move <- function(from,to,q,dx,sigma){
   pr  <- patch2patch(from,to,q,dx,sigma) #prob moving to other patches
   wk  <- myrmultinom(1,t(pr))            #movement vector
   which(wk == 1,arr.ind=T)[,2]           #extract new patch location
} 


####################################################

smooth.na <- function(x,y){   

  #remove missing values
  #x is the index
  #y is a matrix with rows indexed by x
  
  if(!is.matrix(y))y <- matrix(y,ncol=1)

    wy <- which(!is.finite(y),arr.ind =T)
    if(length(wy) == 0)return(cbind(x,y))
    wy <- unique(wy[,1])
      ynew <- y[-wy,]
      xnew <- x[-wy]

    return(cbind(xnew,ynew))
}
####################################################

smooth.ma <- function(y,wt){   

  #moving average filter with weights (w0,w1,...), assumed symmetric

  if(length(wt) > length(y))wt <- wt[1:length(y)]
  nw <- length(wt)
  ny <- length(y)
  w <- c(rev(wt),wt[2:nw])
  ymat <- matrix(NA,ny,length(w))
  kb <- nw
  ke <- ny
  ky <- ny - kb + 1

  kb <- c(nw:-(nw-2)); kb[kb < 1] <- 1
  ke <- c((ny+nw-1): (ny-nw+1)); ke[ke > ny] <- ny
  yb <- rev(kb)
  ye <- rev(ke)

  for(kj in 1:(2*nw-1)){
    ymat[kb[kj]:ke[kj],kj] <- y[yb[kj]:ye[kj]]
  }

  wmat <- matrix(rep(w,ny),ny,length(w),byrow=T)
  wmat[which(is.na(ymat))] <- NA
  wmat <- wmat/rowSums(wmat,na.rm=T)
  newy <- rowSums(ymat*wmat,na.rm=T)
  newy

}



####################################################

distmat <- function(x1,y1,x2,y2){
    xd <- outer(x1,x2,function(x1,x2) (x1 - x2)^2)
    yd <- outer(y1,y2,function(y1,y2) (y1 - y2)^2)
    t(sqrt(xd + yd)) 
}
####################################################

tnorm <- function(n,lo,hi,mu,sig){   

  #normal truncated lo and hi

  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))

  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 

  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == Inf]  <- lo[z == Inf]
  z[z == -Inf] <- hi[z == -Inf]
  z
}

nextTimeMat <- function(mat,last=rep(ncol(mat),nrow(mat)),INC=F){

  # mat is n by time matrix, returns same shifted to left

  til <- cbind(c(1:nrow(mat)),last)

  x <- cbind(mat[,-1],mat[,ncol(mat)])
  x[ til ] <- mat[ til ]

  if(INC){
    inc <- x[ cbind(c(1:nrow(mat)),last-1) ] - x[ cbind(c(1:nrow(mat)),last-2) ]
    x[ til ] <- x[ til ] + inc
  }

  x
}

lastTimeMat <- function(mat,first=rep(1,nrow(mat)),INC=F,minInc=.001){  

  # mat is n by time matrix, returns same shifted to right
  # first repeats first value at the first time for obs i 

  nc  <- ncol(mat)
  tif <- cbind(c(1:nrow(mat)),first)
 
  x <- cbind(mat[,1],mat[,-ncol(mat)])
  x[ tif ] <- mat[ tif ]

  if(INC){   #increment first value

    f2 <- first+2
    f1 <- first+1
    
    f2[f2 > nc] <- nc
    f1[f1 >= f2] <- f2[f1 >= f2] - 1

    inc <- x[ cbind(c(1:nrow(mat)),f2) ] - x[ cbind(c(1:nrow(mat)),f1) ]
    inc[inc < minInc] <- minInc
    wf  <- which(first > 1)
    ff  <- x[ tif[wf,] ] - inc[wf]
    ff[ff < minInc] <- minInc
    x[ tif[wf,] ] <- ff

  }
    
  x
}


diffTimeMat <- function(mat,index,first = 0){  # increment matrix from obs by time matrix, pad last value

  md <- mat - lastTimeMat(mat,first) 
  vc <- md[index]
  list(dmat = md, vec = vc)
}


####################################################

tnorm.mvt <- function(avec,muvec,smat,lo=rep(-Inf,length(avec)),hi=rep(Inf,length(avec)),
                      whichSample=c(1:length(avec)),times=1){   

  # truncated multvariate normal
  # muvec is the vector of means
  # smat is the covariance matrix 
  # whichSample indicates which variables to sample

  if(length(lo) == 1)lo <- rep(lo,length(avec))
  if(length(hi) == 1)hi <- rep(hi,length(avec))

  for(j in 1:times){
   for(k in whichSample){

    tmp <- conditionalMVN(avec,muvec,smat,k)
    muk <- tmp$mu
    sgk <- tmp$vr

    if(length(muk) == 0)next

    avec[k] <- tnorm(1,lo[k],hi[k],muk,sqrt(sgk))
   }
  }
  avec
}


priorMVReg <- function(x,y,sigma){

  #Minka (2001) Bayesian Linear Regression

  # V - error covariance (sigma)
  # m - columns in x
  # d - columns in y
  #

  m <- ncol(x)
  d <- ncol(y)

  xx  <- crossprod(x)
  xy  <- crossprod(x,y)
  ixx <- invMat(xx)

  p1    <- solve(sigma)%*%t(xy)%*%ixx%*%xy
  p2    <- m*d
  alpha <- p2/(sum(diag(p1)) + p2) 

  priorAcov <- kronecker(sigma,ixx/alpha)

  list(alpha = alpha, priorAcov = priorAcov)
}


#######################################################

dmvnormLog <- function(x,mu,sigma){    #multivariate normal 

  #mv normal density

    if (is.vector(x))x <- matrix(x, ncol = length(x))

    distval <- mahalanobis(x, mu, sigma)
    logdet  <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    logretval
}


#############################################################

rwish <- function(df,S){

  z  <- matrix(rnorm(df*nrow(S)),df,nrow(S))%*%chol(S)
  crossprod(z)
}

#####################################################
riwish <- function(v,S){

  solve(rwish(v,solve(S)))
}


pmvnormCond <- function(q=0,xx,mu,smat,whichVar=c(1:nrow(smat))){   
  
  # conditional multvariate normal for Pr of each individual in vector conditioned on others
  # xx - matrix of values, n X p
  # mu - matrix of means, n X p
  # smat is the covariance matrix, p X p
  # whichVar indicates which variables to evaluate
  
  pj <- mu*0
  
  for(k in whichVar){
    tmp <- conditionalMVNVec(xx,mu,smat,k)
    muk <- tmp$mu
    sgk <- tmp$vr
    pj[,k] <- pnorm(q,muk,sd=sqrt(sgk),lower.tail=T)
  }
  pj
}




###################################
conditionalNorm <- function(x,mu,sigMat,cindex){  # bi- or trivariate
  
  if(ncol(x) == 2)return( conditionalBiVarNorm(x,mu,sigMat,cindex) )
  if(ncol(x) == 3)return( conditionalTriVarNorm(x,mu,sigMat,cindex) )
}
  
###########################################
conditionalTriVarNorm <- function(x,mu,sigMat,cindex){
  
  # cindex: conditional for this variable, condition on other variable
  # sigMat is nine columns
  
  nn <- nrow(x)
  nc <- length(cindex)
  r <- ncol(x)
  
  rs <- c(1:r^2)
  rm <- cm <- matrix(rs,r,r)
  tm <- rm
  tm[-cindex,] <- NA
  tm[,cindex]  <- NA
  rm[cindex,] <- rm[,cindex] <- NA
  cm[-cindex,] <- cm[,-cindex] <- NA
  
  notC <- sort(rm[is.finite(rm)])
  isC  <- sort(cm[is.finite(cm)])
  coC  <- sort(tm[is.finite(tm)])
  
  sinvMat <- invertcol3(sigMat)$I
  
  if(length(cindex) == 1){
    ct  <- c(1:r)[-cindex]
    p1  <- sigMat[,coC[1]] * sinvMat[,isC]
    p2  <- sigMat[,coC[2]] * sinvMat[,isC]
    mu1 <- mu[,cindex] + p1*(x[,ct[1]] - mu[,ct[1]]) + p2*(x[,ct[2]] - mu[,ct[2]])
    vr <- sigMat[,isC]
  }
  if(length(cindex) == 2){
    p1  <- sigMat[,coC[1]] * sinvMat[,isC[1]] + sigMat[,coC[2]] * sinvMat[,isC[2]]
    p2  <- sigMat[,coC[1]] * sinvMat[,isC[3]] + sigMat[,coC[2]] * sinvMat[,isC[4]]
    mu1 <- mu[,cindex[1]] + p1*(x[,-cindex] - mu[,-cindex])
    mu2 <- mu[,cindex[2]] + p2*(x[,-cindex] - mu[,-cindex])
    mu1  <- cbind(mu1,mu2)
    vr1 <- sigMat[,isC[1]] - p1*sigMat[,coC[1]]
    vr2 <- sigMat[,isC[2]] - p1*sigMat[,coC[2]]
    vr3 <- sigMat[,isC[3]] - p2*sigMat[,coC[1]]
    vr4 <- sigMat[,isC[4]] - p2*sigMat[,coC[2]]
    vr  <- cbind(vr1,vr2,vr3,vr4)
  }
  
  list(mu = mu1, vr = vr)
}
###########################################
conditionalBiVarNorm <- function(xx,mu,sigMat,cindex){

 # cindex: conditional for this variable, condition on other variable
 # sigMat is four columns

  isC  <- 1
  notC <- 4
  if(cindex == 2){
    isC <- 4
    notC <- 1
  }

  sinv <- 1/sigMat[,notC]

  p1  <- sigMat[,2]*sinv

  mu1 <- mu[,cindex] + p1*(xx[,-cindex] - mu[,-cindex])
  vr1 <- sigMat[,isC] - p1*sigMat[,2]

  list(mu = mu1, vr = vr1)
}
  ##################################3

conditionalMVNcdf <- function(up, xx, mu,sigma,k){  
  
  # xx and mu are n by p matrices for which we want the conditional CDF for k
  # k is index for conditional
  
  testv <- try(chol(sigma[-k,-k]),T)
  if(inherits(testv,'try-error')){
    return( list(mu = numeric(0), vr = numeric(0)) )
  }
  sin <- chol2inv(testv)
   
  p1  <- sigma[k,-k]%*%sin
  mu1 <- mu[,k] + t( p1%*%t(xx[,-k] - mu[,-k]) )
  vr1 <- sigma[k,k] - p1%*%sigma[-k,k]
  
  pnorm(up,mu1,sqrt(vr1))
}

###########################################
conditionalMVNVec <- function(xx, mu, sigma, cdex){ 

  #  xx, mu are n by p matrices
  #  sigma is p by p or n by p (one row per x)
  #  cdex is a vector of elements < p

  n <- nrow(xx)
  tiny <- min(diag(sigma))*.0001
  p   <- ncol(mu)
  if(ncol(xx) != p | ncol(sigma) != p)stop('different lengths in conditionalMVNVec')

  testv <- try(chol(sigma[-cdex,-cdex]),T)
  if(inherits(testv,'try-error')){
      return( list(mu = numeric(0), vr = numeric(0)) )
  }
  sinv <- chol2inv(testv)
  
  p1  <- sigma[cdex,-cdex]%*%sinv
  mu1 <- mu[,cdex] + t( p1%*%t(xx[,-cdex] - mu[,-cdex]) )
  vr1 <- sigma[cdex,cdex] - p1%*%sigma[-cdex,cdex]

  list(mu = mu1, vr = vr1)
}


#######################################

conditionalMVN <- function(xx, mu, sigma, cindex){  

  # x and mu are vectors, cindex is vector index for conditional

  tiny <- min(diag(sigma))*.0001
  nm   <- length(mu)
  if(length(xx) != nm)stop('x and mu different length in conditionalMVN')

  xx <- matrix(xx,nrow=1)
  mu <- matrix(mu,nrow=1)

  testv <- try(chol(sigma[-cindex,-cindex]),T)
  if(inherits(testv,'try-error')){
      return( list(mu = numeric(0), vr = numeric(0)) )
  }

  sin <- chol2inv(testv)
  p1  <- sigma[cindex,-cindex]%*%sin

  mu1 <- mu[cindex] +  p1%*%(xx[-cindex] - mu[-cindex]) 
  vr1 <- sigma[cindex,cindex] - p1%*%sigma[-cindex,cindex]

  list(mu = mu1, vr = vr1)
}
###################

processStates <- function(xchains,y,ADD=F,cols=1){
	
	xci <- apply(xchains,2,quantile,c(.5,.025,.975))
	
	yr <- range(xci)
	
	if(!ADD)plot(xci[1,],type='l',lwd=2,ylim=yr,col=cols)
	lines(xci[1,],type='l',lwd=2,ylim=yr,col=cols)
	for(j in 2:3)lines(xci[j,],lty=2,col=cols)
	points(y,col=cols[1])

        invisible(xci)
}

##########################################################
processPars <- function(xgb,xtrue=numeric(0),CPLOT=F,DPLOT=F,
                        sigOnly = F,burnin=1,xlimits = NULL){  

  #xg      - matrix of gibbs chains
  #xtrue   - true values (simulated data)
  #CPLOT   - if T, plot chains
  #DPLOT   - if T, plot density
  #burnin  - analyze chains > burnin
  #xlimits - xlimits for plot
  #sigOnly - plot only parameters that 95% CI does not include 0
  
  if(!is.matrix(xgb))xgb <- matrix(xgb,ncol=1)
  if(is.null(colnames(xgb)))colnames(xgb) <- paste('V',c(1:ncol(xgb)),sep='-')
  
  if(sigOnly){
    wi   <- grep('intercept',colnames(xgb))      #extract covariates for plotting
    btmp <- xgb
    if(length(wi) > 0){
    	btmp <- xgb[,-wi]
      if(length(xtrue) > 0)xtrue <- xtrue[-wi]
    }

    wq   <- apply(btmp,2,quantile,c(.025,.975),na.rm=T)  #extract parameters != 0
    wq   <- which(wq[1,] < 0 & wq[2,] > 0)
    if(length(wq) > 0){
      xgb  <- btmp[,-wq]
      if(length(xtrue) > 0)xtrue <- xtrue[-wq]
    }
   }

  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  if(burnin > 1){
  	     if(burnin > (nrow(xgb) + 100))stop("burnin too large")
  	     xgb <- xgb[-c(1:burnin),]
  }
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  nc <- ncol(xgb)
  nf <- round(sqrt(nc),0)

  out <- t(rbind(apply(xgb,2,mean,na.rm=T),apply(xgb,2,sd,na.rm=T),
                 apply(xgb,2,quantile,c(.025,.975),na.rm=T)))
  if(!is.null(colnames(xgb)))rownames(out) <- colnames(xgb)
  colnames(out) <- c('estimate','se','0.025','0.975')
  if(length(xtrue) > 0){
    out <- cbind(out,xtrue)
    colnames(out) <- c('estimate','se','0.025','0.975','true value')
  }

  armat <- matrix(0,nc,10)  #for AR model
  
 # for(j in 1:nc)armat[j,] <- ar(xgb[,j],aic=F,order.max = 10)$ar[1:10]

 # if(!is.null(colnames(xgb)))rownames(armat) <- colnames(xgb)
#  colnames(armat) <- paste('AR',c(1:10),sep='-')

  if(CPLOT | DPLOT)par(mfrow=c((nf+1),nf))
  if(CPLOT & DPLOT)par(mfrow=c((nf+1),nc))

  if(CPLOT){
      for(j in 1:nc){
       plot(xgb[,j],type='l')
       abline(h=out[j,],lty=2)
       if(length(xtrue) > 0)abline(h=xtrue[j],col='red')
       abline(h = 0, col='grey',lwd=2)
       title(colnames(xgb)[j])
     }
  }
  xlims <- xlimits
  if(DPLOT){
      for(j in 1:nc){
        xj <- density(xgb[,j])
        if(is.null(xlimits))xlims <- range(xj$x)
        plot(xj$x,xj$y,type='l',xlim=xlims)
        abline(v=out[j,],lty=2)
        if(length(xtrue) > 0)abline(v=xtrue[j],col='red')
        title(colnames(xgb)[j])
     }
  }
  list(summary = signif(out,4)
)

}

#####################
distmat <- function(xt,yt,xs,ys){
    xd <- outer(xt,xs,function(xt,xs) (xt - xs)^2)
    yd <- outer(yt,ys,function(yt,ys) (yt - ys)^2)
    t(sqrt(xd + yd)) 
}
#################################
parSummary <- function(x){

    if(!is.matrix(x)){
      xb <- c(mean(x),sd(x),quantile(x,c(.025,.975)))
    }
    if(is.matrix(x)){
      xb <- cbind(apply(x,2,mean,na.rm=T),apply(x,2,sd,na.rm=T),
            t(apply(x,2,quantile,c(.025,.975),na.rm=T)))
    }
    xb
}


multiLogitStates <- function(b,discStates,contStates,maxD){    # update of continuous states based on discrete states

      tiny <- 1e-20
      nn   <- length(discStates)
 
      discStates[is.na(discStates)] <- 0
      prob <- discStates*0 + 1

      if(maxD == 1)return(log(prob))

      wk <- which(discStates == 1 & is.finite(contStates),arr.ind=T)

      prob[wk] <- invlogt(b[1,],contStates[wk])

      if(maxD > 2){
        for(k in 2:(maxD-1)){
         wk       <- which(discStates == k & is.finite(contStates),arr.ind=T)
         prob[wk] <- invlogt(b[k,],contStates[wk]) - 
                     invlogt(b[(k-1),],contStates[wk])
        }
      }
      prob[discStates == maxD] <- 1 - invlogt(b[(maxD-1),],contStates[discStates == maxD])

    prob[prob < tiny] <- tiny
    log(prob)
}


#############################################

y2zMVlogit <- function(y){     #multivar logit to fractions

  if(!is.matrix(y)){
     z <- invlogit(y)
     return(cbind(z,1 - z))
  }

  zs   <- rowSums(exp(y))
  z1   <- 1/(1 + zs)
  zm   <- exp(y)/ (1 + zs)
  cbind(zm,z1)
}

z2yMVlogit <- function(z){     #fractions to multivar logit

  r <- ncol(z)
  if(r == 2){
    ss <- z[,1]/rowSums(z)
    return(log(ss/(1 - ss) ))
  }

  log(z[,-r]/(1 - rowSums(z[,-r])))
  
}

fisherR2Z <- function(r) {
  
  .5*(log(1 + r) - log(1 - r) )

}

compareCorMatrix <- function(r1,r2,n1,n2,alpha=.05){
  
  # zhou (2007) Psycholog Methods 12:399
  # ci - r1 significantly less than or greater than r2
  # sumDis - total distance
  # sumDisST - in units of standard deviations
  
  f1 <- fisherR2Z(r1)
  f2 <- fisherR2Z(r2)
  
  zscore <- (f1 - f2)/sqrt( 1/(n1 - 3) + 1/(n2 - 3) )
  
  p <- 2*pnorm(-abs(zscore))
  
  zalpha <- qnorm(alpha/2)
  zd     <- zalpha*sqrt(1/(n1 - 3))
  
  l1 <- f1 + zd
  u1 <- f1 - zd
  l2 <- f2 + zd
  u2 <- f2 - zd
  
  l1 <- (exp(2*l1) - 1)/(exp(2*l1) + 1)
  u1 <- (exp(2*u1) - 1)/(exp(2*u1) + 1)
  l2 <- (exp(2*l2) - 1)/(exp(2*l2) + 1)
  u2 <- (exp(2*u2) - 1)/(exp(2*u2) + 1)
  
  dr <- r1 - r2
  
  lo <- dr - sqrt( (r1 - l1)^2 + (u2 - r2)^2 )
  hi <- dr + sqrt( (u1 - r1)^2 + (r2 - l2)^2 )
  
  dr <- r1 - r2
  
  wl <- which( lo < 0 & hi < 0)
  wh <- which( lo > 0 & hi > 0)
  
  # stdevs from zero
  stdev <- (dr - lo)/1.96
  from0 <- dr/stdev
  
  cmat <- matrix(0,nrow(r1),nrow(r1))
  cmat[wl] <- -1
  cmat[wh] <- 1
  
  loTri    <- upper.tri(cmat)
  sumPos   <- length( which(cmat[loTri] == 1) )
  sumNeg   <- length( which(cmat[loTri] == -1) )
  fraction <- (sumPos + sumNeg)/length(cmat[loTri])
  meanDisST <- sum( abs(from0[loTri]) )/length(cmat[loTri])
  meanDis   <- sum( abs(dr[loTri]) )/length(cmat[loTri])
  
  
  list(ci = cmat, pvalue = p, standDist = from0, sumPos = sumPos, sumNeg = sumNeg,
       fraction = fraction, meanDis = meanDis, meanDisST = meanDisST)
}
  

clusterPlot <- function(dcor=NULL,dist=NULL,main=' ',xlab='Species',method='complete',
                        cex=1,ncluster=2,xlim=NULL, colCode = NULL){  
  #dcor is a correlation matrix

  
  if(!is.null(dist))diss <- as.dist( dist )
  if(is.null(dist))diss <- as.dist(cov2Dist(dcor))

  tmp  <- hclust(diss,method)
  ctmp <- cutree(tmp,k=1:ncluster)
  
  wclus <- ctmp[,ncluster]
  clusterCol <- NULL
  
  clusterIndex <- ctmp[,ncluster]
  
  clusterList <- character(0)
  
#  mycols <- mapColors(ncluster)
  
  colF   <- colorRampPalette(c('darkblue','green','orange','brown','red'))
  mycols   <- colF(ncluster)
  
  if(is.null(colCode)){
    colCode <- mycols[ctmp[,ncluster]]
    names(colCode) <- rownames(ctmp)
  }
  
  colLab <- function(n) {
 
    if(is.leaf(n)) {
      a <- attributes(n)
      
      attr(n, "nodePar") <-
        c(a$nodePar, list(col = colCode[n[1]],lab.col = colCode[n[1]]))
    }
    n
  }
  
  tdendro <- as.dendrogram(tmp)
  dL      <- dendrapply(tdendro,colLab)
  
  
#  plot(dL,type='triangle',cex=1.5)
  tmp <- plot(dL,nodePar=list(cex=cex,lab.cex=cex),horiz=T,xlim=xlim)
  title(main)
  
  invisible(list( clusterList = clusterList, colCode = colCode, clusterIndex = clusterIndex) )

}

dendrapply1 <- function (X, FUN, ...) 
{
  FUN <- match.fun(FUN)
  if (!inherits(X, "dendrogram")) 
    stop("'X' is not a dendrogram")
  Napply <- function(d) {
    r <- FUN(d, ...)
    if (!is.leaf(d)) {
      if (!is.list(r)) 
        r <- as.list(r)
      if (length(r) < (n <- length(d))) 
        r[seq_len(n)] <- vector("list", n)
      r[] <- lapply(d, Napply)
    }
    print(r)
    r
  }
  Napply(X)
}



dmvnormAR1 <- function(xx,psi=NULL,const,sinv=NULL,sig=NULL,ev=NULL){  
  
  #const - length nn vector of variance scalars
  #xx    - (x - mu)
  #sinv  - inverse AR1 covariance matrix
  #ev    - eigenvalues of sinv
  #must have either (sinv, ev) or (psi, variance)
  
  nn <- nrow(xx)
  nc <- ncol(xx)
  if(length(const) == 1)const <- rep(const,nn)
  
  if(is.null(sinv)){
    sinv <- invertAR1(nc,psi,sig)
    ev   <- eigen(sinv, only.values = TRUE)$values
  }
  
  s0 <- c(2:(nc-1))
  s1 <- c(1:(nc-2))
  s2 <- c(3:nc)
  
  x1 <- xx[,1]*(xx[,1] - psi*xx[,2])
  xn <- xx[,nc]*(xx[,nc] - psi*xx[,nc-1])
  
  r1 <- 1 + psi^2
  
  xi <- xx[,s0]*(xx[,s0]*r1 - (xx[,s1] + xx[,s2])*psi)
  
  expon <- rowSums( cbind(x1,xi,xn) )/sig
  
  distval <- expon/const
  logdet  <- -nc*log(const) - sum(log(ev))
  
  -(nc * log(2 * pi) + logdet + distval)/2
}

###########################################################


pmvnormAR1_approx <- function(neach,xx,psi,const,sig,lower=-Inf,upper=Inf,stages=ncol(xx)){
  
  # neach - number of samples per row of xx
  # xx    - (x - mu)
  # const - length nn vector of variance scalars
  
  xx <- xx[,stages]    # may not include large sizes, mostly empty
  nc <- length(stages)
  
  nn <- nrow(xx)
  
  oo <- makeAR1(nc,psi,ovar)$varmat
  
  
  mm <- nn*neach
  if(length(lower) == 1)lower <- matrix(lower,mm,nc)
  if(length(upper) == 1)upper <- matrix(upper,mm,nc)
  
  tmp <- myrmvnorm(mm,0,oo)
  ww  <- which(tmp > lower & tmp < upper,arr.ind=T)
  
  tmp <- tmp*0
  tmp[ww] <- 1
  
  tmp <- rowSums(tmp)
  tmp[tmp < nc] <- 0
  tmp[tmp != 0] <- 1
  
  ii <- rep(1:nn,each=neach)
  
  tmp <- byFunctionRcpp(tmp,ii,ii*0+1,matrix(0,nn,1),matrix(0,nn,1),MEAN=F)/neach
  tmp
}
#######################3 

dmvnormZeroMean <- function(xx,smat=NULL,sinv=NULL){          #MVN density for mean 0

  if(!is.matrix(xx))xx <- matrix(xx,1)
  
  if(!is.null(sinv)){
    distval <- diag( xx%*%sinv%*%t(xx) )
    ev      <- eigen(sinv, only.values = TRUE)$values
    logdet  <- -sum(log(ev))
  }
  
  if(is.null(sinv)){
    testv <- try(chol(smat),T)
    if(inherits(testv,'try-error')){
       tiny <- min(abs(xx))/100 + 1e-5
       smat <- smat + diag(diag(smat + tiny))
       testv <- try(chol(smat),T)
    }
    covm <- chol2inv(testv)
    distval <- rowSums((xx %*% covm) * xx)
    ev <- eigen(smat, only.values = TRUE)$values 
    if(min(ev) < 0)ev <- nearPD(smat)$eigenvalues
    logdet   <- sum(log( ev ))
  }

    -(ncol(xx) * log(2 * pi) + logdet + distval)/2
}
########################################

getChainMeanVar <- function(chains){  #mean vector, covariance matrix from chains matrix

  chains <- row2Mat(chains)
  nn     <- nrow(chains)
  if(nn == 1){
    mub    <- chains
    vrb    <- diag(1,ncol(chains))
  }
  if(nn > 1){
     mub   <- colMeans(chains)
     vrb   <- cov(chains)
  }
  list(mu = mub, vr = vrb)
}
#####################################
getPost <- function(b,bchain){     #posterior estimate for marginal likelihood

  #gaussian
  # method of Chib 1995, JASA

 # b <- row2mat(b)

  tmp <- getChainMeanVar(bchain)
  mub <- tmp$mu
  vrb <- tmp$vr
  
  mu  <- matrix(mub,nrow(bchain),ncol(bchain),byrow=T) - bchain
 
  log( mean(exp( dmvnormZeroMean(mu,vrb) )) )

}

#######################################3
biVarMoments <- function(x1,x2,wt=1,PLOT = F, POINTS=F, color=1, pr = .95, ADD=F, lty=1,lwd=1,
                        xylab=c('x','y')){  

  #x1, x2 - vectors for variables, wt is weight

  require(cluster)

  if(length(pr) > 1){
     if(length(lty) == 1)lty = rep(lty,length(pr))
     if(length(lwd) == 1)lwd = rep(lwd,length(pr))
     if(length(color) == 1)color = rep(color,length(pr))
  }
  
  if(length(wt) == 1)wt <- rep(wt,length(x1))

  ww <- which(is.finite(x1) & is.finite(x2))

  x1 <- x1[ww]
  x2 <- x2[ww]
  wt <- wt[ww]

  w1 <- x1*wt
  w2 <- x2*wt
  m1 <- sum(w1)/sum(wt)
  m2 <- sum(w2)/sum(wt)

  v1  <- sum(wt*x1^2)/sum(wt) - m1^2
  v2  <- sum(wt*x2^2)/sum(wt) - m2^2
  c   <- sum(wt*(x1 - m1)*(x2 - m2))/sum(wt)

  for(k in 1:length(pr)){
      tmp <- list(loc = c(m1,m2), cov = matrix(c(v1,c,c,v2),2,2), d2 = qchisq(pr[k],1) )
      tmp <- predict.ellipsoid(tmp)
      if(PLOT & !ADD & k == 1)plot(tmp[,1],tmp[,2],type='l',col=color[k],lty=lty[k],
                              lwd=lwd[k],xlab=xylab[1], ylab=xylab[2])
      if(POINTS)points(x1,x2,cex=.3,col='grey')
      if(ADD | k > 1)lines(tmp[,1],tmp[,2],type='l',col=color[k],lwd=lwd[k],lty=lty[k])
  }
  

  invisible( list(loc = c(m1,m2), 
                    cov = matrix(c(v1,c,c,v2),2,2) ,
                    d2 = qchisq(pr[1],1), ellipse = tmp) )
}
#############################################################
pmake <- function(pars){  

  p   <- pars[1]                   #frequency of a
  f   <- pars[2]                   #inbreeding coefficient 
  paa <- p^2 + f*p*(1 - p)	  #frequency of p.aa
  pab <- 2*p*(1 - p)*(1 - f)	  #frequency of p.ab
  pbb <- (1 - p)^2 + f*p*(1 - p)  #frequency of p.bb
  c(paa,pab,pbb)
}

minf <- function(p){   #minimum inbreeding coefficient

  lof <- rep(-1,length(p))
  lop <- -(1 - p)/p
  hip <- -p/(1 - p)
  lof[p > (1 - p)] <- lop[p > (1 - p)]
  lof[p < (1 - p)] <- hip[p < (1 - p)]
  lof
}
pf.update <- function(p,f){  #log posterior

  multlik(c(p,f)) + dbeta(p,g1,g2,log=T) + 
           dnorm(f,priorF,priorFSD,log=T)
}

multlik <- function(pars){  #multinom likelihood
 
  p    <- pmake(pars)
  pmat <- matrix(rep(p,npop),3,npop)
  sum(y*log(pmat))

}

update_pf <- function(){

  propp  <- tnorm(1,.02,.98,pg,.002) #propose pa
  fl     <- minf(propp)              #lower f limit for pa
  propf  <- tnorm(1,fl,1,fg,.01)     #propose f > fl
  pnew   <- pf.update(propp,propf)
  pnow   <- pf.update(pg,fg)
  atmp   <- acceptMH(pnow,pnew,c(pg,fg),c(propp,propf),BLOCK=T)
  pg     <- atmp$x[1]
  fg     <- atmp$x[2]
  ac     <- atmp$accept

  list(pg = pg, fg = fg, accept = ac)
}

updateVariance <- function(yy,mu,s1=1,s2=1,lo=NULL,hi=NULL){
  
  # if yy is a matrix, one variance for each column
  
  require(pscl)

  tiny <- 1e-10
  
  k <- 1
  res <- (yy - mu)^2
  nn  <- length(yy)
  
  if(is.matrix(yy)){
    k  <- ncol(yy)
    nn <- nrow(yy)
    sr <- colSums(res)
  }else{
    sr <- sum(res)
  }
  

  u1 <- s1 + nn/2
  u2 <- s2 + .5*sr

  if(is.null(lo) & is.null(hi))return( 1/rgamma(k,u1,u2) )
  
  if(is.null(lo))lo <- 0

  rtrunc(k,lo,hi,u1,u2,'igamma')

}


rtrunc <- function(nn,lo,hi,p1,p2,FN){    #truncated using inv dist sampling
   #truncated 2 parameter distribution
	
  require(pscl)
  
  if(length(lo) == 1)lo <- rep(lo,nn)
  if(length(hi) == 1)hi <- rep(hi,nn)

  pf <- paste('p',FN,sep='')
  qf <- paste('q',FN,sep='')
  pf1 <- match.fun(pf)
  qf1 <- match.fun(qf)
  
  if(FN == 'igamma'){
    z1 <- 1 - pgamma(1/lo, p1, p2)
    z2 <- 1 - pgamma(1/hi, p1, p2)
  }
  if(FN != 'igamma'){
    z1 <- pf1(lo,p1,p2)
    z2 <- pf1(hi,p1,p2)
  }
  
  z  <- runif(nn,z1,z2)
  
  r  <- hi
  wr <- which(z > 0)
  
  if(length(wr) == 0)return(hi)
  
  if(FN == 'igamma')r[wr] <- 1/qgamma(1 - z[wr],p1[wr],p2[wr])
  if(FN != 'igamma')r[wr] <- qf1(z[wr],p1[wr],p2[wr])
  r[r < lo] <- lo[r < lo]
  r
}


pasteChar2End <- function(x,char2paste) paste(x,char2paste,sep='') #paste character to end

capFirstLetter <- function(x) {     #capiltalize first letter of every word
   s <- unlist(strsplit(x, " "))
   s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
               sep = "", collapse = " ")
   unlist(strsplit(s, " "))
}

capSubString <- function(x,string) {  #replace substring with caps
   s <- unlist(strsplit(x, string))
   ww <- which(nchar(s) == 0)
   s[ww] <- toupper(x[ww])
   s
}

replaceString <- function(xx,now='_',new=' '){  #replace now string in vector with new
  
  ww <- grep(now,xx)
  for(k in ww){
    s  <- unlist( strsplit(xx[k],now) )
    ss <- s[1]
    for(kk in 2:length(s)) ss <- paste( ss,s[kk],sep=new)
    xx[k] <- ss
  }
  xx
}



acceptMH <- function(p0,p1,x0,x1,BLOCK=F){   #accept for M, M-H
	# if BLOCK, then accept as a block,
	# otherwise, accept individually

  nz          <- length(x0)  #no. to accept
  if(BLOCK)nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a,arr.ind=T)
  
  if(BLOCK & length(keep) > 0)x0 <- x1
  if(!BLOCK)                  x0[keep] <- x1[keep]           
  ac <- length(keep)        

  list(x = x0, accept = ac)
}

predVsObs <- function(o,p,ylim=range(p,na.rm=T),xlab=' ',ylab=' ',colors=rep(1,length(o))){ 
	
  #o  - length n vector of obs or true values
  #p - ng by n matrix of estimates
  
  n <- length(o)
  y <- apply(p,2,quantile,c(.5,.025,.975))

  plot(o,y[1,],ylim=ylim,xlab=xlab,ylab=ylab,col=colors)
  for(j in 1:n)lines(c(o[j],o[j]),y[2:3,j],col=colors[j])
  abline(0,1,lty=2)
  invisible(y)
}



profilePlot <- function(x,y,proValues=quantile(x,seq(.1,.9,by=.2)),ylim=range(y),xlab=' ',ylab=' ',yscale=1){

  plot(x,y,xlab=xlab,ylab=ylab,ylim=ylim,cex=.3,bty='n',col='grey')

  nv <- length(proValues)


  wf <- which(is.finite(x) & is.finite(y))
  x  <- x[wf]
  y  <- y[wf]

  mids <- c(proValues[-nv] + diff(proValues)/2,max(x))
  mids <- c(proValues[1] - diff(proValues[1:2])/2,mids)

  tiny <- 1e-10
  huge <- 1 - tiny

  y[y < tiny] <- tiny
  y[y > huge] <- huge

  fi   <- findInterval(x,mids)
  fi[fi > nv] <- nv

  for(j in 1:nv){

    wj  <- which(fi == j)

    tmp <- hist(y[wj],breaks=quantile(y[wj],seq(0,1,length=20)) ,plot=F)
    yy  <- tmp$density
    yy  <- yy * yscale
    yy[!is.finite(yy)] <- 0
    yy  <- yy  + proValues[j]
    xx  <- tmp$mids

    yy  <- c(proValues[j],yy,proValues[j])
    xx  <- c(xx[1],xx,xx[1])
    
    lines(yy,xx,type='s',lwd=2)
  }
 # points(x,y,cex=.3,col='grey')
 
}


#######################################################
plotObsPred <- function(obs,yMean,ySE=NULL,nbin=NULL,nPerBin=NULL,breaks=NULL,
                        log=F,xlimit=NULL,ylimit=NULL,xlabel='Observed',ylabel='Predicted',
                        ptcol=NULL,boxPerc = .6826895, whiskerPerc = .95,
                        fill=NULL,add=F,box.col='black',wide=NULL,POINTS=T){
  
  aa <- (1 - boxPerc)/2
  boxQuant <- c(aa, 1 - aa )
  aa <- (1 - whiskerPerc)/2
  whiskerQuant <- c(aa, 1 - aa )
  
  if(is.null(ptcol)){
    ptcol <- 'black'
    ptcol <- 'grey'
    if(!is.null(nbin))ptcol <- 'grey'
  }
  if(length(ptcol) == 1)ptcol <- rep(ptcol,length(obs))
  
  if(is.null(xlimit))xlimit <- range(obs)
  if(is.null(ylimit))ylimit <- range(yMean)
  
  xxx <- obs
  yyy <- yMean
  
  if(!POINTS){
    xxx <- 0
    yyy <- 0
  }
  
  if(is.null(xlimit))xlimit <- range(obs)
  
  if(!add){
    if(is.null(ylimit)){
      if(!log)plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel)
      if(log) plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,log='xy')
    }
    if(!is.null(ylimit)){
      if(!log)plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,xlim=xlimit,ylim=ylimit)
      if(log) plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,xlim=xlimit,log='xy',ylim=ylimit)
    }
  }
  if(!is.null(ySE)){
    ylo <- yMean - 1.96*ySE
    yhi <- yMean + 1.96*ySE
    for(i in 1:length(obs))lines(c(obs[i],obs[i]),c(ylo[i],yhi[i]),col='grey',lwd=2)
  }
  
 # if(!is.null(nbin) | !is.null(nPerBin)){
    
  if(is.null(breaks)){
    
    if(is.null(nbin))nbin <- 20
    br   <- range(obs,na.rm=T)
    bins <- seq(br[1],br[2],length=nbin)
    if(log)bins <- 10^seq(log10(br[1]),log10(br[2]),length=nbin)
    
    if(!is.null(nPerBin)){
      nbb <- nPerBin/length(obs)
      nbb <- seq(0,1,by=nbb)
      if(max(nbb) < 1)nbb <- c(nbb,1)
      bins <- quantile(obs,nbb,na.rm=T)
      nbin <- length(bins)
    }
  } else {
    bins <- breaks
    nbin <- length(bins)
  }
    
    if(is.null(wide))wide <- diff(bins)/2
    for(k in 1:(nbin-1)){
      qk <- which(is.finite(yMean) & obs >= bins[k] & obs <= bins[k+1])
      q  <- quantile(yMean[qk],c(.5,whiskerQuant[1],boxQuant[1],
                                 boxQuant[2],whiskerQuant[2]),na.rm=T)
      if(log){
        q[q <= 0] <- .001
      }
      ym <- mean(yMean[qk])
      xx <- mean(bins[k:(k+1)])
      yy <- q[c(2,5)]
      yy[1] <- max(yy[1],ylimit[1],na.rm=T) + .0000001
      yy[2] <- max(yy)
      lines(c(xx,xx),yy,lwd=2,col=box.col)
   #   lines(c(xx-.5*wide[k],xx+.5*wide[k]),q[c(1,1)],lwd=4,col=box.col)
      lines(c(xx-.5*wide[k],xx+.5*wide[k]),c(ym,ym),lwd=2,col=box.col)
      
      yy1 <- q[3]
      yy1 <- max(yy1,ylimit[1],na.rm=T) + .00000001
      yy2 <- max(yy1,q[4])
      rect(xx-.4*wide[k],yy1,xx+.4*wide[k],yy2,col=fill,border=box.col)
    }

  
}


pars2List <- function(pars){    #pars is a named list, filled with values here
  
  for(k in 1:length(pars))pars[[k]] <- get(names(pars)[k])
  pars
}

#######################################################
deviance <- function(y,x,b,s=0,LIKE='norm'){
	
	require(mvtnorm)
	    
    if(LIKE == 'norm')  dv <- dnorm(y,x%*%b,sqrt(s),log=T) 
    if(LIKE == 'pois')  dv <- dpois(y,exp(x%*%b),log=T)
    if(LIKE == 'binom') dv <- dbinom(y,1,invlogit(x%*%b),log=T)
    if(LIKE == 'mvnorm')dv <- dmvnormZeroMean(y - x%*%b,s)
    if(LIKE == 'multinom')dv <- multinomLike(y,x,b)$like
    -2*dv
}
#####################################33
dinvGamma <- function(x,a,b,log=FALSE){
	
	p <- a*log(b) - lgamma(a) - (a+1)*log(x) - b/x
	if(log)return(p)
	exp(p)
}
###########################################
dbivarNormFromCols <- function(yy,mu,S){    #bivariate norm, covariances for each y in columns, log density

  #yy  - n by 2
  #mu - n by 2 or scalar 0
  #S  - n by 4: S11, S12, S21, S22

  if(!is.matrix(yy)){
    yy <- matrix(yy,1)
    mu <- matrix(mu,1)
    S  <- matrix(S,1)
  }

  n <- nrow(yy)
  if(length(mu) == 0)mu <- yy*0

  yy <- yy - mu

  invS <- invertcol2(S)
  z    <- yy[,1]^2*(invS[,1] + invS[,2])+ yy[,2]^2*(invS[,3] + invS[,4])

  ldet <- log(S[,1]*S[,4] - S[,2]*S[,3])

  -(n*log(2*pi) + ldet + z)/2
}

###########################################
dtrivarNormFromCols <- function(yy,mu,S){    #trivariate norm, covariances for each y in columns, log density
  
  #y  - n by 3
  #mu - n by 3 or scalar 0
  #S  - n by 9: S11, S12, S31, S21, S22, ...
  
  n <- nrow(y)
  if(length(mu) == 0)mu <- y*0
  
  yy <- yy - mu
  
  tmp <- invertcol3(S)
  invS <- tmp$I
  D    <- tmp$D
  
  y1 <- cbind( rowSums(yy*invS[,1:3]), rowSums(yy*invS[,4:6]), rowSums(yy*invS[,7:9]) )
  z  <- rowSums(y1*yy)
  
  -(n*log(2*pi) + log(D) + z)/2
  
}


#############################################################

invertcol2 <- function(S){   
    # inverse of matrix supplied and returned as nX4 matrix: var1, cov12, cov12, var2

  dt  <- S[,1]*S[,4] - S[,2]*S[,2]
  I1  <- S[,4]
  I2  <- S[,1]
  I12 <- -S[,2]
  I21 <- -S[,2]

  cbind(I1,I12,I21,I2)/dt
}

#############################################################

invertcol3 <- function(S){   
  # inverse of matrix supplied and returned as nX9 matrix: 
  # var1, cov12, cov13, cov21, var2, cov31, cov31, cov32, var3
  
  dt  <- S[,1]*S[,5]*S[,9] + S[,2]*S[,6]*S[,7] + S[,3]*S[,4]*S[,8] - 
        (S[,3]*S[,5]*S[,7] + S[,1]*S[,6]*S[,8] + S[,2]*S[,4]*S[,9])
    
  I1  <- S[,5]*S[,9] - S[,6]*S[,6]
  I2  <- S[,1]*S[,9] - S[,3]*S[,3]
  I3  <- S[,1]*S[,5] - S[,2]*S[,2]
  I12 <- S[,3]*S[,6] - S[,2]*S[,9]
  I13 <- S[,2]*S[,6] - S[,3]*S[,5]
  I23 <- S[,2]*S[,3] - S[,1]*S[,6]
  
  list(I = cbind(I1,I12,I13,I12,I2,I23,I13,I23,I3)/dt, D = dt)
}



############################################################

initialStatesSS <- function(y){

  require(stats)

  if(!is.matrix(y))y <- matrix(y,ncol=1)
  r <- ncol(y)

  n    <- length(y)
  time <- c(1:n)
  wm   <- which(is.na(y))
  notMiss <- c(1:n)
  
  x <- y

  if(length(wm) > 0){
   notMiss <- notMiss[-wm]
   x[wm]   <- predict(smooth.spline(time[-wm],y[-wm]),time)$y[wm]
  }

  list(x = as.vector(x), notMiss = notMiss, miss = wm)
}


############################################################
condProb <- function(mu,V,y){     #bivariate conditional on y for 2nd variables

  cm <- mu[1] + V[1,2]/V[2,2]*(y - mu[2])
  cV <- V[1,1] - V[1,2]/V[2,2]*V[2,1]

  list(muY = cm, cV = cV)
}
###########################################

multiLogitStates <- function(b,discStates,contStates,maxD){    # continuous based on discrete states

      tiny <- 1e-20
      nn   <- length(discStates)
 
      discStates[is.na(discStates)] <- 0
      prob <- discStates*0 + 1

      if(maxD == 1)return(log(prob))

      wk <- which(discStates == 1 & is.finite(contStates),arr.ind=T)

      prob[wk] <- invlogt(b[1,],contStates[wk])

      if(maxD > 2){
        for(k in 2:(maxD-1)){
         wk <- which(discStates == k & is.finite(contStates),arr.ind=T)
         prob[wk] <- invlogt(b[k,],contStates[wk]) - 
                                  invlogt(b[(k-1),],contStates[wk])
        }
      }
      prob[discStates == maxD] <- 1 - invlogt(b[(maxD-1),],contStates[discStates == maxD])

    prob[prob < tiny] <- tiny
    log(prob)
}

############################################
breaks2pars <- function(cmat,breaks){        #break points and coefficents to ordinal multinomial logit pars

  nd <- length(breaks)
  c0 <- rep(0,nrow(cmat))

  for(k in 1:nd){

    qk <- 0
    if(k > 1){
      for(j in 1:(k-1)){
       ej <- exp(c0[j] + cmat[j,2]*breaks[k])
       qk <- qk + ej/(1 + ej)
      } 
    }
    D     <- -log(1/(.5 + qk) - 1)
    cc0   <- D - cmat[k,2]*breaks[k]
    if(k > 1){
      if(cc0 < cmat[k-1,1]){
        cc0 <- cmat[k-1,1] + .5
        cmat[k,2] <- (D - cc0)/breaks[k]
      }
    }
    cmat[k,1] <- cc0
  }
  cmat
}
#################################################
proposeBreak <- function(cmat,br,brLims){      #propose new breakpoints given limits brLims

  c1 <- cmat[,1]
  c2 <- cmat[,2]

  nd <- nrow(cmat)
  brNew <- tnorm(nd,brLims[1,],brLims[2,],br,.1)

  if(nd == 1){
     c1New <- tnorm(1,1,100,cmat[1],.1)
     c2New <- (log(.5) + log(1 + exp(c1New + cmat[2]*brNew)) - c1New)/brNew
     cnew <- matrix(c(c1New,c2New),1,2)
  }

  if(nd > 1)cnew <- getCmat(cmat,brNew)

  list(cmat = cmat,br = brNew)
}
##########################################
pars2p <- function(cmat,h){               #multinomial logit pars and scale h to Pr for multinom logit

  nd <- nrow(cmat)
  nh <- length(h)

  c1 <- matrix(cmat[,1],nh,nd,byrow=T)
  c2 <- matrix(cmat[,2],nh,nd,byrow=T)
  hh <- matrix(h,nh,nd)
  eh <- exp(c1 + c2*hh)

  theta <- matrix(0,nh,nd+1)
  sumt  <- rep(0,nh)

  for(k in 1:nd){

     tk   <- eh[,k]/(1 + eh[,k]) - sumt
     theta[,k] <- tk 
     sumt <- sumt + tk
  }
  theta[,nd+1] <- 1 - rowSums(theta)

  theta
}


######################################################
plotLogit <- function(lims,breaks,h,cmat){  #plot multinomial ordinal logit

  tmp  <- pars2p(cmat,h)

  plot(h,tmp[,1],type='l',lwd=2,xlab='latent health scale',ylab='Probabilty')

  for(j in 1:ncol(lims)){
 #   polygon(c(lims[,j],rev(lims[,j])),c(0,0,1,1),border=NA,col='grey')
    lines(h,tmp[,j],col=j,lwd=2)
    l1   <- 0
    if(j > 1)l1 <- breaks[j-1]
    midx <- (l1 + breaks[j])/2
    text(midx,.8,j,col=j,cex=1.4)
  }
  lines(h,tmp[,j+1],col=j+1,lwd=2)
  midx <- ( breaks[j] + 100 )/2
  text(midx,.8,j+1,col=j+1,cex=1.4)
  abline(h=.5,lty=2)
}   	
###########################

plotSetup <- function(xtic,ytic,xvals = xtic, yvals = ytic, xlabel=' ',ylabel=' ',
                      endFactor=c(.05,.05),
                      fc = NULL,lcHor=NULL,lcVer=NULL){

  xr <- range(xtic)
  yr <- range(ytic)

  xlimit <- c(xr[1] - diff(xr)*endFactor[1],xr[2] + diff(xr)*endFactor[1])
  ylimit <- c(yr[1] - diff(yr)*endFactor[2],yr[2] + diff(yr)*endFactor[2])

  plot(-100,-100,xlim=xlimit,ylim=ylimit,ylab=ylabel,
       cex=1.2,xlab=xlabel,xaxt='n',yaxt='n')

  if(!is.null(fc))rect(xlimit[1],ylimit[1],xlimit[2],ylimit[2],col=fc,border=NA)

  axis(1,at=xtic,labels=xvals,cex.lab=1.2)
  axis(2,at=ytic,labels=yvals,cex.lab=1.2)

  if(!is.null(lcVer))abline(v=xtic,lwd=3,col=lcVer)
  if(!is.null(lcHor))abline(h=ytic,lwd=3,col=lcHor)
}


############################### set up scale for a map
mapSetup <- function(xlim,ylim,scale){  #scale is m per inch

  px   <- diff(xlim)/scale
  py   <- diff(ylim)/scale
  par(pin=c(px,py))
  
  invisible( c(px,py) )

}


getQuant <- function(data,dim,q){

  quantile(data,dim,q,na.rm=T)
}

cov2Dist <- function(sigma){ #distance induced by covariance
	
	n <- nrow(sigma)
	matrix(diag(sigma),n,n) + matrix(diag(sigma),n,n,byrow=T) - 2*sigma
}


myBarPlot <- function(dataMat,widthFactor=1,xlab=' ',ylab=' ',ylim=NULL,barLabels=NULL,
                      textFactor=1,stack=F){
  
  #dataMat - one row per bar, columns are bar segments
  
  nc <- ncol(dataMat)
  nr <- nrow(dataMat)
  
  if(is.null(ylim))ylim=c(range(rowSums(dataMat,na.rm=T)))
  
  wide <- widthFactor
  
  plot(c(0,0),c(0,0),xlim=c(0,nr),cex=.1,ylim=ylim,xlab = xlab,ylab=ylab)
  
  ystart <- numeric(0)
  
  if(!stack){
    rr    <- range(dataMat)
    space <- diff(rr)*.2
    tmp   <- apply(dataMat,2,max,na.rm=T)
    ystart <- cumsum( tmp + space )
  }
  
  for(j in 1:nr){
    
    xj <- c(j - wide,j + wide)
    
    y1 <- 0
    
    for(k in 1:nc){
      
      y2 <- y1 + dataMat[j,k]
      rect(xj[1],y1,xj[2],y2,col=k,border=k)

      y1 <- y2
      if(!stack)y1 <- ystart[k]
    }
    if(!is.null(barLabels))text(xj[1],y2+.1*diff(range(ylim)),
                                barLabels[j],srt=90,pos=4,cex=textFactor)
  }
  
  invisible(ystart)
}

####################################################

myBoxPlot <- function(mu,limits,ORDER = F,myorder=NULL,xvalues=c(1:length(mu)),
                      boxcol='black',xvals=xtic,yvals=ytic,
                      xlab=' ',ylab=' ',
                      xtic=c(1,length(mu)),ytic=NULL,plabels=NULL,
                      shadelo=NULL,shadehi=NULL,add=F,widthScale=1){

  #mu has median
  #limits - a list with lo,hi as paired columns, inner to outer
  #shadehi - shade above this value
  #shadelo - shade below this value

  n  <- length(mu)
  xi <- xvalues
  or <- c(1:n)
  if(ORDER & is.null(myorder))or <- order(mu,decreasing=T)
  if(!is.null(myorder))or <- myorder
  
  if(length(boxcol) == 1)boxcol <- rep(boxcol,n)

  x  <- limits
  nq <- length(x)

  if(is.null(ytic))ytic <- signif(range(mu),1)

  if(!add){
    plotSetup(xtic,ytic,xvals=xvals,yvals=yvals,xlabel=xlab,ylabel=ylab)
  }

  if(!is.null(shadehi))rect(xtic[1],shadehi,max(xtic),max(ytic),col='bisque1',border=NA)
  if(!is.null(shadelo))rect(xtic[1],min(ytic),max(xtic),shadelo,col='bisque1',border=NA)

  points(xi,mu[or],pch=3,col=boxcol)

 # ji <- c(1:nq)/diff(range(xi))*3

  ji <- .2*diff(xi)
  ji <- c(ji[1],ji)

  jk <- seq(.4,.6,length=nq)

  for(i in 1:n){

    maxi <- ytic[1]

    for(k in 1:nq){

      yi <- x[[k]][or[i],]
      if(is.na(yi[1]))next
      
      ki <- jk[nq - k + 1]*ji[i]

      x1 <- xi[i] - ki*widthScale
      x2 <- xi[i] + ki*widthScale

      rect(x1,yi[1],x2,yi[2],lwd=1,col=boxcol[i],border=boxcol[i])
      if(k == nq)lines(c(xi[i],xi[i]),yi,lwd=1,col=boxcol[i])
      if(max(yi) > maxi)maxi <- max(yi)
    }

    if(maxi > max(ytic))maxi <- max(ytic) - .1*diff(range(ytic))
    if(!is.null(plabels))text(i-.4,maxi,plabels[or[i]],srt=90,pos=4,col=boxcol[i])
  }

  invisible(or)  #return order

}

trim <- function(xx,p){     #p - lower,upper %tile

  xq <- quantile(xx,p,na.rm=T)
  
  xx[xx < xq[1]] <- xq[1]
  xx[xx > xq[2]] <- xq[2]

  xx
}

#################################################
weightAve <- function(x,wt){  #x vector for variable, wt is weight

  w1 <- x*wt
  m1 <- sum(w1)/sum(wt)

  v1  <- sum(wt*x^2)/sum(wt) - m1^2

  c(m1,v1)
}

cov2cor <- function(covmat){  #covariance matrix to correlation matrix

  d    <- nrow(covmat)
  di   <- diag(covmat)
  s    <- matrix(di,d,d)
  covmat/sqrt(s*t(s))
}

cor2cov <- function(sigvec,cormat){ #correlation matrix and variance vector to covariance

  d <- length(sigvec)
  s <- matrix(sigvec,d,d)
  cormat*sqrt(s*t(s))
}

  


heatGrid <- function(x,groupIndex=NULL,legendTick=c(0,1),htScale=1){

  nr   <- nrow(x)
  nc   <- ncol(x)

  ncol <- 100
  colseq <- heat.colors(ncol)

  scale <- seq(min(x,na.rm=T),max(x,na.rm=T),length=ncol)

  x[is.na(x)] <- 0

  ww   <- as.matrix(expand.grid(c(1:nr),c(1:nc)))

  icol <- findInterval(x,scale,all.inside=T)
  coli <- colseq[icol]

  bwide <- rep(1,nrow(ww))
  bhigh <- rep(1/nr/5,nrow(ww))*htScale

  symbols(ww[,2],ww[,1]+1,rectangles=cbind(bwide,bhigh),xlim=c(0,nc+2),ylim=c(0,nr+2),
              fg=coli,bg=coli,inches=F,axes=F,xlab=' ',ylab=' ')
  
  cc <- rownames(x)
  if(!is.null(cc))text(rep(nc+1,nc),1 + c(1:nr),cc,cex=.6)
  rr <- colnames(x)
  if(!is.null(rr))text(1:nc,rep(nr+2,nc),rr,srt=45,cex=.6)

  if(!is.null(groupIndex)){
    groups <- unique(groupIndex)
    ns     <- length(groups)

    for(k in 1:ns){

       wk <- which(groupIndex == groups[k])
       points(rep(1 + .2*k,length(wk)),wk,col='white',pch=7,cex=.5)
    }
  }

  colorLegend(c(.1,.4),nr*c(.9,1),ytic=legendTick,scale=scale,cols=colseq,labside='left')
}

#####################################################
corPlot <- function(cmat,slim=NULL,diag0=F,textSize=1){  #correlation or covariance matrix

 # if(diag0)diag(cmat) <- 0

  cmat[lower.tri(cmat)] <- 0

  d <- nrow(cmat)
  xtext <- rep(c(1,100),d/2)
  if(length(xtext) < d)xtext <- c(xtext,1)

  if(d < 20)xtext <- xtext*0 + 1

  xtext <- xtext*0 + 1

  ncol <- 100
  colseq <- mapColors(ncol)


  if(is.null(slim))slim = range(cmat)
  slim  <- signif(slim,1)
  scale <- seq(slim[1],slim[2],length.out=ncol)

  ww   <- as.matrix(expand.grid(c(1:d),c(1:d)))

  ww  <- ww[ww[,1] <= ww[,2],]
  ww  <- ww[order(ww[,1]),]

  icol <- findInterval(cmat[ww],scale,all.inside=T)
  coli <- colseq[icol]

  if(diag0)coli[ww[,1] == ww[,2]] <- 'white'

  symbols(ww[,1],ww[,2]+1,squares=rep(.9,nrow(ww)),xlim=c(0,d+2),ylim=c(0,d+2),
              fg=coli,bg=coli,inches=F,axes=F,xlab=' ',ylab=' ')

  for(kk in 1:(d+1)){
        ks <- kk - .5
        if(kk <= d)lines(c(ks,ks),c(ks+1,d+1.5),col='grey')
        if(kk > 1) lines(c(.5,ks-1),c(ks,ks),col='grey')
        if(kk <= d)lines(c(ks,ks+1),c(ks+1,ks),col='grey')
  }
  lines(c(.5,ks-1),c(ks+1,ks+1),col='grey')
  text(c(1:d)+.1*xtext,c(1:d)+.4,colnames(cmat),pos=4,srt=-45,cex=textSize)

  colorLegend(c(d+1,d+2),c(0,d/2.2),ytick=c(slim[1],0,slim[2]),
              scale,cols=colseq,labside='left',endLabels=slim)
}
##################################################

zeroOneScale <- function(x) (x - min(x,na.rm=T))/(max(x,na.rm=T) - min(x,na.rm=T))

######################################################

mapContoursUSA <- function(xx,yy,z,zlevs=NULL,yname=NULL){ #must have lon,lat in main program

  require(maps)

  map('state',interior=F,col='grey',xlim=c(-97,-65))
  map('state',boundary=F,col='grey',add=T)

  values2contour(xx,yy,z,zlevs,lwd=1,col='black',add=T)
  if(!is.null(yname))title(yname)
 
}
######################################################

arrowField <- function(xy0,xy1=NULL,directionLong=NULL,angle=20,col=1){

  #xy0, xy1, 2 columns each
  #directionLong in radians, length

  if(is.null(xy1) & is.null(directionLong))stop('either xy1 or directionLong must be given')

  if(!is.null(directionLong)){
     xy1     <- xy0*0
     xy1[,1] <- xy0[,1] + directionLong[,2]*cos(directionLong[,1])
     xy1[,2] <- xy0[,2] + directionLong[,2]*sin(directionLong[,1])
     long    <- sqrt( (xy1[,1] - xy0[,1])^2 + (xy1[,2] - xy0[,2])^2)
  }

  lwide <- 2*zeroOneScale(long)

  arrows(xy1[,1],xy1[,2],xy0[,1],xy0[,2],length=.05*long,angle,col,lwd=lwide)
}

#################################################################3
mapArrowsUSA <- function(xx,yy,direction,long,yname=NULL){

  require(spatial)
  require(maps)

  map('state',interior=F,col='grey',xlim=c(-97,-65))
  map('state',boundary=F,col='grey',add=T)
      

    wf <- which(is.finite(direction) & is.finite(long))

    lo <- xx[wf]
    la <- yy[wf]

    arrowField(cbind(lo,la),directionLong=cbind(direction,long)[wf,],angle=20,col=1)
 
  if(!is.null(yname))title(yname)

  map('state',interior=F,add=T,lwd=4,col='white')
  map('state',interior=F,add=T)

 # map('state',boundary=F,col='white',lwd=3,add=T)
  map('state',boundary=F,col='grey',add=T)
 # title(mapname)

}
#######################################3


mapPointsUSA <- function(xx,yy,...,sym='values',cols='black',whiteOutDist=NULL,
                         fill='none',yname=NULL,scale=NULL,legendScale=NULL,KRIG=F){  

  #scale - symbol size
  #legendScale - range of values for color scale

  require(spatial)
  require(maps)

  map('state',interior=F,col='grey',xlim=c(-97,-65))
  map('state',boundary=F,col='grey',add=T)

  zz  <- list(...)

  if(length(zz) > 1){
    if(length(sym) == 1) sym  <- rep(sym,length(x))
    if(length(cols) == 1)cols <- rep(cols,length(x))
    if(length(fill) == 1)fill <- rep(fill,length(x))
  }

  if('ramp' %in% cols | 'ramp' %in% fill){

    ncol   <- 100
    colseq <- mapColors(ncol)

  }

  for(j in 1:length(zz)){

    xj <- zj <- zz[[j]]
    wf <- which(is.finite(xj))
    xj <- xj[wf]

    if(is.null(legendScale)) colScale <- seq(min(xj,na.rm=T),max(xj,na.rm=T),length=ncol)
    if(!is.null(legendScale))colScale <- seq(legendScale[1],legendScale[2],length=ncol)

    lo <- xx[wf]
    la <- yy[wf]

    if(sym[j] == 'ones'){
       w0 <- which(xj == 1)
       xj <- rep(.1,length(w0))
       ll <- la[w0]
       lo <- lo[w0]
    }

    fj <- cols[j]
    bj <- fill[j]
    if(bj == 'none')bj <- NA

    if(cols[j] == 'ramp' | fill[j] == 'ramp'){
      fj <- bj <- colseq[findInterval(xj,colScale)]
    }

    if(KRIG){
        tmp <- surf.gls(2,expcov,x=lo,y=la,xj,d=.7)
        ps  <- prmat(tmp,xl=min(lo),xu=max(lo),yl=min(la),yu=max(la),n=50)
        xy  <- as.matrix(expand.grid(ps$x,ps$y))
        zk  <- as.vector(ps$z)
        zk[zk < min(colScale)] <- colScale[1] + .1
        zk[zk > max(colScale)] <- max(colScale) - .1
        fj <- bj <- colseq[findInterval(zk,colScale)]
        symbols(xy[,1],xy[,2],squares=rep(.5,nrow(xy)),inches=F,fg=fj,bg=bj,add=T)
        zj <- ps$z

        if(!is.null(whiteOutDist)){   #whiteout areas further than this distance from points
            xy1  <- xy[,1]
            xy2  <- xy[,2]
            tmp  <- distmat(xx,yy,xy1,xy2)
            mind <- apply(tmp,1,min,na.rm=T)
            ww   <- which(mind > whiteOutDist)
            symbols(xy1[ww],xy2[ww],squares=rep(.5,length(ww)),inches=F,fg='white',bg='white',add=T)
        }
    }

    if(!KRIG){
      if(is.null(scale)) xj <- zz[[j]]
      if(!is.null(scale))xj <- rep(scale,length(zz[[j]]))
      symbols(xx,yy,circles=xj,inches=F,fg=fj,bg=bj,add=T)
      zj <- xj
    }

    if(is.null(legendScale))legendScale <- signif(range(zj,na.rm=T),1)
    
    tscale <- legendScale
    ttic   <- tscale
    if(tscale[1] < 0 & tscale[2] > 0)ttic   <- c(tscale[1],0,tscale[2])
    tscale <- seq(tscale[1],tscale[2],length=ncol)

    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    colorLegend(c(-70,-69),c(26,35),ytick=ttic,scale=tscale,colseq,labside='left')
  }

  if(!is.null(yname))title(yname)

  map('state',interior=F,add=T,lwd=4,col='white')
  map('state',interior=F,add=T)

 # map('state',boundary=F,col='white',lwd=3,add=T)
  map('state',boundary=F,col='grey',add=T)
 # title(mapname)

}



mapMask <- function(xx,yy,dx=NULL,dy=NULL,whiteOutDist=1,col='white'){    
  
  #mask parts of map with no data, grid density (dx,dy)
  #add to exising map
  #xy is an expand.grid of unmasked pixels
  
  require(mgcv)
  require(RANN)
  
  rx <- range(xx,na.rm=T)
  ry <- range(yy,na.rm=T)
  
  if(is.null(dx)){
    dx <- diff(rx)/50
    dy <- diff(rx)/50
  }
  
  
  xnew <- seq(rx[1],rx[2],by=dx)
  ynew <- seq(ry[1],ry[2],length=length(xnew))
  xxy   <- as.matrix( expand.grid(xnew,ynew) )
  nr   <- nrow(xxy)
  
  nx <- 0
  wd <- whiteOutDist*1.5
  
  while(nx < min(c(5000,nr/4))){
    tmp <- nn2(cbind(xx,yy),xxy,k = 1)$nn.dists
    xy  <- xxy[tmp > wd,]
    nx  <- nrow(xy)
    print(nx)
    wd  <- wd/1.3
  }
  
  symbols(xy[,1],xy[,2],squares=rep(dx,nx),inches=F,fg=col,bg=col,add=T)
  
}


checkDesign <- function(x,intName='intercept'){  # name of intercept column

  p <- ncol(x)
  
  if(ncol(x) < 3){
    message('no design check, x has 2 columns')
    return( list(VIF = 0, correlation = 1, rank = 2, p = 2) )
  }
    

  if(is.null(colnames(x))){
    colnames(x) <- paste('x',c(1:p),sep='_')
    xrange      <- apply(x,2,range)
    wi          <- which(xrange[1,] == 1 & xrange[2,] == 1)
    if(length(wi) > 0)colnames(x)[wi] <- 'intercept'
  }

  wi <- which(colnames(x) == 'intercept')

  xname <- colnames(x)

  VIF <- rep(NA,p)
  names(VIF) <- xname
  for(k in 1:p){

    if(xname[k] == intName)next

    notk <- xname[xname != xname[k]]

    ykk <- x[,xname[k]]
    xkk <- x[,notk]

    tkk <- summary(lm(ykk ~ xkk))$adj.r.squared
    VIF[k] <- 1/(1 - tkk)
  }

  VIF <- VIF[-wi] 

  corx <- cor(x[,-wi])
  rankx <- qr(x)$rank

  list(VIF = round(VIF,1), correlation = round(corx,2), rank = rankx, p = p)
}


############################################3
points2angle <- function(xy0,xy){             #(x,y) columns in xy from reference (x,y) in xy0, returns angle

  h  <- distmat(xy0[1],xy0[2],xy[,1],xy[,2])
  dx <- xy[,1] - xy0[1]
  dy <- xy[,2] - xy0[2]

  list(theta = atan2(dy,dx), dist = h)
}
#################################################
localSlope <- function(xk,yk,zk,direction=NULL,na.rm=F){ 
  
  #if not null, direction is theta
  # theta in radians, 0 to the right, pi to the left
  # positive theta is up, negative theta is down
  # theta points uphill, aspect points downhill
  # aspect is 0, 390 up (north), 90 right (east), 180 down (south)
  
  if(na.rm){
    ww <- which(is.finite(xk) & is.finite(yk) & is.finite(zk))
    xk <- xk[ww]
    yk <- yk[ww]
    zk <- zk[ww]
  }

  x0 <- mean(xk)
  y0 <- mean(yk)
  z0 <- mean(zk)

  xk <- xk - x0
  yk <- yk - y0
  zk <- zk - z0

  u1 <- u2 <- 1

  X <- cbind(xk,yk)

  if(length(xk) == 2){
    fx <- diff(zk)/diff(xk)
    fy <- diff(zk)/diff(yk)
  }

  if(length(xk) > 3){
    b <- invMat(crossprod(X),NEARPD=T)%*%crossprod(X,zk)
    fx <- b[1,]
    fy <- b[2,]
  }

  if(is.null(direction)){     #maximum derivative
    grade <- sqrt(fx^2 + fy^2)
    theta <- atan2(fy,fx)
  }

  if(!is.null(direction)){    #derivative in direction theta
     theta <- direction
     u1    <- cos(direction)
     u2    <- sin(direction)
     grade <- fx*u1 + fy*u2
  }
  
  deg <-  radians2degrees(pi/2 - theta) + 180
  if(deg > 360)deg <- deg - 360

  list(theta = theta, grade = grade, aspect = deg)
}

radians2degrees <- function(rads){
  
  rads/2/pi*360
}


mapColors <- function(ncol){

    colF   <- colorRampPalette(c('darkblue','blue','lightblue',
                                 'green','lightgreen',
                                 'yellow','orange','red','brown'))
    colF(ncol)
}


heatColors <- function(ncol){
  
  colF   <- colorRampPalette(c('darkblue','blue','brown','red'))
  colF(ncol)
}

darkColors <- function(ncol){
  
  colF   <- colorRampPalette(c('blue','green','red','brown','black'))
  colF(ncol)
}


degrees2radians <- function(degrees){
  
  degrees/360 * 2 * pi
}

#############################################

climGradient <- function(xx,yy,z1,z2,dthreshold=1,mapname=NULL,
                  climGrad=F,climDir=F,vegDir=F,ABS=T,DYDX=T,DTHETA=F){  

# directional gradient in z2 relative to maximum gradient in z1
# climGrad - map of climate gradient derivative
# climDir  - map of climate direction
# vegDir   - vegetation, rather than climate direction
# abs value- when there are multiple species

  require(spatial)

  ncol <- 100
  cols <- mapColors(ncol)

  if(!is.matrix(z2))z2 <- as.matrix(z2)
  r   <- ncol(z2)

  out <- matrix(NA,length(xx),5)
  colnames(out) <- c('theta','grad','thetaV','dtheta','dydx')

  distance <- distmat(xx,yy,xx,yy)

  for(i in 1:length(xx)){

    wi <- unique(which(distance[i,] < dthreshold))
    if(length(wi) < 4)wi <- order(distance[i,])[1:4]
    tmp <- localSlope(xx[wi],yy[wi],z1[wi])
    theta <- tmp$theta
    dx    <- tmp$grade

    wz <- which(z2[wi,] == 0)  #zeros do not occur at this location
    cdir <- theta

    gvec <- vslp <- rep(NA,r)

    for(j in 1:r){
      zj  <- z2[wi,j]
   #   wj  <- which(zj != 0)
   #   if(length(wj) < 4)next
      if(vegDir)cdir <- NULL
   #   tmp <- localSlope(xx[wi[wj]],yy[wi[wj]],zj[wj],direction=cdir)

      vslp[j] <- localSlope(xx[wi],yy[wi],zj)$theta
      tmp  <- localSlope(xx[wi],yy[wi],zj,direction=cdir)
      gvec[j] <- tmp$grade
    }

    dydx <- gvec/dx
    if(ABS)dydx <- abs(dydx)
    vslp    <- mean(vslp,na.rm=T)
    dtheta  <- atan2(sin(theta - vslp),cos(theta - vslp))
    out[i,] <- c(theta,dx,vslp,dtheta,mean(dydx,na.rm=T))
  }

  qGrad <- c(.01,.99)
  qdydx <- c(.01,.99)

 # tmap <- cos( out[,'theta'] )
  tmap <- out[,'theta']
  
  gmap <- zeroOneScale( trim( out[,'grad'],qGrad) )
  
  dtmp <- trim( out[,'dydx'],qdydx)
 # dtmp <- log10(dtmp)
  dtmp[dtmp == -Inf] <- min(dtmp[dtmp != -Inf],na.rm=T)

  dmap <- zeroOneScale( dtmp )
  dmap[out[,'dydx'] == 0] <- min(dmap,na.rm=T)

  mmap <- zeroOneScale( abs(out[,'dtheta']) )


 # dmap <- zeroOneScale( trim( out[,'dydx'],qdydx) )

  mscale <- signif( c(0,pi),2)
  tscale <- signif(c(-pi,pi),2)
  gscale <- signif(quantile(out[,'grad'],qGrad,na.rm=T),1)
  dscale <- signif(quantile(out[,'dydx'],qdydx,na.rm=T),1)

  ttic   <- c(tscale[1],0,tscale[2])
  gtic   <- gscale
  dtic   <- dscale
  mtic   <- mscale

  if(gtic[1] < 0 & gtic[2] > 0)gtic <- c(gtic[1],0,gtic[2])
  if(dtic[1] < 0 & dtic[2] > 0)dtic <- c(dtic[1],0,dtic[2])

  tscale <- seq(tscale[1],tscale[2],length=ncol)
  gscale <- seq(gscale[1],gscale[2],length=ncol)
  dscale <- seq(dscale[1],dscale[2],length=ncol)
  mscale <- seq(mscale[1],mscale[2],length=ncol)

  if(climDir){
    mapPointsUSA(xx,yy,tmap,sym='values',cols='ramp',whiteOutDist=1.4,
                         fill='none',yname=paste(mapname,'theta',sep=' '),KRIG=T)
    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)

    colorLegend(c(-70,-69),c(26,35),ytick=ttic,scale=tscale,cols,labside='left')
  }
  if(climGrad){
    mapPointsUSA(xx,yy,gmap,sym='values',cols='ramp',whiteOutDist=1.4,
                         fill='none',yname=paste(mapname,'Gradient',sep=' '),KRIG=T)
    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    colorLegend(c(-70,-69),c(26,35),ytick=gtic,scale=gscale,cols,labside='left')
  }

  if(DYDX){
    mapnew <- paste(mapname,'dxdy',sep=' ')
 
    mapPointsUSA(xx,yy,dmap,sym='values',cols='ramp',whiteOutDist=1.4,
                         fill='none',yname=mapnew,KRIG=T)
    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    colorLegend(c(-70,-69),c(26,35),ytick=dtic,scale=dscale,cols,labside='left')
  }
  if(DTHETA){
    mapnew <- paste(mapname,'del theta',sep=' ')
 
    mapPointsUSA(xx,yy,-(mmap-1),sym='values',cols='ramp',whiteOutDist=1.4,
                         fill='none',yname=mapnew,KRIG=T)
    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    colorLegend(c(-70,-69),c(26,35),ytick=rev(-mtic),scale=mscale,cols,labside='left')
  }

  invisible(out)
}

################################################3333

associationPlot <- function(z,groupIndex=c(1:nrow(z))){       

  #z - n rows by c columns give abundance for each species
  #groupIndex - length n

  allnames <- sort(unique(groupIndex))
  ns       <- length(allnames)
  nf       <- ncol(z)

  ncol <- 50
  scale <- seq(-1,1,length.out=ncol)
  
  freqY <- numeric(0)

  nr <- round(ns/2,0) + 1
  par(mfrow=c(nr,2),bty='n',mar=c(1,1,3,1))

  for(i in 1:ns){

      yi <- z[groupIndex == allnames[i],]

      yc <- cor(yi)

      corPlot(yc,slim=c(-1,1),diag0=T)
      
      title(allnames[i])
 
  }

  dev.print(device=postscript,file='associationPlot.ps',width=6,horizontal=F)
}
############################################################
colorLegend <- function(x,y,ytick=NULL,scale,cols,labside='right',
                        text.col=NULL,
                        bg=NULL,endLabels=NULL){  
  # x = (x1,x2), y = (y1,y2)
  # bg is color of border

  nn <- length(scale)
  ys <- seq(y[1],y[2],length=nn)

  for(j in 1:(length(scale)-1)){
    rect(x[1],ys[j],x[2],ys[j+1],col=cols[j],border=NA)
  }
  if(!is.null(bg))rect(x[1],y[1],x[2],y[2],border=bg,lwd=3)
  

#    dx <- diff(x)
#    yy <- ys[findInterval(ytick,scale)]
#    lx <- c(x[2],x[2]+dx/3)
#    tx <- x[2]+dx/2.5
 #   tp <- 4
 #   if(labside == 'left'){
 #        lx <- c(x[1]-dx/3,x[1])
 #        tx <- x[1]-dx/2.5
 #        tp <- 2
 #   }
 #   for(j in 1:length(ytick)){
 #      lines(lx,c(yy[j],yy[j]))
 #      if(!is.null(ytick))text(tx,yy[j],ytick[j],pos=tp)
 #   }
  

  if(!is.null(endLabels)){ 
    cx <- cols[c(1,nn)]
    if(!is.null(text.col))cx <- text.col
    if(labside == 'right')text(c(x[2],x[2]),y,endLabels,pos=1,col=cx)
    if(labside == 'left')text(c(x[1],x[1]),y,endLabels,pos=2,col=cx)
  }
}
####################################################
bmultiProp <- function(r,k,b = matrix(0,k,r-1),pBVar=diag(.1,k*(r-1)),
                       loB=NULL,hiB=NULL){  
	
    b <- as.vector(b)
    cvec <- myrmvnorm(1,t(b),pBVar)
    cc   <- matrix(cvec,k,r-1)
    
    if(!is.null(loB)){                     
    	cvec <- tnorm.mvt(b,b,pBVar,loB,hiB,times=1)
      cc   <- matrix(cvec,k,r-1)
    }

  list(cc = cc, cvec = cvec)
}

#################################################################
simX <- function(n,loX,hiX){                #generate design matrix

  k <- length(loX)
  x <- matrix(1,n,k)
  for(j in 1:k)x[,j] <- runif(n,loX[j],hiX[j])
  x
}
####################################################
simY <- function(x,b,LIKE,r = 1,size=rep(1,nrow(x)),sigma = 0,Effort = 1){     #simulate response

  u <- x%*%b

  if(LIKE == 'norm')return( rnorm(length(u),u,sqrt(sigma)) )
  if(LIKE == 'pois'){
  	 u <- exp(u)*Effort
  	 return( rpois(nrow(x),u) )
  } 
  if(LIKE == 'binom'){
  	 u <- invlogit(u)
  	 return( rbinom(nrow(x),1,u) )
  }
  if(LIKE == 'multinom'){
     zs <- rowSums(exp(u))
     z1 <- 1/(1 + zs)
     zm <- exp(u)/ (1 + zs)
     u  <- cbind(zm,z1)
     return( myrmultinom(size,u) )
  }
  if(LIKE == 'mvnorm'){
    u <- myrmvnorm(n,u,sigma)
    return( u )
  }
  if(LIKE == 'mvnorm-multinom'){
    u <- myrmvnorm(n,u,sigma)
    zs <- rowSums(exp(u))
    z1 <- 1/(1 + zs)
    zm <- exp(u)/ (1 + zs)
    u  <- cbind(zm,z1)
    return( myrmultinom(size,u) )
  }
  numeric(0)
}
#########################################
simMVData <- function(LIKE,k,r){
	
   xnames <- c('intercept',paste('x',c(1:(k-1)),sep='-') )
   p      <- length(xnames)
   x      <- matrix(runif(n*k,-4,4),n,k)
   x[,1]  <- 1
   colnames(x) <- xnames
   
   b           <- matrix(runif(k*ns,-2,1),k,r)
   specnames   <- paste('S',c(1:r),sep='-')
   colnames(b) <- specnames
   rownames(b) <- xnames
   sigma <- diag(runif(r,0,.1/r),r)
   rownames(sigma) <- specnames
   colnames(sigma) <- specnames

   mu <- x %*% b
   y     <- myrmvnorm(n,mu,sigma)
   sigma <- crossprod(y - x%*%b)/(n-k)
   y     <- myrmvnorm(n,mu,sigma)
   colnames(y) <- specnames
   
   if(LIKE == 'mvnorm-Pois')    z <- matrix(rpois(length(y),exp(y)),n,r,byrow=F)
   if(LIKE == 'mvnorm-multinom'){
   	
   	  p1 <- matrix(rbeta(n*r,.2,.2),n,r)
   	  p1[cbind(c(1:n),sample(c(1:r),n,replace=T))] <- 10
   	  p1 <- p1/matrix(rowSums(p1),n,r)
   	  zz <- myrmultinom(10,p1)
   	  
   	  b  <- bInitMNomLogit(x,zz,size)
   	
   	  sigma <- sigma[1:(r-1),1:(r-1)]
     y  <- myrmvnorm(n, x %*% b,sigma)
   	  zs <- rowSums(exp(y))
     z1 <- 1/(1 + zs)
     zm <- exp(y)/ (1 + zs)
     y  <- cbind(zm,z1)
   	  z  <- myrmultinom(size,y)
   	  plot(y,jitter(z))
   	}
   
   list(x = x, b = b, y = y, z = z, s = sigma, specnames = specnames)
}

gibbsLoop <- function(LIKE,ng,x,y,b,priorB=b*0,priorVB=diag(1000,length(b)),
                      loB=NULL,hiB=NULL,
                      sigma = 0,z = numeric(0),burnin=1){
                      	
  tiny <- 1e-6
  
  k  <- ncol(x)
  n  <- nrow(x)
  kx <- length(b)
  if(!is.matrix(b))b <- as.matrix(b)
  
  priorIVB <- invMat(priorVB)
  
  bgibbs <- matrix(NA,ng,kx)
  colnames(bgibbs) <- paste('b',c(1:kx),sep='-')
  sgibbs <- matrix(NA,ng,max(1,length(sigma)))
  
  r <- 1                       #responses
  if(is.matrix(y))r <- ncol(y)
  
  if(sigma[1] == 0)sgibbs <- numeric(0) #variances
  
  sg <- sigma
  
  if(is.matrix(sigma)){                  #Wishart prior
    sg        <- prior.W
    colnames(sgibbs) <- outer(rownames(sigma),rownames(sigma),paste,sep='_') 
    colnames(bgibbs) <- as.vector(outer(colnames(x),colnames(b),paste,sep='_'))
  }
  
  pBVar <- invMat(crossprod(x))

  pred <- pred2 <- rep(0,nrow(x)*r)

  bg <- b
  yg <- y
  y1 <- yg
  
  if(LIKE == 'pois') bg <- pBVar%*%crossprod(x,log(y + .1))
  if(LIKE == 'binom'){
  	yl <- y
  	yl[yl == 0] <- .1
  	yl[yl == 1] <- .9
  	bg <- pBVar%*%crossprod(x,logit(yl))
  }
  
  if(LIKE == "mvnorm-multinom"){
  	y1   <- prob2Logit(yg)
  }
  if(LIKE == 'multinom'){
  	size <- rowSums(y)
  	pBVar   <- diag(.0001,k*(r-1))
  	priorB  <- b*0
  	priorVB <- diag(100,k*(r-1))
  	priorIVB <- solve(priorVB)
  }
  
  dev <- 0
  np  <- 0

  for(g in 1:ng){

    bg <- bUpdateGibbs(LIKE,x,y1,bg,priorB,priorIVB,loB,hiB,sigma=sg,pBVar)
    
    if(sigma[1] > 0 & length(sigma) == 1){
    	 sg <- updateSigma(y1,x%*%bg)
    }
    
    if(length(grep('mvnorm',LIKE)) > 0){
    	 sinv <- wishsamp(x,y1,bg)
    	 sg   <- solve(sinv)
    }
    
    if(LIKE == "mvnorm-multinom"){
    	y1 <- ysampMvnormMultinom(x,y1,z,bg,sg)$y
    	yg <- logit2Prob(y1)
    }
    
    dev <- dev + sum(deviance(y1,x,bg,sg,LIKE))

    if(g > burnin){
      py    <- as.vector(simY(x,bg,LIKE,r,size,sg))
      pred  <- pred + py
      pred2 <- pred2 + py^2
      np    <- np + 1
    }

    bgibbs[g,] <- bg
    if(length(sgibbs) > 0)sgibbs[g,] <- sg
    
    if(g %in% c(100,200,500,1000)){
    	pBVar <- .5*cov(bgibbs[20:g,]) + diag(tiny,ncol(bgibbs))
    }
  }

  ymean <- pred/np
  yse   <- sqrt(pred2/np - ymean^2)
  
  if(r > 1){
  	ymean <- matrix(ymean,n,r)
  	yse   <- matrix(yse,n,r)
  }
  bmean <- colMeans(bgibbs)
  bmean <- matrix(bmean,nrow(bg),ncol(bg))
  smean <- numeric(0)
  if(length(sg) > 1) smean <- matrix(colMeans(sgibbs),nrow(sg),ncol(sg))
  if(length(sg) == 1)smean <- mean(sgibbs)
  
  meanDev <- dev/ng
  
  pd  <- meanDev - sum(deviance(y1,x,bmean,smean,LIKE))
  dic <- 2*pd + meanDev

  list(bgibbs = bgibbs,sgibbs = sgibbs, ymean = ymean, yse = yse, dic = dic)
}
####################################################
simY <- function(x,b,LIKE,r = 1,size=rep(1,nrow(x)),sigma = 0,Effort = 1){     #simulate response

  u <- x%*%b

  if(LIKE == 'norm')return( rnorm(length(u),u,sqrt(sigma)) )
  if(LIKE == 'pois'){
  	 u <- exp(u)*Effort
  	 return( rpois(nrow(x),u) )
  } 
  if(LIKE == 'binom'){
  	 u <- invlogit(u)
  	 return( rbinom(nrow(x),1,u) )
  }
  if(LIKE == 'multinom'){
     zs <- apply(exp(u),1,sum)
     z1 <- 1/(1 + zs)
     zm <- exp(u)/ (1 + zs)
     u  <- cbind(zm,z1)
     return( myrmultinom(size,u) )
  }
  if(LIKE == 'mvnorm'){
    u <- myrmvnorm(n,u,sigma)
    return( u )
  }
  if(LIKE == 'mvnorm-multinom'){
    u <- myrmvnorm(n,u,sigma)
    zs <- apply(exp(u),1,sum)
    z1 <- 1/(1 + zs)
    zm <- exp(u)/ (1 + zs)
    u  <- cbind(zm,z1)
    return( myrmultinom(size,u) )
  }
  numeric(0)
}
##########################################################

sensIntercept <- function(bgibbs,xnames,ynames){     

  #sensitivity coeffs for multinomial regression
  #bgibbs - multinomial logit par chains
  #xnames - length-k vector of names for input x 
  #ynames - length-(r-1) vector of names for response y

  r <- length(ynames)
  k <- length(xnames)

  fullnames <- as.vector(outer(xnames,ynames[-r],paste,sep='-')) #columns of bgibbs

  sbgibbs <- numeric(0)

  speccol <- matrix(unlist(strsplit(fullnames,'-')),ncol=2,byrow=T)[,2]

  for(h in 1:(r-1)){
    wh1 <- which(speccol == ynames[h])
    wh2 <- wh1[ !wh1 %in% c(grep('health',fullnames),grep('gap',fullnames)) ]
    wh0 <- wh1[ !wh1 %in% wh2 ]
    fmean  <- rowMeans(bgibbs[,wh2])
    int    <- sweep(bgibbs[,wh2],1,fmean) 
    int    <- cbind(fmean,bgibbs[,wh0],int)
    colnames(int)[1] <- paste('int',ynames[h],sep='-')
    sbgibbs <- cbind(sbgibbs,int)
  }
  sbgibbs
}
#######################################################
multinomLike <- function(y,x,b){  #log likelihood multinomial logit
	
  tiny <- 1e-20
  huge <- 1 - tiny

     z <- x%*%b
     zs    <- rowSums(exp(z))
     z1    <- 1/(1 + zs)
     zm    <- exp(z)/ (1 + zs)
     z2  <- cbind(zm,z1)
     z2[z2 < tiny] <- tiny
     z2[z2 > huge] <- huge
     list(like = y*log(z2), theta = z2)
}
####################################################
bmultiProp <- function(r,k,b = matrix(0,k,r-1),pBVar=diag(.1,k*(r-1)),loB=NULL,hiB=NULL){  
	
    bvec <- as.vector(b)
    cvec <- myrmvnorm(1,t(bvec),pBVar)
    cc   <- matrix(cvec,nrow(b),ncol(b))
    
    if(!is.null(loB)){                         #if lob and hib available, use tnorm.mvt
    	cvec <- tnorm.mvt(bvec,bvec,pBVar,loB,hiB,times=4)
      cc   <- matrix(cvec,nrow(b),ncol(b))
    }

  list(cc = cc, cvec = cvec)
}

####################################################
bUpdateMNom <- function(x,y,b,priorB,priorIVB,
              pBVar=diag(.1,nrow(b)*(ncol(b)-1)),loB=NULL,hiB=NULL){

  bvec <- as.vector(b)
  priorVB <- solve(priorIVB)

  tmp  <- bmultiProp(ncol(b)+1,nrow(b),b,pBVar,loB,hiB)
  cc   <- tmp$cc
  cvec <- tmp$cvec
  
  pnow <- multinomLike(y,x,b)$like
  pnew <- multinomLike(y,x,cc)$like

  pnow <- sum(pnow) + mydmvnorm(bvec,as.vector(priorB),priorVB,log=T)
  pnew <- sum(pnew) + mydmvnorm(cvec,as.vector(priorB),priorVB,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)b <- cc
  b
}


updateWishart <- function(yy,predy,priorS,priorSdf,INVERSEONLY=F){

  n  <- nrow(yy)
  ss <- crossprod(yy - predy) + priorS*priorSdf
  df <- n + priorSdf

  sinv <- rwish(df,invMat(ss,NEARPD=T))
  if(INVERSEONLY)return(sinv)
  sig  <- invMat(sinv)
  list(sigma = sig, sinv = sinv)
}

getScWishMat <- function(SS,df,delta,priorO){
  
  di    <- diag(1/diag(delta))
  tt    <- priorO + di%*%SS%*%di
  IO    <- rwish(df,solve(tt)) 
  O     <- solve(IO)
  sigma <- delta%*%O%*%delta
  
  list(omega = O, omegaInv = IO, sigma = sigma)
}
  
updateScInvWish <- function(yz,predy,delta=diag(1,ncol(yz)),
                            priorO=diag(1,ncol(yz)),priorOdf=(1+ncol(yz)),
                            priorDmu=rep(1,ncol(yz)),priorDvar=rep(10,ncol(yz)),
                            varLo=rep(1e-8,ncol(yz)),varHi=rep(10000,ncol(yz)) ){ 
  
  #varBound is permissible range of variances
  
  ww <- which(diag(delta) < 1e-5)
  if(length(ww) > 0)delta <- delta*10
  
#  ww <- which(diag(delta) > varBound[2])
#  if(length(ww) > 0)diag(delta)[ww] <- varBound[2]
  
  
  k  <- ncol(yz)
  n  <- nrow(yz)
  SS <- crossprod(yz - predy)
  df <- n + priorOdf
  
 # del <- diag(exp(priorDmu))  #prior
  
  tmp <- getScWishMat(SS,df,delta,priorO)
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


updateInvWish <- function(yy,predy,v=2,sigma,aa,A=rep(1000,ncol(yy))){   
  #Huang and Wand, Bayesian Analy (2013)
  
  k  <- ncol(yy)
  nn <- nrow(yy)
  SS <- crossprod(yy - predy)
  df <- nn + v + k
  
  sinv <- solve(sigma)
  
  aa <- 1/rgamma(k,(v + k)/2,v*diag(sinv) + 1/A^2)
  
  sinv  <- rwish(df,solve(SS + 2*v*diag(aa))) 
  smat  <- solve(sinv)
  
  list(sigma = smat, sinv = sinv, aa = aa)
}


################################################
bUpdateNorm <- function(xx,yy,b,
                priorB=matrix(ncol(xx)*0,ncol(xx)),priorIVB=diag(1/100,ncol(xx)),
                loB=NULL,hiB=NULL,sigma){

  V <- invMat( crossprod(xx)/sigma + priorIVB )
  v <- crossprod(xx,yy)/sigma + priorIVB%*%priorB
  if( is.null(loB) & is.null(hiB) )return( t( myrmvnorm(1,t(V%*%v),V) ) )
  
  tnorm.mvt(V%*%v,V%*%v,V,loB,hiB)
}
 

###########################################33
bUpdateGibbs <- function(LIKE,x,y,b,priorB,priorIVB,
                         loB=NULL,hiB=NULL,sigma = 0,pBVar=0){

  if(LIKE == 'norm')    return( bUpdateNorm(x,y,b,priorB,priorIVB,loB,hiB,sigma) )
  if(LIKE == 'multinom')return( bUpdateMNom(x,y,b,priorB,priorIVB,pBVar,loB,hiB) )
  if(LIKE %in% c('mvnorm','mvnorm-multinom'))return( bUpdateMVNorm(x,y,b,sigma) )

  b <- matrix(b,length(b),1)
  if( is.null(loB)) c <- t(myrmvnorm(1,t(b),pBVar))   #proposal
  if(!is.null(loB)) c <- tnorm.mvt(b,b,pBVar,loB,hiB,times=1)

  znow <- x%*%b
  znew <- x%*%c

  if(LIKE == 'pois'){
     pnow <- dpois(y,exp(znow),log=T)
     pnew <- dpois(y,exp(znew),log=T)
  }
  if(LIKE == 'binom'){
     pnow <- dbinom(y,1,invlogit(znow),log=T)
     pnew <- dbinom(y,1,invlogit(znew),log=T)
  }

  priorVB <- solve(priorIVB)
  
  pnow <- sum(pnow) + mydmvnorm(t(b),priorB,priorVB,log=T)
  pnew <- sum(pnew) + mydmvnorm(t(c),priorB,priorVB,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)b <- c
  b
}
####################################################
updateSigma <- function(y,mu,s1=1,s2=1){

  u1 <- s1 + length(y)/2
  u2 <- s2 + .5*sum( (y - mu)^2 )
  1/rgamma(1,u1,u2)
}

######################################
wishsamp <- function(x,yy,b){   #sample from Inv Wishart

   r    <- ncol(b)
   scp  <- crossprod((yy - x %*% b))
   vmat <- solve(scp + prior.W*prior.WDF)
   v2   <- ns + prior.WDF
   stmp <- myrmvnorm(v2,matrix(0,v2,r),vmat)
   crossprod(stmp)
}
###################################
bUpdateMVNorm <- function(x,yy,b,sigma,alpha=0,lo=NULL,hi=NULL,X=NULL){  # update b's for mvnorm
	
	require(mvtnorm)

  r      <- ncol(yy)
  k      <- ncol(x)

  if(is.null(X))X <- crossprod(x)
  
  bigv   <- solve(X)      #multivariate scale invariate prior (Minka)
  smallv <- crossprod(x,yy)
                 
  mu     <- bigv%*%smallv/(alpha + 1)
  vaa    <- kronecker(sigma,bigv/(alpha + 1))   
  if(is.null(lo)) bg <- matrix( rmvnorm(1,as.vector(mu),sigma=vaa) ,k,r,byrow=F)
  if(!is.null(lo))bg <- matrix(tnorm.mvt(as.vector(mu),as.vector(mu),vaa,lo,hi,times=1),k,r)
  
  colnames(bg) <- colnames(yy)
  rownames(bg) <- colnames(x)
  bg
}

#################################

allomConvert <- function(xx,specCodes,allomFile,codeColumn='species',
                         int='htInt',slope='htSlope',
                         defaultSpec=NULL, invert=F){
  
  # x          - vector or matrix of values to convert (e.g., diam)
  # specCodes  - vector equal to nrow(x) 
  # codeColumn - column in allomFile containing specCodes
  # defaultSpec- value to use if specCodes missing from allomFile (e.g., 'other')
  #              there is a row in allomFile with defaultValue in codeColumn
  
  if(!is.matrix(xx))xx <- matrix(xx,ncol=1)
  if(!is.matrix(specCodes))specCodes <- matrix(specCodes,ncol=1)
  
  coeff  <- read.table(allomFile,header=T)
  smatch <- match(specCodes,coeff[,codeColumn])
  
  ws <- which(coeff[,'species'] == 'other')
  
  wmiss <- which(is.na(smatch))
  if(length(wmiss) > 0){
    ss <- sort( unique( specCodes[wmiss] ) )
    sc <- 'missing from allometry:'
    for(k in 1:length(ss))sc <- paste(sc,ss[k],sep=', ')
    message( sc )
    smatch[wmiss] <- which(coeff[,'species'] == 'other')
  }
  
  c1 <- coeff[smatch,int]
  c2 <- coeff[smatch,slope]
  c1[is.na(c1)] <- coeff[ws,int]
  c2[is.na(c2)] <- coeff[ws,slope]
  
  tmp <- 10^( c1 + c2*log10(xx) )
  if(invert)tmp <- 10^( (log10(xx) - c1 )/ c2 )
  tmp
}



diam2mass <- function(d,alloB,alloL){  #allo has c(int,slope) on log10 scale

  b <- 10^alloB[,1] *d^alloB[,2]
  l <- 10^alloL[,1] *d^alloL[,2]
  l[is.finite(l) & l > b] <- b[is.finite(l) & l > b]*.9
  cbind(b,l)

}

stemMass2diam <- function(stem,leaf,alloB){

  10^( (log10(stem + leaf) - alloB[,1])/alloB[,2] )
} 

mass2diam <- function(mass,allo){

  10^((log10(mass) - allo[,1])/allo[,2])
}


#################################
ysampMvnormMultinom <- function(x,y,z=NULL,b,sigma){  

  #sample y on MVlogit scale
  #z - counts, same dimensions as y

  r  <- ncol(y)
  k  <- nrow(b)
  n  <- nrow(y)
  
  propy <- matrix(rnorm(length(y),y,.01),n,r,byrow=F)
  
  pnow <- rep(0,n)
  pnew <- rep(0,n)
  
  for(i in 1:n){
    pnow[i] <- mydmvnorm(y[i,],(x %*% b)[i,], sigma,log=T)
    pnew[i] <- mydmvnorm(propy[i,],(x %*% b)[i,], sigma,log=T)
   }

  zs    <- rowSums(exp(y))
  z1    <- 1/(1 + zs)
  zm    <- exp(y)/ (1 + zs)
  znow  <- cbind(zm,z1)

  zs    <- rowSums(exp(propy))
  z1    <- 1/(1 + zs)
  zm    <- exp(propy)/ (1 + zs)
  znew  <- cbind(zm,z1)

  if(!is.null(z)){
    pnow <- pnow + z*log(znow)
    pnew <- pnew + z*log(znew)
  }

  pnow <- rowSums(pnow)
  pnew <- rowSums(pnew)

  a  <- exp(pnew - pnow)
  zz <- runif(length(a),0,1)
  y[zz < a,] <- propy[zz < a,]
  accept <- length(zz[zz < a])

  list(y = y, a = accept)
}

####################################################
getSens <- function(xv,bchain){

  #names(xv) are repeated for each class in bgibbs
  #ncol(bgibbs) = (r - 1)*length(xv)

  kk <- length(xv)

  nsim <- 2000
  sens <- matrix(NA,nsim,ncol(bchain))
  colnames(sens) <- colnames(bchain)

  for(j in 1:nsim){

    gj  <- sample(ng,1)
    bgj <- matrix(bchain[gj,],kk,r-1)
    sens[j,] <- as.vector( multinomSens(bgj,xvals=xv) )
  }
  sens
}

###############################################
plotSens <- function(sens,ytic,xnames,ynames){  #plot multinomial sensitivities

# sens - output from getSens

  colF <- colorRampPalette(c('darkblue','blue','green','yellow','orange','red'))

  if('int' %in% xnames)xnames <- xnames[xnames != 'int']

  r <- length(ynames)
  k <- length(xnames)

  fullnames <- as.vector(outer(xnames,ynames[-r],paste,sep='-')) #columns of bgibbs

  tmp <- processPars(sens)$summary 

  par(mfrow=c(1,1),bty='n')

  xtic <- c(0:k)+.5

  plotSetup(xtic,ytic,xvals=rep(' ',length(xtic)),ylabel='Sensitivity')

  text(c(1:k),ytic[1],xnames,srt=90,cex=1.2)

  specindex <- matrix(unlist(strsplit(fullnames,'-')),ncol=2,byrow=T)[,2]

  jseq <- seq(-1,1,length=r)*.3
  cols <- colF(r-1)

  for(j in 1:(r-1)){
	
	tj <- tmp[grep(ynames[j],rownames(tmp)),]
        w0 <- grep('int',rownames(tj))
        if(length(w0) > 0)tj <- tj[-w0,]
	if(length(tj) == 0)next
	for(jk in 1:nrow(tj)){
		lines( c(jk,jk)+jseq[j],tj[jk,2:3],lwd=4,col=cols[j])
		points( jk+jseq[j],tj[jk,1],col=cols[j],pch=3)
	}
  }	

  legend('topleft',ynames[-r],text.col=cols[-r],cex=.9,bty='n')
  dev.print(device=postscript,file='multinomModel.ps',width=6,horizontal=F)

}

#######################################################
multinomSens <- function(bb,x=NULL,xvals=NULL){  #sensitivity coefficients multinomial logit
	
  tiny <- 1e-20
  huge <- 1 - tiny
  kk   <- nrow(bb)
  rr   <- ncol(bb)

  if(is.null(xvals))xvals <- matrix(colMeans(x),1)

  z     <- xvals%*%bb
  zs    <- rowSums(exp(z))
  zm    <- exp(z)/ (1 + zs)
  theta <- matrix(zm,kk,rr,byrow=T)
  bb*theta*(1 - theta)
}
####################################################

prob2Logit <- function(y){     #fractions to multivar logit      
  
  log(y[,-r]/(1 - rowSums(y[,-r])))
  
}


par2sens <- function(pchain,vname,xrangej,reverseSign){
  
  #parChain is a MCMC matrix, colnames are 
  
  wc <- grep(vname,colnames(pchain))   #all
  wi <- grep('X',colnames(pchain))     #interactions
  wm <- wc[!wc %in% wi]                #main effect
  xi <- matrix( unlist( strsplit(colnames(pchain)[wi],'X') ) , ncol=2,byrow=T)[,1]
  xx <- pchain[,wm]
  if(vname %in% reverseSign)xx <- -xx
  for(ii in 1:length(xi)){
    ix <- pchain[,xi[ii]] * .5
    if(xi[ii] %in% reverseSign)ix <- -ix
    xx <- xx + ix
  }
  
  #incorporate range
  xx <- xx/(xrangej[2,vname] - xrangej[1,vname])
  xx
}
##########################################


plotDensColumns <- function(var2plot=NULL,chain2plot,allChains=chainList,chainLength=1,
                            vnames=xnames,MAINEFFECT=T,xtic=NULL,labs=NULL,
                            htFraction=.5,textSize=1,ORD=NULL,COLCODE=NULL,
                            reverseSign=character(0),omitSpecs=character(0),FILL=F){
  
  #if !MAINEFFECT then interactions added to main effect
  
  quadCols <- numeric(0)

  wvar2 <- character(0)
  wvar <- var2plot
  wd   <- grep( paste(var2plot,'-',sep='') ,vnames)  # a climate variable with months
  if(length(wd) > 0)wvar <- vnames[wd[1]]          # if so, new name
  
  reverseM <- F
  if(wvar %in% reverseSign)reverseM <- T
  hasQuad <- hasInt <- F
  isInt <- F

  xx <- grep('X',vnames)           #some are interactions
  x2 <- grep('2',vnames)           #some are quadratic
  xi <- grep('X',var2plot)         #variable is an interaction

  wint <- grep(var2plot,vnames[xx])  #variable to plot has interaction terms
  w2   <- grep(var2plot,vnames[x2])  #variable to plot has quadratic terms

  if(length(wint) > 0)hasInt  <- T
  if(length(w2) > 0)  hasQuad <- T
  if(length(xi) > 0)  isInt   <- T
  if(isInt)hasInt <- hasQuad <- F

  
  xtAll <- ytAll <- numeric(0)
  ymax  <- numeric(0)
  xrange <- numeric(0)
  
  modeBySpec <- numeric(0)
  allSpecs   <- character(0)
  
  for(jj in 1:chainLength){
    
    ww <- grep('Chain',names(allChains))
    if(length(ww) > 0){
      chainj <- allChains
    } else {
      chainj <- allChains[[jj]]
    }
    
    namej  <- names(chainj)
    namej  <- unlist( strsplit(namej,'Chain') )
    
    wc     <- which(namej == chain2plot)
    mat    <- chainj[[wc]]
    cnames <- colnames(mat)
    cm     <- cnames[grep(wvar,cnames)]   #names for main effects
    
    if(length(cm) == 0){
      xtAll <- append(xtAll,list(numeric(0)))
      ytAll <- append(ytAll,list(numeric(0)))
      modeBySpec <- append(modeBySpec,list(numeric(0)))
      next
    }
    
    if(hasQuad){
      wvar2 <- paste(var2plot,'2',sep='')
      ww    <- sort( grep(wvar,cnames) )
      cm    <- cnames[ww]
    }
    mat    <- mat[,cm]
    cm     <- colnames(mat)
    
    if(length(omitSpecs) > 0){
      woc <- numeric(0)
      for(kk in 1:length(omitSpecs)){
        woc <- c(woc,grep(omitSpecs[kk],colnames(mat)))
      }
      mat <- mat[,-woc]
      cm  <- colnames(mat)
    }
    
    intCols  <- sort( unique( c( grep( paste(wvar,'X',sep=''),cm), grep( paste('X',wvar,sep=''),cm ) ) ) )
    if(length(wvar2) > 0)quadCols <- grep( wvar2,cm)
    mainCols <- grep(wvar,cm)
    mainCols <- mainCols[!mainCols %in% c(intCols,quadCols)]
    matM    <- mat[,mainCols]  # only main effects
    
    # if(reverseM)matM <- -matM
    
    if(isInt){
      xxn <- unlist(strsplit(var2plot,'X')) 
    }
    
    if(!MAINEFFECT & hasInt){                                     # full effect: add interactions
      wi <- xx[wint]
      for(j in wi){
        xi <- unlist( strsplit(vnames[j],'X') )
        viname <- xi[xi != wvar] 
        
        ii <- 1
        if(viname %in% reverseSign)ii <- -1 
        icol  <- grep(vnames[j],cm)
        if(length(icol) > 0)matM <- matM + mat[,icol]*.5*ii      # assume interaction term is at mean value (.5)
      }
    }
    if(hasQuad){
      wi <- x2[w2]
      for(j in wi){
        qcol  <- grep(vnames[j],cm)
        matM <- matM + mat[,qcol]*2*.5
      }
    }
    
    mat <- matM
    cm <- colnames(mat)
    
  #  if(is.null(labs))
    labs  <- matrix( unlist( strsplit(cm,'_') ),ncol=2,byrow=2)[,2]
    nj    <- length(labs)
    
    rr <- range(matM,na.rm=T)
    if(!is.finite(rr[1]) & !is.finite(rr[2])){
      xtAll <- append(xtAll,list(numeric(0)))
      ytAll <- append(ytAll,list(numeric(0)))
      modeBySpec <- append(modeBySpec,list( numeric(0) ))
      next
    }
    
    tmp    <- chains2density(matM,labs=labs,reverseM=reverseM)   
    xt     <- tmp$x
    yt     <- tmp$y
    ymax   <- max( c(ymax,quantile(yt,.8)) ) 
    xrange <- range( c(xrange,tmp$xrange))
    
    wmax <-  apply(yt,1,which.max )
    xmax <- xt[ cbind(1:nrow(xt),wmax) ]
    names(xmax) <- names(wmax)
    modeBySpec <- append(modeBySpec,list( xmax ))
    
  #  rownames(xt) <- rownames(yt) <- paste('r',jj,rownames(xt),sep='-')
    xtAll <- append(xtAll,list(xt))
    ytAll <- append(ytAll,list(yt))
    
    allSpecs <- sort(unique(c(allSpecs,names(xmax))))
  }
  
  
  minMaxSpec <- matrix(NA,length(allSpecs),2)
  rownames(minMaxSpec) <- allSpecs
  
  for(s in 1:length(allSpecs)){
    
    minmax <- c(Inf,-Inf)
    
    for(j in 1:chainLength){
      
      wsj <- which( names(modeBySpec[[j]]) == allSpecs[s] )
      if(length(wsj) == 0)next
      
      if( minmax[1] > modeBySpec[[j]][wsj] )minmax[1] <- modeBySpec[[j]][wsj]
      if( minmax[2] < modeBySpec[[j]][wsj] )minmax[2] <- modeBySpec[[j]][wsj]
    }
    minMaxSpec[s,] <- minmax
  }
    
  
  if(is.null(xtic)){
    sc     <- diff(xrange)
    x1     <- xrange[1] - sc/3
    x2     <- xrange[2] + sc/3
    xtic   <- seq(xrange[1],xrange[2],length=4)
  }
  
  if(is.null(ORD)){
    ordSpec <- order( minMaxSpec[,2] )   #order by minimum value
    if(min(xtic) <= 0 & abs(xtic[1]) > max(abs(xtic[-1])))ordSpec <- order( minMaxSpec[,1] )
    ORD <- allSpecs[ordSpec]
  }
  
  if( length(xtAll[[1]]) == 0 ){
    xtmp <- ytmp <- numeric(0)
    for(j in 2:length(xtAll)){
      xtmp <- append(xtmp,list(xtAll[[j]]))
      ytmp <- append(ytmp,list(ytAll[[j]]))
    }
    xtAll <- xtmp
    ytAll <- ytmp
    COLCODE <- COLCODE[-1]
  }
                     
  
  tmp <- plotPosteriorOrder(xx=xtAll[[1]],yy=ytAll[[1]],x1=xtAll,y1=ytAll,
                            xtic=xtic,
                            ymax=ymax,textSize=textSize,ORD=ORD,COLCODE=COLCODE,
                            MAIN=paste(chain2plot,var2plot),htFraction=htFraction,
                            FILL=FILL)

  #  tmp <- plotPosteriorOrder(xx=xt,yy=yt,xtic=xtic,
#                          ymax=ymax,textSize=textSize,ORD=ORD,COLCODE=COLCODE,
#                          MAIN=paste(chain2plot,var2plot),htFraction=htFraction,FILL=FILL)
  invisible(tmp)

}

chains2density <- function(chainMat,labs=NULL,reverseM=F,varName=NULL){
  
  #assumes column names are varName or 'something_varname'
  
  #chainMat - MCMC output [samples,chains]
  
  chNames <- colnames(chainMat)
  
  if(!is.null(varName)){
    wc <- grep(varName,colnames(chainMat))
    if(length(wc) == 0)stop('varName not found in colnames(chainMat)')
    
    ww <- grep('_',colnames(chainMat))
    if(length(ww) > 0){
      vnameAll <- matrix( unlist(strsplit(colnames(chainMat),'_')),ncol=2,byrow=T)[,2]
      wc <- which(vnameAll == varName)
    }
    chainMat <- chainMat[,wc]
    if(!is.matrix(chainMat))chainMat <- matrix(chainMat,ncol=1)
    colnames(chainMat) <- chNames[wc]
  }
  
  nj <- ncol(chainMat)
  nd <- 512
  
  clab <- colnames(chainMat)
  if(is.null(labs) & !is.null(clab))labs <- clab
  if(is.null(labs) & is.null(clab)) labs <- paste('v',c(1:nj),sep='-')
  
  xt <- yt <- matrix(NA,nj,nd)
  rownames(xt) <- rownames(yt) <- labs
  
  xrange <- signif(range(chainMat),2)
  
  for(j in 1:nj){
    
 #   lj  <- labs[j]
    xj  <- chainMat[,j]
    tmp <- density(xj,n = nd, cut=0)
    xt[j,]  <- tmp$x
    yt[j,]  <- tmp$y
    
  }
  yymax <- max(yt,na.rm=T)
  
  if(reverseM){
    xt <- -t( apply(xt,1,rev) )
    yt <- t( apply(yt,1,rev) )
  }
  
  list(x = xt, y = yt, xrange = xrange, ymax = yymax, chainMat = chainMat)
}

#######################################
plotPosteriorOrder <- function(xx,yy, x1 = NULL, y1 = NULL, 
                               xtic=NULL,ymax=NULL,xlab=' ',
                               ALTLABEL=F,textSize=1,ORD=NULL,
                               COLCODE=NULL,
                               MAIN=NULL,htFraction=1,
                               X5=F,trimx=0,CI=.95,FILL=F){    

  #  yy - each row is a density
  #  xx - each row is the value for the density in y
  #  x1, y1 - lists of additional densities to be plotted with x,y
  #  names are rownames for y
  #  COLCODE  - names used to match with y
  #  htFraction < 1 makes overlap with next density
  #  X5 - plot 5x exaggeration
  #  trimx/2 is fraction of tail to trimx
  
  print(ymax)
  
  wsord   <- apply(yy,1,which.max)
  word    <- order( xx[cbind(1:nrow(xx),wsord)] )
  specOrd <- rownames(yy)[word]
  
  if(length(y1) > 0){
    
    wxx <- numeric(0)
    
    for(j in 1:length(y1)){
      
      xj      <- x1[[j]]
      yj      <- y1[[j]]
      
      wsx   <- matrix(xj[apply(yj,1,which.max)],nrow=1)
      colnames(wsx) <- rownames(xj)
      
      wxx <- appendMatrix(wxx,wsx,SORT=T)
    }
    
    minBySpec <- apply(wxx,2,min,na.rm=T)
    word      <- order(minBySpec)
    specOrd   <- names(minBySpec)[word]
      
  }
  
  specOrd <- rev(unique(specOrd))
  yvalw   <- seq(0,ymax,length=3)
  nc      <- length(specOrd)
  
  xrw     <- signif(range(xx),1)
  yrw     <- round(htFraction*quantile(yy,.95),-1)
  if(yrw == 0)yrw <- signif(quantile(yy,.98),1)

  
  
  if(!is.null(ORD)){
    specOrd <- ORD
    word    <- match(specOrd,rownames(yy))
  }
  
  tail1 <- (1 - CI)/2
  tail2 <- 1 - tail1
  
  postInt <- matrix(0,length(specOrd),3)
  rownames(postInt) <- specOrd
  colnames(postInt) <- c(tail1,tail2,'neg/pos')

 # if(!is.null(x1)){
 #   nc      <- nc + nrow(x1)
 #   wnord   <- apply(y1,1,which.max)
 #   xrw     <- range(c(xx,x1))
 #   yrw     <- round(htFraction*max(rbind(yy,y1)),-2) 
 #   word    <- order( c(xx[cbind(1:nrow(xx),wsord)],x1[cbind(1:nrow(x1),wnord)]) )
 #   specOrd <- c(rownames(yy),rownames(y1))[word]
 # }

 # totalHt <- yrw*nc
  
#  yrw     <- htFraction*ymax
  totalHt <- yrw*nc

  postext <- 4

  if(is.null(ymax))ymax <- yrw


  if(!is.null(COLCODE)){
  #  colSpecs <- COLCODE[specOrd]
    colSpecs <- COLCODE[word]
    if(length(x1) > 0){
      repSpecs <- COLCODE[1:length(x1)]
      colSpecs <- rep(COLCODE[1],length(x1))
    }
  }
  if(is.null(COLCODE)){
    colS     <-  colorRampPalette( c("darkblue","orange") )
    colSpecs <- colS(max(1,length(y1)))
    repSpecs <- colS(max(1,length(y1)))
  }
print(colSpecs)
  ytic <- seq(0,totalHt,by=yrw)
  if(length(ytic) < (nc-1))ytic <- c(ytic,max(ytic) + diff(ytic)[1])

  if(is.null(xtic)){
    xinc <- signif(diff(xrw)/5,1)
    xtic <- signif( seq(xrw[1],xrw[2],by=xinc),1)
  }
  jj   <- rev(ytic)[-1]

################## color box if != 0
  plotSetup(xtic,ytic,xlabel=xlab,yvals=rep(' ',length(ytic)),
            ylabel='Density',lcVer='grey')
  
  
  if(length(y1) > 1){
    postInt <- cbind(postInt, matrix(0,nrow(postInt),(length(y1)-1)) )
  }
  
  for(j in 1:length(specOrd)){

    w0 <- which(rownames(yy) == specOrd[j])
    w1 <- numeric(0)
    
    for(jk in 1:length(y1)){
      ww <- which(rownames(y1[[jk]]) == specOrd[j])
      if(length(ww) == 0)ww <- 0
      w1 <- c(w1,ww)
    }
    
    if(length(c(w0,w1)) == 0)next
    
    if(length(w0) > 0){
      
      dx <- diff(xx[w0,])[1]
      cj <- cumsum(dx*yy[w0,])
      if(cj[1] > 0)cj[1] <- 0
      if(cj[length(cj)] < 1)cj[length(cj)] <- 1    #cdf
      rj <- xx[w0, findInterval(c(trimx/2,1 - trimx/2),cj) ]   #determine x limits
      yj <- jj[j] + yy[w0,]
      yj[yj > ymax] <- ymax
      
      postInt[j,1:2] <- xx[w0,findInterval( c(tail1,tail2),cj)]
      if(postInt[j,1] > 0) postInt[j,3] <- 1
      if(postInt[j,2] < 0) postInt[j,3] <- -1
    }
    
    if(length(y1) > 1){
      
      for(jk in 2:length(y1)){

        if(w1[jk] == 0)next
        
        xj <- x1[[jk]][w1[jk],]
        dx <- diff(xj)[1]
        cjj <- cumsum(dx*y1[[jk]][w1[jk],])
        if(cjj[1] > 0)cjj[1] <- 0
        if(cjj[length(cjj)] < 1)cjj[length(cjj)] <- 1    #cdf
        rjj <- xj[ findInterval(c(tail1,tail2),cjj) ] 
        if(rjj[1] > 0) postInt[j,2 + jk] <- 1
        if(rjj[2] < 0) postInt[j,2 + jk] <- -1
      }               
    }                    
    
    anyNeg <- anyPos <- F
    if(max(postInt[j,-c(1:2)]) > 0)anyPos <- T
    if(min(postInt[j,-c(1:2)]) < 0)anyNeg <- T

    jmax <- totalHt
    if(j > 1)jmax <- jj[j-1]
    if(anyPos)rect(max(c(0,xtic[1])),jj[j],max(xtic),jmax,col='wheat',border='wheat')
    if(anyNeg)rect(xtic[1],jj[j],min(c(0,max(xtic))),jmax,col='wheat',border='wheat')
    if(anyNeg & anyPos)rect(xtic[1],jj[j],max(xtic),jmax,col='wheat',border='wheat')
  }

  abline(v=0,col='grey',lwd=3)
###########################################

  ldens <- 20

  for(j in 1:length(specOrd)){

    w0 <- which(rownames(yy) == specOrd[j])
    w1 <- numeric(0)
    
    for(jk in 1:length(y1)){
      ww <- which(rownames(y1[[jk]]) == specOrd[j])
      if(length(ww) == 0)ww <- 0
      w1 <- c(w1,ww)
    }

    if(length(c(w0,w1)) == 0)next
  
    xl <- numeric(0)

    if(length(w0) > 0){

       xl <- range(xx[w0,])
       yj <- yy[w0,]
       yj[1] <- yj[length(yj)] <- 0
       if(X5){
         yk <- 5*yj
         yk[yk > ymax] <- ymax
         yk <- jj[j] + yk
       }
       yj[yj > yrw] <- yrw
       yj <- jj[j] + yj

       lines(xx[w0,],yj,lwd=2,col=colSpecs[j])

       if(FILL){
         polygon(c(xx[w0,],xx[w0,1]),c(yj,yj[1]),
                 border=repSpecs[1],col=colSpecs[j])
       }
       if(X5)polygon(c(xx[w0,],xx[w0,1]),c(yk,yk[1]),border=colSpecs[j],
              col=colSpecs[j],density=ldens)
    }
    
    if(length(w1) > 0){
      
      for(jk in 1:length(y1)){
        
        if(w1[jk] == 0)next
        x2 <- range(x1[[jk]][w1[jk],])
        
    #    if(diff(x2) < .1){
          
          yj <- y1[[jk]][w1[jk],]
          xj <- x1[[jk]][w1[jk],]
          yj[yj > yrw] <- yrw
          yj <- jj[j] + yj
        
     #     lines(xj,yj,lwd=5,col='white')
          lines(xj,yj,lwd=2,col=repSpecs[jk])
          if(FILL){
            polygon( c(xj[1],xj,xj[length(xj)], xj[1]) , c(jj[j],yj,jj[j],jj[j]),
                     border=repSpecs[jk],col=repSpecs[jk])
          }
          xl <- range(c(xl,x2))
  #      }
      }
    }
  #  lines(xl,c(jj[j],jj[j]),col=colSpecs[j],lty=2,lwd=2)

    yt <- htFraction*max(yy[w0,])
    xt <- xx[w0,which.max(yy[w0,]) - 12]
    xxl <- max(xl)
    
    leftDist <- min(xl) - min(xtic)
    rightDist <- max(xtic) - max(xl)
    
    LEFT <- T
    
    if(leftDist < 0)LEFT <- F
    if(leftDist > 0 & leftDist > rightDist)LEFT <- T
    if(rightDist > 0 & rightDist > leftDist)LEFT <- F
    

    if(!ALTLABEL){
      if(LEFT){
        xxl <- min(xl)
        postext <- 2
      }
      if(!LEFT){
        postext <- 4
        xxl <- max(xl)
      }
    }

    if(xxl > max(xtic) & !ALTLABEL){
        xxl <- min(xl)
        postext <- 2
    }

   if(postext == 4)xxl <- max(xl)
   if(postext == 2)xxl <- min(xl)

   text(xxl,jj[j],specOrd[j],col=colSpecs[j],pos=postext,cex=textSize)

   if(ALTLABEL & postext == 4){
     postext <- 2
   } else {postext <- 4}

  }

  ys2 <- yvalw[length(yvalw)]
  
  if(!is.null(MAIN))title(MAIN)
  
  
  tk    <- postInt[,-c(1:2)]
  
  if(is.matrix(tk)){
    sump  <- which(tk == 1,arr.ind=T)[,2]
    sumn  <- which(tk == -1,arr.ind=T)[,2]
    sump  <- round(byIndex(sump*0+1,sump,sum,coerce=T)/nrow(tk),3)*100
    sumn  <- round(byIndex(sumn*0+1,sumn,sum,coerce=T)/nrow(tk),3)*100
  }
  if(!is.matrix(tk)){
    sump  <- which(tk == 1,arr.ind=T)
    sumn  <- which(tk == -1,arr.ind=T)
    sump  <- length(sump)/length(tk)*100
    sumn  <- length(sumn)/length(tk)*100
  }
  
  sump <- paste(sump,'%',sep='')
  sumn <- paste(sumn,'%',sep='')
  
  list(specOrd=specOrd, postInt = postInt, repSpecs = repSpecs, percPos = sump, percNeg = sumn)

}

getScoreNorm <- function(x,mu,xvar){  #Gneiting and Raftery's proper scoring rule

  #outcome x, prediction mean variance (mu, xvar)

  - ( (x - mu)^2)/xvar - log(xvar)

}

getScoreBinary <- function(x,p){    #binary outcome x (0,1), probability p

  (log(p))^x*(log(1 - p))^(1 - x)

}
advection <- function(u,vel,dt,dx){
	
	#finite difference left to right with velocity vel
	#Crank-Nicholson, http://farside.ph.utexas.edu/teaching/329/lectures/node92.html
	#u - current concentration
	
	n <- length(u)
	const <- vel*dt/dx
	w <- rep(0,n)
	
	b <- rep(1,n)
   a <- b*(.25*const)
   c <- -a
   
   for(j in 2:(n-1))w[j] <- u[j] - .25*const*(u[j+1] - u[j-1])
   
   w <- triDiagonal(b,a,c,w)
   w[1] <- w[n] <- 0
   w[w < 0] <- 0
   w
}
   
triDiagonal <- function(a,b,c,y){
	
  #  solves Ax = y for tridiagonal matrix A
  #  a       main diagonal 
  #  b       upper diagonal 
  #  c       lower diagonal
  #  y       right-hand side vector

  n <- length(a)

  #   factorization
  for (i in 1:(n-1)){
     b[i]   <- b[i] / a[i]
     a[i+1] <- a[i+1] - c[i] * b[i]
  }

  #   forward substitution
  y[1] <- y[1]/a[1]
  for(i in (2:n))y[i] <- ( y[i] - c[i-1] * y[i-1] ) / a[i]

  #   back substitution
  for(i in (n-1):1)y[i] <- y[i] - y[i+1] * b[i]
  
  y
}

mcmcSetup <- function(chains=NULL,sums=NULL,chainVals = NULL, 
                      sumVals = NULL,ng=1000){

  chainList <- chainSum <- numeric(0)

  if(!is.null(chainVals))for(k in 1:length(chains))assign(chains[k],chainVals[[k]])
  if(!is.null(sumVals))  for(k in 1:length(sums))  assign(sums[k],sumVals[[k]])

  if(!is.null(chains)){

      ne <- length(chains)
      for(k in 1:ne){

         kv   <- get(chains[k])
         vn   <- paste(chains[k],'Chain',sep='')
         kmat <- matrix(NA,ng,length(kv))

         gnames <- c(1:length(kv))

         if(!is.matrix(kv) & !is.null(names(kv)))gnames = names(kv)

         if(is.matrix(kv)){
            if( is.null(rownames(kv)) ) rownames(kv) <- c(1:nrow(kv))
            if( is.null(colnames(kv)) ) colnames(kv) <- c(1:ncol(kv))
            gnames <- outer(rownames(kv),colnames(kv),paste,sep='_')
         }
         colnames(kmat) <- gnames
         assign(vn,kmat)
         chainList <- append(chainList,list(get(vn)))
         names(chainList)[k] <- vn
     }

   }
         
   if(!is.null(sums)){
      
     ns <- length(sums)
     for(k in 1:ns){

        kv <- get(sums[k])
        v1 <- paste(sums[k],'Sum',sep='')
        v2 <- paste(sums[k],'Sum2',sep='')
        kg <- rep(0,length(kv))
        assign(v1,kg)
        assign(v2,kg)
        chainSum <- append(chainSum,list(get(v1)))
        names(chainSum)[length(chainSum)] <- v1
        chainSum <- append(chainSum,list(get(v2)))
        names(chainSum)[length(chainSum)] <- v2
     }
   }
   list(chainList = chainList, chainSum = chainSum)
}

loadChains <- function(g,burnin=1,chains=chains,sums=NULL, pars, sumVals,
                     chainList, chainSum){


  for(k in 1:length(chains)){
     chainList[[k]][g,] <- pars[[ which(names(pars) == chains[k]) ]]
  }

   if(g >  burnin & !is.null(sums)){
     for(k in 1:length(sums)){

       w1 <- which(names(chainSum) == paste(sums[k],'Sum',sep='') )
       w2 <- which(names(chainSum) == paste(sums[k],'Sum2',sep='') )

       ws <- which(names(sumVals) == sums[k])

       chainSum[[w1]] <- chainSum[[w1]] + sumVals[[ws]]
       chainSum[[w2]] <- chainSum[[w2]] + sumVals[[ws]]^2
     }
   }
   list(chainList = chainList, chainSum = chainSum)
}

processMCMC <- function( chains=NULL, sums = NULL, burnin=1, ng=nrow(chainList[[1]]), PLOTS=F,
                         outfolder=character(0),outfile='mcmc',chainList = chainList, chainSum = chainSum,
                         chainlength= 50000){

  CPLOT <- T
  file  <- outfile

  if(length(outfolder) > 0){
    if(!file.exists(outfolder))dir.create(outfolder)
    file <- paste(outfolder,'/',outfile,sep='')
  }

  postSummary <- postMuSe <- chainOut <- parMeans <- parSds <- numeric(0)

  kk <- 0

  if(!is.null(chains)){

    for(k in 1:length(chains)){

      if(PLOTS)CPLOT <- T
      x  <- col2Mat( chainList[[k]][burnin:ng,] )
      vk <- apply(x,2,var)
      kx <- which(vk > 0)
      if(length(kx) == 0)next

      kk <- kk + 1
 #     x  <- col2Mat(x[,kx])

      nx <- ncol(x)

      if(nx > 12)CPLOT <- F

      tmp <- processPars(x,rep(0,nx),CPLOT=CPLOT,DPLOT=F,
                  burnin=burnin)$summary[,c('estimate','se','0.025','0.975')]

      tmp <- row2Mat(tmp)
      rownames(tmp) <- paste(chains[k],rownames(tmp),sep='_')
      postSummary <- rowBind(postSummary,tmp,rownames(tmp))

      parMeans <- append(parMeans,list(tmp[,1]))
      parSds   <- append(parSds,list(tmp[,2]))
      names(parMeans)[kk] <- names(parSds)[kk] <- chains[k]

      if(nx > 12){

          i1 <- 1
          i2 <- 12
          while(i1 < ncol(x)){
             ii <- i1:i2
             processPars(x[,ii],rep(0,length(ii)),CPLOT=PLOTS)
             if(PLOTS)dev.copy2pdf(file=paste(file,i1,'_',i2,chains[k],'.pdf',sep=''))
             i1 <- i1 + 12
             i2 <- i2 + 12
             if(i2 > ncol(x))i2 <- ncol(x)
          }
       }

      if(nrow(x) > chainlength){
        thin <- round( seq(1,nrow(x),length=chainlength), 0)
        x <- as.matrix(x[thin,])
      }
      chainOut <- append(chainOut,list(x))

      if(PLOTS & CPLOT){
        dev.copy2pdf(file=paste(file,chains[k],'.pdf',sep=''))
      }
    }
  }

  if(!is.null(sums)){

    ntot <- ng - burnin

    for(k in 1:length(sums)){

       w1 <- which(names(chainSum) == paste(sums[k],'Sum',sep='') )
       w2 <- which(names(chainSum) == paste(sums[k],'Sum2',sep='') )

       mu <- chainSum[[w1]]/ntot
       se <- sqrt(chainSum[[w2]]/ntot - mu^2)

        postMuSe <- append(postMuSe,list(cbind(mu,se)))
     }
    names(postMuSe) <- sums
  }

  list(postSummary = postSummary, postMuSe = postMuSe, parMeans = parMeans, parSds = parSds, chainOut = chainOut)
}


priorInteractionIndex <- function(vnames){  #has all variable names, including interactions
  
  if('years' %in% vnames)vnames <- vnames[vnames != 'years']
  vn  <- vnames
  ivn <- grep('X',vnames)
  
  m <- numeric(0)
  
  if(length(ivn) == 0)return( list(isIntPrior = m, noIntPrior = m, intPrior = m) )
  
  aprior <- allPrior[allPrior %in% vn]
  
  vn <- vn[-ivn]
  ivn <- vnames[ivn]
  
  isIntPrior <- matrix(0,length(aprior),length(ivn))    #indicator: prior on main effects in interaction
  colnames(isIntPrior) <- ivn
  rownames(isIntPrior) <- aprior
  
  tmp <- matrix( unlist(strsplit(ivn,'X')), nrow=2)
  
  for(m in 1:ncol(tmp)){
    mm <- match(tmp[,m],aprior)
    wm <- which(is.finite(mm))
    isIntPrior[mm[wm],m] <- 1
  }
  
  cc <- length(vn) + which(colSums(isIntPrior) == 0)    # neither has prior
  
  noIntPrior <- c(1:length(vn),cc)                    # main effects, interactions without priors
  intPrior   <- c(1:length(vnames))                      # interactions with priors
  if(length(noIntPrior) > 0)intPrior   <- intPrior[-noIntPrior]
  
  list(isIntPrior = isIntPrior, noIntPrior = noIntPrior, intPrior = intPrior)
}



updateIntBeta <- function(xx,yy,bb,lo,hi,sigma,isIntPrior,intPrior,noIntPrior){  
  
  #allPrior in main program
  
  if(is.null(lo))lo <- bb*0 - 100
  if(is.null(hi))hi <- bb*0 + 100

  if(!is.matrix(lo))lo <- as.matrix(lo)
  if(!is.matrix(hi))hi <- as.matrix(hi)
  
  wi <- grep('X',rownames(bb))
  
  aprior <- rownames(isIntPrior)
  
  #main effects conditional on interaction terms
  
  yz <- yy - as.matrix(xx[,intPrior])%*%bb[intPrior,]   #from interactions
  
  bb[noIntPrior,]  <- bUpdateMVN_Rcpp(xx[,noIntPrior],yz,b=bb[noIntPrior,],
                                      lo=lo[noIntPrior,],hi=hi[noIntPrior,],sigma)
  
  nk <- ncol(bb)
  minB <- -bb[aprior,]      #min bound
  if(!is.matrix(minB))minB <- matrix(minB,nrow=length(aprior),ncol=nk)
  
  ss <- sample(length(wi))
  for(k in ss){

    wk <- wi[k]
    yz <- yy - xx[,-wk]%*%bb[-wk,]
    
    mink <- minB*isIntPrior[,k]
    mink[isIntPrior[,k] == 0] <- -500
    mm   <- apply(mink,2,which.max)
                    
    lok <- mink[cbind(mm,1:nk)]
    
    bb[wk,] <- bUpdateMVN_Rcpp(xx[,wk],yz,b=bb[wk,],lo=lok,hi=hi[wk,],sigma)
    
    minB <- minB - matrix(bb[wk,],nrow(minB),nk,byrow=T)*isIntPrior[,k]
  }
  bb
}

updateBetaMVN <- function(xx,yy,bb,lo=NULL,hi=NULL,sigma,delta,priorO=diag(1,nrow(sigma)),
                          priorDmu=rep(.1,nrow(sigma)),priorDvar=rep(1,nrow(sigma)),INT=F,
                          isIntPrior=NULL,intPrior=NULL,
                          noIntPrior=NULL,varLo=rep(1e-8,ncol(yy)),varHi=rep(10000,ncol(yy)),
                          SIGMA=T){
  #bb must have rownames

  ss <- dd <- sigma*0
  
  if(INT)bb <- updateIntBeta(xx,yy,bb=bb,lo=lo,hi=hi,sigma,
                               isIntPrior=isIntPrior,intPrior=intPrior,
                               noIntPrior=noIntPrior)
  
  if(!INT)bb  <- bUpdateMVN_Rcpp(xx,yy,b=bb,lo=lo,hi=hi,sigma)
  
  if(SIGMA){
 #   tmp   <- updateScInvWish(yy,xx%*%bb,delta=delta,
 #                            priorDmu=priorDmu,priorDvar=priorDvar,
 #                            varLo=rep(1e-8,ncol(yy)),varHi=rep(10000,ncol(yy)) )
    
    tmp   <- updateScInvWish(yy,xx%*%bb,delta=delta,
                             priorDmu=priorDmu,priorDvar=priorDvar,
                             varLo=varLo,varHi=varHi )
    
    ss <- tmp$sigma
    dd <- tmp$delta
  }
  
  list(beta = bb, sigma = ss, delta = dd)
  
} 






qprob <- function(p){  #evaluate likelihood for pathogen model: multinomial

    theta <- p[1]
    phi   <- p[2]
    s0    <- p[3]
    s1    <- p[4]
    q00 <- (1 - theta)*(1 - s0) + theta*(1 - phi)*(1 - s1)
    q01 <- (1 - theta)*s0 + theta*(1 - phi)*s1
    q10 <- theta*phi*(1 - s1)
    q11 <- theta*phi*s1

    c(q00,q01,q10,q11)
}

like_multinom <- function(pars){
   q <- qprob(pars)
  -sum(nvec*log(q))
}

outDataFrame <- function(...){  # a list of variable names coerced to dataframe (e.g. outDataFrame('a','b'))

  x  <- list(...)

  nx <- length(x)
  nc <- 0
  nr <- 0
  xn <- character(0)
  

  for(j in 1:nx){
    xj <- get(x[[j]])
    if(!is.null(dim(xj)) )xj <- as.matrix(xj)
    if( is.null(dim(xj)) )xj <- matrix(xj,ncol=1)
    nj <- nrow(xj)
    if(is.null(rownames(xj))){
      if(ncol(xj) == 1)colnames(xj) <- x[[j]]
      if(ncol(xj) > 1 )colnames(xj) <- paste(x[[j]],c(1:nj),sep='-')
   }
   xj <- t(xj)

    nr <- max(nr,ncol(xj))
    nc <- nc + nrow(xj)
    xn <- c(xn,rownames(xj))
  }

  out <- data.frame(matrix(NA,nr,nc))
  colnames(out) <- xn

  j1 <- 1

  for(j in 1:nx){
     xj <- get(x[[j]])
     if(is.matrix(xj)) xj <- t(xj)
     if(!is.matrix(xj))xj <- as.matrix(xj,ncol=1)
     j2 <- j1 + ncol(xj) - 1
     out[1:nrow(xj),j1:j2] <- xj
     j1 <- j2 + 1
  }
  format(out,justify='right')
}
    
write2file <- function(x,file,append=F,row.names=T){  #x is the name of an object (in quotes)

  if(append){
    xx <- as.data.frame( matrix(c(' ',x),2,1) )
    write.table(xx,file,row.names=row.names,col.names=F,append=append,quote=F)
  }

  z <- get(x)
  rnames <- F

 # if(is.null(rownames(z)))rnames <- F

  if(!is.null(rownames(z)) & row.names){
    rnames <- rownames(z)
    rownames(z) <- NULL
    nr <- nchar(rnames)
    nc <- max(nr)
    spaces <- rnames
    nl <- length(spaces)
    while(nl > 0){
      spaces[nchar(spaces) < nc] <- paste(spaces[nchar(spaces) < nc],' ',sep='')
      nl <- length( which(nchar(spaces) < nc) )
    }
    rnames <- spaces
  }
   

 # if(!is.data.frame(z))z <- data.frame(z)

  z <- format(z,justify='right')
  write.table(z,file,append=append,quote=F,row.names=rnames)
}

updateSSRW <- function(states,y,missing,tg,sg){        #state-space random walk 
	#update continuous states, random walk
	#missing times, obs y, obs error tg, process error sg

  for(t in 1:nt){

    VI <- 0
    v  <- 0

    if(!t %in% missing){          #observations
      v  <- y[t]/tg
      VI <- 1/tg
    }

    if(t < nt){              #t+1 term excluded for last 
      v  <- v + states[t+1]/sg 
      VI <- VI + 1/sg
    }

   if(t > 1){                #t-1 term excluded for 1st 
      v  <- v + states[t-1]/sg
      VI <- VI + 1/sg
   }

   V     <- 1/VI
   states[t] <- rnorm(1,V*v,sqrt(V))
  }
  states
}

updateSSB <- function(states,sg,priorB,priorIVB){  #intercept SS model
	
	V <- 1/( (length(y) - 1)/sg + priorIVB)
	v <- (states[length(states)] - states[1])/sg + priorB*priorIVB
	rnorm(1,V*v,sqrt(V))
}


updateSSparsMet <- function(b,priorB,priorVB,
                            lo=rep(-Inf,length(b)),hi=rep(Inf,length(b)),
                            x,sigma,dt=1,propVar){
                            	
  require(mvtnorm)
  
  k <- length(b)
  if(k == 1)pb <- tnorm(1,lo,hi,b,propVar)
  if(k > 1) pb <- tnorm.mvt(b,b,propVar,lo,hi) 
  
  xnow <- xg[-nt] + fx(xg[-nt],b)*dt     #predicted x
  xnew <- xg[-nt] + fx(xg[-nt],pb)*dt
	
  pnow <- sum(dnorm(xg[-1],xnow,sqrt(sigma*dt),log=T)) +
          dmvnorm(t(b),priorB,priorVB,log=T)
  pnew <- sum(dnorm(xg[-1],xnew,sqrt(sigma*dt),log=T)) +
	      dmvnorm(t(pb),priorB,priorVB,log=T)
  atmp <- acceptMH(pnow,pnew,b,pb,BLOCK=T)
  b   <- atmp$x

  list(b = b, aa = atmp$accept)

}
updateSSstatesMet <- function(x,y,b,sigma,tau,lo=-Inf,hi=Inf,notMiss=c(1:length(x)),
                              dt=1,propSd=.1){   #propose/accept x as a block
  nt   <- length(x)
  xp   <- tnorm(nt,lo,hi,x,rexp(nt,1/propSd) )       #proposed x
  
  xnew <- xp[-nt] + fx(xp[-nt],b)*dt      #mean proposed
  xnow <- x[-nt]  + fx(x[-nt],b)*dt       #mean current

  pnow <- sum(dnorm(x[-1],xnow,sqrt(sigma*dt),log=T)) +
	       sum(dnorm(y[notMiss],x[notMiss],sqrt(tau),log=T))
	       
  pnew <- sum(dnorm(xp[-1],xnew,sqrt(sigma*dt),log=T)) +
	       sum(dnorm(y[notMiss],xp[notMiss],sqrt(tau),log=T))
	       
  tmp <- acceptMH(pnow,pnew,x[-1],xp[-1],BLOCK=T)
  x[-1]   <- tmp$x

  list(x = x, aa = tmp$accept)
}

diagXcovar <- function(diagmat,covmat){

  #  diagmat %*% sigma %*% t(diagmat) 
  #  diagmat - rows are diagonal of a k by k diagonal matrix
  #  covmat  - covariance matrix (symm, pos-def)
  #  returns nn by k^2 matrix, each row are elements of a k by k matrix
  
  nn <- nrow(diagmat)
  k  <- nrow(covmat)
  lm <- matrix(0,nn,k^2)

  km <- matrix(1:(k^2),k,k)

  ii <- km[!upper.tri(km)]
  jj <- t(km)[!upper.tri(km)]

  kk <- 0
  for(i in 1:k){
    for(j in i:k){
      kk <- kk + 1
      lm[,ii[kk]] <- lm[,jj[kk]] <- diagmat[,i]*diagmat[,j]*covmat[i,j]
    }
  }
  lm
}

histWeight <- function(xx,wt,bins,removeZeros=T,useMids=F){

  # xx   - data to be binned
  # wt   - wt assigned to each x
  # bins - partition for x

  #  example: wt is 1/(plotarea), then results are no/area

  xx <- as.vector(xx)
  wt <- as.vector(wt)

  nb <- length(bins)
  db   <- diff(bins)

  if(useMids){

    mids <- bins[-nb] + db/2
    bins <- c(mids,bins[nb]+db[nb-1]/2)
  }
  
  if(!removeZeros)xx[xx < bins[1]] <- bins[1]
  
  x1 <- findInterval(xx,bins) 
  if(removeZeros)x1[xx == 0] <- NA


  freq <- by(as.vector(wt),as.vector(x1),sum,na.rm=T) 

  wb <- as.numeric( names(freq) )
  whist <- bins*0
  whist[wb] <- unlist(freq)
  
  db <- c(db,db[length(db)])
  dens <- whist/sum(whist)

  list(wtHist = whist, bins = bins, dens=dens)
}
  

histByYr <- function(xx,yr,xbin=1,breaks=NULL,dens=F,dt=1,weight=xx*0 + 1,deathyr=NULL){

  # xx has observations as rows, years as columns
  # yr is the vector of years corresponding to columns
  # change - if true, return change in xx distribution
  # dens   - if true, return density
  # deathyr    - length nn vector of death years

  ny <- length(yr)
  nn <- nrow(xx)
  
  if(is.null(deathyr))deathyr <- rep(NA,nn)
  
  dt <- 1

  maxx   <- round( max(xx,na.rm=T) ,0) + xbin
  if(is.null(breaks))breaks <- seq(0,maxx,by=xbin)
  mx     <- length(breaks) 

  times  <- c(1:round(ny/dt,0) )
  nt     <- length(times)
  
  gMat <- matrix(0,mx,nt)

  nh <- length(breaks)
  
  xx[xx == 0] <- NA
  
  ymat <- matrix(yr,nrow(xx),length(yr),byrow=T)
  incr <- t( apply(xx,1,diff)/apply(ymat,1,diff) )
  if(nrow(incr) == 1 & nrow(xx) > 1)incr <- t(incr)
  
  liveMat <- ymat*(xx*0 + 1)
  wlive <- which(is.finite(incr),arr.ind=T)
  wlNo  <- which(is.finite(xx),arr.ind=T)
  
  incr[incr < 0] <- 0
  
  inow  <- findInterval(xx[wlive],breaks)
  iyr   <- match(ymat[wlive],yr)
  
  inowNo  <- findInterval(xx[wlNo],breaks)
  iyrNo   <- match(ymat[wlNo],yr)
  
  groMu <- round(byFunctionRcpp(incr[wlive],inow,iyr,gMat*0,gMat*0,MEAN=T),3)
  groVr <- signif(byFunction(incr[wlive],inow,iyr,gMat*0,var),3)
  groNo <- byFunctionRcpp(xx[wlNo]*0+1,inowNo,iyrNo,gMat*0,gMat*0,MEAN=F)
  
  wt <- byFunctionRcpp(weight[wlNo],inowNo,iyrNo,gMat*0,gMat*0,MEAN=T)
  
  dist  <- groNo*wt
  ddist <- groNo[,-1] - groNo[,-nt]
  
  deadTable <- groNo*0
  
  wd        <- which(is.finite(deathyr))
  
  if(length(wd) > 0){
    
    wdead     <- cbind( wd, match(deathyr[wd],yr) - 1 )
    
    diamd     <- xx[wdead]
    wdd       <- which(!is.finite(diamd))
    if(length(wdd) > 0){
      dd <- xx[wdead[wdd,1],]
      if(!is.matrix(dd))dd <- matrix(dd,nrow=1)
      dd <- apply(dd,1,max,na.rm=T)

      diamd[wdd] <- dd
      missyr       <- findInterval( deathyr[wd[wdd]] ,yr) + 1 
      missyr[missyr > nt] <- nt
      missyr[missyr == 0] <- 1
      wdead[wdd,2] <- missyr
      wi <- which(!is.finite(diamd))
      if(length(wi) > 0){
        diamd <- diamd[-wi]
        wdead <- wdead[-wi,]
        wd    <- wd[-wi]
        if(!is.matrix(wdead))wdead <- matrix(wdead,1)
      }
    }
    
    idead     <- findInterval(diamd,breaks)
    
    if(length(wd) == 1){
      cc <- cbind(idead,wdead[2])
      deadTable[cc] <- 1
    }
    if(length(wd) > 1)deadTable <- byFunctionRcpp(rep(1,length(wd)),idead,wdead[,2],gMat*0,gMat*0,MEAN=F)
  }

  dsurv     <- groNo - deadTable
  dsurv[dsurv < 0] <- 0              #check this
  
  if(!is.matrix(ddist))ddist <- matrix(ddist,ncol=1)


  groNo[is.na(groNo)] <- 0
  w0 <- which(groNo == 0,arr.ind=T)
  groMu[w0] <- groVr[w0] <- NA
  
  dsurv[,nt] <- dist[,nt]
  
  colnames(dist) <- yr
  colnames(ddist) <- yr[-ny]

  list(dist = dist, change = ddist, size = breaks+.5*dt,
       surv = dsurv, groMu  = groMu, groVr = groVr, groNo = groNo)
}


getMonthNames  <- function(){

  mm <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
  invisible(mm)
}


subETOPO5 <- function(lon,lat){

  #lon is negative degrees W, positive degrees E
  #lat is positive N
  
  require(geomapdata)
  

  data(ETOPO5)

  ilon <- nrow(ETOPO5)/360
  ilat <- ncol(ETOPO5)/180
  
  x <- lon
  y <- lat
  
  
  xx <- 360 + lon
  
#  lon[lon <= 0] <- -lon[lon <= 0]  + 180  #(0, 360) scale
  lat[lat > 0] <- 90 - lat[lat > 0]
  lat[lat < 0] <- -lat[lat < 0] + 90

#  rlon <- lon*ilon
  rlon <- xx*ilon
  rlat <- lat*ilat

  print(rlon)
  print(rlat)

  z <- ETOPO5[rlon[1]:rlon[2],]
  z <- z[,rlat[1]:rlat[2]]

  x <- seq(x[1],x[2],length=nrow(z))
  y <- seq(y[1],y[2],length=ncol(z))

  list(x = x, y = y, z = z)

}

long2UTMZone <- function(lon){     #assumes west of prime meridian, (-180, 0)
  floor((lon + 180)/6 %% 60) + 1
}

lonLat2UTM <- function(lon,lat){
  
  #xy has 2 columns, long, lat
  #long is longitude, used to find UTM zone
  
  require(PBSmapping)
  
  xy <- cbind(lon,lat)
  
  if(!is.matrix(xy))xy <- matrix(xy,1,2)
  colnames(xy) <- c('X','Y')
  
  zone <- long2UTMZone(xy[,1])
  
  zz <- sort(unique(zone))
  nz <- length(zz)
  
  ll <- xy*0
  
  for(j in zz){
    xyj <- xy[zone == j,]
    if(!is.matrix(xyj))xyj <- matrix(xyj,1)
    colnames(xyj) <- c('X','Y')
    attr(xyj,'projection') <- 'LL'
   
    ll[zone == j,] <- as.matrix(convUL(xyj,km=F))
  }
  ll
}

UTM2latlon <- function(xy,long=-80){
  
  #xy has 2 columns, 'X', 'Y' for UTMx, UTMy
  #long is longitude, used to find UTM zone
  
  require(PBSmapping)
  
  colnames(xy) <- c('X','Y')
  
  zone <- long2UTMZone(long)
  if(length(zone) == 1)zone <- rep(zone,nrow(xy))
  
  zi <- unique(zone)
  utm <- xy*0
  
  for(i in zi){
    wi  <- which(zone == i)
    xyi <-  matrix( xy[wi,], ncol=2 )
    colnames(xyi) <- c('X','Y')
    attr(xyi,'zone') <- i
    attr(xyi,'projection') <- 'UTM'
    utm[wi,] <- as.matrix( convUL(xyi,km=F) )
  }
  utm
}

distributionMapLittle <- function(spec, dataLoc="/nfs/clark/clark.unix/fia/little/" ){
  #returns shape file
  
  require(maptools)
  require(classInt)
  require(foreign)
  
  treeCodes <- read.table('/nfs/clark/clark.unix/allocationmodel/datafiles/treeCodesDuke.txt',
                          header=T)
  fiaCode <- as.character(treeCodes[ match(spec,treeCodes[,'code']),'fiaCode'] )

  fold <- tolower(spec)
  ff <- list.files( paste(dataLoc,fiaCode,sep='') )
  ff <- paste(dataLoc,fiaCode,ff,sep='/')

  frame <- read.dbf( ff[grep('dbf',ff)] )
 # attr  <- read.dta( ff[grep('dta',ff)] )
   shap  <- readShapePoly( ff[grep('shp',ff)] )
  shap
}

dayLength <- function(JD,lat){
  
  nl  <- length(lat)
  nj  <- length(JD)
  latNew <- lat/180
  
  JDmat <- JD
  
  if(nl > 1 & nj > 1){
    latNew <- matrix(latNew,nl,nj)
    JDmat  <- matrix(JD,nl,nj,byrow=T)
  }
  
  p <- asin( .39795*cos(.2163108 + 2*atan(.9671396*tan(.0086*(JDmat - 186)))) )
  tmp <- 24 - (24/pi)*acos( (sin(.8333*pi/180) + sin(latNew*pi)*sin(p)) /
                       (cos(latNew*pi)*cos(p)))
  tmp
}
monthlyPET <- function(yi,mi,tempMatrix,precMatrix,lat){
  
  # m  = nyr*12 (serial months)
  # mi = month index (1:12)
  # yi = yr index    (e.g., 1990, 1991)
  # returns n X month*nyr matrix and n X nyr matrix 
  # tempMatrix - n X m location by month
  # precMatrix - n X m location by month
  # yi     - length-m year index
  # mi     - length-m month index
  
  if(!is.matrix(tempMatrix))tempMatrix <- matrix(tempMatrix,1)
  if(!is.matrix(precMatrix))precMatrix <- matrix(precMatrix,1)
  
  n <- nrow(tempMatrix)
  m <- ncol(tempMatrix)
  yrvec <- min(yi):max(yi)
  yseq  <- yrvec - yrvec[1] + 1
  nyr   <- length(yrvec)
  ymat  <- matrix(yseq,n,nyr,byrow=T)
  
  yindex <- match(yi,yrvec)
  
  JD <- c(1:365)
  
  dl <- dayLength(JD,lat)
  mo <- daysSinceDate(1,1,yrvec[1],c(2:12),1,yrvec[1],nineteen=F)
  di <- findInterval(JD,mo) + 1
  
  yvec <- as.vector( matrix(yindex,n,m,byrow=T) )
  ivec <- as.vector( matrix(c(1:n),n,m) )
  
  # heat index by year
  tmp1   <- tempMatrix
  tmp1[tmp1 < 0] <- 0
  HI    <- byFunctionRcpp(as.vector((tmp1/5)^1.514),ivec,yvec,ymat*0,ymat*0,MEAN=F)
  HI    <- HI[yvec]
  alpha <- (6.75*1e-7)*HI^3 - (7.71*1e-5)*HI^2 + (1.792*1e-2)*HI + .49239
  
  tmp <- aggregateSequence(di,dl,action='mean')  #mean monthly daylength
  dlMon <- tmp[[1]]
  dlLen <- tmp[[2]]
  dlLen <- dlLen[,mi]
  
  PET     <- 16*dlLen/12*(10*tmp1/HI)^alpha
  list(PET = PET, daylength = dlLen)
}

monthlyPHr <- function(yi,mi,tempMatrix=tempj,precMatrix=precj,lat=lonlatj[2],tempThresh=4){
  
  tmp <- monthlyPET(yi,mi,tempMatrix,precMatrix,lat)
  PET <- tmp$PET
  daylen <- tmp$daylength
  
  pExcess <- precMatrix - PET
  
  dayFraction <- daylen
  #  pHours  <- dayFraction*pExcess
  #  pHours[pHours < 0] <- 0
  
  posE <- pExcess        #surplus
  posE[posE > 0] <- 1
  posE[posE <= 0] <- 0
  posD <- 1 - posE       #deficit
  
  DD    <- tempMatrix
  DD[DD < tempThresh] <- 0   #threshold for DD

  degreeHrs      <- DD*dayFraction
  degreeHrExcess <- degreeHrs*pExcess
  degreeHrPos    <- degreeHrs*posE        #degreeMonthHrs
  degreeHrNeg    <- degreeHrs*posD
  
  tmp <- aggregateSequence(yi,degreeHrExcess,action='sum')  #total annual daylength
  degHrExYr <- tmp$data
  
  tmp <- aggregateSequence(yi,degreeHrPos ,action='sum')  #total annual daylength
  degHrPosYr <- tmp$data
  
  tmp <- aggregateSequence(yi,degreeHrNeg ,action='sum')  #total annual daylength
  degHrNegYr <- tmp$data
  
  tmp <- aggregateSequence(mi,degreeHrExcess ,action='sum')  #total monthly daylength
  degHrExMo <- t(tmp$data)
  
  tmp <- aggregateSequence(mi, degreeHrPos ,action='sum')  #total monthly daylength
  degHrPosMo <- tmp$data
  
  list(degHrExYr = degHrExYr, degHrPosYr = degHrPosYr, degHrNegYr = degHrNegYr,
       degHrExMo = degHrExMo, degHrPosMo = degHrPosMo,
       degreeHrExcess = degreeHrExcess, degreeHrPos = degreeHrPos, degreeHrNeg = degreeHrNeg)
}

##############################################################################3
plotstart <- function(plotfile,REMOTE=F){
  
  # initiates plot
  # postscript plotfile ends in .ps
  # pdf plotfile ends in .pdf
  
  
  cc  <-NULL
  PDF <- F
  if(length(grep('.ps',plotfile)) > 0) cc <- 1
  if(length(grep('.pdf',plotfile)) > 0)cc <- 2
  
  if(is.null(cc)){
    warning('plotfile must end in .ps or .pdf')
    return()
  }
  if(cc == 2)PDF <- T
  
  if(!REMOTE)graphics.off()
  if(REMOTE & !PDF)postscript(file=plotfile,width=6, height=9,horizontal=FALSE)
}

####################################################

plotend <- function(plotfile,REMOTE=F){
  
  # terminates plot
  # postscript plotfile ends in .ps
  # pdf plotfile ends in .pdf
  
  cc  <- NULL
  PDF <- F
  if(length(grep('.ps',plotfile)) > 0) cc <- 1
  if(length(grep('.pdf',plotfile)) > 0)cc <- 2
  
  if(is.null(cc)){
    warning('plotfile must end in .ps or .pdf')
    return()
  }
  if(cc == 2)PDF <- T
  
  
  if (REMOTE)dev.off()
  if(!REMOTE)dev.print(device=postscript,file=plotfile,width=7, 
                       height=10, horizontal=FALSE)
  if(REMOTE & PDF)dev.copy2pdf(file=plotfile)

  
}

speciesCountSummary <- function(y,aggregate=F,TYPE='continuous'){
  
  #y is a sample by species matrix of counts
  
  S <- ncol(y)
  n <- nrow(y)
  
  ii <- list(spec = as.vector(matrix(c(1:S),n,S,byrow=T)), ss = as.vector(y))
  classBySpec <- imat <- byIndex(as.vector(y)*0+1,ii,sum)
  
  classAll <- table(y)
  
  lowestClass <- matrix(1:ncol(classBySpec),S,ncol(classBySpec),byrow=T)
  
  tmp <- classBySpec
  tmp[tmp > 0] <- 1
  tmp <- tmp*lowestClass  
  tmp[tmp == 0] <- NA
  lowestClass <- apply(tmp,1,min,na.rm=T)
  highestClass <- apply(tmp,1,max,na.rm=T)
  
  ii <- list(spec = as.vector(matrix(c(1:S),n,S,byrow=T)), ss = as.vector(y))
  imat <- byIndex(as.vector(y)*0+1,ii,sum)
  imat[imat > 1] <- 1
  isum <- rowSums(imat)
  novary <- which(isum < 2)

  if(length(novary) > 0 & !TYPE %in% c('continuous', 'composition') ){
    warning( paste('novariation in species: ',names(isum)[novary],
                                      sep='') )
  }
                                                                          
  #aggregate spp that don't vary
  
  if(aggregate){
    novary <- 100
    while(length(novary) > 0){
      ii <- list(spec = as.vector(matrix(c(1:S),n,S,byrow=T)), ss = as.vector(y))
      classBySpec <- imat <- byIndex(as.vector(y)*0+1,ii,sum)
      imat[imat > 1] <- 1
      isum <- rowSums(imat)
      novary <- which(isum < 2)
      if(length(novary) == 1){
        y <- y[,-novary]
      }
      if(length(novary) > 1){
        other <- rowSums(y[,novary])
        y <- cbind(y[,-novary],other)
      }
      S <- ncol(y)             
      snames <- colnames(y)
    }
  }
  print(S)
  
  list(y = y, classBySpec = classBySpec, classAll = classAll, 
       lowestClass = lowestClass, highestClass = highestClass)
}


UJSDM_updateBetaNoPrior <- function(wix,ixx,sg,...){
  
  myrmvnorm(1,as.vector(wix),kronecker(sg,ixx))
}

UJSDM_updateBetaPrior <- function(wix,ixx,sg,alpha,loBeta,hiBeta){
  
  tmp <- tnorm.mvtRcpp(avec=as.vector(alpha),muvec=as.vector(wix),smat=kronecker(sg,ixx),
                       lo=loBeta,hi=hiBeta)
  tmp[!is.finite(tmp)] <- alpha[!is.finite(tmp)]
  tmp
}

UJSDM_specBreakI <- function(z){
  
  S <- ncol(z)
  n <- nrow(z)
  i1     <- cbind(rep(c(1:S),each=n),as.vector(z))  
  i2     <- i1
  i2[,2] <- i2[,2] + 1
  
  list(i1 = i1, i2 = i2)
}

UJSDM_setupDISC <- function(y,z,breaks,maxBreaks,TYPE){
  
  n <- nrow(y)
  S <- ncol(y)
  
  if(is.null(breaks))breaks <- c(-Inf,seq(0,(max(y)-1)),Inf)
  if(length(breaks) > 50){
    warning('breaks created')
    breaks <- round( seq(min(y),(max(y) - 1),length=maxBreaks), 0)
    breaks <- c(-Inf,breaks,Inf)
  }
  
  z      <- matrix( findInterval(y-.001,breaks),n,S)
  
  ncut <- length(breaks)
  cuts <- matrix(breaks,S,ncut,byrow=T)
  
  w <- UJSDM_Y2W(y,breaks,cuts,lohi=10)  
  
  tmp <- speciesCountSummary(z,aggregate=F,TYPE=TYPE)
  if(min(z) == 0)z <- z + 1
  classBySpec  <- tmp$classBySpec
  classAll     <- tmp$classAll
  lowestClass  <- tmp$lowestClass
  highestClass <- tmp$highestClass
  
  list(w = w, z = z, breaks = breaks, cuts = cuts, ncut = ncut, 
       classBySpec = classBySpec, classAll = classAll, 
       lowestClass = lowestClass, highestClass = highestClass)
}

UJSDM_setupORD <- function(y,z,breaks,maxBreaks,TYPE){
  
  n <- nrow(y)
  S <- ncol(y)
  
  if(min(y) == 0)y <- y + 1    #zeros are first class
  breaks <- c(-Inf,c(0:(max(y)-2)),Inf)
  z <- y
  
  ncut <- length(breaks)
  cuts <- matrix(breaks,S,ncut,byrow=T)
  
  w <- UJSDM_Y2W(y,breaks,cuts,lohi=1)  
  
  tmp <- speciesCountSummary(z,aggregate=F,TYPE=TYPE)
  if(min(z) == 0)z <- z + 1
  classBySpec  <- tmp$classBySpec
  classAll     <- tmp$classAll
  lowestClass  <- tmp$lowestClass
  highestClass <- tmp$highestClass
  
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
  
  list(w = w, z = z, breaks = breaks, cuts = cuts, ncut = ncut, 
       classBySpec = classBySpec, classAll = classAll, 
       lowestClass = lowestClass, highestClass = highestClass,
       fixTheta = fixTheta, infTheta = infTheta)
}

UJSDM_setupCOMP <- function(y,z,breaks,TYPE){
  
  n <- nrow(y)
  S <- ncol(y)
  
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
  
  
  list(w = w, z = z, breaks = breaks, cuts = cuts, ncut = ncut, 
       classBySpec = classBySpec)
}

UJSDM_setupCONT <- function(y,breaks,TYPE){
  
  n <- nrow(y)
  S <- ncol(y)
  w <- z <- y
  
  loSpec <- hiSpec <- numeric(0)
  
  maxy <- max(y)
  
  #censored values set to 2nd to last cut
  if(!is.null(breaks) & !is.matrix(breaks)){
    breaks <- matrix(breaks,n,length(breaks),byrow=T)
  }
  if(is.null(breaks)){
    breaks <- matrix(c(-Inf,0,Inf),n,3,byrow=T)
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
  
  z[y == 0] <- 1            #classes defined for discrete/continuous
  z[y > 0]  <- 2
  z[y < 0]  <- 2                       #negative values must be predictors
  if(ncol(breaks) > 3){                #censored class
    z3 <- matrix(breaks[,3],n,S,byrow=T)
    z[y > z3] <- 3
  }
  
  tmp <- speciesCountSummary(z,aggregate=F)
  classBySpec <- tmp$classBySpec
  
  list(w = w, z = z, breaks = breaks, cuts = cuts, ncut = ncut, 
       classBySpec = classBySpec,wy0 = wy0,wy1 = wy1, wy3 = wy3, 
       loSpec = loSpec, hiSpec = hiSpec)
}



UJSDM_gibbs <- function(x,y,breaks=NULL,holdoutN=0,holdoutIndex=numeric(0),
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
  # holdoutN   - number to predict out-of-sample
  # holdoutIndex - which to predict out-of-sample
  
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
  
  if(ncol(x) > 2){
    tmp <- checkDesign(x)
    print(tmp)
    if(tmp$rank < tmp$p)stop( 'x not full rank' )
  }
  maxy <- max(y)
  
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
    xx     <- matrix( unlist(strsplit(betaPrior$posBeta,'_')),ncol=2,byrow=T )
    loBeta[ cbind( match(xx[,2],xnames),match(xx[,1],snames) ) ] <- 0
    loBeta <- as.vector(loBeta)
    
    hiBeta <- matrix(Inf,q,S)
    xx     <- matrix( unlist(strsplit(betaPrior$negBeta,'_')),ncol=2,byrow=T )
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
    
    UJSDM_updateW <- UJSDM_updateWord
    w             <- UJSDM_Y2W(y,breaks,cuts,lohi=1)  
    
    tmp <- speciesCountSummary(z,aggregate=F,TYPE=TYPE)
    classBySpec <- tmp$classBySpec
  }  
  
  if(TYPE == 'ordinal'){
    
    UJSDM_updateW <- UJSDM_updateWord
    ZEROINFL <- F
    
    tmp <- UJSDM_setupORD(y,z,breaks,maxBreaks,TYPE)
               w <- tmp$w
               z <- tmp$z
          breaks <- tmp$breaks
            cuts <- tmp$cuts
            ncut <- tmp$ncut
     classBySpec <- tmp$classBySpec
        classAll <- tmp$classAll
     lowestClass <- tmp$lowestClass
    highestClass <- tmp$highestClass
        fixTheta <- tmp$fixTheta
        infTheta <- tmp$infTheta
  }
  
  if(TYPE == 'discrete'){
    
    UJSDM_updateW <- UJSDM_updateWdiscrete
    
    tmp <- UJSDM_setupDISC(y,z,breaks,maxBreaks,TYPE)
               w <- tmp$w
               z <- tmp$z
          breaks <- tmp$breaks
            cuts <- tmp$cuts
            ncut <- tmp$ncut
     classBySpec <- tmp$classBySpec
        classAll <- tmp$classAll
     lowestClass <- tmp$lowestClass
    highestClass <- tmp$highestClass
  }
  
  if(TYPE == 'composition'){
    
    UJSDM_updateW <- UJSDM_updateWcomp
    
    tmp <- UJSDM_setupCOMP(y,z,breaks,TYPE)
              w <- tmp$w
              z <- tmp$z
         breaks <- tmp$breaks
           cuts <- tmp$cuts
           ncut <- tmp$ncut
    classBySpec <- tmp$classBySpec
  }
  
  if(TYPE == 'continuous'){
    
    UJSDM_updateW <- UJSDM_updateWcont
    
    tmp <-  UJSDM_setupCONT(y,breaks,TYPE)
              w <- tmp$w
              z <- tmp$z
         breaks <- tmp$breaks
           cuts <- tmp$cuts
           ncut <- tmp$ncut
    classBySpec <- tmp$classBySpec
            wy0 <- tmp$wy0
            wy1 <- tmp$wy1
            wy3 <- tmp$wy3
         loSpec <- tmp$loSpec
         hiSpec <- tmp$hiSpec
  }
  
  tmp <- UJSDM_specBreakI(z)
  i1  <- tmp$i1  
  i2  <- tmp$i2
  
  #observed zeros
  y0 <- y
  y0[y > 0] <- 1
  y0 <- 1 - y0
  
  b <- y*0
  thetag  <- 1                              # initial Pr fail to observe when present
  wabs   <- which(w <= 0,arr.ind=T)         # if not ZEROINF, absent where not seen
  wabs   <- wabs[wabs[,2] %in% compCols,]
  wpres  <- numeric(0)                       # present but not observed
  
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
  
  if(length(holdoutIndex) > 0)holdoutN <- length(holdoutIndex)
  
  HOLDOUT <- F
  allSamples <- c(1:n)
  if(holdoutN > 0){
    HOLDOUT <- T
    if(length(holdoutIndex) == 0)holdoutIndex <- sort( sample(n,holdoutN) )
    allSamples <- allSamples[-holdoutIndex]
  }
  
  nall <- length(allSamples)
    
  
  clim <- Inf
  
  cnames <- paste('C',1:ncut,sep='-')
  
  XX  <- crossprod(x[allSamples,])
  IXX <- solve(XX)
  WX  <- crossprod(x[allSamples,],w[allSamples,])
  WIX <- IXX%*%WX
  
  tg       <- cutg <- cuts
  sg       <- nugget*100 + diag(.1,S)
  sigmaDf  <- nall - q + S - 1
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
  
  
  priorXIV <- diag(1e-5,ncol(x))
  priorX       <- colMeans(x)
  predx <- predx2    <- x*0
  
  ypred  <- ypred2 <- sumb <- y*0
  sumDev <- 0   #for DIC
  sMean  <- sg*0
  ntot   <- 0
  predx  <- predx2 <- x*0
  
  pbar <- txtProgressBar(min=1,max=ng,style=1)
  
  for(g in 1:ng){
    
    WX  <- crossprod(x[allSamples,],w[allSamples,])                                #covariance
    WIX <- IXX%*%WX
    SS  <- crossprod(w) - t(WX)%*%WIX + nugget
    SI  <- solve(SS)
    sinv <- rwish(sigmaDf,SI)
    sg  <- solve( sinv )
    
    alpha <- bg <- matrix( updateBeta(WIX,IXX,sg,alpha,loBeta,hiBeta),q,S )
    muw   <- x%*%bg
    
    if(TYPE == 'ordinal' | TYPE == 'presenceAbsence') bg    <- UJSDM_A2B(alpha,sg)  
    
    if(TYPE == 'ordinal'){
      tg    <- UJSDM_updateTheta(w,z,tg,fixTheta,infTheta,lohi=100)
      cutg  <- UJSDM_theta2cuts(tg,sg)                 # correlation scale
      cgibbs[g,] <- cutg[,-c(1,2,ncut)]
    }
    
    tmp <- UJSDM_updateW(w,muw,sigma = sg,y = y,
                         b = b,theta = tg,loHi = 100,
                         i1 = i1, i2 = i2,
                         wabs = wabs, wfill = wfill, wy1 = wy1, wy3 = wy3, 
                         noZeros = noZeros,
                         loSpec = loSpec, hiSpec = hiSpec,compCols=compCols)
    w  <- tmp$w
    yp <- tmp$yp
    
    
    if(ZEROINFL){
      thetag <- updateZero(thetaPrior,y0,b)
      b      <- UJSDM_sampleB(muw,w,y,sg,thetag)    # latent indicator of absence
      wabs   <- which(b == 1,arr.ind=T)             # absent -- required for sampleW
      wpres  <- which(b == 0 & y == 0,arr.ind=T)    # present but not observed
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
      ntot   <- ntot + 1
      ypred  <- ypred + yp
      ypred2 <- ypred2 + yp^2
      sumDev <- sumDev + sum(deviance(w,x,alpha,s=sg,LIKE='mvnorm'))
      sMean  <- sMean + sg
      
      tmp <- predictY2X_linear(x,yy=w,bg,sg, 
                               priorIV = priorXIV, 
                               priorX=priorX)
      predx  <- predx + tmp
      predx2 <- predx2 + tmp^2
    }
    
    setTxtProgressBar(pbar,g)
  }
  
  sMean <- sMean/ntot
  
  chains <- list( rgibbs = rgibbs, lgibbs = lgibbs, 
                  sgibbs = sgibbs, bgibbs = bgibbs) 
  
  xpredMu <- predx/ntot
  xpredSd <- sqrt(predx2/ntot - xpredMu^2)
  
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
  
  score <- mean( getScoreNorm(y,yMu,ySd^2),na.rm=T )  # gaussian w
  
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
  
  colnames(bMu)   <- colnames(bSe) <-
    rownames(rMu) <- colnames(rMu) <- rownames(rSe) <- 
    colnames(rSe) <- snames
  rownames(bMu)   <- rownames(bSe) <- xnames
  
  absProb <- sumb/ng
  
  list(classBySpec = classBySpec, TYPE = TYPE, absProb = absProb,
       chains = chains, breaks = breaks, x = x, y = y,
       cutMu = cMu, cutSe = cSe, betaMu = bMu, betaSe = bSe, 
       corMu = rMu, corSe = rSe, sigMu = sMu, sigSe = sSe,
       wcorMu = lMu, wcorSe = lSe, betaCov = betaCov, xpredMu = xpredMu,xpredSd = xpredSd,
       yMu = yMu, ySd = ySd, prAbs = prAbs, compCols = compCols,
       DIC = DIC, score = score, holdoutIndex = holdoutIndex)
}


UJSDMtraits_gibbs <- function(x,y,breaks=NULL,holdoutN=0,holdoutIndex=numeric(0),
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
    
    UJSDM_updateW <- UJSDM_updateWord
    w             <- UJSDM_Y2W(y,breaks,cuts,lohi=1)  
    
    tmp <- speciesCountSummary(z,aggregate=F,TYPE=TYPE)
    classBySpec <- tmp$classBySpec
  }  
  
  if(TYPE == 'ordinal'){
    
    UJSDM_updateW <- UJSDM_updateWord
    ZEROINFL <- F
    
    tmp <- UJSDM_setupORD(y,z,breaks,maxBreaks,TYPE)
    w <- tmp$w
    z <- tmp$z
    breaks <- tmp$breaks
    cuts <- tmp$cuts
    ncut <- tmp$ncut
    classBySpec <- tmp$classBySpec
    classAll <- tmp$classAll
    lowestClass <- tmp$lowestClass
    highestClass <- tmp$highestClass
    fixTheta <- tmp$fixTheta
    infTheta <- tmp$infTheta
  }
  
  if(TYPE == 'discrete'){
    
    UJSDM_updateW <- UJSDM_updateWdiscrete
    
    tmp <- UJSDM_setupDISC(y,z,breaks,maxBreaks,TYPE)
    w <- tmp$w
    z <- tmp$z
    breaks <- tmp$breaks
    cuts <- tmp$cuts
    ncut <- tmp$ncut
    classBySpec <- tmp$classBySpec
    classAll <- tmp$classAll
    lowestClass <- tmp$lowestClass
    highestClass <- tmp$highestClass
  }
  
  if(TYPE == 'composition'){
    
    UJSDM_updateW <- UJSDM_updateWcomp
    
    tmp <- UJSDM_setupCOMP(y,z,breaks,TYPE)
    w <- tmp$w
    z <- tmp$z
    breaks <- tmp$breaks
    cuts <- tmp$cuts
    ncut <- tmp$ncut
    classBySpec <- tmp$classBySpec
  }
  
  if(TYPE == 'continuous'){
    
    UJSDM_updateW <- UJSDM_updateWcont
    
    tmp <-  UJSDM_setupCONT(y,breaks,TYPE)
    w <- tmp$w
    z <- tmp$z
    breaks <- tmp$breaks
    cuts <- tmp$cuts
    ncut <- tmp$ncut
    classBySpec <- tmp$classBySpec
    wy0 <- tmp$wy0
    wy1 <- tmp$wy1
    wy3 <- tmp$wy3
    loSpec <- tmp$loSpec
    hiSpec <- tmp$hiSpec
  }
  
  tmp <- UJSDM_specBreakI(z)
  i1  <- tmp$i1  
  i2  <- tmp$i2
  
  #observed zeros
  y0 <- y
  y0[y > 0] <- 1
  y0 <- 1 - y0
  
  b <- y*0
  thetag  <- 1                              # initial Pr fail to observe when present
  wabs   <- which(w <= 0,arr.ind=T)         # if not ZEROINF, absent where not seen
  wabs   <- wabs[wabs[,2] %in% compCols,]
  wpres  <- numeric(0)                       # present but not observed
  
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
  

  if(length(holdoutIndex) > 0)holdoutN <- length(holdoutIndex)
  
  HOLDOUT <- F
  allSamples <- c(1:n)
  if(holdoutN > 0){
    HOLDOUT <- T
    if(length(holdoutIndex) == 0)holdoutIndex <- sort( sample(n,holdoutN) )
    allSamples <- allSamples[-holdoutIndex]
  }
  nall <- length(allSamples)
  
  clim <- Inf
  
  cnames <- paste('C',1:ncut,sep='-')
  
  XX  <- crossprod(x[allSamples,])
  IXX <- solve(XX)
  WX  <- crossprod(x[allSamples,],w[allSamples,])
  WIX <- IXX%*%WX
  
  tg       <- cutg <- cuts
  sg       <- nugget*100 + diag(.1,S)
  sigmaDf  <- nall - q + S - 1
  alpha    <- matrix(0,q,S)
  
  alpha <- matrix( updateBeta(WIX,IXX,sg,alpha,loBeta,hiBeta),q,S )
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
  lgibbs <- matrix(0,ng,M*M)
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
  
  priorXIV <- diag(1e-5,ncol(x))
  priorX       <- colMeans(x)
  predx <- predx2    <- x*0
  
  #temporary
  
  agibbs <- matrix(0,ng,M*q)
  mgibbs <- matrix(0,ng,M*M)
  tpred <- tpred2 <- traitMuAll*0
  
  colnames(agibbs) <- as.vector( t(outer(tnames,xnames,paste,sep='_')) )
  colnames(mgibbs) <- as.vector( t(outer(tnames,tnames,paste,sep='_')) )
  
  pbar <- txtProgressBar(min=1,max=ng,style=1)
  
  for(g in 1:ng){
    
    WX  <- crossprod(x[allSamples,],w[allSamples,])             #covariance
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
    
    tmp <- UJSDM_updateW(w,muw,sigma = sg,y = y,
                         b = b,theta = tg,loHi = 100,
                         i1 = i1, i2 = i2,
                         wabs = wabs, wfill = wfill, wy1 = wy1, wy3 = wy3, 
                         noZeros = noZeros,
                         loSpec = loSpec, hiSpec = hiSpec,compCols=compCols)
    w  <- tmp$w
    yp <- tmp$yp
    
    
    Atrait <- bg%*%t(traitMat)
    Strait <- traitMat%*%sg%*%t(traitMat)
    Ttrait <- yp%*%t(traitMat)
  
    
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
    lgibbs[g,] <- cor(Ttrait)
    sgibbs[g,] <- sg
    
    agibbs[g,] <- Atrait
    mgibbs[g,] <- Strait
    
    if(g > burnin){
      ntot   <- ntot + 1
      ypred  <- ypred + yp
      ypred2 <- ypred2 + yp^2
      sumDev <- sumDev + sum(deviance(w,x,alpha,s=sg,LIKE='mvnorm'))
      sMean  <- sMean + sg
      
      
      tpred  <- tpred + Ttrait
      tpred2 <- tpred2 + Ttrait^2
      
   #   yy <- w
   #   yy[yy < 0] <- 0
   #   yy <- yy%*%t(traitMat)
      
      tmp <- predictY2X_linear(x,yy=w,bg,sg, 
                               priorIV = priorXIV, 
                               priorX=priorX)
      predx  <- predx + tmp
      predx2 <- predx2 + tmp^2
      
    }
    
    setTxtProgressBar(pbar,g)
  }
  
  sMean <- sMean/ntot
  
  chains <- list( rgibbs = rgibbs, lgibbs = lgibbs, 
                  sgibbs = sgibbs, bgibbs = bgibbs,
                  agibbs = agibbs, mgibbs = mgibbs) 
  
  xpredMu <- predx/ntot
  xpredSd <- sqrt(predx2/ntot - xpredMu^2)
  
  tmp <- processPars(agibbs)$summary
  atMu <- matrix(tmp[,'estimate'],q,M)
  atSe <- matrix(tmp[,'se'],q,M)
  
  tmp <- processPars(mgibbs)$summary
  mMu <- matrix(tmp[,'estimate'],M,M)
  mSe <- matrix(tmp[,'se'],M,M)
  rownames(mMu) <- colnames(mMu) <- tnames
  
  tMu <- tpred/ntot
  tSd <- sqrt(tpred2/ntot - tMu^2)
  
  
  tmp <- processPars(bgibbs)$summary
  bMu <- matrix(tmp[,'estimate'],q,S)
  bSe <- matrix(tmp[,'se'],q,S)
  # aMu <- bMu
  
  yMu <- ypred/ntot
  ySd <- sqrt(ypred2/ntot - yMu^2)
  cMu <- cuts
  cSe <- numeric(0)
  
  # note: on latent w scale
  meanDev <- sumDev/ntot
  pd  <- meanDev - sum(deviance(yMu,x,bMu,sMean,LIKE='mvnorm'))
  DIC <- 2*pd + meanDev
  
  score <- mean( getScoreNorm(y,yMu,ySd^2),na.rm=T )  # gaussian w
  
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
  lMu <- matrix(tmp[,'estimate'],M,M)
  lSe <- matrix(tmp[,'se'],M,M)
  
  tmp <- processPars(sgibbs)$summary
  sMu <- matrix(tmp[,'estimate'],S,S)
  sSe <- matrix(tmp[,'se'],S,S)
  
  
  betaCov <- t(bMu)%*%cov(x)%*%bMu
  
  prAbs <- numeric(0)
  if(ZEROINFL)prAbs <- sumb/ntot
  
  colnames(bMu)   <- colnames(bSe) <-
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
       DIC = DIC, score = score, xpredMu = xpredMu,xpredSd = xpredSd,
       atMu = atMu, mMu = mMu, mSe = mSe, tMu = tMu, tSd = tSd, 
       holdoutIndex = holdoutIndex)
}

predictY2X_linear <- function(xx,yy,bb,ss, 
                              priorIV = diag(1e-10,ncol(xx)), 
                              priorX=matrix(0,ncol(xx)),
                              predCols=c(2:ncol(xx))){
  
  #inverse prediction for multivariate linear in x
  
  nn <- nrow(yy)
  notPred <- c(1:ncol(xx))[-predCols]

 # xb <- xx%*%bb
  
  bn <- matrix(bb[notPred,],length(notPred))
  bp <- matrix(bb[predCols,],length(predCols))
  
  yi <- yy - xx[,notPred]%*%bn
  pp <- ncol(xx) - 1
  
  bs <- bp%*%ss
  
  V <- solve( bs%*%t(bp) + priorIV[predCols,predCols] )
  v <- yi%*%t(bs) + matrix(priorIV[predCols,predCols]%*%priorX[predCols],nn,pp,byrow=T)
  mu <- v%*%V
  xn <- myrmvnorm(nn,mu,V)
  xx[,predCols] <- xn
  xx
}



UJSDM_w2y <- function(w,sigma){
  
  S <- nrow(sigma)
  n <- nrow(w)
  w/matrix( sqrt(diag(sigma)),n,S,byrow=T)
}

imputX_MVN <- function(xx,yy,beta,wmiss,sinv,xprior=0){
  
  wcol <- unique(wmiss[,2])
  S    <- nrow(sinv)
  
  for(j in wcol){
    
    wj <- wmiss[wmiss[,2] == j,]
    if(!is.matrix(wj))wj <- matrix(wj,1,2)
    wr <- wj[,1]
    wc <- wj[1,2]
    xp <- xprior[wmiss[,2] == j]
    
    z <- yy[wr,] - xx[wr,-wc]%*%beta[-wc,]
    V <- 1/( beta[wc,]%*%sinv%*%beta[wc,] + 1 )
    v <-  rowSums( z * matrix( beta[wc,]%*%sinv,length(wr),S,byrow=T) + xp/1)
    
    xx[wj] <- rnorm(length(wr),v*V,sqrt(V))
  }
  xx
}
  
  
UJSDM_updateWcont <- function(w,muw,sigma,y,
                              b,wabs,wfill,noZeros,
                              theta,wy1,wy3,loSpec,hiSpec,...){
  
  S  <- nrow(sigma)
  n  <- nrow(w)
  
  minSigma <- diag(sigma)*.00001
  
  svec <- unique(c(loSpec,hiSpec))
  
  WHI <- F
  if(length(wy3) > 0)WHI <- T
  
  for(s in svec){
    
    w1 <- which(wy1[,2] == s)
    
    tmp <- conditionalMVNVecRcpp(w,muw,sigma,cdex=s)
    
    mu <- tmp$mu
    vr <- tmp$vr
    if(vr < minSigma[s])vr <- minSigma[s]
    lo <- hi <- rep(0,n)
    lo[wy1[w1,1]] <- -Inf
    ws <- wy1[w1,1]
    if(length(ws) == 0)next
    if(WHI){
      w3 <- which(wy3[,2] == s)
      ws <- c(wy1[w1,1],wy3[w3,1])
      #     lo[wy3[w3,1]] <- theta[1,3]
      lo[wy3[w3,1]] <- theta[3]
      hi[wy3[w3,1]] <- Inf
    }
    w[ws,s] <- tnorm(length(ws),lo[ws],hi[ws],mu[ws],sqrt(vr))
  }

  yp <- w
  yp[yp < 0] <- 0
  
  list(w = w, yp = yp)
}


UJSDM_updateWdiscrete <- function(w,xb,sigma,y,
                                  b,wabs,wfill,noZeros,
                                  theta,i1,i2,loHi= 100,...){
  
  S    <- nrow(sigma)
  n    <- nrow(w)
  ncut <- ncol(theta)
  lo <- matrix( theta[i1],n,S)
  hi <- matrix( theta[i2],n,S)
  mt <- max(theta[,-ncut]) + loHi
  
  lo[lo < -loHi] <- -loHi
  hi[hi > mt]  <- mt
  
  for(s in 1:S){
    tmp <- conditionalMVNVec(w,xb,sigma,cdex=s)
    mu <- tmp$mu
    vr <- tmp$vr
    w[,s] <- tnorm(n,lo[,s],hi[,s],mu,sqrt(vr))
  }

  yp <- w
  yp[yp < 0] <- 0
  
  list(w = w, yp = yp)
}

UJSDM_updateWord <- function(w,xb,sigma,y,
                             b,wabs,wfill,noZeros,
                             theta,i1,i2,loHi= 100,...){
  
  S    <- nrow(sigma)
  n    <- nrow(w)
  ncut <- ncol(theta)
  lo <- matrix( theta[i1],n,S)
  hi <- matrix( theta[i2],n,S)
  mt <- max(theta[,-ncut]) + loHi
  
  lo[lo < -loHi] <- -loHi
  hi[hi > mt]  <- mt
  
  for(s in 1:S){
    tmp <- conditionalMVNVec(w,xb,sigma,cdex=s)
    mu <- tmp$mu
    vr <- tmp$vr
    w[,s] <- tnorm(n,lo[,s],hi[,s],mu,sqrt(vr))
  }
  
  yp <- w
 # yp[yp < 0] <- 0
  yp <- UJSDM_w2y(yp,sigma)
  
  list(w = w, yp = yp)
}



UJSDM_theta2cuts <- function(theta,sigma){
  ncut <- ncol(theta)
  S    <- nrow(sigma)
  theta/matrix( sqrt(diag(sigma)),S,ncut)
}

UJSDM_updateTheta <- function(w,z,theta,fixTheta,infTheta,lohi=20){ 
  
  #cut points for ordinal data
  tiny <- 1e-4
  ncut  <- ncol(theta)
  nc    <- ncut - 2
  
  theta[,1]    <- theta[,2] - lohi
  theta[,ncut] <- theta[,ncut-1] + lohi
  
  for(k in 2:nc){
    
    w0 <- which(z == k,arr.ind=T)
    w1 <- which(z == (k+1),arr.ind=T)
    
    lo <- byIndex(w[w0],w0[,2],max)
    hi <- byIndex(w[w1],w1[,2],min)
    
    ww <- intersect(names(lo),names(hi))
    
    lot <- lo[ww]
    hit <- hi[ww]
    
    hit[hit < lot] <- lot[hit < lot] + tiny
    
    theta[as.numeric(ww),k+1] <- runif(length(ww),lot,hit)
  }
  theta[,1]    <- -Inf
  theta[,2]    <- 0
  theta[,ncut] <- Inf
  
  if(length(fixTheta) > 0)theta[fixTheta] <- fixTheta[,2] - 2
  if(length(infTheta) > 0)theta[infTheta] <- Inf
  
  theta
}

UJSDM_A2B <- function(alpha,sigma){
  S <- nrow(sigma)
  q <- nrow(alpha)
  alpha/matrix( sqrt(diag(sigma)),q,S,byrow=T)
}


UJSDM_B2A <- function(beta,sigma){
  S <- nrow(sigma)
  q <- nrow(beta)
  beta*matrix( sqrt(diag(sigma)),q,S,byrow=T)
}

UJSDM_Y2W <- function(yy,breaks,cuty,lohi=10){
  
  #lohi - width for lowest and highest
  
  ncut <- ncol(cuty)
  
  if(min(yy) > 0)yy <- yy - 1
  
  cuty[,1]   <- cuty[,2] - lohi
  cuty[,ncut] <- cuty[,ncut-1] + lohi
  
  w <- yy*0
  
  for(k in 1:(ncut-1)){
 #   wk <- which(yy == k,arr.ind=T)
    wk <- which(yy > breaks[k] & yy <= breaks[k+1],arr.ind=T)
    lo <- cuty[wk[,2],k]
    hi <- cuty[wk[,2],k+1]
    w[wk] <- runif(nrow(wk),lo,hi)
  }
  w
}

UJSDM_Z2Y <- function(z,cuts){
  
  y  <- z*0 + 1
  nc <- ncol(cuts) - 1
  for(k in 2:nc)y[z > cuts[,k]] <- k
  y
}


shadeInterval <- function(xvalues,loHi,col='grey',PLOT=T,add=T,xlab=' ',ylab=' '){
  
  #draw shaded interval
  
  tmp <- smooth.na(xvalues,loHi)
  xvalues <- tmp[,1]
  loHi <- tmp[,-1]
  
  xbound <- c(xvalues,rev(xvalues))
  ybound <- c(loHi[,1],rev(loHi[,2]))
  if(!add)plot(xvalues,loHi[,1]*0,cex=.01,ylim=c(range(loHi,na.rm=T)),
               xlab=xlab,ylab=ylab)
  if(PLOT)polygon(xbound,ybound,border=NA,col=col)
  
  invisible(cbind(xbound,ybound))
  
}


UJSDM_simData <- function(n=1000,S=10,q=5,cuts=c(-Inf,0,Inf),nmiss=0,TYPE){
  
  # zero class for multivariate data
  # TYPE  = 'ordinal', 'discrete', 'composition', 'continuous'
  # nmiss = number of missing values in x
  
  x <- matrix( runif(n*q,0,3), n,q)  #design matrix (insert your own predictors)
  x[,1] <- 1
  
  rb   <- c(-2,2)
  ss    <- diag(.02,S)    
  if(TYPE == 'composition'){
    rb <- c(-.1,.3)
    ss <- diag(.001,S)
  }
  if(TYPE == 'ordinal'){
    ss <- diag(1,S)
    rb   <- c(-.4,2)
  }
  if(TYPE == 'discrete'){
    ss <- diag(1,S)
    rb   <- c(-3,5)
  }
  if(TYPE == 'continuous'){
    ss <- diag(.05,S)
    rb   <- c(-1,10)
  }
  if(TYPE == 'presenceAbsence'){
    ss <- diag(1,S)
    rb   <- c(-.5,1)
  }
  
  beta <- matrix( runif(S*q,rb[1],rb[2]),q,S ) #coefficient matrix: predictors by species
  
  sigma <- cov( myrmvnorm(5,0,ss) )        # covariance of sample
  if(TYPE == 'ordinal' | TYPE == 'presenceAbsence')sigma <- cov2cor(sigma)
  z     <- x%*%beta + myrmvnorm(n,0,sigma) # simulate observations
  
  
  xnames <- paste('X',1:q,sep='-')
  snames <- paste('S',1:S,sep='-')   
  w <- ybreaks <- NULL
  ybreaks <- cuts <- c(-Inf,0,Inf)
  
  CONT <- F
  
  if(TYPE == 'presenceAbsence'){
    ybreaks < c(-Inf,0,Inf)
    ncut <- 3
    cuts  <- matrix( ybreaks,S,ncut,byrow=T)
    rownames(cuts) <- snames

    y <- UJSDM_Z2Y(z,cuts) - 1         
    w <- z
    z <- y + 1
  }
  
  if(TYPE == 'ordinal'){
    maxy <- floor(max(z))
    ncut <- maxy + 2
    cuts  <- matrix( c(-Inf, seq(0,(maxy-1),length=(ncut-2)) ,Inf),S,ncut,byrow=T)
    rownames(cuts) <- snames
    ybreaks <- cuts[1,]
    y <- UJSDM_Z2Y(z,cuts)         
    w <- z
    z <- y
  }
  
  if(TYPE == 'discrete'){
    w       <- z
    ym      <- round(max(z),0) - 1
    ybreaks <- c(-Inf,seq(0,ym),Inf)  
    cuts    <- ybreaks
    y <- findInterval(z,ybreaks)
    #y <- ybreaks[y+1]
    y <- matrix(y,n,S) - 1    #include zeros
    colnames(beta) <- snames
    z <- y + 1
    
    #approximately true parameter values
 #   beta  <- solve(crossprod(x))%*%crossprod(x,w - 1) 
 #   sigma <- crossprod(w - x%*%beta)/(n - q)
    colnames(w) <- snames
  }
  
  if(TYPE == 'composition'){
    w <- z
    z[z < 0] <- 0
    ss <- rowSums(z)
    ww <- which(ss == 0)
    if(length(ww) > 0){
      z <- z[-ww,]
      w <- w[-ww,]
      x <- x[-ww,]
      n <- nrow(x)
    }
    z <- z/matrix(rowSums(z),n,S)
    w[z > 0] <- z[z > 0]
    y <- z
    
    #approximately true parameter values
    beta  <- solve(crossprod(x))%*%crossprod(x,w)
    sigma <- crossprod(w - x%*%beta)/(n - q)
    colnames(w) <- snames
    CONT  <- T
    ybreaks <- cuts <- c(-Inf,0,1)
  }
  
  if(TYPE == 'continuous'){
    y <- z
    ybreaks <- cuts <- c(-Inf,0,max(y) - 1,Inf)
    y[y < ybreaks[2]] <- 0
    y[y > ybreaks[3]] <- ybreaks[3]
    CONT    <- T
  }
  ncut    <- length(ybreaks)
  
  if(nmiss > 0){
    x[ sample(length(x),nmiss) ] <- NA
    x[,1] <- 1
    wmiss <- which(is.na(x),arr.ind=T)
    nmiss <- nrow(wmiss)
  }
  
  colnames(z) <- colnames(y) <- colnames(beta) <- 
    rownames(sigma) <- colnames(sigma) <- snames
  colnames(x) <- rownames(beta) <- xnames
  
  list(x = x, y = y, w = w, z = z, beta = beta, sigma = sigma, 
       cuts = cuts, ybreaks = ybreaks)
}
  
UJSDM_updateZeroSample <- function(thetaPrior,y0,b){  
  
  # detection for zero inflation
  # theta - [not observed|present] = [y == 0|b == 0]
  S <- ncol(b)
  
  u1 <- rowSums( (1 - b)*y0 ) + thetaPrior[,1]  # present, not observed
  u2 <- rowSums( 1 - y0 ) + thetaPrior[,2]        # present, observed
  
  matrix( rbeta(n,u1,u2), n,S )
}

UJSDM_updateZeroSpecies <- function(thetaPrior,y0,b){  
  
  # detection for zero inflation
  # theta - [not observed|present] = [y == 0|b == 0]
  S <- ncol(b)
  
  u1 <- colSums( (1 - b)*y0 ) + thetaPrior[,1]  # present, not observed
  u2 <- colSums( 1 - y0 ) + thetaPrior[,2]      # present, observed
  
  matrix( rbeta(S,u1,u2), n,S,byrow=T )
}
UJSDM_updateZeroTG <- function(thetaPrior1,thetaPrior2,b,wpres){  
  
  # detection for zero inflation
  # theta - [not observed|present] = [y == 0|b == 0]
  S <- ncol(b)
  
  yy <- y
  yy[y > 0] <- 1
  
  u1 <- (1 - b)*(1 - yy) + thetaPrior1  # present, not observed
  u2 <- (1 - b)*yy + thetaPrior2        # present, observed
  
  matrix( rbeta(n*S,u1,u2), n,S)
  
}

UJSDM_sampleB <- function(xb,w,y,sg,theta){ #probability absent
  
  tiny <- 1e-10
  n <- nrow(w)
  S <- ncol(w)
  
  # theta - Pr present, but missed due to sampling
  
  pjs <- exp( pmvnormCondRcpp(q=0,xx=w,mu=xb,smat=sg,whichVar=c(1:S)) ) # Pr < 0
  pjs[pjs == 0] <- tiny
  pb  <- pjs/(pjs + (1 - pjs)*theta)
  b   <- matrix(rbinom(n*S,1,pb),n,S)
  b[y > 0] <- 0
  b
}

UJSDM_updateWcomp <- function(w,muw,sigma,y,
                              b,wabs,wfill,noZeros,compCols,...){  #latent states 
  
  S  <- ncol(y)
  
 # bcomp <- b[,compCols]
  
  lo <- y*0 - 10000
  hi <- y*0 + 10000
 # lo[,compCols] <- 0
 # hi[,compCols] <- 1
  
  lc <- lo[,compCols]*0
  hc <- hi[,compCols]*0 + 1
  
  lc[b[,compCols] == 1] <- -1
  hc[b[,compCols] == 1] <- 0
  
  lo[,compCols] <- lc
  hi[,compCols] <- hc
  
 # tmp <- tnorm.mvtMatrixRcpp(avec=w, muvec=muw, smat=sigma,lo=lo,hi=hi)
  
  for(s in 1:S){
    tmp <- conditionalMVNVec(w,muw,sigma,cdex=s)
    mu <- tmp$mu
    vr <- tmp$vr
    w[,s] <- tnorm(n,lo[,s],hi[,s],mu,sqrt(vr))
  }
  
  wz <- w
  wz[b == 0 & y > 0] <- y[b == 0 & y > 0]
  wz[wz < 0] <- 0

  
 # ww[wfill] <- w[wfill]   #absent or present but not observed
 # wz <- ww*(1 - b)
  sz <- rowSums(wz[,compCols])
  wz[,compCols] <- wz[,compCols]/sz
  wz[w < 0]      <- w[w < 0]
  
  yp <- w
  yp[yp < 0] <- 0
  sz <- rowSums(yp[,compCols])
  yp[,compCols] <- yp[,compCols]/sz
  
 # wz[noZeros,] <- y[noZeros,]
  
  list(w = wz, yp = yp)
}


oceanValues <- function(minz=0,minc=0,xx,yy,zz,lonLat,climVec){
  
  # adds low values to climVec so ocean is not shaded on map
  # minz - minimum elevation to define ocean
  # minc - minimum value for climate variable
  
  ww <- which(zz < minz,arr.ind=T)  #ocean
  
  if(length(ww) == 0) return( list(clim = climVec, lonLat = lonLat) )
  
  if(nrow(ww) > 5000)ww <- ww[sample(nrow(ww),5000),]
  
  zl     <- cbind(xx[ww[,1]],yy[ww[,2]])
  colnames(zl) <- colnames(lonLat)
  lonLat <- rbind(lonLat,zl)
  
  clim   <- rep(minc,nrow(ww))
  clim   <- c(climVec,clim)
  
  list(clim = clim, lonLat = lonLat, ocean = zl)
}

plotMCMC <- function(out, outfile=character(0), trueValues = NULL, sdScale=F,
                     ylog=F,ncolPlot=2){
  
  # output     - list generated by UJSDM_gibbs
  # outfile    - name for output, can include path
  # trueValues - 'beta', 'y', 'cor', 'cuts'
  # sdScale    - plot posteriors on standardized scale
  
  palette('default')
  
  chainNames <- names(out$chains)
  
  bgibbs <- out$chains$bgibbs
  rgibbs <- out$chains$rgibbs
  cgibbs <- out$chains$cgibbs
  specs  <- rownames(out$classBySpec)
  DIC    <- out$DIC
  score  <- out$score
  compCols <- out$compCols
  x      <- out$x
  y      <- out$y
  TYPE   <- out$TYPE
  
  cutMu  <- out$cutMu
  cutSe  <- out$cutSe
  betaMu <- out$betaMu
  betaSe <- out$betaSe
  corMu  <- out$corMu
  corSe  <- out$corSe
  wcorMu <- out$wcorMu
  wcorSe <- out$wcorSe
  yMu    <- out$yMu
  ySd    <- out$ySd
  

  
  ncut   <- ncol(out$classBySpec) + 1
  S      <- ncol(y)
  snames <- colnames(y)
  xnames <- colnames(x)
  
  xSd <- sqrt( diag(cov(x)) )
  ySd <- sqrt( diag(cov(y)) )
  
  outfile <- paste(out$TYPE,outfile,sep='_')
  
  bCoeffTable <- processPars(bgibbs,sigOnly=T)
  rCoeffTable <- processPars(rgibbs)
  sigBeta     <- rownames(bCoeffTable[[1]])
  
  TB <- F
  if('beta' %in% names(trueValues))TB <- T
  
  if(TB){
    graphics.off()
    beta <- trueValues$beta
    bCoeffTable <- processPars(bgibbs,xtrue=as.vector(beta) )
    predVsObs(as.vector(beta),bgibbs,xlab='true',ylab='estimated',ylim=range(beta))
    title('beta')
    dev.copy2pdf(file=paste(outfile,'betaTrue.pdf',sep='_')) 
  }
  
  if('y' %in% names(trueValues)){
    graphics.off()
    yy <- trueValues$y
    if(TYPE == 'presenceAbsence'){
      plotObsPred(jitter(yy),jitter(yMu),xlab='true',ylab='predicted',
                  breaks=c(-.5,.5,1.5),wide=c(.04,.04))
    } else {
      plot(yy,yMu,xlab='true',ylab='predicted',cex=.2)
    }
    title('observations')
    abline(0,1,lty=2)
    dev.copy2pdf(file=paste(outfile,'yTrue.pdf',sep='_')) 
  }
  
  graphics.off()
  if(TYPE == 'presenceAbsence'){
    plotObsPred(jitter(y),jitter(yMu),xlabel='observed',ylabel='predicted',
                breaks=c(-.5,.5,1.5),wide=c(.04,.04))
    abline(0,1,lty=2)
  } else {
    par(mfrow=c(1,2))
    plotObsPred(y[,compCols],yMu[,compCols],nbin=10,
                ylimit=c(0,max(yMu[,compCols])),xlabel='Observed',ylabel='Predicted')
    abline(0,1,lty=2)
    title('compCols')
    if(length(compCols) < ncol(yMu)){
      ww <- c(1:S)[-compCols]
      plot(c(0,0),c(0,0),xlim=c(0,1),ylim=c(0,1),xlab='Observed (fraction)',
           ylab='Predicted (fraction)',cex=.01)
      for(j in 1:length(ww)){
        y1 <- y[,ww[j]]/max(y[,ww[j]])
        y2 <- yMu[,ww[j]]/max(y[,ww[j]])
        points(y1,y2,cex=.2,col=j)
      }
      abline(0,1,lty=2)
      title('observations')
    }
  }
  
  dev.copy2pdf(file=paste(outfile,'predictY.pdf',sep='_')) 
  
  
  if('cor' %in% names(trueValues)){
    graphics.off()
    cc <- trueValues$cor
    predVsObs(as.vector(cc),rgibbs,xlab='true',ylab='estimated')
    title('residual species correlation matrix')
    dev.copy2pdf(file=paste(outfile,'residCorTrue.pdf',sep='_'))            
  }
  
  graphics.off()
  if('lgibbs' %in% chainNames){
    lgibbs <- out$chains$lgibbs
    cc <- cor( output$y )
    #  predVsObs(as.vector(cc),lgibbs,xlab='true',ylab='estimated')
    predVsObs(as.vector(cc),lgibbs)
    title('species correlation: obs vs latent scale')
    dev.copy2pdf(file=paste(outfile,'corrObsLatent.pdf',sep='_'))   
  }
  
#  if('mgibbs' %in% chainNames){
#    mgibbs <- out$chains$mgibbs
  
  
  if(output$TYPE == 'ordinal' & 'cuts' %in% names(trueValues)){
    graphics.off()
    cc <- trueValues$cuts
    ncut <- ncol(cc)
    cvec <- as.vector(cc[,-c(1,2,ncut)])
    ww   <- which(is.finite(cgibbs[1,]))
    predVsObs(cvec[ww],cgibbs[,ww],xlab='true',ylab='estimated')
    title('breaks')
    dev.copy2pdf(file=paste(outfile,'breaksTrue.pdf',sep='_'))            
  }
  
  
  rmspe <- sqrt( mean( (y[,compCols] - yMu[,compCols])^2 ) )
  
  fit <- signif( c(DIC,score,rmspe), 5)
  names(fit) <- c('DIC','score','rmspe')
  
  # beta coefficients
  betaSigGibbs <- bgibbs[,sigBeta]
  
  if(sdScale){
    for(j in 2:ncol(x)){
      ww <- grep(xnames[j],colnames(betaSigGibbs))
      betaSigGibbs[,ww] <- betaSigGibbs[,ww]*xSd[j]
    }
    for(j in 1:ncol(y)){
      ww <- grep(snames[j],colnames(betaSigGibbs))
      betaSigGibbs[,ww] <- betaSigGibbs[,ww]/ySd[j]
    }
  }
      
  
  
  colF   <- colorRampPalette(c('black','brown','orange'))
  cols <- colF(S)
  
  
  tmp <- matrix( unlist( strsplit(colnames(betaSigGibbs),'_') ),ncol=2,byrow=T)
  snam <- tmp[,1]
  vnam <- tmp[,2]
  
  vnames <- sort(unique(vnam))
  vnames <- vnames[vnames != 'intercept']
  
  nrr <- ceiling(length(vnames)/ncolPlot)
  
  graphics.off()
  par(mfrow=c(nrr,ncolPlot),bty='n')
  
  for(k in 1:length(vnames)){
    
    tname <- vnames[k]
    
    tmp <- chains2density(betaSigGibbs,varName=vnames[k])
    xt  <- tmp$x
    yt  <- tmp$y
    xrange <- tmp$xrange
    mr     <- round( log10(diff(xrange)),0 )
    yymax  <- tmp$ymax
    chains <- tmp$chainMat
    
    snamek <- matrix( unlist(strsplit(colnames(chains),'_')),ncol=2,byrow=T)[,1]
    
    nn <- nrow(chains)
    plot(0,0,xlim=c(0,(1.4*nn)),ylim=range(chains),xlab='Iteration'
         ,ylab='Parameter values',cex=.01)
    title(tname)
    
    for(j in 1:ncol(chains)){
      lines(chains[,j],col=cols[j])
      if(ncol(chains) < 25)text(nn,chains[nn,j],snamek[j],col=cols[j],pos=4)
    }
    abline(h=0,lwd=4,col='white')
    abline(h=0,lty=2)
  }
  
  dev.copy2pdf(file=paste(outfile,'_betaChains.pdf',sep=''))
  
    
  
  graphics.off()
  par(mfrow=c(nrr,ncolPlot),bty='n')
  
  for(k in 1:length(vnames)){
    
    tname <- vnames[k]
    
    tmp <- chains2density(betaSigGibbs,varName=vnames[k])
    xt  <- tmp$x
    yt  <- tmp$y
    xrange <- tmp$xrange
    mr     <- round( log10(diff(xrange)),0 )
    yymax  <- tmp$ymax
    chains <- tmp$chainMat
    
    nn <- nrow(chains)
    
    snam <- matrix( unlist( strsplit(rownames(xt),'_') ),ncol=2,byrow=T)[,1]
    rownames(xt) <- rownames(yt) <- snam
    
    sc     <- diff(xrange)
    x1     <- xrange[1] - sc/3
    x2     <- xrange[2] + sc/2
    
    ends <- round(c(x1,x2),-mr+1)
    xtic   <- seq(ends[1],ends[2],length=5)
    
    if(ends[1] < 0 & ends[2] > 0){
      nseq <- seq(ends[1],0,length=3)
      pseq <-  seq(0,ends[2],by=diff(nseq)[1])[-1]
      xtic <- c(nseq,pseq)
    }

    if(!ylog)plot(10,10,xlim=range(xtic),ylim=c(0,1.8*max(yt)),
         xlab='Standardized parameter value',ylab='Density',cex=.01)
    
    if(ylog)plot(10,10,xlim=range(xtic),ylim=c(.00001,2*max(yt)),log='y',
                  xlab='Standardized parameter value',ylab='Density',cex=.01)
    title(tname)
    
    for(j in 1:nrow(xt)){
      
      cj <- cumsum(yt[j,])
      cj <- cj/max(cj)
      
      ci <- xt[j, findInterval(c(.02,.98),cj) ]
      
      label <- rownames(xt)[j]
  
      wm <- which.max(yt[j,])
      lines(xt[j,],yt[j,],lwd=2,col=cols[j])
      lines(range(xt[j,]),c(0,0),lwd=2,col=cols[j])
      text(xt[j,wm],1.1*yt[j,wm],label,srt=55,pos=4,col=cols[j],cex=.8)
    }
    abline(v=0,lty=2,lwd=2,col='grey')
  }
  dev.copy2pdf(file=paste(outfile,'_betaPost.pdf',sep=''))
  
  invisible( list(betaEstimates = bCoeffTable, corEstimates = rCoeffTable,
                  fit = fit) )
  
}

myMBA.surf <- function (xyz, no.X, no.Y, n = 1, m = 1, h = 8, extend = FALSE, 
                    sp = FALSE, ...) {
  
  require(MBA)
  
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  if (!extend) {
    hpts <- chull(xyz[, c(1, 2)])
    hpts <- c(hpts, hpts[1])
  }
  else {
    hpts <- NULL
  }
  x.min.dom <- min(xyz[, 1])
  x.max.dom <- max(xyz[, 1])
  y.min.dom <- min(xyz[, 2])
  y.max.dom <- max(xyz[, 2])
  if ("b.box" %in% elip.args) {
    b.box <- list(...)$b.box
    if (b.box[1] < x.min.dom) {
      x.min.dom <- b.box[1]
    }
    else {
      warning("b.box[1] set to min x value")
    }
    if (b.box[2] > x.max.dom) {
      x.max.dom <- b.box[2]
    }
    else {
      warning("b.box[2] set to max x value")
    }
    if (b.box[3] < y.min.dom) {
      y.min.dom <- b.box[3]
    }
    else {
      warning("b.box[3] set to min y value")
    }
    if (b.box[4] > y.max.dom) {
      y.max.dom <- b.box[4]
    }
    else {
      warning("b.box[4] set to max y value")
    }
  }
  xyz <- as.matrix(xyz)
  storage.mode(xyz) <- "double"
  storage.mode(x.min.dom) <- "double"
  storage.mode(x.max.dom) <- "double"
  storage.mode(y.min.dom) <- "double"
  storage.mode(y.max.dom) <- "double"
  out <- .Call("MBASurf", xyz, as.integer(no.X), as.integer(no.Y), 
               as.integer(m), as.integer(n), as.integer(h), as.integer(extend), 
               as.integer(hpts), 
               x.min.dom, x.max.dom, y.min.dom, y.max.dom)
  if (sp) {
    xy <- expand.grid(out[["x"]], out[["y"]])
    grid <- data.frame(z = matrix(out[["z"]], 
                                  as.integer(length(out[["x"]]) * length(out[["y"]])), 1), 
                       x = xy[, 1], y = xy[, 2])
    coordinates(grid) = ~x + y
    gridded(grid) <- TRUE
  } else {
    grid <- out[c("x", "y", "z")]
  }
  out <- list()
  out$xyz.est <- grid
    out$no.X <- no.X
    out$no.Y <- no.Y
    out$n <- n
    out$m <- m
    out$h <- h
    out$extend <- extend
    out$sp <- sp
    out$b.box <- c(x.min.dom, x.max.dom, y.min.dom, y.max.dom)
    out
}

getCoverageNorm <- function(meanVec,sdVec,values,fraction){
  
  #coverage for values contained in interval fraction
  #  for mean and sd
  
  q  <- -qnorm( (1 - fraction)/2 )
  lo <- meanVec - q*sdVec
  hi <- meanVec + q*sdVec
  ii <- which(values >= lo & values <= hi)
  length(ii)/length(lo)
}


fitSummary <- function(yy,mu,stdev,coverage=.95){
  
  wf <- which(stdev > 0 & is.finite(stdev))
  
  ww <- which(stdev == 0 | !is.finite(stdev))
  stdev[ww] <- min(stdev,na.rm=T)/100
  
  rmse  <-  sqrt( mean((yy - mu)^2) )
  score <-  mean( getScoreNorm(yy[wf],mu[wf],stdev[wf]^2),na.rm=T )
  cover <- getCoverageNorm(mu,stdev,yy,coverage)
  mvar  <- quantile(stdev^2,.5)
  
  fitTable <- matrix(c(rmse,score,mvar,cover),ncol=1)
  rownames(fitTable) <- c('rmse','score','median var',
                          paste(100*coverage,'% coverage',sep=''))
  signif(fitTable,3)
}

updateWishartNoPrior <- function(x,y,df,IXX=NULL){
  
  if(is.null(IXX))IXX  <- solve( crossprod(x) )
  WX  <- crossprod(x,y)
  WIX <- IXX%*%WX
  
  SS  <- crossprod(y) - t(WX)%*%WIX
  SI  <- solve(SS)
  solve( rwish(df,SI) )
}
