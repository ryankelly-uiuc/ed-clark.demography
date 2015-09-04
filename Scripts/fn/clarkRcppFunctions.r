
library(inline)
library(RcppArmadillo)

##################################################

srccvec <- '
using namespace Rcpp;
using namespace arma;

uvec cdex = as<uvec>(cdex_fR);
uvec gdex = as<uvec>(gdex_fR);

mat xx = as<mat>(xx_fR);
mat mu = as<mat>(mu_fR);
mat sigma = as<mat>(sigma_fR);


mat sinv = inv(sigma.submat(gdex,gdex));
mat p1 = sigma.submat(cdex, gdex) * sinv;
mat mu1 = mu.cols(cdex) + trans(p1 * trans(xx.cols(gdex) - mu.cols(gdex)));
mat vr1 = sigma.submat(cdex, cdex) - p1 * sigma.submat(gdex,cdex);

List z; z["mu"] = mu1; z["vr"] = vr1;
return(z);

'

conditionalMVNVec_cpp <- cxxfunction( signature(cdex_fR = "numeric", gdex_fR = "numeric", 
                                                xx_fR = "numeric", mu_fR = "numeric", 
                                                sigma_fR = "numeric" ) , 
                                                body = srccvec, plugin = "RcppArmadillo" )


###################### pmvnormCond #################

srcpc <- '
using namespace Rcpp;
using namespace arma;

double q = as<double>(q_fR);

mat xx = as<mat>(xx_fR);
mat mu = as<mat>(mu_fR);
mat sigma = as<mat>(sigma_fR);

uvec cdex = as<uvec>(cdex_fR);
umat gdex = as<umat>(gdex_fR);

int n = mu.n_rows;
int p = mu.n_cols;
int nSelVars = gdex.n_cols;
int nGivenVars = gdex.n_rows;

mat sinv(nGivenVars, nGivenVars);
mat p1(1, nGivenVars);
mat mu1(n,1);
mat vr1(1,1);
mat ps(n, nSelVars);
mat aux_row(1,p);
mat aux_col(p,1);

int i,j;


for(j = 0; j < nSelVars; j++){
sinv = inv(sigma.submat(gdex.col(j),gdex.col(j)));
aux_row = sigma.row(cdex(j));
p1 = aux_row.cols(gdex.col(j)) * sinv;

mu1 = mu.col(cdex(j)) + trans(p1 * trans(xx.cols(gdex.col(j)) - mu.cols(gdex.col(j))));
aux_col = sigma.rows(gdex.col(j));
vr1 = sigma(cdex(j), cdex(j)) - p1 * aux_col.col(cdex(j));

for(i = 0; i < n; i++)
ps(i,j) = R::pnorm(q, mu1(i,0), sqrt(vr1(0,0)),1,1);  
}


return(wrap(ps));

'

pmvnormCond_cpp <- cxxfunction( signature(q_fR = "numeric", cdex_fR = "numeric", 
                                          gdex_fR = "numeric", xx_fR = "numeric", 
                                          mu_fR = "numeric", sigma_fR = "numeric" ) , 
                                          body = srcpc, plugin = "RcppArmadillo" )


comp_outindex <- function(index_in, total) return((1:total)[-index_in]-1)

pmvnormCond_Rcpp <- function(q=0,xx,mu,smat,whichVar=c(1:nrow(smat))){  
  
  # log conditional MVN PDF, Pr(xx < q|mu,smat)
 
  pj <- mu*0  
  gindexes      <- sapply(whichVar, comp_outindex, total = ncol(xx))
  pj[,whichVar] <- pmvnormCond_cpp(q, (whichVar-1) , gindexes, xx, mu, smat)
  pj  
}

  
############################### tnorm, mvtnorm

Inc <- '

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
  using namespace arma;

double tnorm(double lo, double hi, double mu, double sig){
  double q1, q2, z;
  

  q1 = Rf_pnorm5(lo,mu,sig,1,0);
  q2 = Rf_pnorm5(hi,mu,sig,1,0);
  z = q1 + unif_rand()*(q2-q1);
  z = Rf_qnorm5(z, mu, sig, 1, 0);

  if(!is_finite(z)){
    if(z > hi)
      z = lo;
    else
      z = hi;
  }
  
  return(z);

}

void conditionalMVN(rowvec &x, rowvec &mu, mat &sigma, int cindex, vec &mAs, double tiny, int nm, umat &idxALLm){
  
  rowvec p1(nm-1);
  mat sin(nm-1, nm-1);
  uvec cid(1);
  mat m1(1,1);
  mat s1(1,1);
  
  
  cid(0) = cindex;

  uvec idx = idxALLm.col(cindex);

  sin = inv_sympd(sigma.submat(idx, idx));

  p1 = trans(sigma.submat(idx, cid)) * sin;

  m1 = mu[cindex] + dot(p1, (x.elem(idx) - mu.elem(idx)));
  s1 = sigma(cindex,cindex) - dot(p1, sigma.submat(cid, idx)) ;

  mAs[0] = m1(0,0);
  mAs[1] = s1(0,0);
  if(mAs[1] < 0) mAs[1] = tiny;  
}

' 

Src1 <- '
rowvec avec = as<rowvec>(avec_fR), muvec = as<rowvec>(muvec_fR);
uvec idxALL = as<uvec>(idxALL_fR), whichSample = as<uvec>(whichSample_fR);
mat smat = as<mat>(smat_fR);
vec lo = as<vec>(lo_fR), hi = as<vec>(hi_fR);
int cindex, times = as<int>(times_fR);
vec mAs(2);
int nm = smat.n_rows;
unsigned int i,k;
double tiny = min(smat.diag())*.0001;

umat idxALLm(nm-1, nm);
for(int j=0; j< nm; j++)
  idxALLm.col(j) = idxALL.elem( find(idxALL != j) );

for(i = 0; i < times ; i++)
  for(k = 0; k < whichSample.n_elem; k++){
    cindex = whichSample[k]-1;
    conditionalMVN(avec, muvec, smat, cindex, mAs, tiny, nm, idxALLm);
    
    avec[cindex] = tnorm(lo[cindex], hi[cindex], mAs[0], sqrt(mAs[1]));
  }

return wrap(avec);
'

Src2 <- '
rowvec avec = as<rowvec>(avec_fR), muvec = as<rowvec>(muvec_fR);
uvec idxALL = as<uvec>(idxALL_fR), whichSample = as<uvec>(whichSample_fR);
mat smat = as<mat>(smat_fR);
vec lo = as<vec>(lo_fR), hi = as<vec>(hi_fR);
int cindex, times = as<int>(times_fR);
vec mAs(2);
int nm = smat.n_rows;
int nr = times;
unsigned int i,k;
double tiny = min(smat.diag())*.0001;

mat A(nr, nm); A.fill(NA_REAL);

umat idxALLm(nm-1, nm);
for(int j=0; j< nm; j++)
  idxALLm.col(j) = idxALL.elem( find(idxALL != j) );

for(i = 0; i < times ; i++){
  for(k = 0; k < whichSample.n_elem; k++){
    cindex = whichSample[k]-1;
    conditionalMVN(avec, muvec, smat, cindex, mAs, tiny, nm, idxALLm);
    
    avec[cindex] = tnorm(lo[cindex], hi[cindex], mAs[0], sqrt(mAs[1]));
    A(i,cindex) = avec(cindex);
  }
}

return wrap(A);
'



tnorm.mvt_cpp <- cxxfunction(signature(avec_fR = "numeric", muvec_fR = "numeric", 
                                       smat_fR = "numeric", lo_fR = "numeric",
                                       hi_fR = "numeric", whichSample_fR = "numeric", 
                                       times_fR = "numeric", idxALL_fR ="numeric"), 
                                       includes = Inc, body = Src2, plugin = "RcppArmadillo")


tnorm.mvtRcpp <- function(avec, muvec, smat, lo=rep(-Inf,length(muvec)), hi=rep(Inf,length(muvec)),
                           whichSample = c(1:length(muvec)), times = 1){
  
  if(max(whichSample) > length(muvec))stop('whichSample outside length(muvec)')
  
  if(length(lo) == 1)lo <- rep(lo,length(muvec))
  if(length(hi) == 1)hi <- rep(hi,length(muvec))
  
  r <- matrix(avec,times,length(avec),byrow=T)
  
  a <- tnorm.mvt_cpp(avec, muvec, smat, lo, hi, whichSample, times, 0:(length(avec)-1))  
  
  r[,whichSample] <- a[,whichSample]
  r
  
}

############for matrix

Src3 <- '
uvec idxALL = as<uvec>(idxALL_fR), whichSample = as<uvec>(whichSample_fR);
mat smat = as<mat>(smat_fR), avec = as<mat>(avec_fR), muvec = as<mat>(muvec_fR), lo = as<mat>(lo_fR), hi = as<mat>(hi_fR);
int cindex;
rowvec av;
rowvec mv;
vec mAs(2);
int nm = smat.n_rows;
int nr = muvec.n_rows;
int i,k;
double tiny = min(smat.diag())*.0001;

mat A(nr, nm); A.fill(NA_REAL);

umat idxALLm(nm-1, nm);
for(int j=0; j< nm; j++)
  idxALLm.col(j) = idxALL.elem( find(idxALL != j) );

for(i = 0; i < nr ; i++){
  for(k = 0; k < whichSample.n_elem; k++){
    cindex = whichSample[k]-1;
    av = avec.row(i);
    mv = muvec.row(i);
    conditionalMVN(av, mv, smat, cindex, mAs, tiny, nm, idxALLm);
    
    avec(i,cindex) = tnorm(lo(i,cindex), hi(i,cindex), mAs[0], sqrt(mAs[1]));
    A(i,cindex) = avec(i,cindex);
  }
}

return wrap(A);
'


tnorm.mvtMatrix_cpp <- cxxfunction(signature(avec_fR = "numeric", muvec_fR = "numeric", 
                                       smat_fR = "numeric", lo_fR = "numeric",
                                       hi_fR = "numeric", whichSample_fR = "numeric", 
                                       idxALL_fR ="numeric"), 
                                       includes = Inc, body = Src3, plugin = "RcppArmadillo")

tnorm.mvtMatrixRcpp <- function(avec, muvec, smat, 
                                lo=matrix(-1000,nrow(muvec),ncol(muvec)), 
                                hi=matrix(1000,nrow(muvec),ncol(muvec)),
                                whichSample = c(1:nrow(smat))){
  
  #lo, hi must be same dimensions as muvec,avec
  
  lo[lo < -1000] <- -1000
  hi[hi > 1000]  <- 1000
  
  if(max(whichSample) > length(muvec))stop('whichSample outside length(muvec)')

  r <- avec

  a <- tnorm.mvtMatrix_cpp(avec, muvec, smat, lo, hi, whichSample, 0:(nrow(smat)-1))  

  r[,whichSample] <- a[,whichSample]
  r

}
############################### lagGrowth
src <- ' 
using namespace Rcpp;
using namespace arma;

int nr = as<int>(nr_fR);
int nc = as<int>(nc_fR);
int tl = as<int>(tl_fR) - 1;
double minDinc = as<double>(minDinc_fR);
double maxDinc = as<double>(maxDinc_fR);
vec den = as<vec>(den_fR);
vec mscal = as<vec>(mscal_fR);

ivec firstTime = as<ivec>(firstTime_fR);
ivec lastTime = as<ivec>(lastTime_fR);

mat dt0 = as<mat>(dt0_fR);
mat A(nr, nc); A.fill(NA_REAL);

int i, j, k, maxk, jMf;
double s;

for(i = 0; i < nr; i++){
  for(j = firstTime[i]; j <= lastTime[i]; j++){
    s = 0;
    jMf = j-firstTime[i];
    maxk = (jMf < tl ? jMf : tl);
    for(k = maxk; k >=0; k--)
      s +=   mscal(k) * dt0(i, j-k-1);
    
    s /= den(maxk);

    A(i,j-1) = s;

    if(s > maxDinc)
      A(i,j-1) = maxDinc;   
     
    if(s < minDinc)
      A(i,j-1) = minDinc;
      
  }
}  
return wrap(A);

'

lagGrowth_cpp <- cxxfunction(signature(nr_fR = 'numeric', nc_fR = 'numeric', tl_fR = 'numeric',
                             minDinc_fR ='numeric', maxDinc_fR ='numeric',
                             den_fR = 'numeric', mscal_fR = 'numeric', firstTime_fR ='numeric',
                             lastTime_fR ='numeric',dt0_fR = 'numeric'), body = src, 
                             plugin = "RcppArmadillo")


lagGrowthRcpp <- function(dimat, firstTime, lastTime,timelag, discount){  #requires vectors firstTime, lastTime
  
  nr  <- ntree
  nc  <- nyr
  tmp <- lagGrowth_cpp(nr, nc, timelag, minDinc, maxDinc, den, mscal, firstTime, lastTime, dimat)
  tmp
}
##################################seedprob
gsp <- ' 
using namespace Rcpp;
using namespace arma;

int nr = as<int>(nr_fR);
int nc = as<int>(nc_fR);
int nt = as<int>(nt_fR);

mat fecmt = as<mat>(fecmt_fR);
mat kern = as<mat>(kern_fR);
mat A(nr, nt); A.fill(NA_REAL);

int i, j, k;
double s;

for(k = 0; k < nt; k++){

  for(j = 0; j < nr; j++){
    s = 0;
    for(i = 0; i < nc; i++){

      s += kern(j,i) * fecmt(i,k);

    } 
    A(j,k) = s;
  }
}  
return wrap(A);

'

getSeedProb_cpp <- cxxfunction(signature(nr_fR = 'numeric', nc_fR = 'numeric', nt_fR = 'numeric',
                                         fecmt_fR = 'numeric',kern_fR = 'numeric'),
                               body = gsp, plugin = "RcppArmadillo")


##################################byFunction
bsp <- ' 
using namespace Rcpp;
using namespace arma;

int nm = as<int>(nm_fR);
mat frommat = as<mat>(frommat_fR);

mat totmat = as<mat>(totmat_fR);
mat summat = as<mat>(summat_fR);

unsigned int k, i, j;
double s;

for(k = 0; k < nm; k++){

  i = frommat(k,0) - 1;
  j = frommat(k,1) - 1;
  s = frommat(k,2);
  totmat(i,j) = totmat(i,j) + 1;
  summat(i,j) = summat(i,j) + s;
}

List z; z["total"] = totmat; z["sum"] = summat;

return wrap(z);

'

byFunction_cpp <- cxxfunction(signature(nm_fR = 'numeric', frommat_fR = 'numeric',
                                        totmat_fR = 'numeric',summat_fR = 'numeric'),
                              body = bsp, plugin = "RcppArmadillo")

byFunctionRcpp <- function(x, i, j, summat=matrix(0,max(i),max(j)), totmat=summat, MEAN=T){  #
  
  nn <- length(x)
  if( nn != length(i) | nn != length(j) )stop('vectors unequal in byFunctionRcpp')
  if( nrow(summat) < max(i) | ncol(summat) < max(j) )stop('matrix too small')
  
  ww <- which(is.na(x))
  if(length(ww) > 0){
    x <- x[-ww]
    i <- i[-ww]
    j <- j[-ww]
  }
  
  frommat <- cbind(i,j,x)
  
  nm  <- nrow(frommat)
  tmp <- byFunction_cpp(nm, frommat, totmat, summat)
  
  if(!MEAN)return(tmp$sum)
  mu <- tmp$sum/tmp$total
  mu[is.na(mu)] <- 0
  mu
}

#getGamMatrix from FIAfunctions.r, written by Joao

src_loop_GM <- '
using namespace Rcpp;
using namespace arma;

mat kdist = as<mat>(kdist_fr);
mat dxMat = as<mat>(dxMat_fR);
uvec index = as<uvec>(index_fR);
mat mmat = as<mat>(mmat_fR);
mat sdMat = as<mat>(sdMat_fR);
mat ptot = as<mat>(ptot_fR);
mat p0 = as<mat>(p0_fR);
mat gg = as<mat>(gg_fR);
mat ggam = as<mat>(ggam_fR);
mat all_powered = as<mat>(all_powered_fR);

vec kdex = as<vec>(kdex_fR);
uvec indexKoriginal = as<uvec>(indexKoriginal_fR);
uvec indexIN = indexKoriginal;
uvec indexIN2 = indexKoriginal;

int kstep = as<int>(kstep_fR);
int nbreak = as<int>(nbreak_fR);
int sizeindexIN = indexKoriginal.n_elem;

int j, k, l, m, sizeindexIN2, iRow, iCol;
int nr_kdist = kdist.n_rows;

uvec to(sizeindexIN);
uvec from(sizeindexIN);

for(k = 1; k <= kstep; k++){
l = 0; m = 0;
for(j = 0; j < sizeindexIN; j++){    
if(kdex(indexIN(j)) > k){
to(l) = indexKoriginal(indexIN(j));
indexIN(l) = indexIN(j);
l++;      
}    
if(kdex(indexIN2(j)) <= (nbreak-k)){
from(m) = indexKoriginal(indexIN2(j));
indexIN2(m) = indexIN2(j);
m++;      
}    
}
sizeindexIN = l;
sizeindexIN2 = m;
indexIN.resize(sizeindexIN);  
indexIN2.resize(sizeindexIN2);
to.resize(sizeindexIN);
from.resize(sizeindexIN2);

kdist.cols(to) += dxMat.submat(index, to);

for(iRow = 0; iRow < nr_kdist ; iRow++)
for(iCol = 0; iCol < sizeindexIN; iCol++)
p0(iRow,to(iCol)) = R::pnorm(kdist(iRow,to(iCol)), mmat(iRow,from(iCol)), sdMat(iRow,from(iCol)),1,0) - ptot(iRow, to(iCol));

gg.cols(to) += ggam.cols(from) % p0.cols(to) % all_powered.cols(from);
ptot.cols(to) += p0.cols(to);

}


return(wrap(gg));'



loop_in_getGamMatrix_cpp <- cxxfunction( signature(kdex_fR = "numeric", indexKoriginal_fR = "numeric", 
                                                   kstep_fR = "numeric", nbreak_fR = "numeric", 
                                                   kdist_fr = "numeric", index_fR = "numeric", 
                                                   dxMat_fR = "numeric", mmat_fR = "numeric",
                                                   sdMat_fR="numeric", p0_fR = "numeric",
                                                   ptot_fR = "numeric", gg_fR = "numeric", 
                                                   ggam_fR = "numeric",all_powered_fR = "numeric") ,
                                         body = src_loop_GM, plugin = "RcppArmadillo" )


