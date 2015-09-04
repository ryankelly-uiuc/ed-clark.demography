byFunctionRK = function(x, i, j, summat=matrix(0,max(i),max(j)), totmat=summat, MEAN=T){  #
  nn = length(x)
  if( nn != length(i) | nn != length(j) )stop('vectors unequal in byFunctionRK')
  if( nrow(summat) < max(i) | ncol(summat) < max(j) )stop('matrix too small')
  
  ww = which(is.na(x))
  if(length(ww) > 0){
    x = x[-ww]
    i = i[-ww]
    j = j[-ww]
  }
  
#   frommat = cbind(i,j,x)
  
  nx  = length(x)
#   tmp = byFunction_cpp(nm, frommat, totmat, summat)

  for(q in 1:nx) {
    summat[ i[q], j[q] ] = summat[ i[q], j[q] ] + x[q]
    totmat[ i[q], j[q] ] = totmat[ i[q], j[q] ] + 1
  }
  
  if(MEAN) {
    mu = summat/totmat
    mu[is.na(mu)] = 0
    return(mu)
  } else {
    return(summat)
  }
}