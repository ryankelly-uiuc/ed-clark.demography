rm(list=ls())
setwd('/Work/Research/Macrosystems/Clarkmodel')

load('./Data/allData.Rdata')


# ----- Check out variable dims
  allvars = ls()
  n.vars = length(allvars)
  dims = matrix(NA, nrow=n.vars,ncol=3)
  i=1
  for(i in 1:n.vars) {
    x = get(allvars[i])
    if(!is.null(dim(x))) {
      if(length(dim(x))==2) {
        dims[i,1:2] = dim(x)
      } else if(length(dim(x))==3) {
        dims[i,1:3] = dim(x)
      }
    } else {
      dims[i,1]=length(x)
    }
  }
  dims = data.frame(dims)
    rownames(dims) = allvars
    names(dims) = paste("dim", 1:3, sep=".")

    bigdims = dims[which(apply(dims,1,sum, na.rm=T)>1),]
    bigdims



# ----- Find Macrosystems Demography Plots
  tmp <- read.table('./Data/plotDataDuke.txt',header=T)[,'plotname']
  demPlots <- as.character(tmp)

  idem <- numeric(0)
  for(j in 1:length(demPlots)){
    idem <- c(idem,grep(demPlots[j],rownames(growMeanEstAll)))
  }



# ----- Overwrite big vars with only the demplot data
  nall <- nrow(ytotDataAll)          #no. plots
  plotmats = which(dims[,1]==nall)
  i=1
  for(i in 1:length(plotmats)) {
    x = get(allvars[plotmats[i]])
    if(!is.null(dim(x))) {
      assign(allvars[plotmats[i]], x[idem,])
    } else {
      assign(allvars[plotmats[i]], x[idem])
    }
  }


# ----- Get rid of some big variables
  # List variables by size
  sort( sapply(ls(),function(x){object.size(get(x))}))
  sum( sapply(ls(),function(x){object.size(get(x))}))

  rm(chainsAll) # Will need to reload the big file if I want all the chains


# ----- Save stripped down version
  save.image('./Data/allData_demPlots.Rdata')


# ----- TEST
  rm(list=ls())
  load('./Data/allData_demPlots.Rdata')

