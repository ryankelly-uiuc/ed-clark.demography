rm(list=ls())

# library(ade4)
# library(fpc)

setwd('/Work/Research/Macrosystems/Clarkmodel')

source('./Scripts/fn/clarkRcppFunctions.r')
# source('./Scripts/fn/demFunctions.r')
# source('./Scripts/fn/clarkFunctions.R')
# source('./Scripts/fn/FIAfunctions.r')
source('./Scripts/fn/spMap.r')


dat.file = './Data/allData_demPlots.Rdata'

breaks.new = c(0,10,50,100)
site  = "DF_HW"

vars      = c("growMeanDatAll","growMeanDatAll","growSdDatAll","survDatAll","survEstAll")

# These are by species only, no size class divs. Not handling for now
# vars      = c("reproDatAll","reproEstAll","phiEstAll")

# ------------------------------------------------------
  load(dat.file)

# --- indices
  ii = matrix( unlist( strsplit(colnames(ytotDataAll),'-') ),ncol=2,byrow=T)
  bk = as.numeric(ii[,2])

  specs  = unique(ii[,1])
  breaks = unique(bk)
  sdex   = match(ii[,1],specs)   #species index columns
  kdex   = match(bk,breaks)      #size species columns

  nspec  = length(specs)
  nb     = length(breaks)
  nsite  = nrow(ytotDataAll)
  nvar   = length(vars)

  j      = grep(site, rownames(ytotDataAll))

# --- Aggregate by size class / PFT
  nb.new     = length(breaks.new)
  kdex.new   = findInterval( kdex, breaks.new )


# --- Map species to PFT
  pfts = spMap(specs)
  pdex = pfts[sdex]


# --- Aggrgate variables
vars.out = paste(site, vars, sep=".")
v=1
for(v in 1:nvar) {
  y.j = get(vars[v])[j,]
  assign(vars.out[v], 
          byFunctionRcpp(y.j, pdex, kdex.new, matrix(0, 19, nb.new), MEAN=T))
}


  


####### demography: growMeanDatAll,reproEstAll,phiEstAll,
####### growSdDatAll,survDatAll,survEstAll
#presTreeNCAll, presTreeAll,p50 = presPopAll



# byFunctionRCPP example: Sum over all size classes
#   nall = nrow(ytotDataAll)          #no. plots
#   nsnb = ncol(ytotDataAll)          # Sp X # Size class
#   jj   = as.vector( matrix( c(1:nall),nall,nsnb) )      #plot index
#   ss   = as.vector( matrix( sdex,nall,nsnb,byrow=T) )   #species index
# 
#   yall = byFunctionRcpp(as.vector(ytotDataAll),jj,ss,
#                          matrix(0,nall,nspec),matrix(0,nall,nspec),
#                          MEAN=T)




