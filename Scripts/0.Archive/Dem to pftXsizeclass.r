rm(list=ls())

library(ade4)
library(fpc)

setwd('/Work/Research/Macrosystems/Clarkmodel')

source('./Scripts/fn/clarkRcppFunctions.r')
source('./Scripts/fn/demFunctions.r')
source('./Scripts/fn/clarkFunctions.R')
source('./Scripts/fn/FIAfunctions.r')
source('./Scripts/fn/spMap.r')


load('./Data/allData_demPlots.Rdata')
regSections <- c('SW','NE','SE','MW')
nsecs <- length(regSections)


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


# byFunctionRCPP example: Sum over all size classes
# yall <- byFunctionRcpp(as.vector(ytotDataAll),jj,ss,
#                        matrix(0,nall,nspec),matrix(0,nall,nspec),
#                        MEAN=T)


# Aggregate by size class / PFT
breaks.new = c(0,10,50,100)
nb.new     = length(breaks.new)
kdex.new   = findInterval( kdex, breaks.new )


# Choose variable y to aggregate
y     = ytotDataAll

site  = "DF_HW"
j     = grep(site, rownames(ytotDataAll))


y.j = y[j,]


# Map species
source('./Scripts/fn/spMap.r')
pfts = spMap(specs)

pdex = pfts[sdex]

pftXsize.j = byFunctionRcpp(y.j, pdex, kdex.new, matrix(0, 19, nb.new), MEAN=T)
  

# spXsize.j = byFunctionRcpp(y.j, sdex, kdex.new,
#                        matrix(0, nspec, nb.new), MEAN=T)

pfts[which(apply(spXsize.j, 1, sum) >0)]


# Check byFunctionRcpp is doing what I think
# timer=Sys.time()
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
# Sys.time()-timer
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















# List of PFTs to consider
  pft.names = c(
    "temperate.Late_Hardwood",        # ED PFT 11
    "temperate.Early_Hardwood",       # ED PFT 9
#     "temperate.Northern_Pine",        # ED PFT 6
    "temperate.Late_Conifer",         # ED PFT 8
    "temperate.Southern_Pine",        # ED PFT 7
    "temperate.North_Mid_Hardwood"    # ED PFT 10
  )
  
  pft.nums = c(11,9,8,7,10)

# BETY settings on BU pecan server
  db.params = list(
    driver   = "PostgreSQL",
    host     = "128.197.230.32",
    user     = "bety",
    password = "bety",
    dbname   = "bety",
    write    = "true" )

# --- Get PFT-species mapping
  pfts = pftSpeciesList(db.params, pft.names)
