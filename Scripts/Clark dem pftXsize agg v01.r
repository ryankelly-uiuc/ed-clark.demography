rm(list=ls())

setwd(  '/Work/Research/Macrosystems/Clarkmodel' )
load(   './Data/allData_demPlots.Rdata'          )
source( './Scripts/fn/clarkDem_pftXsize.r'       )

breaks.new = c(0,10,50,100)
site       = "DF_HW"
vars       = c("growMeanEstAll","growMeanDatAll","growSdDatAll","survDatAll","survEstAll")
  # These are by species only, no size class divs. Not handling for now
  # vars      = c("reproDatAll","reproEstAll","phiEstAll")

pftlist.file = '/Work/Research/Macrosystems/Clarkmodel/Data/pftlist.txt'

# ---------------------------
# --- Aggrgate variables
nvar = length(vars)
vars.out = paste0(site, ".", vars, ".agg")
v=1
for(v in 1:nvar) {
  X = get(vars[v])
  assign(vars.out[v], clarkDem_pftXsize(X, site, breaks.new))
}



