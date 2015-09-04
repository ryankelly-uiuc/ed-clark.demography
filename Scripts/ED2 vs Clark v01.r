rm(list=ls())
source("~/edanalysis/scripts/ED2.analy_fn.data_v02.r")
source("~/edanalysis/scripts/ED2.analy_fn.co2py_v04.r")
source("~/edanalysis/scripts/ED2.analy_fn.plot_v01.r")
source("~/edanalysis/scripts/ED2.analy_fn.varlookup_v01.r")
source("~/ed2-clark_compare/Scripts/clarkDem_pftXsize_fn.r")

# ----------------
# --- ED2 Output Params
  ed.path     = "~/edoutputs/analy/"
  ed.res.flag = '-E-'
  ed.prefix   = "u18.DF-HW.99-07.tiny.reprodoff.ddmort0.8.mortsch0"
  ed.tind     = 49:60   # which samples in time to take. 

# --- Clark Data Params
  clark.data.file = '~/ed2-clark_compare/Data/allData_demPlots.Rdata'
  pftlist.file    = '~/ed2-clark_compare/Data//pftlist.txt'
  site            = "DF_HW"

# --- Size classes
  breaks.new = seq(0,100,10)
  pfts       = 7:11   # Will aggregate all, but only these get plotted

# --- Variable Names
#   ed.varnames = c(
#     "AGB_CO", "DBH", "BA_CO", "HITE", 
# 
#     "PFT", "NPLANT", "LAI_CO", "MMEAN_TRANSP_CO", 
# 
#     "MMEAN_GPP_CO", "MMEAN_NPP_CO", "CBR_BAR", 
#   
#     "MMEAN_MORT_RATE_CO",
#   
#     "BALIVE", "BDEAD", "BLEAF", "BROOT", 
#     "BSAPWOODA", "BSAPWOODB", "BSEEDS_CO", "BSTORAGE", 
#   
#     "DAGB_DT", "DDBH_DT", "DBA_DT", 
#     "DLNAGB_DT", "DLNDBH_DT", "DLNBA_DT", 
# 
#     "MMEAN_LEAF_RESP_CO", "MMEAN_ROOT_RESP_CO", "MMEAN_GROWTH_RESP_CO", "MMEAN_STORAGE_RESP_CO", 
# 
#     "MMEAN_NPPLEAF_CO", "MMEAN_NPPCROOT_CO", "MMEAN_NPPFROOT_CO", 
#     "MMEAN_NPPSAPWOOD_CO", "MMEAN_NPPWOOD_CO", "MMEAN_NPPSEEDS_CO"
#   )
    # Some vars to consider in the future?
    # "MMEAN_RH_PA", "MMEAN_CB_CO", "PATCH_COUNT", "MMEAN_LEAF_GSW_CO", "MMEAN_SOIL_TEMP_PA", "MMEAN_SOIL_WATER_PA", "MMEAN_FAST_SOIL_C", "MMEAN_SLOW_SOIL_C", "MMEAN_STRUCT_SOIL_C", "MMEAN_STRUCT_SOIL_L"


#   clark.varnames    = c("growMeanEstAll","growMeanDatAll","growSdDatAll","survDatAll","survEstAll")
    # These are by species only, no size class divs. Not handling for now
    # vars      = c("reproDatAll","reproEstAll","phiEstAll")

  ed.varnames    = c('NPLANT'     , 'DDBH_DT',        'DDBH_DT',        
                     'MMEAN_MORT_RATE_CO', 'MMEAN_MORT_RATE_CO')
  clark.varnames = c('ytotDataAll', 'growMeanDatAll', 'growMeanEstAll', 
                     'survDatAll',         'survEstAll')
  ed.trans.fn    = list(
    function(x) {x*100*100},
    function(x) {x},
    function(x) {x},
    function(x) {1-x},
    function(x) {1-x}
    )



# --- Save options
  save.fig = T
  save.dat = F




# -------------------------------------------------------------------------------------- #

  # ----- Read in data
    ed.dat = ed.read(ed.path, ed.prefix, ed.res.flag)
    load(clark.data.file)


  # ----- Output files
    if(save.fig) {
      dir.create("/fs/data2/rykelly/ed2-clark_compare/Figs/", recursive=T, showWarnings=F)
      fig.file = file.path("/fs/data2/rykelly/ed2-clark_compare/Figs/",paste0(ed.prefix,"_vsJC.pdf"))
      pdf(fig.file, height=8, width=8)
    }
#     if(save.dat) {
#       dir.create("~/edanalysis/output/dat", recursive=T, showWarnings=F)
#       dat.file = file.path("~/edanalysis/output/dat",paste0(prefix,"_COdat.rds"))
#     }


n.var = length(ed.varnames)

for(var.num in 1:n.var) {
  ed.varname    = ed.varnames[var.num]
  clark.varname = clark.varnames[var.num]


  # ----- Aggregate Clark variables
    clark.X = get(clark.varname)
    clark.X.agg = clarkDem_pftXsize(clark.X, site, breaks.new)


  # --- Aggregate ED variable
    ed.pft    = ed.var(ed.dat, 'PFT')
    ed.dbh    = ed.var(ed.dat, 'DBH')
    ed.X      = ed.var(ed.dat, ed.varname)
    
    if(ed.varname=="MMEAN_MORT_RATE_CO") {
      tmp = list()
      for(i in 1:length(ed.X$time)) {
          tmp[[i]] = ed.X$co[[1]][[i]]+ed.X$co[[2]][[i]]+ed.X$co[[3]][[i]]+
                     ed.X$co[[4]][[i]]+ed.X$co[[5]][[i]]
      }
      ed.X$co = tmp
    }
    
    ed.X.agg = ed.cobin(ed.X, ed.pft, ed.dbh, mean, 1:19, TRUE, breaks.new)$xbin
    
    ed.X.agg = ed.trans.fn[[var.num]](ed.X.agg)

    # Average over specified time period
    ed.X.agg = apply(ed.X.agg[ed.tind,,], c(2,3), mean)
    
    

 plot(clark.X.agg[pfts,], ed.X.agg[pfts,], col=rep(pftcols[pfts], length(breaks.new)), pch=19,
  ylab=ed.varname, xlab=clark.varname)
  abline(0,1)
#   abline(lm(as.vector(ed.X.agg[pfts,])~as.vector(clark.X.agg[pfts,])))
}

if(save.fig) dev.off()