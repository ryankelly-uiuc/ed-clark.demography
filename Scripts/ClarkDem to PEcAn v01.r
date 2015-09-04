rm(list=ls())
source("~/edanalysis/scripts/fn/ED2.analy_fn.data_v02.r")
source("~/edanalysis/scripts/fn/ED2.analy_fn.co2py_v04.r")
source("~/edanalysis/scripts/fn/ED2.analy_fn.plot_v01.r")
source("~/edanalysis/scripts/fn/ED2.analy_fn.varlookup_v01.r")
source("~/ed2-clark_compare/Scripts/fn/clarkDem_pftXsize.r")

# ----------------
# --- Clark Data Params
  clark.data.file = '~/ed2-clark_compare/Data/allData_demPlots.2015.08.05.Rdata'
  pftlist.file    = '~/ed2-clark_compare/Data//pftlist.txt'
  clark.site            = "DF_HW"

  clark.varnames = c('ytotDataAll', 
                     'growMeanDatAll', 'growMeanEstAll', 
                     'survDatAll',     'survEstAll',
                     
                     'phiEstAll')

  clark.units = c('count?',
                  'cm/yr?', 'cm/yr?',
                  'proportion surviving?', 'proportion surviving?',
                  'recruits per basal area per yr')

  out.varnames = gsub('All', '', clark.varnames)



# --- Save options
  save.dat = T
  save.dir = '~/ed2-clark_compare/Data/site_demog'

  store.db     = T
  overwrite.db = T

  pecan.siteid     = 755
  pecan.mimetype   = 'application/x-rds'
  pecan.formatname = 'Clark Model Demography Data'
  pecan.startdate  = '1999-01-01 00:00:00' # Making up for now
  pecan.enddate    = '2007-12-31 11:59:59' # Making up for now
  
  db.settings = list(
    user='bety', 
    password='bety', 
    host='psql-pecan', 
    dbname='bety', 
    driver='PostgreSQL', 
    write=TRUE
  )
  


# -------------------------------------------------------------------------------------- #
# --- Read in data
  load(clark.data.file)


# --- Extract data for this site
  ## Can place in loop later to do multiple sites
  # Setup
  out = list()
  n.var = length(clark.varnames)

  # Loop over variables
  for(i in 1:n.var) {
    # Select the variable
    clark.varname = clark.varnames[i]
    clark.X = get(clark.varname)

    # Find row corresponding to site
    site.ind = grep(clark.site, rownames(clark.X))

    # Add new empty element to output list
    out[[i]] = list()
      names(out)[i] = out.varnames[i]

    # Store clark data and any units provided
    out[[i]]$data = clark.X[site.ind, ]
    out[[i]]$units = clark.units[i]
  }

  # Save output list
  if(save.dat) {
    dir.create(save.dir, recursive=T, showWarnings=F)
    save.file = paste0(clark.site, '.rds')
    saveRDS(out, file.path(save.dir, save.file))

    # Store record in BETY
    if(store.db) {
      library(PEcAn.DB)
      library(RPostgreSQL)
      con = db.open(db.settings)
      
      if(overwrite.db) {
        check = dbfile.input.check(pecan.siteid, pecan.startdate, pecan.enddate, pecan.mimetype, pecan.formatname, parentid=NA, con=con, hostname=fqdn())
        if(nrow(check) > 0) {
          db.query(paste0('DELETE FROM dbfiles WHERE id=', check$id), con)
        }
      }

      dbfile.input.insert(normalizePath(save.dir), save.file, 
        pecan.siteid, pecan.startdate, pecan.enddate, pecan.mimetype, pecan.formatname, 
        parentid=NA, con=con, hostname=fqdn())
    
      db.close(con)
    }
  }
