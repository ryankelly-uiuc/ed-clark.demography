spMap = function(spp, showWarnings=F) {
  if(!showWarnings) warnstate = getOption("warn"); options(warn=-1)
  
  spp = sapply(spp, sp.lookup)
  pft = sapply(spp, pft.lookup)
  
  if(!showWarnings) options(warn=warnstate)

  return(pft)
}

sp.lookup = function(sp) {
  if(sp=='acerBarb') sp = 'ACBA3'  else 
  if(sp=='acerRubr') sp = 'ACRU'   else 
  if(sp=='carpCaro') sp = 'CACA18' else 
  if(sp=='caryGlab') sp = 'CARYA'  else # Carya glabra (pignut). All Carya already mapped are mid-hardwood. 
  if(sp=='caryOvat') sp = 'CAOV2'  else 
  if(sp=='caryTome') sp = 'CARYA'  else # Carya alba (mockernut) not in a PFT
  if(sp=='caryUNKN') sp = 'CARYA'  else 
  if(sp=='cercCana') sp = 'QURU'   else # Cercis canadensis (redbud) ==mid successional (USFS tree index) 
  if(sp=='chioVirg') sp = 'QURU'   else # Chionanthus virginicus ==no idea 
  if(sp=='cornFlor') sp = 'COFL2'  else 
  if(sp=='diosVirg') sp = 'DIVI5'  else 
  if(sp=='fraxAmer') sp = 'FRAM2'  else 
  if(sp=='fraxPenn') sp = 'FRPE'  else 
  if(sp=='ilexDeci') sp = 'ILAM'   else # Ilex decidua is present in all stages. Can colonize early (USFS) 
  if(sp=='ilexOpac') sp = 'ILMO'   else # Very shade tolerant (USFS) else 
  if(sp=='ilexVert') sp = 'ILAM'   else # Ilex verticillata req. full-partial sun (mtcubacenter.org)
  if(sp=='juniVirg') sp = 'JUVI'   else 
  if(sp=='liquStyr') sp = 'LIST2'  else 
  if(sp=='liriTuli') sp = 'LITU'   else 
  if(sp=='moruRubr') sp = 'MORU2'  else 
  if(sp=='nyssSylv') sp = 'QURU'   else # Typically mid-canopy
  if(sp=='ostrVirg') sp = 'OSVI'   else 
  if(sp=='pinuTaed') sp = 'PITA'   else 
  if(sp=='pinuVirg') sp = 'PIVI2'  else 
  if(sp=='prunSero') sp = 'PRSE2'  else 
  if(sp=='querAlba') sp = 'QUAL'   else # Note all oaks mapped already are mid-hardwoods
  if(sp=='querFalc') sp = 'QURU'   else # Southern red oak; mid- to low shade tolerance (USFS)
  if(sp=='querPhel') sp = 'QURU'   else # Willow oak (note USFS says shade intolerant--should be early-hardwood?) 
  if(sp=='querStel') sp = 'QURU'   else # Post oak (note USFS says shade intolerant--should be early-hardwood?)
  if(sp=='querVelu') sp = 'QURU'   else 
  if(sp=='ulmuAlat') sp = 'ULMUS'  else # Classing as mid- like all other Ulmus currently mapped else 
  if(sp=='ulmuAmer') sp = 'ULAM'   else 
  if(sp=='ulmuRubr') sp = 'ULRU'   else
  if(sp=='ulmuUNKN') sp = 'ULMUS'  else

  if(sp=='other')    sp = 'OTHER'  else {
    warning( paste0("Species ", sp, " not found! Set to NA.") )
    sp = NA
  }

  return(sp)
}

pft.lookup = function(sp, num=TRUE, tab='/Work/Research/Macrosystems/Clarkmodel/Data/pftlist.txt', pft.other=18, pft.na=19) {
  dat = read.table(tab, header=T, sep="\t")
  if(is.na(sp)) {
    pft = pft.na
  } else if(sp=="OTHER") {
    pft = pft.other
  } else {
    ind = which(dat$Symbol == sp)
    if(length(ind)==0) {
      warning( paste0("No PFT found for species ", sp, "! Set to NA.") )
      pft = NA
    } else if(length(ind)>1) {
      warning( paste0("Multiple PFT matches found for species ", sp, "! Set to NA.") )
      pft = NA
    } else {

      pft = dat$pft[ ind ]
  
      if(num) {
        if(pft=="temperate.Northern_Pine")      pft = 6 else
        if(pft=="temperate.Southern_Pine")      pft = 7 else
        if(pft=="temperate.Late_Conifer")       pft = 8 else
        if(pft=="temperate.Early_Hardwood")     pft = 9 else
        if(pft=="temperate.North_Mid_Hardwood") pft = 10 else
        if(pft=="temperate.Late_Hardwood")      pft = 11
      }
    }
  }
  return(pft)
}