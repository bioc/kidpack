##------------------------------------------------------------
## 1. average over duplicate spots
## 2. average over duplicate hybs (where available)
##------------------------------------------------------------
library(Biobase)

makeeset = function(intdupcorthresh, mediancvthresh, 
  bgmethod = "bgfilt", useref = TRUE, verbose=TRUE) {

  if(verbose) 
    cat("makeeset:\nbgmethod =", bgmethod, "\nuseref   =", useref, 
        "\nintdupcorthresh =", intdupcorthresh, 
        "\nmediancvthresh =", mediancvthresh, "\n")

  load(paste("ny175", bgmethod, "rda", sep="."))
  load(paste("qc",    bgmethod, "rda", sep="."))
  load("pat.rda")
  load("hybanno.rda")
  load("cloneanno.rda")

  dat  = get(paste("ny175", bgmethod, sep="."))@h
  nslides = ncol(dat)/2
  nprobe = 4224

  if (useref) {
    ## subtract red from green
    kdat = dat[, 1:nslides] - dat[, (nslides+1):ncol(dat)]
  } else {
    kdat = dat[, 1:nslides]
  }
  
  ## average over duplicate spots per clone:
  kdat = (kdat[cloneanno$spot1,] + kdat[cloneanno$spot2,])/2

  ## "good" results: intdupcorthresh = 0.75, mediancvthresh = 0.7
  
  hyb.sel = intdupcor > intdupcorthresh & propneg < 0.7 & mediancv < mediancvthresh

  ## if (method == "subtractbg" & !subtractred)
  ##  hyb.sel = (intdupcorgreen > 0.75 & propneggreen < 0.6 & mediancvgreen < 1.1)

  cat("using", length(which(hyb.sel)), "of", length(hyb.sel), "hybs: ")
  
  ## average over replicate hybs
  npat  = length(unique(hybanno$patientid[hyb.sel]))
  exprs = matrix(nrow=nprobe, ncol=npat)
  for (i in 1:npat) {
    hybs = which(hyb.sel & hybanno$patientid == unique(hybanno$patientid[which(hyb.sel)])[i])
    if (length(hybs)==1) {
      exprs[,i] = kdat[,hybs]
    } else {
      exprs[,i] = rowMeans(kdat[,hybs])
    }
  }
  cat(ncol(exprs), "patients.\n")

  ## the patient table - only those with at least one good hyb
  pat.sel = (pat$patientid %in% hybanno$patientid[which(hyb.sel)])
  
  ##------------------------------------------------------------
  ## Constructor exprSet
  ##------------------------------------------------------------
  rownames(pat) = pat$patientid
  eset = new('exprSet', 
    exprs     = exprs,
    phenoData = new('phenoData',
      pData     = pat[pat.sel,],
      varLabels = as.list(colnames(pat))),
    notes = paste(length(which(hyb.sel))))

  return(eset)  
} 


##  filenam =  paste("eset", length(which(hyb.sel)), ".rda", sep="")
##    save(eset, file=filenam)
