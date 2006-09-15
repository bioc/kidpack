library(Biobase)
library(multtest)
library(globaltest)
thresh = 0.005
sink(file=paste("Results2/test", thresh, ".txt", sep=""))
#sink(file="test3.txt")
source("mtest.R")
if (!exists("eset"))
  load("eset74.rda")
if (!exists("cloneanno"))
  load("cloneanno.rda")
write.clonelists = F
do.lm = F
do.cox = T
do.Ftest = T
do.cctests = T
save.results = T
estimate.FDR = T
source("write.clonelist.R")
pat = eset@phenoData@pData
npat = nrow(pat)
x11(height=, width=8)
par(mfrow=c(3,4))
SubtypeNum = rep(0, nrow(pat))
SubtypeNum[which(pat$Subtype=="pRCC")] = 1
SubtypeNum[which(pat$Subtype=="chRCC")] = 2
print(table(pat$Subtype))
set.seed(28764)

comparisons1 = c("ccRCC", "chRCC", "pRCC", "newcluster", "pRCC.vs.ccRCC", "ccRCCstage", "ccRCCgrade", "ccRCCprogress", "ccRCCgain5q", "ccm")
comparisons2 = c("no.of.net.changes") #, "tumor.size.cm")
comparisons = c("F-test", comparisons1, comparisons2, "Cox")
if (estimate.FDR & !exists("FDR")) {
  nrep = 500
  nfp = matrix(nrow=nrep, ncol=length(comparisons))
  enfp = FDR = rep(NA, length(comparisons))
  colnames(nfp) = names(enfp) = names(FDR) = comparisons
}

if (do.Ftest) {
  f = ftestn(eset@exprs, as.factor(pat$Subtype))
  hist(f$pvalue, breaks=100, main=paste("F-test for subtype"))
  cat("F-test for subtype: Bonferroni-adjusted minimal pvalue:", min(1, min(f$pvalue)*nrow(eset@exprs)), "; ", length(which(f$pvalue<thresh)), "clones with p <", thresh, "\n")
  if (estimate.FDR) {
    for (rep in 1:nrep) {
      if (rep %% 100 ==0)
        cat(rep, "\t")
      perm = sample(ncol(eset@exprs))
      fperm = ftestn(eset@exprs[,perm], as.factor(pat$Subtype))
      nfp[rep,"F-test"] = length(which(fperm$pvalue < thresh))
    }
    enfp["F-test"] = mean(nfp[,"F-test"])*(2*length(which(f$pvalue>0.5)))/nrow(eset@exprs)
    FDR["F-test"] = min(1, enfp["F-test"]/length(which(f$pvalue<thresh)))
    cat("estimated FDR:", FDR["F-test"], "\n")
  }
  if (write.clonelists) {
    table = cbind(I(f$pvalue), cloneanno)
    colnames(table)[1] = "p-value"
    write.clonelist(table[order(f$pvalue)[1:length(which(f$pvalue<thresh))], ], filename = "Results2/Ftest.subtype", title = "F-test for subtype", sortby="p-value")
  }
  cat("hi\n")
}
###################################################################
#different tests within the clear cell tumors:
###################################################################

if (do.cctests) {
cccoeff = ccp = matrix(nrow=nrow(eset@exprs), ncol= length(comparisons1)+ length(comparisons2))
colnames(cccoeff) = colnames(ccp) = c(comparisons1, comparisons2)
for (comp in comparisons1) {
  if (comp == "ccRCC") {
    patients = 1:npat
    classlabels = rep(0, npat)
    classlabels[which(pat$Subtype=="ccRCC")] = 1
  }
  if (comp == "chRCC") {
    patients = 1:npat
    classlabels = rep(0, npat)
    classlabels[which(pat$Subtype=="chRCC")] = 1
  }
  if (comp == "pRCC") {
    patients = 1:npat
    classlabels = rep(0, npat)
    classlabels[which(pat$Subtype=="pRCC")] = 1
  }
  if (comp == "pRCC.vs.ccRCC") {
    patients = which(pat$Subtype=="ccRCC" | pat$Subtype=="pRCC")
    classlabels = SubtypeNum[patients]
  }
  if (comp == "chRCC.vs.ccRCC") {
    patients = which(pat$Subtype=="ccRCC" | pat$Subtype=="chRCC")
    classlabels = 0.5*SubtypeNum[patients]
  }
  if (comp == "chRCC.vs.pRCC") {
    patients = which(pat$Subtype=="pRCC" | pat$Subtype=="chRCC")
    classlabels = SubtypeNum[patients]-1
  }
  if (comp == "ccRCCstage") {
    patients = which(pat$Subtype=="ccRCC")
    classlabels = rep(0, length(patients))
    classlabels[which(pat$clinical.stage[patients] > 2)] = 1
  }
  if (comp == "ccRCCgrade") {
    patients = which(pat$Subtype=="ccRCC"& !(is.na(pat$grading)))
    classlabels = rep(0, length(patients))
    classlabels[which(pat$grading[patients] > 1)] = 1
  }
  if (comp == "ccRCCprogress") {
    patients = which(pat$Subtype=="ccRCC" & !(is.na(pat$progress)))
    classlabels = pat$progress[patients]
  }
  if (comp == "ccRCCgain5q") {
    patients = which(pat$Subtype=="ccRCC" & !is.na(pat[,"5q"]))
    classlabels = rep(0, length(patients))
    classlabels[which(pat[patients,"5q"] > 2)] = 1
  }
  if (comp == "ccm") {
    patients = which(pat$Subtype=="ccRCC" & !is.na(pat[,"m"]))
    #patients = which(pat$Subtype=="ccRCC" & !is.na(pat[,"progress"]))
    classlabels = rep(0, length(patients))
    #classlabels[which(pat[patients,"rf.survival"] == 0 & pat[patients,"progress"])] = 1
    classlabels[which(pat[patients,"m"] == 1)] = 1
  }
  if (comp == "surv24") {
    patients = which(pat$Subtype=="ccRCC" & !is.na(pat[,"survival.time"]) & (pat[,"died"]==1 | pat[,"survival.time"] > 24))
    classlabels = rep(0, length(patients))
    classlabels[which(pat[patients,"died"] ==1 & pat[patients,"survival.time"] <= 24)] = 1 
  }
  if (comp == "newcluster") {
    group1 = which(pat$id.neu %in% c(40,31,63,39,70,46,74,30,59,69,27,72,71,58,53,56,65,48,23,62,52,37))
    group2 = which(pat$id.neu %in% c(64,54,45,34,32,57,41,36,47,60,68,28,66,24,42,50,33,51,55,38,67,43,73,26))
    patients = c(group1, group2)
    stopifnot(all(pat$Subtype[patients]=="ccRCC"))
    classlabels = c(rep(0, length(group1)), rep(1, length(group2)))
  }
  cccoeff[,comp] = apply(eset@exprs[,patients[which(classlabels==1)]],1,mean) - apply(eset@exprs[,patients[which(classlabels==0)]],1,mean)
  t = mt.teststat(eset@exprs[,patients], classlabels, test="t.equalvar")
  ccp[,comp] = 2*pmin(pt(t, df=length(patients)-2), 1-pt(t, df=length(patients)-2))
  hist(ccp[,comp], breaks=100, main=comp)
  cat(comp, length(patients)-sum(classlabels), "vs.", sum(classlabels), "patients: Bonferroni-adjusted minimal p-value:", min(1, min(ccp[,comp])*nrow(eset@exprs)), "; ", length(which(ccp[,comp]<thresh)), "clones with p <", thresh, "\n")
  z= globaltest(eset@exprs[,patients], classlabels,permutation=T)
  cat("global p-value (globaltest): p=", z@p.value, "\n")
  if (estimate.FDR) {
    for (rep in 1:nrep) {
      if (rep %% 100 ==0)
        cat(rep, "\t")
      perm = sample(length(patients))
      tperm = mt.teststat(eset@exprs[,patients[perm]], classlabels)
      pperm = 2*pmin(pt(tperm, df=length(patients)-2), 1-pt(tperm, df=length(patients)-2))
      nfp[rep,comp] = length(which(pperm < thresh))
    }
    enfp[comp] = mean(nfp[,comp])*(2*length(which(ccp[,comp]>0.5)))/nrow(eset@exprs)
    FDR[comp] =  min(1, enfp[comp]/length(which(ccp[,comp]<thresh)))
    cat("estimated FDR:", FDR[comp], "\n")
  }
  if (write.clonelists) {
    table = cbind(I(ccp[,comp]), I(cccoeff[,comp]), cloneanno)
    colnames(table)[1] = "p-value"
    colnames(table)[2] = "coefficient"
    write.clonelist(table[order(ccp[,comp])[1:length(which(ccp[,comp]<thresh))], ], filename = paste("Results2/", comp, sep=""), title = comp, sortby="p-value")
  }
}

#tumor size and no of changes:
options(contrasts=c("contr.sum", "contr.poly")) #constraint on coefficients is s.th. they sum up to 0
for (comp in comparisons2) {#no.of.net.changes
  patients = which(pat$Subtype=="ccRCC" & !is.na(pat[,comp]))
  gl = globaltest(eset@exprs[,patients], pat[patients,comp], model="linear", permutation=T)
  for (j in 1:nrow(eset@exprs)) {
    if (j %% 200 ==0)
      cat(j, "\t")
    z <- lm(pat[patients,comp] ~ eset@exprs[j,patients])
    ccp[j, comp] = summary(z)$coefficients[2,4]
    cccoeff[j, comp] = summary(z)$coefficients[2,1]
  }
  cat("\n", comp, "Bonferroni-adjusted minimal p-value:", min(1, min(ccp[,comp])*nrow(eset@exprs)), "; ", length(which(ccp[,comp]<thresh)), "clones with p < ", thresh, "\n")
  cat("global p-value (globaltest): p=", gl@p.value, "\n")
  if (estimate.FDR) {
    for (rep in 1:nrep) {
      pperm = rep(NA, nrow(eset@exprs))
      if (rep %% 10 ==0)
        cat(rep, "\t")
      perm = sample(length(patients))
      for (j in 1:nrow(eset@exprs)) {
        zperm = lm(pat[patients,comp] ~ eset@exprs[j,patients[perm]])
        pperm[j] = summary(zperm)$coefficients[2,4]
      }
      nfp[rep,comp] = length(which(pperm < thresh))
    }
    enfp[comp] = mean(nfp[,comp])*(2*length(which(ccp[,comp]>0.5)))/nrow(eset@exprs)
    FDR[comp] =  min(1, enfp[comp]/length(which(ccp[,comp]<thresh)))
    cat("estimated FDR:", min(1, enfp[comp]/length(which(ccp[,comp]<thresh))), "\n")
  }
  if (write.clonelists) {
    table = cbind(I(ccp[,comp]), I(cccoeff[,comp]), cloneanno)
    colnames(table)[1] = "p-value"
    colnames(table)[2] = "coefficient"
    write.clonelist(table[order(ccp[,comp])[1:length(which(ccp[,comp]<thresh))], ], filename = paste("Results2/cc.", comp, sep=""), title = comp, sortby="p-value")
  }
}
}

if (do.cox) {
  library(survival)
  patients = which(pat$Subtype=="ccRCC" & !is.na(pat[,"survival.time"]))
  sur = Surv(pat$survival.time[patients], pat$died[patients])
  pcox = rep(1, nrow(eset@exprs))
  coeffcox = rep(0, nrow(eset@exprs))
  for (j in 1:nrow(eset@exprs)) {
    if (j %% 200 ==0)
      cat(j, "\t")
    #if (!(j %in% c(112, 799, 980, 1192, 1246, 1535, 1627, 2473, 2657, 2895, 3043))) {
      cp = coxph(sur ~ eset@exprs[j,patients])
      pcox[j] = 1 - pchisq(2*(cp$loglik[2] - cp$loglik[1]), 1)
      coeffcox[j] = cp$coefficients
    #}
  }
  cat("survival (ccRCC), Bonferroni-adjusted minimal p-value:", min(1, min(pcox*nrow(eset@exprs))), "; ", length(which(pcox<thresh)), "clones with p <", thresh, "\n")
  hist(pcox, breaks=100, main="survival (Cox model)")
  if (estimate.FDR) {
    comp = "Cox"
    for (rep in 1:nrep) {
      pperm = rep(NA, nrow(eset@exprs))
      if (rep %% 10 ==0)
        cat(rep, "\t")
      perm = sample(length(patients))
      for (j in 1:nrow(eset@exprs)) {
        cp = coxph(sur ~ eset@exprs[j,patients[perm]])
        pperm[j] = 1 - pchisq(2*(cp$loglik[2] - cp$loglik[1]), 1)
      }
      nfp[rep,comp] = length(which(pperm < thresh))
    }
    enfp[comp] = mean(nfp[,comp])*(2*length(which(pcox>0.5)))/nrow(eset@exprs)
    FDR[comp] =  min(1, enfp[comp]/length(which(pcox<thresh)))
    cat("estimated FDR:", FDR[comp], "\n")
   } 
  if (write.clonelists) {
    comp = "ccRCCsurvival"
    table = cbind(I(pcox), I(coeffcox), cloneanno)
    colnames(table)[1] = "p-value"
    colnames(table)[2] = "coefficient"
    write.clonelist(table[order(pcox)[1:length(which(pcox<thresh))], ], filename = paste("Results2/", comp, sep=""), title = paste(comp, "(Cox model)"), sortby="p-value")
  }
}

if (save.results) {
  if (do.cctests) 
    save(cccoeff, ccp, file = "Results2/cctest.results.rda")
  if (do.Ftest)
    save(f, file = "Results2/subtype.results.rda")
  if (do.cox)
    save(pcox, coeffcox, file = "Results2/cox.results.RData")
  if (estimate.FDR)
  save(nfp, enfp, FDR, file = paste("Results2/test.FDR.results", thresh, ".RData", sep=""))
}
sink()

if (0) {
  interest = which(p1<thresh)
  patlist1 = which(pat$Subtype=="pRCC")
  patlist2 = intersect(grep("+5",pat$clonal.net.changes), which(pat$Subtype=="ccRCC")) 
  patlist3 = setdiff(which(pat$Subtype=="ccRCC"), patlist2)
  patients = c(patlist1, patlist2, patlist3)
  subdata = eset@exprs[interest, patients]
  x11()
  par(mfrow=c(5,6))
  for (j in 1:length(interest)) {
    plot(subdata[j,], col=c(rep("Black", length(patlist1)), rep("Red", length(patlist2)),rep("Blue", length(patlist3))))
  }
}

if (do.lm) {
  phenotypes = c("lmcc", "lmch", "lmpap") 
  pa = coeff = matrix(nrow = nrow(eset@exprs), ncol=length(phenotypes))
  colnames(pa) = colnames(coeff) = phenotypes
  options(contrasts=c("contr.sum", "contr.poly")) #constraint on coefficients is s.th. they sum up to 0
  asubtype <- rep("apRCC", nrow(pat))
  asubtype[which(pat$Subtype=="ccRCC")] = "ccRCC"
  asubtype[which(pat$Subtype=="chRCC")] = "chRCC"
  for (j in 1:nrow(eset@exprs)) {
    if (j %% 200 ==0)
      cat(j, "\t")
    z <- lm(eset@exprs[j,] ~ as.factor(pat$Subtype))
    pa[j, "lmcc"] = summary(z)$coefficients[2,4]
    pa[j, "lmch"] = summary(z)$coefficients[3,4]
    coeff[j, "lmcc"] = summary(z)$coefficients[2,1]
    coeff[j, "lmch"] = summary(z)$coefficients[3,1]
    z2 <- lm(eset@exprs[j,] ~ as.factor(asubtype))
    pa[j, "lmpap"] = summary(z2)$coefficients[2,4]
    coeff[j, "lmpap"] = summary(z2)$coefficients[2,1]
  }
  
  hist(pa[, "lmcc"], breaks=100, main=paste("lm clear cell"))
  cat("\nlm clear cell, Bonferroni-adjusted minimal p-value:", min(pa[, "lmcc"])*nrow(eset@exprs), "; ", length(which(pa[, "lmcc"]<thresh)), "clones with p <", thresh, "\n")

  hist(pa[, "lmch"], breaks=100, main=paste("lm chromophobe"))
  cat("lm chromophobe, Bonferroni-adjusted minimal p-value:", min(pa[, "lmch"])*nrow(eset@exprs), "; ", length(which(pa[, "lmch"]<thresh)), "clones with p <", thresh, "\n")

  hist(pa[, "lmpap"], breaks=100, main=paste("lm papillary"))
  cat("lm papillary, Bonferroni-adjusted minimal p-value:", min(pa[, "lmpap"])*nrow(eset@exprs), "; ", length(which(pa[, "lmpap"]<thresh)), "clones with p <", thresh, "\n")
  if (write.clonelists) {
    for (comp in phenotypes) {
      table = cbind(I(pa[,comp]), I(coeff[,comp]), cloneanno)
      colnames(table)[1] = "p-value"
      colnames(table)[2] = "coefficient"
      write.clonelist(table[order(pa[,comp])[1:length(which(pa[,comp]<thresh))], ], filename = paste("Results2/", comp, sep=""), title = paste(comp, "(linear model)"), sortby="p-value")
    }
  }
  if (save.results) {
    save(pa, coeff, file = "Results2/subtypeLM.rda")
  }
}

