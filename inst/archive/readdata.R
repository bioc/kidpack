#Read data:
#1. Patient data from ii) the general Excel table Genereal-NZK*
#   and i) a special table . In the case of conflicting entries,
#   ii) is preferred.
#2. List of hybs from the corresponding directory
#3. The intensity data files
#
# Output: pat.rda, hybanno.rda, qua.rda

#--------------------------------------------------
# filename, path parameters and other constants
#--------------------------------------------------

patientfile1 = "SampleData(Kidney2).txt"
patientfile2 = "Table1-General.RData"
#patientfile2 = "/home/heydebre/cyto_dir/General-NZK-04-06-03.RData" #old
#patientfile3 = "/home/heydebre/cyto_dir/barplot.data4.cnside.RData" #data on most frequent aberrations
hybdir      = "Echip"
read.intensities = F
#------------------------------------------------------------
# read patient data file
#------------------------------------------------------------
pat <- read.table(file.path(hybdir,patientfile1), header=TRUE, skip=4, 
                  sep="\t", fill=FALSE,
                  colClasses=c("numeric","character",   ##   lfdnr, patientid
                               "numeric","numeric",     ##   RNA.c., µg.RNA
                               "numeric","numeric",     ##   µl.RNA, µl.H2O
                               "character","character", ##   Dye, Subtype
                               "character","character", ##   clonal.net.changes, grading
                               "numeric"))              ##   clinical.stage
colnames(pat) = c("lfdnr", "patientid", "RNA.c.", "µg.RNA", "µl.RNA", "µl.H2O", "Dye",
          "Subtype", "clonal.net.changes", "grading", "clinical.stage")
#because the column names have changed in the text file!

cat("Read", nrow(pat), "rows with", ncol(pat), "columns from", patientfile1, "\n")
pat$patientid = gsub(" ", "", pat$patientid)

#some corrections concerning the E-Nr.
pat$patientid[which(pat$patientid=="98-U09256")] = "98-09256"
pat$patientid[which(pat$patientid=="98-U04607")] = "98-04607"
pat$patientid[which(pat$patientid=="98-U04349")] = "98-04349"
pat$patientid[which(pat$patientid=="98-U00100")] = "98-00100"
pat$patientid[which(pat$patientid=="97-09688")]  = "97-U09688"

rownames(pat) = pat$patientid
## plausibility checks
stopifnot(all(pat$lfdnr==1:nrow(pat)))
stopifnot(all(!duplicated(pat$patientid)))
stopifnot(all(pat$Dye=="Cy3"))
stopifnot(all(sort(unique(pat$Subtype))==c("ccRCC", "chRCC", "pRCC")))

if(!is.numeric(pat$grading)){
  pat$grading <- gsub(" ", "", pat$grading)
  pat$grading <- gsub("G", "", pat$grading)
  pat$grading[which(pat$grading=="1-2")]  <- "2"
  pat$grading[which(pat$grading=="2-3")]  <- "3"
  pat$grading[which(pat$grading=="")]  <- NA
  irregular.grade = which(!(pat$grading %in% c("1","2","3",NA)))
  cat("irregular grade entries replaced by NA:", pat$grading[irregular.grade],"\n")
  pat$grading[irregular.grade] = NA
}
pat$grading <- as.numeric(pat$grading)
stopifnot(all(pat$grading %in% c(1:3) | is.na(pat$grading)))
stopifnot(all(pat$clinical.stage %in% 1:4 | is.na(pat$clinical.stage)))
pat$clonal.net.changes[which(pat$clonal.net.changes=="?")] = NA
pat= cbind(pat, rep(F, nrow(pat)))
colnames(pat)[ncol(pat)] = "in.Table1"

## add information from the General-table:
load(patientfile2)
x[grep("klarz", x$histo),"histo"] = "ccRCC"
x[grep("chromo", x$histo),"histo"] = "chRCC"
x[grep("pap", x$histo),"histo"] = "pRCC"

xpat = character(0)
for (p in 1:nrow(pat)) {
  #replace conflicting entries by those from the General-table:
  if (any(rownames(x) == rownames(pat)[p])) {
    pat$in.Table1[p] = T
    newrow = rownames(pat)[p]
    stopifnot(length(newrow)==1)
    xpat = c(xpat, newrow)
    if (!identical(pat$Subtype[p], x[newrow,"histo"])) {
      cat(rownames(pat)[p], ": subytpe:", x[newrow,"histo"], "instead of", pat$Subtype[p], "\n")
      pat$Subtype[p] = x[newrow,"histo"]
    }
    if (!identical(pat$clinical.stage[p], x[newrow,"clinical.stage"])) {
      cat(rownames(pat)[p], ": clinical stage:", x[newrow,"clinical.stage"], "instead of", pat$clinical.stage[p], "\n")
      pat$clinical.stage[p] = x[newrow,"clinical.stage"]
    }
    if (!identical(pat$grading[p], x[newrow,"grading"])) {
      cat(rownames(pat)[p], ": grade:", x[newrow,"grading"], "instead of", pat$grading[p], "\n")
      pat$grading[p] = x[newrow,"grading"]
    }
    if (is.na(x[newrow,"histo"])) {#I checked the other patient's subtype manually with check.pat.data.R
      cat(rownames(pat)[p], ": Subtype NA instead of",pat$Subtype[p],"\n")  
      pat$Subtype[p] = NA
    }
    pat$clonal.net.changes[p] = gsub(" ", "", pat$clonal.net.changes[p])
    pat$clonal.net.changes[p] = gsub("\\\"", "", pat$clonal.net.changes[p])
    if (!identical(pat$clonal.net.changes[p], x[newrow,"clonal.netto.changes"]) & !(x[newrow,"clonal.netto.changes"]=="" & is.na(pat$clonal.net.changes[p]))) {
      cat(pat$patientid[p], ": clonal net changes:", as.character(x[newrow,"clonal.netto.changes"]), "instead of", pat$clonal.net.changes[p], "\n")
      pat$clonal.net.changes[p] = as.character(x[newrow,"clonal.netto.changes"])
    }
  }
}
#add additional variables
columns = c("organ", "t", "n", "m", "progress", "rf.survival", "died", "survival.time.in.months", "age", "sex", "tumor.size.cm", "no.of.net.changes", "id.neu")
pat = cbind(pat, matrix(NA, nrow=nrow(pat), ncol=length(columns)))
colnames(pat)[(ncol(pat)-length(columns)+1):ncol(pat)] = columns
for (column in columns) 
  pat[xpat, columns] =  x[xpat, columns]
  #pat[xpat, columns] =  x[xpat, columns]
colnames(pat)[which(colnames(pat)=="survival.time.in.months")] = "survival.time"


#chromosomal aberrations from sidelines (copy number data):
tmp = matrix(NA, nrow = nrow(pat), ncol=ncol(patcnarmside))
colnames(tmp) = colnames(patcnarmside)
for (pt in 1:nrow(pat)) {
  if (any(rownames(patcnarmside)==rownames(pat)[pt]))
    tmp[pt,] = patcnarmside[rownames(pat)[pt],]
}
pat = cbind(pat, tmp)
save(pat,     file="pat.rda")

stop()
#------------------------------------------------------------
# construct hybanno table with patient names and filenames
#------------------------------------------------------------
hybfiles = sort(dir(hybdir, pattern="[0-9][0-9]*.txt$"))
sts      = strsplit(hybfiles, "_")
hybanno  = data.frame(filename  = I(hybfiles),
                      patientid = I(sapply(sts, function(x) return(x[1]))),
                      slideid   = I(sub(".txt", "", sapply(sts, function(x) return(x[2])))))
hybanno$patientid[which(hybanno$patientid=="98-U09256")] = "98-09256"
hybanno$patientid[which(hybanno$patientid=="98-U04607")] = "98-04607"
hybanno$patientid[which(hybanno$patientid=="98-U04349")] = "98-04349"
hybanno$patientid[which(hybanno$patientid=="98-U00100")] = "98-00100"
hybanno$patientid[which(hybanno$patientid=="97-09688")] = "97-U09688"

stopifnot(all(hybanno$patientid %in% pat$patientid))
save(hybanno, file="hybanno.rda")

cat( "table(table(hybanno$patientid))\n")
print(table(table(hybanno$patientid)))

#------------------------------------------------------------
# read intensity files
#------------------------------------------------------------
if (read.intensities) {
  col.spotanno = c("Spot.labels")
  col.qua      = c("VOL...Levels.x.mm21","Bkgd4", "VOL...Levels.x.mm2", "Bkgd") ## Data, Control
  col.diff     = "Diff..sVOL....Ctrl...Data"
  col.numeric  = c(col.qua, col.diff)
  
  for (h in 1:nrow(hybanno)) {
    fn  = file.path(hybdir, hybanno$filename[h])
    cat(h, fn, "\n", sep="\t")
    dat = read.table(fn, sep="\t", header=TRUE, skip=1, as.is=TRUE)
    if(h==1){
      qua      = array(NA, dim=c(nrow(dat), 4, nrow(hybanno)))
      spotanno = dat[,col.spotanno]
    } else {
      stopifnot(all(dat[,col.spotanno]==spotanno))
    }
    stopifnot(all(sapply(dat[,col.numeric], is.numeric)))
    for (i in 1:length(col.qua)) {
      qua[,i,h] = dat[, col.qua[i]]
    }
    
    ## check for consistency
    ch1 = qua[,1,h] - qua[,2,h]
    ch2 = qua[,3,h] - qua[,4,h]
    ch1[ch1<0] = 0
    ch2[ch2<0] = 0
    stopifnot(all(abs(ch2 - ch1 - dat[, col.diff]) < 1e-2))
  }

  dimnames(qua) = list(1:8704, c("fg.green","bg.green", "fg.red", "bg.red"), 1:175)
  save(qua,      file="qua.rda")
}

## save(spotanno, file="spotanno.Rdata") #is now made in makeeset.holger.R
