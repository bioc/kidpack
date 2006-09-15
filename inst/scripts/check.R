compdf <- function (df1, df2, name) {
  cat("***", name, "***\n")
  for(i in 1:ncol(df1)) {
    bad <- which(df1[,i]!=df2[,i])
    if(length(bad)>0)
      cat(name, i, colnames(df1)[i], "\n", df1[bad,i], "\n", df2[bad,i], "\n")
  }
}

load("hybanno.rda")
hybanno.new <- hybanno
load("~/Kidney3/hybanno.rda")
compdf(hybanno.new, hybanno, "hybanno")

load("spotanno.rda")
spotanno.new <- spotanno
load("~/Kidney3/spotanno.rda")
compdf(spotanno.new, spotanno, "spotanno")

load("cloneanno.rda")
cloneanno.new <- cloneanno
load("~/Kidney3/cloneanno.rda")
compdf(cloneanno.new, cloneanno, "cloneanno")

load("qua.rda")
qua.new <- qua
load("~/Kidney3/qua.rda")
stopifnot(identical(qua, qua.new))

load("eset.rda")
eset.new <- eset
load("~/Kidney3/eset74.rda")
stopifnot(identical(exprs(eset), exprs(eset.new)))
compdf(pData(eset.new), pData(eset), "eset")

