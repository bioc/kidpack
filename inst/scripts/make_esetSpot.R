## *** This script ran exactly ONE to convert the serialized old style
##     (pre R 2.4) esetSpot object (file) into one that is up-to-date.
##
##                       wh. 9.9.2007
##
library("Biobase")

## assume we are in kidpack/inst/script and .rda files come and go to
##   kidpack/data
datadir = function(x) file.path("..", "..", "data", x)

load(datadir("eset.rda"))

## the old, pre R2.4 object, as of from package version 1.4.3:
load(datadir("esetSpot.rda"))

stopifnot(nrow(exprs(eset))    ==4224L,
          nrow(exprs(esetSpot))==8448L)

ex = exprs(esetSpot)
e1 = (ex[1:4224,] +  ex[(1:4224)+4224,])/2
e2 =  exprs(eset)
stopifnot(identical(colnames(e1), colnames(e2)),
          max(abs(e1-e2))<0.02)

rownames(ex) = paste(seq_len(nrow(ex)), rownames(ex), sep=" ")
esetSpotNew = new("ExpressionSet",
  phenoData = new("AnnotatedDataFrame", data=pData(esetSpot)),
  exprs=ex)


esetSpot = esetSpotNew
save(esetSpot, file=datadir("esetSpot.rda"))
