##------------------------------------------------------------
## normalize the set of hybs selcted by Holger (1 per patient,
## see selected_chiphybs.txt)
## 3 different methods
## - "default"
## - "fgonly"
## - "bgfilt" (actually used, smoothed bg)
##------------------------------------------------------------
library(vsn)

if(!exists("qua"))      load("qua.rda")
if(!exists("spotanno")) load("spotanno.rda")
if(!exists("hybanno")) load("hybanno.rda")
## subtbg: use background as is
if(FALSE) {
  bqua  = cbind(qua[,"fg.green",]-qua[,"bg.green",], qua[,"fg.red",]-qua[,"bg.red",])
  ny175.subtbg = vsn(bqua)
  save(ny175.subtbg, file="ny175.rda")
}

## fgonly: do not use background
if(FALSE) {
  ny175.fgonly = vsn(cbind(qua[,"fg.green",], qua[,"fg.red",]))
  save(ny175.fgonly, file="ny175.fgonly.rda")
}

## bgfilt: use running-median-filtered background
nrcol = 17
nrrow = 16
stopifnot(all(sort(unique(spotanno$Column)) == 1:nrcol))
stopifnot(all(sort(unique(spotanno$Row))    == 1:nrrow))

## physical 2d coordinates on the slide
nrbl = 4
spotanno$x = spotanno$Column + ((spotanno$Block-1)  %% nrbl) * nrcol
spotanno$y = spotanno$Row    + ((spotanno$Block-1) %/% nrbl) * nrrow

## make sure x and y make sense
tmp = read.table("Echip/00-P09206_E44-1.txt", sep='\t', header=T, skip=1, as.is=T)
graphics.off(); x11(width=7, height=7); par(mfrow=c(2,2))
plot(tmp$Pos.X...mm,  spotanno$x, pch=".")
plot(tmp$Pos.Y...mm,  spotanno$y, pch=".")
plot(tmp$Pos.X...mm2, spotanno$x, pch=".")
plot(tmp$Pos.Y...mm3, spotanno$y, pch=".")

epsilon = 1.5^2  ## neighboorhood size
nrslide = dim(qua)[3]
bqua    = matrix(NA, nrow=nrow(spotanno), ncol=2*nrslide)
for (k in 1:nrow(spotanno)) {
  nb <- (spotanno$x - spotanno$x[k])^2 + (spotanno$y - spotanno$y[k])^2 <= epsilon
  #cat(sprintf("%5d:%2d ", as.integer(k), as.integer(length(which(nb)))))
  for(h in 1:dim(qua)[3]) {
    bqua[k, h        ] = qua[k, "fg.green", h] - median(qua[nb, "bg.green", h])
    bqua[k, h+nrslide] = qua[k, "fg.red",   h] - median(qua[nb, "bg.red",   h])
  }
}

hyb.sel = rep(F, nrslide)
holger.tab = read.table("selected_chiphybs.txt", header=T, sep="\t")
hyb.sel[holger.tab[,1]] = T
hyb.sel[which(hybanno$patientid=="92-26315")] = F #unknown subtype
hyb.sel[which(hybanno$patientid=="01-U04275")] = F #is a metastasis
ny74.bgfilt = vsn(bqua[, c(which(hyb.sel), which(hyb.sel)+nrslide)])
save(ny74.bgfilt, hyb.sel, file="ny74.bgfilt.rda")
