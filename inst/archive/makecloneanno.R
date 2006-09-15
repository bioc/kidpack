
  spotanno = read.table("Echip/Kidney-2-E(GM2).txt", sep="\t", skip=0, header=T, as.is=T)
  colnames(spotanno)[9] = "vendor"
  colnames(spotanno)[which(colnames(spotanno)=="Name")] = "description"
  spotanno = cbind(spotanno, rep(NA, nrow(spotanno)))

  ## bring into the same order as the intensity files:
  int.index = numeric(nrow(spotanno))
  for (spot in 1:nrow(spotanno)) 
    int.index[spot] = (spot-1)%/%(17*16*4)*(17*16*4) + (spotanno$Row[spot]-1)*17*4 + (spotanno$Block[spot]-1)%%4*17 + spotanno$Column[spot]
  ## the number of completed "block rows" + the number of completed rows in the new "block row" + the current columnn
  ## compute the inverse permutation:
  new = numeric(nrow(spotanno))
  new[int.index] = 1:8704
  spotanno = spotanno[new,]

  ## add a column for the cloneanno index:
  colnames(spotanno)[ncol(spotanno)] = "probe"
  cloneanno = data.frame(plate=numeric(nprobe),
                         SrcRow=numeric(nprobe),SrcCol=numeric(nprobe),
                         imageid=numeric(nprobe),
                         AccNumber=I(character(nprobe)),
                         spot1=numeric(nprobe),
                         spot2=numeric(nprobe),
                         description=I(character(nprobe)),
                         vendor=I(character(nprobe)))
  z=0
  for (j in unique(spotanno$SrcPlt)) {
    for (k in unique(spotanno$SrcRow)) {
      for (l in unique(spotanno$SrcCol)) {
        #256 spots with spotanno$SrcPlt=0 are "Blank"
        spots = which(spotanno$SrcPlt == j & spotanno$SrcRow == k & spotanno$SrcCol == l)
        stopifnot(length(spots) %in% c(0,2) | spotanno$SrcPlt[spots[1]]==0)
        if (length(spots) == 2) {
          #includes 24 Kanamycin and 34 empty spots
          z=z+1
          spotanno$probe[spots] = z 
          cloneanno$plate[z] = spotanno$SrcPlt[spots[1]]
          cloneanno$SrcRow[z] = spotanno$SrcRow[spots[1]]
          cloneanno$SrcCol[z] = spotanno$SrcCol[spots[1]]
          cloneanno$imageid[z] = spotanno$ImageID[spots[1]]
          cloneanno$AccNumber[z] = spotanno$AccNumber[spots[1]]
          cloneanno$description[z] = spotanno$Name[spots[1]]
          cloneanno$vendor[z] = spotanno$vendor[spots[1]]
          cloneanno$spot1[z] = spots[1]
          cloneanno$spot2[z] = spots[2]
        }
      }
    }
  }
  save(spotanno, file="spotanno.rda")
  save(cloneanno, file="cloneanno.rda")

