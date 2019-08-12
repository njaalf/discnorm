pc_freq <- function(Y1, Y2) {#from lavaan source
  max.y1 <- max(Y1, na.rm=TRUE); max.y2 <- max(Y2, na.rm=TRUE)
  bin <- Y1 - 1L; bin <- bin + max.y1 * (Y2 - 1L); bin <- bin[!is.na(bin)]
  if (length(bin)) bin <- bin + 1L
  array(tabulate(bin, nbins = max.y1*max.y2), dim=c(max.y1, max.y2))
}