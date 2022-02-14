## C3 and C4 crash when computing probabilities.
pops <- c("DH", "F4", "F4DH", "BC3", "BC3DH", "BC1S2", "BC1S2DH", "C3",
          "C3DH", "C3S4", "C3S4DH", "C4", "C4DH", "C4S3", "C4S3DH")

## input file doesn't accept relative paths.
popDirs <- paste0("C:/Projects/R_packages/statgenIBD/testScripts/simQTL/pop", pops)

## Simulate populations.
for (popDir in popDirs) {
  statgenIBD::simQTL(inputfile = paste0(popDir, "/sim_script.txt"),
                     dir_name = popDir)
}

## Calculate IBDs.
popIBDs <- lapply(X = seq_along(pops), FUN = function(i) {
  cat(pops[i], "\n")
  calcIBD(popType = pops[i],
          markerFile = paste0(popDirs[i], "/cross.loc"),
          mapFile = paste0(popDirs[i], "/mapfile.map"))
})

## Plot pedigrees.
pdf("C:/Projects/R_packages/statgenIBD/testScripts/simQTL/pedigrees.pdf")
for (popIBD in popIBDs) {
  plot(popIBD, plotType = "pedigree", title = popIBD$popType)
}
dev.off()






