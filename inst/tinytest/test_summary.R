### Test calcIBD summary.

## Define file locations.
SxMloc <- system.file("extdata/SxM", "SxM_geno.txt", package = "statgenIBD")
SxMmap <- system.file("extdata/SxM", "SxM_map.txt", package = "statgenIBD")

SxMIBD <- calcIBD(poptype = "DH", locfile = SxMloc, mapfile = SxMmap)

SxMSumm <- capture.output(summary(SxMIBD))
expect_true(any(grepl("population type:  DH", SxMSumm)))
expect_true(any(grepl("Number of evaluation points:  116", SxMSumm)))
expect_true(any(grepl("Number of individuals:  150", SxMSumm)))
expect_true(any(grepl("Parents:  Morex Steptoe", SxMSumm)))
