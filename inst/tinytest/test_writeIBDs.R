### Test writeIBDs function.

## Define file locations.
SxMloc <- system.file("extdata/SxM", "SxM_geno.txt", package = "statgenIBD")
SxMmap <- system.file("extdata/SxM", "SxM_map.txt", package = "statgenIBD")

## Compute IBDs
SxMIBD <- calcIBD(popType = "DH", markerFile = SxMloc, mapFile = SxMmap)

## Check that input checks are working correctly.
expect_error(writeIBDs(SxMIBD, outFile = 1),
             "outFile should be a single character string ending in .ibd")
expect_error(writeIBDs(SxMIBD, outFile = "a/b.ibd"),
             "No permission to write to")

## Define tempt output file.
ibdOut <- tempfile(fileext = ".ibd")

expect_error(writeIBDs(SxMIBD, outFile = ibdOut, decimals = "a"),
             "decimals should be a single numerical value")
expect_error(writeIBDs(SxMIBD, outFile = ibdOut, minProb = 0.75),
             "minProb should be a numerical value between 0 and")

## There is very little that can actually be tested.
## Just checking silent execution.
expect_silent(writeIBDs(SxMIBD, outFile = ibdOut))