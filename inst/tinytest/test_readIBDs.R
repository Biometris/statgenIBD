### Test readIBDs function.

## Define input file.
ibdFile <- system.file("extdata/SxM", "SxM_IBDs.ibd", package = "statgenIBD")

## Checks for correct input.
expect_error(readIBDs(file = 1), "File path should be a single character string.")
expect_error(readIBDs(file = "a/b.txt"), "File not found.")

## Check for successful read and returned structure.
expect_silent(ibds <- readIBDs(file = ibdFile))
expect_inherits(ibds, "IBDprob")
expect_inherits(ibds$markers, "array")
expect_equal(dim(ibds$markers), c(116, 150, 2))
expect_equal(length(ibds$parents), 2)
