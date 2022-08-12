### Test readIBDs function.

## Define input file.
ibdFile <- system.file("extdata/SxM", "SxM_IBDs.ibd", package = "statgenIBD")

## Checks for correct input.
expect_error(readIBDs(infile = 1),
             "infile path should be a single character string.")
expect_error(readIBDs(infile = "a/b.txt"),
             "infile not found.")

## Check for successful read and returned structure.
expect_silent(ibds <- readIBDs(infile = ibdFile))
expect_inherits(ibds, "IBDprob")
expect_inherits(ibds$markers, "array")
expect_equal(dim(ibds$markers), c(150, 116, 2))
expect_equal(length(ibds$parents), 2)
