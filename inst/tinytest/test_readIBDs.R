### Test readIBDs function.

## Define input file.
ibdFile <- system.file("extdata/SxM", "SxM_IBDs.txt", package = "statgenIBD")
mapFile <- system.file("extdata/SxM", "SxM_map.txt", package = "statgenIBD")
map <- read.delim(mapFile, header = FALSE)
rownames(map) <- map$V1
map <- map[,-1]
colnames(map) <- c("chr", "pos")

## Checks for correct input.
expect_error(readIBDs(infile = 1),
             "should be a character string indicating a readable .txt file")
expect_error(readIBDs(infile = "a/b.txt"),
             "should be a character string indicating a readable .txt file")

## Check for successful read and returned structure.
expect_silent(ibds <- readIBDs(infile = ibdFile, map = map))
expect_inherits(ibds, "IBDprob")
expect_inherits(ibds$markers, "array")
expect_equal(dim(ibds$markers), c(150, 116, 2))
expect_equal(length(ibds$parents), 2)
