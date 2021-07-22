### Test calcIBD function

## Define file locations.
SxMloc <- system.file("extdata/SxM", "SxM_geno.txt", package = "statgenIBD")
SxMmap <- system.file("extdata/SxM", "SxM_map.txt", package = "statgenIBD")

## Check the input check are working correctly.

expect_error(calcIBD(poptype = "tst", locfile = SxMloc, mapfile = SxMmap),
             "unknown type tst")
expect_error(calcIBD(poptype = "DH", locfile = "tst", mapfile = SxMmap),
             "Cannot open file tst")
expect_error(calcIBD(poptype = "DH", locfile = SxMloc, mapfile = "tst"),
             "Cannot read file tst")

## Check that the output structure is correct.
expect_silent(SxMIBD <- calcIBD(poptype = "DH", locfile = SxMloc,
                                mapfile = SxMmap))

expect_inherits(SxMIBD, "calcIBD")
expect_equal(names(SxMIBD), c("map", "markers", "poptype", "multiCross"))
expect_inherits(SxMIBD$map, "data.frame")
expect_inherits(SxMIBD$markers, "array")
expect_equal(dim(SxMIBD$markers), c(116, 150, 2))
expect_inherits(SxMIBD$poptype, "character")
expect_inherits(SxMIBD$multiCross, "logical")

## Check that output content is correct.

expect_equal_to_reference(SxMIBD$map, "SxMIBD_map")
expect_equal_to_reference(SxMIBD$markers, "SxMIBD_markers", tolerance = 10e-6)
expect_equal(SxMIBD$poptype, "DH")
expect_false(SxMIBD$multiCross)


