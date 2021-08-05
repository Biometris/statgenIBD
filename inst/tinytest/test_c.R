### Test IBDprob concatenation.

## Define file locations.
ABloc <- system.file("extdata/multipop", "AxB.txt", package = "statgenIBD")
ACloc <- system.file("extdata/multipop", "AxC.txt", package = "statgenIBD")
ABCmap <- system.file("extdata/multipop", "mapfile.txt", package = "statgenIBD")

## IBD calculations for two populations separately.
AB <- calcIBD(poptype = "F4DH", locfile = ABloc, mapfile = ABCmap)
AC <- calcIBD(poptype = "F4DH", locfile = ACloc, mapfile = ABCmap)

## Alternative calculations for AC to test input checks.
ACalt1 <- calcIBD(poptype = "F4", locfile = ACloc, mapfile = ABCmap)
ACalt2 <- calcIBD(poptype = "F4DH", locfile = ACloc, mapfile = ABCmap,
                  evaldist = 5)

## Check that input checks are working correctly.
expect_error(c(AB, "tst"),
             "All inputs should be of class IBDprob")
expect_error(c(AB, ACalt1),
             "All inputs should have the same population type")
expect_error(c(AB, ACalt2),
             "All inputs should have the same map")

## c with single argument should return input.
expect_equal(c(AB), AB)

## Multiple outputs should be combined correctly.
ABC <- c(AB, AC)

expect_inherits(ABC, "IBDprob")
expect_equal(ABC$map, AB$map)
expect_equal_to_reference(ABC$markers, "ABC_markers")
expect_equal(ABC$parents, c("A", "B", "C"))
expect_equal(ABC$poptype, AB$poptype)
expect_true(ABC$multicross)

## Check that genoCross attribute is added correctly.
genoCross <- attr(x = ABC, which = "genoCross")

expect_inherits(genoCross, "data.frame")
expect_equal(genoCross[["cross"]],
             rep(c("cross1", "cross2"), times = c(100, 80)))
expect_equal(genoCross[["geno"]],
             c(dimnames(AB$markers)[[2]], dimnames(AC$markers)[[2]]))






