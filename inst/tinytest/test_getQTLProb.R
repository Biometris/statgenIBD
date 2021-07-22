### Test getQTLProb function..

## Define file locations.
ABloc <- system.file("extdata/multipop", "AxB.loc", package = "statgenIBD")
ACloc <- system.file("extdata/multipop", "AxC.loc", package = "statgenIBD")
ABCmap <- system.file("extdata/multipop", "mapfile.map", package = "statgenIBD")

## IBD calculations for two populations separately.
AB <- calcIBD(poptype = "F4DH", locfile = ABloc, mapfile = ABCmap)
AC <- calcIBD(poptype = "F4DH", locfile = ACloc, mapfile = ABCmap)
ABC <- c(AB, AC)

## Check that input checks are working correctly.
expect_error(getQTLProb(IBD = "tst", QTLS = "M1_1"),
             "should be an object of class calcIBD")
expect_error(getQTLProb(IBD = AB, QTLS = 1),
             "QTLS should be a character vector")
expect_error(getQTLProb(IBD = AB, QTLS = "M1"),
             "The following QTLs are not in AB: M1")

## Check that output is correct.
AB_M1_1 <- getQTLProb(IBD = AB, QTLS = "M1_1")

expect_inherits(AB_M1_1, "data.frame")
expect_equal_to_reference(AB_M1_1, "AB_M1_1")

AB_M1_1_M3_3 <- getQTLProb(IBD = AB, QTLS = c("M1_1", "M3_3"))

expect_equal(AB_M1_1, AB_M1_1_M3_3[, 1:3])

## For multicross there should be an extra column cross.

ABC_M1_1 <- getQTLProb(IBD = ABC, QTLS = "M1_1")
expect_equal(colnames(ABC_M1_1),
             c("cross", "geno", "M1_1_pA", "M1_1_pB", "M1_1_pC"))

