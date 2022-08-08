library(statgenIBD)

popDir <- "C:/Projects/R_packages/statgenIBD/testScripts/simQTL/popBC1S2-2/"


QTLdf <- data.frame(name = c("qtl1", "qtl2"), chr = c(1, 1),
                    pos = c(30, 75), add = c(1, 1), dom = c(1, 0))

fnddf1 <- data.frame(name = c("BLUE", "YELLOW"),
                    e1 = c("+", "-"), e2 = c("-", "+"))

epiMat <- matrix(c(0,1,1,0), nrow = 2)

fnddf2 <- data.frame(name = c("BLUE", "YELLOW"))

# Two QTLs
statgenIBD::simQTL(inbFnd = fnddf1,
                   popType = "BC1S2",
                   nInd = 13,
                   dir_name = popDir,
                   start_seed = -1234,
                   chrLength = c(100.0, 120),
                   nlocChr = c(11.0, 14),
                   QTLPos = QTLdf,
                   epiInt = epiMat,
                   nr_alleles = 3)

# No QTLs
statgenIBD::simQTL(inbFnd = fnddf2,
                   popType = "BC1S2",
                   nInd = 13,
                   dir_name = popDir,
                   start_seed = -1234,
                   chrLength = c(100.0, 120),
                   nlocChr = c(11.0, 14),
                   nr_alleles = 3)


statgenIBD::simQTL(inbFnd = fnddf2,
                   popType = "BC1S2",
                   nInd = 13,
                   dir_name = popDir,
                   start_seed = -1234,
                   chrLength = c(100.0, 120),
                   nlocChr = c(11.0, 14),
                   nr_alleles = 0)
