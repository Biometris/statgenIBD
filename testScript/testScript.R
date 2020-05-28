library(statgenIBD)

# example biparental cross:
setwd("C:/Projects/R packages/calcIbd/SimExamples/popF4")

res1 <- calcIBD(poptype     = "F4",
                locfile     = "cross.loc",
                mapfile     = "mapfile.map",
                #evalposfile = "eval.txt",
                evaldist = 3
)
print(res1$markers, digits=4)
# sum of prob should be equal to one:
all.equal(rowSums(res1$markers), rep(ncol(res1$markers), nrow(res1$markers)))

getQTL(res1, "M1_1")
getQTL(res1, c("M1_1", "M1_6"))

# example three-way cross:
setwd("C:/Projects/R packages/calcIbd/SimExamples/popC3S4DH")

res2 <- calcIBD(poptype     = "C3S4DH",
                locfile     = "cross.loc",
                mapfile     = "mapfile.map",
                evalposfile = NULL
)
print(res2$markers, digits=4)
# sum of prob should be equal to one:
all.equal(rowSums(res2$markers), rep(ncol(res2$markers), nrow(res2$markers)))

getQTL(res2, "M1_1")
getQTL(res2, c("M1_1", "M1_6"))

setwd("C:/Projects/R packages/calcIbd/SimExamples/popC4S3")

res3 <- calcIBD(poptype     = "C4S3",
                locfile     = "cross.loc",
                mapfile     = "mapfile.map",
                evalposfile = "eval.txt")

print(res3$markers, digits=4)
# sum of prob should be equal to one:
all.equal(rowSums(res3$markers), rep(ncol(res3$markers), nrow(res3$markers)))

## Doesn't work here, all QTL's are named __EVALPOS
getQTL(res3, "M1_1")
getQTL(res3, c("M1_1", "M1_6"))

