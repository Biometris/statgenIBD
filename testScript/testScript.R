library(statgenIBD)

# example biparental cross:
setwd("C:/Projects/R packages/calcIbd/SimExamples/popF4")

res1 <- calcIBD(poptype     = "F4",
                locfile     = "cross.loc",
                mapfile     = "mapfile.map",
                evalposfile = "eval.txt",
                evaldist = 3
)
print(res1$markers, digits=4)
# sum of prob should be equal to one:
all.equal(rowSums(res1$markers), rep(ncol(res1$markers), nrow(res1$markers)))

getQTLProb(res1, "M1_1")
getQTLProb(res1, c("M1_1", "M1_6"))

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

getQTLProb(res2, "M1_1")
getQTLProb(res2, c("M1_1", "M1_6"))

setwd("C:/Projects/R packages/calcIbd/SimExamples/popC4S3")

res3 <- calcIBD(poptype     = "C4S3",
                locfile     = "cross.loc",
                mapfile     = "mapfile.map",
                evalposfile = "eval.txt")

print(res3$markers, digits=4)
# sum of prob should be equal to one:
all.equal(rowSums(res3$markers), rep(ncol(res3$markers), nrow(res3$markers)))

getQTLProb(res3, "EVAL_1_0")
getQTLProb(res3, c("EVAL_1_0", "EVAL_1_10"))

