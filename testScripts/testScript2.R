library(statgenIBD)

setwd(r"(C:\Projects\R packages\calcIbd\SimExamples\multipop\)")

res1 <- statgenIBD::calcIBD(poptype = "F4DH",
                            locfile     = "AxB.loc",
                            mapfile     = "mapfile.map",
                            evalposfile = "eval.txt")

res2 <- statgenIBD::calcIBD(poptype = "F4DH",
                            locfile     = "AxC.loc",
                            mapfile     = "mapfile.map",
                            evalposfile = "eval.txt")

res12 <- c(res1, res2)

df <- getQTLProb(res12, c("EVAL_1_25", "EVAL_2_75", "EVAL_3_75"))
head(df)
tail(df)


res1 <- statgenIBD::calcIBD(poptype = "F4DH",
                            locfile     = "AxB.loc",
                            mapfile     = "mapfile.map",
                            evalposfile = "eval2.txt")




microbenchmark::microbenchmark(
res1 <- statgenIBD::calcIBD(poptype = "F4DH",
                            locfile     = "AxB.loc",
                            mapfile     = "mapfile.map",
                            evalposfile = "eval.txt"))


Unit: milliseconds
expr
res1 <- statgenIBD::calcIBD(poptype = "F4DH", locfile = "AxB.loc",      mapfile = "mapfile.map", evalposfile = "eval.txt")
min       lq     mean   median       uq      max neval
565.6722 572.7497 595.4176 574.6099 578.5152 1529.842   100

