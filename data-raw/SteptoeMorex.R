## Load phenotypic data.
phenoFile <- system.file("extdata", "SxM_pheno.csv", package = "statgenIBD")
SteptoeMorexPheno <- read.csv(phenoFile)

usethis::use_data(SteptoeMorexPheno, overwrite = TRUE)
