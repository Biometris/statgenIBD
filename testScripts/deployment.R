## Check with winbuilder - develop only.
devtools::check_win_devel()

## Check with rhub.
rhub::check_for_cran(env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "false"))

## Rebuild readme.
devtools::build_readme()

## Submit to CRAN
devtools::release()


## Build site.
# pkgdown::clean_site()
# pkgdown::build_site()

detach("package:statgenIBD", unload = TRUE)
covr::gitlab()


SxMIBD_Ext <- calcIBD(poptype = "DH",
                      locfile = system.file("extdata/SxM", "SxM_geno.txt",
                                            package = "statgenIBD"),
                      mapfile = system.file("extdata/SxM", "SxM_map.txt",
                                            package = "statgenIBD"),
                      evaldist = 5,
                      grid = FALSE)
plot(SxMIBD_Ext, genotype = "dh016")
