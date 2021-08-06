## Check with winbuilder - develop only.
devtools::check_win_devel()

## Check with rhub.
rhub::check_for_cran(env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "false"))

## Rebuild readme.
devtools::build_readme()

## Submit to CRAN
devtools::release()


## Build site.
pkgdown::clean_site()
pkgdown::build_site()

detach("package:statgenIBD", unload = TRUE)
covr::gitlab()


