## Check with winbuilder - develop only.
devtools::check_win_devel()

## Check with rhub.
rhub::check_for_cran()

## Check with rhub.
rhub::check_for_cran(platforms = c("ubuntu-gcc-release",
                                   "fedora-clang-devel"),
                     path = "C:/Projects/R_packages/statgenIBD/")

## Rebuild readme.
devtools::build_readme()

## Submit to CRAN
devtools::release()


## Build site for local check.
pkgdown::clean_site()
pkgdown::build_site()

## Code coverage - local.
detach("package:statgenIBD", unload = TRUE)
covr::gitlab()


