## Check with winbuilder - develop only.
devtools::check_win_devel()
## Check with mac builder
devtools::check_mac_release()

rhub::platforms()

## Check with rhub.
rhub::check_for_cran(path = "C:/Projects/R_packages/statgenIBD/")

## Check with rhub.
rhub::check_for_cran(platforms = c("ubuntu-gcc-release",
                                   "fedora-clang-devel",
                                   "debian-clang-devel",
                                   "debian-gcc-devel"),
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
library(statgenIBD)

## Check reverse dependencies. (takes about 45 minutes).
Sys.setenv(R_BIOC_VERSION = "3.20")
revdepcheck::revdep_check()
