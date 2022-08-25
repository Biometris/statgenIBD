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

## Check reverse dependencies.
# First build source package and put it in testScripts/revDep folder.
tools::check_packages_in_dir("testScripts/revDep", reverse = list())
