## Patch release

- Patch release to fix problems with CRAN check for gcc11

----

## Test environments

* local Windows 10 install, R 4.1
* winbuilder (release)
* Ubuntu (on github actions, devel and release)
* macOS (on github actions, release)
* R-hub (devel and release)

----

## R CMD check results

There were no ERRORs or WARNINGs.

There were 2 NOTES:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Bart-Jan van Rossum <bart-jan.vanrossum@wur.nl>'

Found the following (possibly) invalid URLs:
  URL: https://www.jstor.org/stable/29713
    From: DESCRIPTION
    Status: 403
    Message: Forbidden

This link works fine when opened in a browser.

* checking installed package size ... NOTE
  installed size is  9.1Mb
  sub-directories of 1Mb or more:
    libs   8.5Mb

Installed size depends on platform. This only occurs at some of the tested platforms.    
