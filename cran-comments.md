## Release with some additions and bug fixes.

- A new plot option is added and a few minor bugs are fixed.

----

## Test environments

* local Windows 10 install, R 4.1.2
* winbuilder (release)
* macbuilder (release)
* Ubuntu (on github actions, devel and release)
* macOS (on github actions, release)
* R-hub (devel and release)

----

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE

Found the following (possibly) invalid URLs:
  URL: https://www.jstor.org/stable/29713
    From: DESCRIPTION
    Status: 403
    Message: Forbidden

This link works fine when opened in a browser.
