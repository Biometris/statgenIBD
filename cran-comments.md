## Minor release

- Some extra functionality and an extra vignette are added.

----

## Test environments

* local Windows 10 install, R 4.2.1
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
  URL: https://doi.org/10.1073/pnas.1100465108
    From: inst/doc/IBDCalculations.html
    Status: 503
    Message: Service Unavailable
    
  URL: https://www.jstor.org/stable/29713
    From: DESCRIPTION
          man/statgenIBD-package.Rd
          inst/doc/IBDCalculations.html
    Status: 403
    Message: Forbidden
    
These links works fine when opened in a browser.
    
