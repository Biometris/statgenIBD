## Release with 'bug' fixes.

- Though not really a bug fix, after the latest release we noticed a major bottleneck in the code that slowed down calculations a lot. Calculations are now sped up by about a factor 3.

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
    
Found the following (possibly) invalid DOIs:
  DOI: 10.1073/pnas.1100465108
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503    

The DOI is correct. The corresponding link can be opened in a browser.
