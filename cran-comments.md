## Initial release

- Initial package release, updated DESCRIPTION following remarks from previous submission

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
   
New submission

Possibly mis-spelled words in DESCRIPTION:
  Biometris (10:18)
  IBD (3:23, 4:50)
  biparental (5:5)

These are all spelled correctly.  
   
Found the following (possibly) invalid URLs:
  URL: https://www.jstor.org/stable/29713
    From: DESCRIPTION
    Message: Forbidden
    Status: 403
  
This link is fine when opened from a browser.

* checking installed package size ... NOTE
  installed size is  9.1Mb
  sub-directories of 1Mb or more:
    libs   8.5Mb

Installed size depends on platform. This not only occurs at some of the tested platforms.

  
