# statgenIBD 1.0.4

* Increased speed of IBD calculations by improving an inefficient part of the algorithm.

# statgenIBD 1.0.3

* It is now possible to plot the pedigree of an object of class `IBDprob` by specifying `plotType = "pedigree"`.
* A small bug in the calculation of probabilities for populations C3 and C4 is fixed.

# statgenIBD 1.0.2

* It is now possible to specify chromosomes in character format, i.e. 1a, 1b, 1c as is common in certain crops.
* The plot function for calculated IBD probabilities now has an extra option `plotType` that can have values "singleGeno", for the original plot, and "allGeno" for a new plot showing probabilities for all genotypes in a single plot.
* A bug that in some cases could lead to misnaming of the markers when combining multiple imputed populations is fixed.

# statgenIBD 1.0.1

* A bug where written Flapjack files could not be read is fixed.

# statgenIBD 1.0.0

* Initial CRAN version
