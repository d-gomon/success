# success 1.1.0

* Added references to published article in The R Journal
* Removed print statements from bernoulli_ARL_* functions 

# success 1.0.1

* Added bernoulli_ARL, bernoulli_RL_cdf functions which allow to calculate ARL and RL distribution functions using 2 different approaches 
  * Some tests put into place so that these functions work as expected
* Fixed figure display in README on CRAN
* Added some references
* funnel plot functions changed to display unit labels, as well as prediction intervals by using shaded ribbons instead of lines.
* Added tests for most functions in the package
* Added CITATION file


# success 1.0.0

* bernoulli_cusum() no longer re-calculates failure probability at every time point. Major speed-up
* Added references to published articles in cgr_cusum() and funnel_plot().
* Added baseline, cumulative and inverse cumulative Weibull hazards (see ?haz_weib)
* Updated README file with references.
* Interactive plot now correctly displays downward CUSUM charts.
* Added Vignette describing the package.



# success 0.1.1

* Added a `NEWS.md` file to track changes to the package.
* Fixed lower CGR_CUSUM to calculate correct values
* Changed the way in which upper and lower BK- and CGR-CUSUMs are visualised.
  * Both charts now recalculate the highest and lowest value at failure times, aiding interpretation of the charts, but increasing computation time 2 fold for CGR-CUSUMs
* Reduced computation time of CGR-CUSUM approximately 4 fold (depends on size of data) using Rfast package. 
  * Future: Rfast binary_search() is not optimal for this problem. Either write an adjusted binary_search fitting for the problem manually or look for other implementations.
* Created a hexagon for the package (can be found in /hexagon)
* Added `breast` data set based on real cancer trial.
