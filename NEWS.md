# success 0.1.2

* Added a `NEWS.md` file to track changes to the package.
* Fixed lower CGR_CUSUM to calculate correct values
* Changed the way in which upper and lower BK- and CGR-CUSUMs are visualised.
  * Both charts now recalculate the highest and lowest value at failure times, aiding interpretation of the charts, but increasing computation time 2 fold for CGR-CUSUMs
* Reduced computation time of CGR-CUSUM approximately 4 fold (depends on size of data) using Rfast package. 
  * Future: Rfast binary_search() is not optimal for this problem. Either write an adjusted binary_search fitting for the problem manually or look for other implementations.
* Created a hexagon for the package (can be found in /hexagon)
