Tools to describe spatial structures and test for the presence of spatial correlation through univariate structure functions (correlograms, variograms).

## Scripts

*./scripts*
* *batch_correlations.R* -- batch assessment of the correlation between two spatial processes for a set of raster files; relies on the modified.ttest function from the SpatialPack package.
* *geostats_points_function.R* -- the function performs preliminary data analysis and calculates the parameters necessary for the assesmant of semivariance and autocorrelation, namely sample size, lag, number of distance classes. It also creates four directional profiles and a random sample.
* *variograms.R* -- using variog4 function from geoR library the script calculates 5 variograms (4 directional + 1 omnidirectional).
* *correlograms.R* -- using eco.correlog function from EcoGenetics library the script calculates 5 correlograms (4 directional + 1 omnidirectional). Both Moran's I and Geary's c autocorrelations can be computed.


## Usage

```
source(geostats_points_function.R)
geostats.points.function(raster.data = "raster.tif") #sample.size = 1000, dist.of.interest = NULL, lag.of.interest = NULL
source(variograms.R)
source(correlograms.R)
```

---

*./data*

the catalogue contains some output data examples, produced during the testing

---

Read more in this [presentation](https://www.slideshare.net/DariaSvidzinska/an-approach-to-the-rapid-univariate-analysis-of-spatial-pattern-of-landscape-gradients)

*./2017_vnz-conf_presentation*

all the presentation-related stuff


