#-------------------------------------------------------------------------------
# UTF-8
# 
# A simple function for a batch assessment of the correlation between
# two spatial processes for a set of raster files.
# It relies on the modified.ttest function from the SpatialPack package.
# 
# Arguments:
# 
# data.path - a path to the directory with the rasters (in GeoTIFF format
# and projected CRS).
# The rasters are expected to have the same extent and resolution.
# In case the resolutions of the rasters are different, pixel coordinates of
# the raster with the smaller resolution (for example, 90 m) are used to sample
# the raster with the larger resolution (for example, 30 m).
# 
# n.class - a single number giving the number of cells for Moran's index.
# The default is 13. If this argument is NULL Sturges' formula is used.
# 
# Returns:
# RDS-file with the list of the results (just for backup).
# 
# CSV-file with the sample correlation coefficient (corr) and additional
# parameters (consult the modified.ttest help for naming conventions
# and meaning) and the overall number of data points (N) for each pair of
# the rasters.
#-------------------------------------------------------------------------------


batch.correlations <- function(data.path = "path to the directory", n.class = 13){
    
    # load libraries
    library(raster)
    library(SpatialPack)
    
    # get a list of GeoTIFF rasters from the working directory and
    # all possible unique pairwise combinations
    setwd(data.path)
    files <- list.files(pattern = "\\.tif$")
    xy.df <- expand.grid(files, files, stringsAsFactors = F)
    xy.df <- xy.df[!(xy.df[, 1] == xy.df[, 2]), ]
    xy.df <- data.frame(t(apply(xy.df, 1, sort)), stringsAsFactors = F)
    xy.df <- xy.df[!duplicated(xy.df),]
    # add columns for resolutions
    res.names          <- c("res.max", "res.min")
    xy.df[, res.names] <- NA
    
    # check rasters' resolutions, if they are equal then NA, else name of a raster
    for (i in 1:nrow(xy.df)){
        res.x <- res(raster(xy.df[i, 1]))[1]
        res.y <- res(raster(xy.df[i, 2]))[2]
        if (res.x > res.y) {
            xy.df[i, 3] <- xy.df[i, 1]
            xy.df[i, 4] <- xy.df[i, 2]
        } else if (res.x < res.y) {
            xy.df[i, 3] <- xy.df[i, 2]
            xy.df[i, 4] <- xy.df[i, 1]
        }
    }
    
    # compute pairwise correlations
    correlations.list <- list()
    for (i in 1:nrow(xy.df)){
        if (is.na(xy.df[i, 3])) {
            raster1 <- raster(xy.df[i, 1])
            raster2 <- raster(xy.df[i, 2])
            coords  <- as.data.frame(raster1, xy = T, na.rm = T)[, 1:2]
            var1    <- as.data.frame(raster1, na.rm = T)[, 1]
            var2    <- extract(raster2, coords) # in case rasters have the same extent and resolution,
            #but different number and location of NA cells
        } else {
            raster1 <- raster(xy.df[i, 3])
            raster2 <- raster(xy.df[i, 4])
            coords  <- as.data.frame(raster1, xy = T, na.rm = T)[, 1:2]
            var1    <- as.data.frame(raster1, na.rm = T)[, 1]
            var2    <- extract(raster2, coords)
        }
        correlations.list[[i]] <- modified.ttest(x = var1, y = var2,
                                                 coords = coords,
                                                 nclass = n.class)
        message <- paste("Processed pair", i, "out of", nrow(xy.df), "pairs")
        print(message)
        saveRDS(correlations.list, "correlations.rds")
    }
    
    # create a data frame to store the results
    correlations.df <- data.frame(variable1 = xy.df[, 1], variable2 = xy.df[, 2])
    params.names    <- c("corr", "ESS", "Fstat", "dof", "p.value", "N")
    correlations.df[, params.names] <- NA
    
    # write the results into the dataframe
    for (i in 1:length(correlations.list)) {
        correlations.df[i,3] <- correlations.list[[i]]$corr
        correlations.df[i,4] <- correlations.list[[i]]$ESS
        correlations.df[i,5] <- correlations.list[[i]]$Fstat
        correlations.df[i,6] <- correlations.list[[i]]$dof
        correlations.df[i,7] <- correlations.list[[i]]$p.value
        correlations.df[i,8] <- correlations.list[[i]]$dims[1]
    }
    
    # export output to the csv-file
    write.table(correlations.df, file = "correlations.csv", row.names = F,
                sep = ",", fileEncoding = "UTF-8")
    print("Done!")
}