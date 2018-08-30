#-------------------------------------------------------------------------------
# UTF-8
# 
# This function performs preliminary data analysis and calculates the parameters
# necessary for the assesment of univariate structure functions, namely
# semivariance and covarance.
#
# Arguments:
#       - mandatory:
#                     - raster.data: the raster to be analyzed
#                     - sample.size: number of pixels to be sampled. Usually
#                       sample from 1000 (the default) to 3000 points is good
#                       enough to produce reliable results. Too large samples
#                       require significant computing resources and cannot be
#                       digested by GeoR and EcoGenetics libraries       
#       - optional:
#                     - dist.of.interest: distance of interest to limit the
#                       maximum distance calculations (in the measurment units
#                       of the raster projection, usually meters)
#                     - lag.of.interest: lag of interest that defines the
#                       increment of distances up to maximim distance (in the
#                       measurment units of raster projection, usually meters).
#                       While you can limit the distance of interest, defining
#                       the lag is not recommended because this virtually
#                       deprives the analysis of its sense
#
# Returns:
#       - ESRI Shapefile:
#                     - point shapefile with a random sample of the defined size
#       - RDS-files:
#                     - these files store intermediary outputs and data objects
#                       necessary  for further analysis, namely: sample geodata,
#                       information on distance bins (size, centers, breaks, sum),
#                       number of distance classes
#                       
#-------------------------------------------------------------------------------

geostats.points.function <- function(raster.data = "path to the raster",
                                     sample.size = 1000,
                                     dist.of.interest = NULL,
                                     lag.of.interest = NULL){
        
        
# 0. LOAD LIBRARIES-------------------------------------------------------------
        library(raster)  # work with rasters and sampling
        library(sp)      # spatial points and lines
        library(rgdal)   # writing outputs to ESRI shapes
        library(geoR)    # create geoR object


# 1. IMPORT RASTER AND DEFINE ITS DIAGONAL--------------------------------------
        gradient.raster <- raster(raster.data)
        # obtain x coordinates from extent
        x.min <- xmin(gradient.raster)
        x.max <- xmax(gradient.raster)
        # obtain y coordinates from extent
        y.min <- ymin(gradient.raster)
        y.max <- ymax(gradient.raster)
        # combine coordintes into corner points and calculate diagonal
        corner.nw <- c(x.min, y.max)
        corner.se <- c(x.max, y.min)
        diagonal.points   <- rbind(corner.nw, corner.se)
        diagonal.dist     <- dist(diagonal.points)[1]
              
        
# 2. CREATE A RANDOM SAMPLE-----------------------------------------------------
        random.sample.spdf <- sampleRandom(gradient.raster, size = sample.size,
                                           cells = T, xy = T, sp = T, na.rm = T)
        writeOGR(obj = random.sample.spdf, dsn ="./data", layer = "sample_points",
                 driver = "ESRI Shapefile", overwrite_layer = T)

#3. DEFINE OPTIMUM NUMBER OF CLASSES AND BINS ----------------------------------

        # calculate maximum possible number of classes
        valid.pixels <- cellStats(!is.na(gradient.raster), 'sum')
        n.pairs      <- (valid.pixels*(valid.pixels-1))/2 # the number of pairs in the distance matrix
        n.classes    <- round(1 + 3.322*log10(n.pairs)) # Sturges' rule
        dist.of.interest <- ifelse(is.null(dist.of.interest) == T,
                                   diagonal.dist/3, dist.of.interest)
        
        if (is.null(lag.of.interest) == T) {
                bin.size <- ifelse(test = ((dist.of.interest/n.classes)%%5) <= 2.5,
                                   yes = (dist.of.interest/n.classes) - ((dist.of.interest/n.classes)%%5),
                                   no = (dist.of.interest/n.classes) + (5 - ((dist.of.interest/n.classes)%%5)))
                
                # calculate bins' breaks
                bins.sum <- 0
                bins.breaks <- 0
                for (i in 1:n.classes){
                        bins.sum <- bins.sum + bin.size
                        bins.breaks <- append(bins.breaks, bins.sum)
                        if ((bins.sum + bin.size) > dist.of.interest) break
                }
                
        } else {
                bin.size <- lag.of.interest
                bins.sum <- 0
                bins.breaks <- 0
                for (i in 1:n.classes){
                        bins.sum <- bins.sum + bin.size
                        bins.breaks <- append(bins.breaks, bins.sum)
                        if ((bins.sum + bin.size) > dist.of.interest) break
                }
        }      
        bins.centers <- bins.breaks[2:length(bins.breaks)] - (bin.size/2)
        
        # create and save data for further analysis
        sample.geodata <- as.geodata(random.sample.spdf, data.col = 4)
        saveRDS(sample.geodata, "./data/sample_geodata.rds")
        saveRDS(bins.centers, "./data/bins_centers.rds")
        saveRDS(bins.breaks, "./data/bins_breaks.rds")
        saveRDS(bins.sum, "./data/bins_sum.rds")
        saveRDS(n.classes, "./data/n_classes.rds")
        saveRDS(bin.size, "./data/bin_size.rds")
        
        print("Processing is completed.")        
}
                    