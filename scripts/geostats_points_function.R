#-------------------------------------------------------------------------------
# UTF-8
# 
# This function performs preliminary data analysis and calculates the parameters
# necessary for the assesmant of univariate structure functions, namely
# semivariance and autocorrelation.
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
#       - ESRI Shapefiles:
#                     - 4 point shapefiles for 4 directional profiles (0, 45, 90,
#                       and 135 degrees respectively)
#                     - point shapefile with a random sample of the defined size
#       - CSV file:
#                     - the data sampled along the directional profiles
#                       structured in the so-called long form
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
        
        
#------------------------------------------------------------------------------#
# 0. LOAD LIBRARIES
#-------------------------------------------------------------------------------
        library(raster)  # work with rasters and sampling
        library(sp)      # spatial points and lines
        library(inlmisc) # sampling points along transects
        library(rgdal)   # writing outputs to ESRI shapes
        library(geoR)    # create geoR object
#-------------------------------------------------------------------------------

#------------------------------------------------------------------------------#
# 1. IMPORT RASTER, BUILD PROFILES AND SAMPLE DATA ALONG THEM
#-------------------------------------------------------------------------------
        gradient.raster <- raster(raster.data)
        # obtain information on projection
        proj.string <- projection(gradient.raster, asText = TRUE)
        cell.size <- res(gradient.raster)[1]
        
        # obtain x coordinates from extent
        x.min <- xmin(gradient.raster)
        x.max <- xmax(gradient.raster)
        x.cen <- x.min + ((x.max - x.min)/2)
        
        # obtain y coordinates from extent
        y.min <- ymin(gradient.raster)
        y.max <- ymax(gradient.raster)
        y.cen <- y.min + ((y.max - y.min)/2)
        
        # combine coordintes into corner points
        corner.nw <- c(x.min, y.max)
        corner.ne <- c(x.max, y.max)
        corner.se <- c(x.max, y.min)
        corner.sw <- c(x.min, y.min)
        # combine coordinates into center points
        center.n <- c(x.cen, y.max)
        center.s <- c(x.cen, y.min)
        center.e <- c(x.max, y.cen)
        center.w <- c(x.min, y.cen)
        
        # combine points into the start-end points
        points.w_e <- rbind(center.w, center.e) #0 degrees
        points.nw_se <- rbind(corner.nw, corner.se) #45 degrees
        points.n_s <- rbind(center.n, center.s) #90 degrees
        points.ne_sw <- rbind(corner.ne, corner.sw) #135 degrees
        
        # merge points and their names into a single list
        points.list <- list(points.w_e, points.nw_se, points.n_s, points.ne_sw)
        names(points.list)  <- c('deg0_w-e', 'deg45_nw-se', 'deg90_n-s', 'deg135_ne-sw')
        directions <- c(0, 45, 90, 135)
        
        
        sp <- NULL
        transect <- list()
        trans.points <- list()
        profiles.df <- data.frame(direction = character(0),
                                  distance = numeric(0),
                                  value = numeric(0))
        obs.n <- NULL
        for (i in 1:length(points.list)){
                sp <- SpatialPoints(points.list[i], proj4string = CRS(proj.string))
                transect[[i]] <- as(sp, "SpatialLines")
                trans.points[[i]] <- remove.duplicates(ExtractAlongTransect(transect[[i]], gradient.raster)[[1]])
                writeOGR(obj = trans.points[[i]], dsn ="./data",
                         layer = paste0("profile_points_", names(points.list[i])),
                         driver = "ESRI Shapefile", overwrite_layer=TRUE)
                obs.n <- nrow(trans.points[[i]]@data)
                temp.df <- data.frame(id = 1:obs.n)
                temp.df <- cbind(temp.df,
                                 rep(directions[i], obs.n),
                                 trans.points[[i]]@data[[1]],
                                 trans.points[[i]]@data[[2]])
                profiles.df <- rbind(profiles.df, temp.df[, -1])
        }
        colnames(profiles.df) <- c('direction', 'distance', 'value')
        write.table(profiles.df, file = "./data/profiles.csv", row.names=FALSE,
                    sep=",", fileEncoding = "UTF-8")      
        
#-------------------------------------------------------------------------------

#------------------------------------------------------------------------------#
# 2. CREATE RANDOM SAMPLE                       
#-------------------------------------------------------------------------------
        random.sample.spdf <- sampleRandom(gradient.raster, size = sample.size,
                                           cells = TRUE, xy = TRUE, sp = TRUE,
                                           na.rm = TRUE)
        writeOGR(obj = random.sample.spdf, dsn ="./data", layer = "sample_points",
                 driver = "ESRI Shapefile", overwrite_layer=TRUE)
#-------------------------------------------------------------------------------
        
#------------------------------------------------------------------------------#
# 3. DEFINE OPTIMUM NUMBER OF CLASSES AND BINS                                 #
#-------------------------------------------------------------------------------
        # calculate maximum possible number of classes
        valid.pixels <- cellStats(!is.na(gradient.raster), 'sum')
        n.pairs <- (valid.pixels*(valid.pixels-1))/2 # the number of pairs in the distance matrix
        n.classes <- round(1 + 3.322*log10(n.pairs)) # Sturges' rule
        dist.of.interest <- ifelse(is.null(dist.of.interest)==TRUE,
                                   max(profiles.df$distance)/2,
                                   dist.of.interest)
        
        if (is.null(lag.of.interest) == TRUE) {
                bin.size <- ifelse(test = ((dist.of.interest/n.classes)%%5) <= 2.5,
                                   yes = (dist.of.interest/n.classes) - ((dist.of.interest/n.classes)%%5),
                                   no = (dist.of.interest/n.classes) + (5 - ((dist.of.interest/n.classes)%%5)))
                
                # calculate bins' breaks
                bins.sum <- 0
                bins.breaks <- 0
                for (i in 1:n.classes){
                        bins.sum <- bins.sum + bin.size
                        bins.breaks <- append(bins.breaks, bins.sum)
                        if (bins.sum >= dist.of.interest) break
                }
        } else {
                bin.size <- lag.of.interest
                bins.sum <- 0
                bins.breaks <- 0
                for (i in 1:n.classes){
                        bins.sum <- bins.sum + bin.size
                        bins.breaks <- append(bins.breaks, bins.sum)
                        if (bins.sum >= dist.of.interest) break
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
                    