#-------------------------------------------------------------------------------
# UTF-8
# 
# This function performs preliminary data analysis and calculates the parameters
# necessary for the assesmant of univariate structure functions, namely
# semivariance and covarance.
#
# Arguments:
#       - mandatory:
#                     - raster.data: the raster to be analyzed
#                     - alpha: significance level alpha (e.g. 0.05)
#                       corresponding to a certain p-value (e.g. 0.95) or
#                       confidence level (e.g. 95%). Should take values of 0.05
#                       (default), 0.01, or 0.001.          
#       - optional:
#                     - r.matrix: radius of circular matrix (in the raster
#                       projection measurment units, usually meters) to calculate
#                       Moran's I index of spatial autocorrelation
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
                                     alpha = 0.05, r.matrix = NULL,
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
# 1. INTRODUCE HELPING FUNCTIONS
#-------------------------------------------------------------------------------
        # the function to remove NA's from SpatialDataFrame Object
        # from here http://tinyurl.com/nas-omit
        # x - SpatialDataFrame Object
        # margin - remove rows (1) or columns (2)
        sp.na.omit <- function(x, margin=1) {
                if (!inherits(x, "SpatialPointsDataFrame") & !inherits(x, "SpatialPolygonsDataFrame")) 
                        stop("MUST BE sp SpatialPointsDataFrame OR SpatialPolygonsDataFrame CLASS OBJECT") 
                na.index <- unique(as.data.frame(which(is.na(x@data),arr.ind=TRUE))[,margin])
                if(margin == 1) {  
                        cat("DELETING ROWS: ", na.index, "\n") 
                        return( x[-na.index,]  ) 
                }
                if(margin == 2) {  
                        cat("DELETING COLUMNS: ", na.index, "\n") 
                        return( x[,-na.index]  ) 
                }
        }
        
        # function to make a circular weights matrix of given radius and resolution
        # from here http://tinyurl.com/circ-matrix
        # Important: radius must be an even multiple of resolution
        make.circ.matrix<-function(radius, res){
                circ.matrix<-matrix(0, nrow=1+(2*radius/res), ncol=1+(2*radius/res))
                dimnames(circ.matrix)[[1]]<-seq(-radius, radius, by=res)
                dimnames(circ.matrix)[[2]]<-seq(-radius, radius, by=res)
                sweeper<-function(mat){
                        for(row in 1:nrow(mat)){
                                for(col in 1:ncol(mat)){
                                        dist<-sqrt((as.numeric(dimnames(mat)[[1]])[row])^2 +
                                                           (as.numeric(dimnames(mat)[[1]])[col])^2)
                                        if(dist<=radius) {mat[row, col]<-1}
                                }
                        }
                        return(mat)
                }
                out<-sweeper(circ.matrix)
                center <- ((radius/res)+1)
                out[center, center] <- 0
                return(out)
        }
        
        # sample size function
        # so far the 3 most common alpha values are possible
        # 0.05 (default), 0.01, and 0.001
        sample.size = function(raster, alpha=0.05, matrix=FALSE) {
                if (alpha==0.05){
                        z.value <- 1.96
                } else if (alpha==0.01){
                        z.value <- 2.58
                } else if (alpha==0.001){
                        z.value <- 3.09
                }
                sd <- cellStats(raster, stat = 'sd', na.rm=TRUE)
                mean <- cellStats(raster, stat = 'mean', na.rm=TRUE)
                pop <- cellStats(!is.na(gradient.raster), 'sum')
                if (matrix == FALSE){
                        ss1 <- (z.value^2*sd^2)/(mean*alpha)^2
                        ss2 <- ss1/(1+(ss1/pop))
                        sample <- list((round(ss2, 0)))
                        names(sample) <- 'sample size'
                } else if (matrix == TRUE) {
                        ss1 <- (z.value^2*sd^2)/(mean*alpha)^2
                        ss2 <- ss1/(1+(ss1/pop))
                        moran.cor <- Moran(raster, w = weight.matrix)
                        moran.cor <- ifelse(moran.cor<0, 0, moran.cor)
                        adj.ss <- ss2*((1-moran.cor)/(1+moran.cor))
                        ss <- (round(ss2, 0))
                        adj.ss <- round(adj.ss, 0)
                        sample <- list(adj.ss, ss)
                        names(sample) <- c('adjust_sample', 'full_sample') 
                }
                return(sample)
                
        }
#-------------------------------------------------------------------------------

        
#------------------------------------------------------------------------------#
# 2. IMPORT RASTER, BUILD PROFILES AND SAMPLE DATA ALONG THEM
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
                if (anyNA(trans.points[[i]]@data[,2]) == TRUE){
                        trans.points[[i]] <- sp.na.omit(trans.points[[i]])
                }
                # correct distance origin to 0
                trans.points[[i]]@data[[1]] <- trans.points[[i]]@data[[1]] - min(trans.points[[i]]@data[[1]])
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
# 3. DEFINE OPTIMAL SAMPLE SIZE AND CREATE RANDOM SAMPLE                       
#-------------------------------------------------------------------------------
        # create matrix to measure global autocorrelation
        matrix <- ifelse(is.null(r.matrix)==TRUE, FALSE, TRUE)
        if (matrix==TRUE){
                weight.distance <- r.matrix
                weight.radius <- (as.integer(weight.distance/cell.size))*cell.size
                weight.matrix <- make.circ.matrix(radius = weight.radius, res = cell.size)
        }
        
        sample <- sample.size(gradient.raster, alpha = alpha, matrix = matrix)
        
        # create random sample
        random.sample.spdf <- sampleRandom(gradient.raster, size = sample[[1]],
                                           cells = TRUE, xy = TRUE, sp = TRUE,
                                           na.rm = TRUE)
        writeOGR(obj = random.sample.spdf, dsn ="./data", layer = "sample_points",
                 driver = "ESRI Shapefile", overwrite_layer=TRUE)
#-------------------------------------------------------------------------------
        
#------------------------------------------------------------------------------#
# 4. DEFINE OPTIMUM NUMBER OF CLASSES AND BINS                                 #
#-------------------------------------------------------------------------------
        # calculate number of classes from the sample size
        n.pairs <- (sample[[1]]*(sample[[1]]-1))/2 # the number of pairs in the distance matrix
        n.classes <- as.integer(1 + 3.322*log10(n.pairs)) # Sturges' rule
        dist.of.interest <- ifelse(is.null(dist.of.interest)==TRUE,
                                   max(profiles.df$distance)/2,
                                   dist.of.interest)
        
        if (is.null(lag.of.interest) == TRUE) {
                bin.size <- ifelse(test = ((dist.of.interest/n.classes)%%10) <= 5,
                                   yes = (dist.of.interest/n.classes) - ((dist.of.interest/n.classes)%%10),
                                   no = (dist.of.interest/n.classes) + (10 - ((dist.of.interest/n.classes)%%10)))
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
        saveRDS(sample, "./data/sample.rds")
        saveRDS(sample.geodata, "./data/sample_geodata.rds")
        saveRDS(bins.centers, "./data/bins_centers.rds")
        saveRDS(bins.breaks, "./data/bins_breaks.rds")
        saveRDS(bins.sum, "./data/bins_sum.rds")
        saveRDS(n.classes, "./data/n_classes.rds")
        saveRDS(bin.size, "./data/bin_size.rds")
        
        print("Processing is completed.")        
}
                    