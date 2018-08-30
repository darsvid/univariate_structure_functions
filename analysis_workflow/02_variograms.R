#------------------------------------------------------------------------------#
# UTF-8
# 
# CALCULATE VARIOGRAMS                                                      
#------------------------------------------------------------------------------#

library(geoR) # variograms

# load all necessary data from previous steps
sample.geodata <- readRDS("./data/sample_geodata.rds")
bins.centers   <- readRDS("./data/bins_centers.rds")
bins.breaks    <- readRDS("./data/bins_breaks.rds")
bins.sum       <- readRDS("./data/bins_sum.rds")
n.classes      <- length(bins.centers)



# create an empty dataframe to store the results
variograms.df <- data.frame(direction = numeric(0), bin.center = numeric(0),
                            v = numeric(0), n = numeric(0), sd = numeric(0),
                            bin.limit = numeric(0))
col.names     <- names(variograms.df)


directions.names <- c(0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5, 180, 360)
directions.deg   <- directions.names[-10]
sample.variograms <- list()

# calculate directional variograms
for (i in 1:length(directions.deg)){
    sample.variograms[[i]] <- variog(sample.geodata, sample.geodata$coords,
                                     sample.geodata$data, uvec = bins.centers,
                                     breaks = bins.breaks, max.dist = bins.sum,
                                     pairs.min = 50, direction = directions.deg[i],
                                     tolerance = 11.25, unit.angle = "degree")
}


# calculate omnidirectional variogram
sample.variograms[[10]]   <- variog(sample.geodata, sample.geodata$coords,
                                    sample.geodata$data, uvec = bins.centers,
                                    breaks = bins.breaks, max.dist = bins.sum,
                                    pairs.min = 50, direction = "omnidirectional",
                                    tolerance = 11.25, unit.angle = "degree")

# write the results into the dataframe
for (i in 1:length(sample.variograms)) {
        var.length <- length(sample.variograms[[i]][[1]])
        temp.df    <- data.frame(id = 1:n.classes)
        if (var.length < n.classes){
            dif     <- n.classes - var.length
            temp.df <- cbind(temp.df,
                             rep(directions.names[i], n.classes),
                             c(sample.variograms[[i]][[1]], rep(NA, dif)),
                             c(sample.variograms[[i]][[2]], rep(NA, dif)),
                             c(sample.variograms[[i]][[3]], rep(NA, dif)),
                             c(sample.variograms[[i]][[4]], rep(NA, dif)),
                             c(sample.variograms[[i]][[5]][2:(var.length + 1)],
                               rep(NA,dif)))
        } else {
            temp.df <- cbind(temp.df,
                             rep(directions.names[i], n.classes),
                             sample.variograms[[i]][[1]],
                             sample.variograms[[i]][[2]],
                             sample.variograms[[i]][[3]],
                             sample.variograms[[i]][[4]],
                             sample.variograms[[i]][[5]][2:length(bins.breaks)])
        }
        variograms.df <- rbind(variograms.df, temp.df[, -1])
}
colnames(variograms.df) <- col.names

#standardize the variograms by the sample variance
variance <- var(sample.geodata$data)
v.stand  <- variograms.df$v/variance
variograms.df <- cbind(variograms.df, v.stand)

# export output to the csv-file
write.table(variograms.df, file = "./data/variograms.csv", row.names = F,
            sep = ",", fileEncoding = "UTF-8")