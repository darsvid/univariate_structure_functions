#------------------------------------------------------------------------------#
# UTF-8
# 
# CALCULATE VARIOGRAMS                                                      
#------------------------------------------------------------------------------#

library(geoR) # variograms

# load all necessary data from previous steps
sample.geodata <- readRDS("./data/sample_geodata.rds")
bins.centers <- readRDS("./data/bins_centers.rds")
bins.breaks <- readRDS("./data/bins_breaks.rds")
bins.sum <- readRDS("./data/bins_sum.rds")
n.classes <- length(bins.centers)


# calculate variograms
sample.variograms <- variog4(sample.geodata, sample.geodata$coords,
                             sample.geodata$data, uvec=bins.centers,
                             breaks=bins.breaks, max.dist=bins.sum, pairs.min=8)

# create an empty dataframe to store the results
variograms.df <- data.frame(direction = character(0), bin.center = numeric(0),
                            v = numeric(0), n = numeric(0), sd = numeric(0),
                            bin.limit = numeric(0))
directions.names <- c('0', '45', '90', '135', 'omnidir')
col.names <- names(variograms.df)

# write the results into the dataframe
for (i in 1:length(sample.variograms)) {
        temp.df <- data.frame(id = 1:n.classes)
        temp.df <- cbind(temp.df,
                         rep(directions.names[i], n.classes),
                         sample.variograms[[i]][[1]],
                         sample.variograms[[i]][[2]],
                         sample.variograms[[i]][[3]],
                         sample.variograms[[i]][[4]],
                         sample.variograms[[i]][[5]][2:length(bins.breaks)])
        variograms.df <- rbind(variograms.df, temp.df[, -1])
}
colnames(variograms.df) <- col.names

# export output to the csv-file
write.table(variograms.df, file = "./data/variograms.csv", row.names=FALSE,
            sep=",", fileEncoding = "UTF-8")