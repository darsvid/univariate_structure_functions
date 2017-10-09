#------------------------------------------------------------------------------#
# UTF-8
# 
# CALCULATE   CORRELOGRAMS                                                  
#------------------------------------------------------------------------------#

library(EcoGenetics) # correlograms

# load all necessary data from previous steps
sample.geodata <- readRDS("./data/sample_geodata.rds")
bins.centers <- readRDS("./data/bins_centers.rds")
bins.breaks <- readRDS("./data/bins_breaks.rds")
bins.sum <- readRDS("./data/bins_sum.rds")
bin.size <- readRDS("./data/bin_size.rds")
n.classes <- length(bins.centers)

# calculate directional correlograms
directions.deg <- c(0, 45, 90, 135)
sample.correlograms <- list()
for (i in 1:length(directions.deg)){
        # change method='C' if you want to calculate Geary's C autocorrelation
        # increase the nsim parameter if you want to increase the precision and
        # reliability of p-values.
        # please consult this thread https://tinyurl.com/p-permutation to define
        # an appropriate number of permutations
        # Note: the higher nsim the longer the time of the calculation,
        # especially for Geary's C autocorrelation
        sample.correlograms[[i]] <- eco.correlog(Z = sample.geodata$data,
                                                 XY = sample.geodata$coords,
                                                 method = 'I', int = bin.size,
                                                 smax = bins.sum, nsim = 500,
                                                 angle = directions.deg[i],
                                                 alpha = 0.01)
}
# calculate omnidirectional correlogram
# change method='C' if you want to calculate Geary's C autocorrelation increase
# the nsim parameter if you want to increase the precision and reliability of
# p-values.
# please consult this thread https://tinyurl.com/p-permutation to define
# an appropriate number of permutations
# Note: the higher nsim the longer the time of the calculation, especially for
# Geary's C autocorrelation
sample.correlograms[[5]] <- eco.correlog(Z = sample.geodata$data,
                                         XY = sample.geodata$coords,
                                         method = 'I', int = bin.size,
                                         smax = bins.sum, nsim = 500,
                                         alpha = 0.01)


# create an empty dataframe to store the results
correlograms.df <- data.frame(direction = character(0), d.mean = numeric(0),
                              breaks = numeric(0), obs = numeric(0),
                              p.val = numeric(0), size = numeric(0))
directions.names <- c('0', '45', '90', '135', 'omnidir')
col.names <- names(correlograms.df)

# write the results into the dataframe
for (i in 1:length(sample.correlograms)) {
        temp.df <- data.frame(id = 1:n.classes)
        temp.df <- cbind(temp.df,
                         rep(directions.names[i], n.classes),
                         sample.correlograms[[i]]@OUT[[1]]$d.mean,
                         sample.correlograms[[i]]@BREAKS[2:length(bins.breaks)],
                         sample.correlograms[[i]]@OUT[[1]]$obs,
                         sample.correlograms[[i]]@OUT[[1]]$p.val,
                         sample.correlograms[[i]]@OUT[[1]]$size)
        correlograms.df <- rbind(correlograms.df, temp.df[, -1])
}
colnames(correlograms.df) <- col.names

# export output to the csv-file
write.table(correlograms.df, file = "./data/correlograms.csv", row.names=FALSE,
            sep=",", fileEncoding = "UTF-8")