#------------------------------------------------------------------------------#
# UTF-8
# 
# CALCULATE   CORRELOGRAMS                                                  
#------------------------------------------------------------------------------#

library(EcoGenetics) # correlograms

# load all necessary data from previous steps
sample.geodata <- readRDS("./data/sample_geodata.rds")
bins.centers   <- readRDS("./data/bins_centers.rds")
bins.breaks    <- readRDS("./data/bins_breaks.rds")
bins.sum       <- readRDS("./data/bins_sum.rds")
bin.size       <- readRDS("./data/bin_size.rds")
n.classes      <- length(bins.centers)

# calculate directional correlograms
directions.names <- c(0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5, 180, 360)
directions.deg   <- directions.names[-10]

# this operation tends to fail sometimes, so we save intermediary results for
# each step and force the cycle to start from where it was interrupted
if (file.exists("./data/correlograms.rds")){
  sample.correlograms <- readRDS("./data/correlograms.rds")
} else {
  sample.correlograms <- list()
}
n <- length(sample.correlograms) + 1

for (i in n:length(directions.deg)){
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
                                                 smax = bins.sum, nsim = 5000,
                                                 angle = directions.deg[i],
                                                 alpha = 0.01)
        saveRDS(sample.correlograms, "./data/correlograms.rds")
        
        
}
# calculate omnidirectional correlogram
sample.correlograms[[10]]        <- eco.correlog(Z = sample.geodata$data,
                                                 XY = sample.geodata$coords,
                                                 method = 'I', int = bin.size,
                                                 smax = bins.sum, nsim = 5000,
                                                 alpha = 0.01)
saveRDS(sample.correlograms, "./data/correlograms.rds")

# create an empty dataframe to store the results
correlograms.df <- data.frame(direction = numeric(0), d.mean = numeric(0),
                              breaks = numeric(0), obs = numeric(0),
                              p.val = numeric(0), size = numeric(0))

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
write.table(correlograms.df, file = "./data/correlograms_moran.csv", row.names = F,
            sep = ",", fileEncoding = "UTF-8")