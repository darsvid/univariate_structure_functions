library(corrplot)

# read the correlations data
corr.df <- read.csv("correlations_100.csv", stringsAsFactors = F)

# create the matrix of correlations
vars <- unique(unlist(corr.df[1:2]))
corr.df[1:2] <- lapply(corr.df[1:2], factor, levels = vars)
mt.cor <- matrix(0, nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
mt.cor[cbind(match(corr.df$variable1, rownames(mt.cor)), 
         match(corr.df$variable2, colnames(mt.cor)))] <- corr.df$corr
vars.names <- c("easterness", "EVI", "local dif. from mean", "local mean", "local std. deviation", "northerness", "slope")
rownames(mt.cor) <- vars.names
colnames(mt.cor) <- vars.names
corr.matrix <- mt.cor + t(mt.cor)

# create the matrix of p-values
mt.p <- matrix(0, nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
mt.p[cbind(match(corr.df$variable1, rownames(mt.p)), 
         match(corr.df$variable2, colnames(mt.p)))] <- corr.df$p.value
rownames(mt.p) <- vars.names
colnames(mt.p) <- vars.names
pval.matrix <- mt.p + t(mt.p)

# visualize correlations
col <- colorRampPalette(c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"))
png("correlations_plot.png", res = 300, height = 10, width = 10, units = "cm", pointsize = 6)
corr.plot <- corrplot(corr.matrix, method="color", col=col(10),  addgrid.col = "grey",
                      type="upper", addCoef.col = "black", number.digits = 3, # Add coefficient of correlation
                      tl.col="black", tl.srt=45, #Text label color and rotation
                    # Combine with significance
                      p.mat = pval.matrix, sig.level = 0.01, insig = "blank", 
                    # hide correlation coefficient on the principal diagonal
                      diag = FALSE)
dev.off()