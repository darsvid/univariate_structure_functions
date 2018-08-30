library(ggplot2)
library(latex2exp)

# SEMIVARIANCE PLOTS------------------------------------------------------------
# 
# prep the data ----------------------------------------------------------------

# read and reorganize the data
var.files <- list.files(pattern = "\\_variograms.csv$")
var.df <- data.frame(variable = character(0), direction = numeric(0),
                     bin.center = numeric(0), v.stand = numeric(0),
                     v.omni = numeric(0))
for (i in 1:length(var.files)){
    df       <- read.csv(var.files[i], stringsAsFactors = F)
    dir.df   <- subset(df, direction <= 180)
    omni.df  <- subset(df, direction == 360)
    v.omni   <- rep(omni.df$v.stand, 9)
    name     <- gsub(pattern = "(^.*)(\\_variograms.csv$)", "\\1", var.files[i])
    variable <- rep(name, nrow(dir.df))
    temp.df  <- cbind(variable, dir.df[, c(1, 2, 7)], v.omni)
    var.df   <- rbind(var.df, temp.df)
    
}

# redefine variables names and levels for proper naming
levels(var.df$variable) <- c("local dif. from mean", "easterness", "EVI",
                             "local mean", "northerness", "slope")
var.df$variable <- factor(var.df$variable,
                          levels = c("local mean", "slope", "local dif. from mean",
                                     "northerness", "easterness", "EVI"))
palette <- c('#440154', '#443a83', '#31688e', '#20908d',
             '#35b779', '#8fd744', '#fde725')

# semivariance - 6 vars---------------------------------------------------------

gv6 <- ggplot(data = var.df, aes(direction, bin.center/1000)) +
    stat_contour(aes(x = direction + 180, z = v.omni, colour = ..level..),
                 binwidth = 0.1, size = 0.6) +
    stat_contour(aes(z = v.stand, colour =..level..), binwidth=0.1, size = 0.6) +
    scale_colour_gradientn(colours = palette, name = TeX("$\\frac{\\hat{\\gamma}(h)}{s^2}$"), breaks = seq(0, 1.5, 0.2)) +
    scale_x_continuous(breaks = seq(45, 360, 45), limits = c(0, 360)) +
    scale_y_continuous(breaks = seq(0, 30,5)) +
    coord_polar(theta = "x", start = 1.5*pi, direction = -1) +
    labs( x = "direction", y = "distance, km") +
    theme(panel.background = element_rect(fill = "white", colour = "white",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.3, linetype = 'solid',
                                          colour = "grey80"),
          panel.grid.minor = element_line(size = 0.3, linetype = 'solid',
                                          colour = "grey80"),    
          plot.margin = rep(unit(1,'pt'), 4)) +
         #strip.text = element_text(size=16,face="bold"),
         #plot.title = element_text(size=16,face="bold"),
         #axis.text=element_text(size=14),
         #axis.title=element_text(size=14)) +
    facet_wrap(~variable, ncol = 2)
gv6
ggsave("semivariance_pplot.png", gv6, "png", width = 13.48, height = 18, units = "cm", dpi = 500)



# AUTOCORRELATION PLOTS---------------------------------------------------------
#
# prep the data-----------------------------------------------------------------

cor.files <- list.files(pattern = "\\_correlograms_moran.csv$")
cor.df <- data.frame(variable = character(0), direction = numeric(0),
                     d.mean = numeric(0), cor.dir = numeric(0),
                     p.dir = numeric(0), cor.omni = numeric(0),
                     p.omni = numeric(0))

for (i in 1:length(cor.files)){
    df       <- read.csv(cor.files[i], stringsAsFactors = F)
    dir.df   <- subset(df, direction <= 180)
    omni.df  <- subset(df, direction == 360)
    cor.omni <- rep(omni.df$obs, 9)
    p.omni   <- rep(omni.df$p.val, 9)
    name     <- gsub(pattern = "(^.*)(\\_correlograms_moran.csv$)", "\\1",
                     cor.files[i])
    variable <- rep(name, nrow(dir.df))
    temp.df  <- cbind(variable, dir.df[, c(1, 2, 4, 5)], cor.omni, p.omni)
    cor.df   <- rbind(cor.df, temp.df)
    
}

# redefine variables names and levels for proper naming
levels(cor.df$variable) <- c("local dif. from mean", "easterness", "EVI", "local mean",
                             "northerness", "slope")
cor.df$variable <- factor(cor.df$variable,
                          levels = c("local mean", "slope", "local dif. from mean",
                                     "northerness", "easterness", "EVI"))

# autocorrelation - 6 vars------------------------------------------------------
gc6 <- ggplot(data = cor.df, aes(direction, d.mean / 1000)) +
    stat_contour(aes(x = direction + 180, z = cor.omni, colour = ..level..),
                 binwidth = 0.1, size = 0.6) +
    stat_contour(aes(z = obs, colour =..level..), binwidth = 0.1, size = 0.6) +
    scale_colour_gradient2(high = '#d53e6b', mid = '#e3e3e3', low = '#4a6fe4',
                           name = TeX("$I_{Moran}$"), breaks = seq(-1, 1, 0.2)) +
    geom_point(data = subset(cor.df, p.val > 0.01),
               aes(x = direction, y = d.mean / 1000), size = 0.01, inherit.aes = F) +
    geom_point(data = subset(cor.df, p.omni > 0.01),
               aes(x = direction + 180, y = d.mean / 1000), size = 0.01, inherit.aes = F) +
    scale_x_continuous(breaks = seq(45, 360, 45), limits = c(0, 360)) +
    scale_y_continuous(breaks = seq(0, 30,5)) +
    coord_polar(theta = "x", start = 1.5*pi, direction = -1) +
    labs( y = "distance, km")+
    theme(panel.background = element_rect(fill = "white", colour = "white",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.3, linetype = 'solid',
                                          colour = "grey80"),
          panel.grid.minor = element_line(size = 0.3, linetype = 'solid',
                                          colour = "grey80"),
          plot.margin = rep(unit(1,'pt'), 4)) +
    facet_wrap(~variable, ncol = 2)
gc6
ggsave("autocorrelation_pplot.png", gc6, "png", width = 13.48, height = 18, units = "cm", dpi = 500)