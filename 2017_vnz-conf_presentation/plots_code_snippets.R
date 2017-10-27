library(ggplot2)
library(dplyr)

## PROFILES---------------------------------------------------------------------
profiles.df <- read.csv('./profiles.csv')

# profiles plot
profiles.plot <- ggplot(profiles.df, aes(x = distance, y = value)) +
        #geom_line(aes(colour = value), size = 0.8) +
        geom_line(size = 0.8) +
        # mean ratio = 10
        # std dev ratio = 100
        # dev from mean ratio = 500
        # evi ratio = 3000
        coord_fixed(ratio = 3000) +
        labs(x = 'distance, m', y = 'gradient, units') +
        #remove option limits if you do not want to cut off profiles
        scale_x_continuous(expand = c(0, 0), limits = c(0, 14000),
                           breaks=seq(0, 14000, 2000)) +
        #scale_colour_gradientn(colours = terrain.colors(12)) +
        facet_wrap(~direction, ncol = 1) +
        theme(#legend.position = c(0.96, 0.915),
              plot.margin = rep(unit(0,'null'),4),
              strip.text = element_text(size=16,face="bold"),
              plot.title = element_text(size=16,face="bold"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=16,face="bold"))
profiles.plot
ggsave('profiles_plot.png', dpi = 300)


## SEMIVARIANCE-----------------------------------------------------------------
semivariance.df <- read.csv('./variograms.csv',
                            stringsAsFactors = FALSE)
dir.semivar.df <- subset(semivariance.df,
                         direction %in% c('0', '45', '90', '135')) %>%
        mutate(direction = as.numeric(direction))
omni.semivar.df <- subset(semivariance.df, direction == 'omnidir')       

variograms.plot <- ggplot(dir.semivar.df, aes(bin.center, v), color='black',
                             size = 0.5) +
        annotate(geom='line', x=omni.semivar.df$bin.center, y=omni.semivar.df$v,
                 color = 'grey40', size = 0.8) +
        geom_line(size = 0.8) +
        # mean ratio = 1
        # std dev ratio = 300
        # dev from mean ratio = 500000
        # evi ratio = 75000
        coord_fixed(ratio = 75000) +
        labs(x = 'distance, m', y = 'semivariance') +
        scale_x_continuous(expand=c(0,0), limits = c(0, 7000),
                           breaks=seq(0, 7000, 1000)) +
        facet_wrap(~direction, ncol = 1) +
        theme(plot.margin = rep(unit(0,'null'),4),
              strip.text = element_text(size=16,face="bold"),
              plot.title = element_text(size=16,face="bold"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=16,face="bold"))
variograms.plot
ggsave('variograms_plot.png', dpi = 300)

## AUTOCORRELATION I------------------------------------------------------------
autocorr.df <- read.csv('./correlograms_I.csv',
                            stringsAsFactors = FALSE)
dir.autocorr.df <- subset(autocorr.df,
                         direction %in% c('0', '45', '90', '135')) %>%
        mutate(direction = as.integer(direction))
omni.autocorr.df <- subset(autocorr.df, direction == 'omnidir')

correlograms.moran.plot <- ggplot(dir.autocorr.df, aes(d.mean, obs), color='black') +
        annotate(geom = 'segment', x = 0, xend = 7000, y = 0, yend = 0,
                 colour = 'blue', alpha = 0.5) +
        annotate(geom='line', x=omni.autocorr.df$d.mean, y=omni.autocorr.df$obs,
                 color = 'grey50', size = 0.8) +
        geom_line(size = 0.3) +
        geom_point(aes(fill = ifelse(p.val<=0.01, "<= 0.01", "> 0.01")),
                   shape = 21, color = 'black', size = 2) +
        scale_fill_manual(values = c("black", "white"), name = 'p-value') +
        # mean ratio = 1000
        # std dev ratio = 1000
        coord_fixed(ratio = 1000) +
        labs(x = 'distance, m', y = "Moran's I") +
        scale_x_continuous(expand=c(0,0), limits = c(0, 7000),
                           breaks=seq(0, 7000, 1000)) +
        scale_y_continuous(limits = c(-1, 1)) +
        facet_wrap(~direction, ncol = 1) +
        theme(legend.position = c(0.90, 0.95),
              legend.title=element_text(size=14),
              legend.text=element_text(size=14),
              plot.margin = rep(unit(0,'null'),4),
              strip.text = element_text(size=16,face="bold"),
              plot.title = element_text(size=16,face="bold"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=16,face="bold"))
correlograms.moran.plot
ggsave('correlograms_moran_plot.png', dpi = 300)

## AUTOCORRELATION C------------------------------------------------------------
autocorrc.df <- read.csv('./correlograms_C.csv',
                        stringsAsFactors = FALSE)

dir.autocorrc.df <- subset(autocorrc.df,
                          direction %in% c('0', '45', '90', '135')) %>%
        mutate(direction = as.integer(direction))
omni.autocorrc.df <- subset(autocorrc.df, direction == 'omnidir')

correlograms.geary.plot <- ggplot(dir.autocorrc.df,
                                  aes(d.mean, obs),
                                  color='black') +
        annotate(geom = 'segment', x = 0, xend = 7000, y = 1, yend = 1,
                 colour = 'blue', alpha = 0.5) +
        annotate(geom='line', x=omni.autocorrc.df$d.mean, y=omni.autocorrc.df$obs,
                 color = 'grey50', size = 0.8) +
        geom_line(size = 0.3) +
        geom_point(aes(fill = ifelse(p.val<=0.01, "<= 0.01", "> 0.01")),
                   shape = 21, color = 'black', size = 2) +
        scale_fill_manual(values = c("black", "white"), name = 'p-value') +
        coord_fixed(ratio = 1000) +
        labs(x = 'distance, m', y = "Geary's c") +
        scale_x_continuous(expand=c(0,0), limits = c(0, 7000),
                           breaks=seq(0, 7000, 1000)) +
        scale_y_continuous(limits = c(0, 1.5)) +
        facet_wrap(~direction, ncol = 1) +
        theme(legend.position = c(0.91, 0.85),
              legend.title=element_text(size=14),
              legend.text=element_text(size=14),
              plot.margin = rep(unit(0,'null'),4),
              strip.text = element_text(size=16,face="bold"),
              plot.title = element_text(size=16,face="bold"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=16,face="bold"))
correlograms.geary.plot
ggsave('correlograms_geary_plot.png', dpi = 300)

