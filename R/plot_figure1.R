# script to generate Figure 1 from Srikrishnan et al (2018)

# set user controlled parameters
year <- 2070

# load libraries
library(plyr)
library(extRemes)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(scales)

# find the index of the target year
year.ind <- year - 1850 + 1

# set directories
plot.dir <- 'figures'
data.dir <- 'data'
out.dir <- 'output'

filename.SOWs <- file.path(out.dir, 'SL-SOWs.rds')

SOWs <- readRDS(filename.SOWs)

# compute 100-year return levels
rl.100 <- lapply(SOWs, quantile, prob=c(.99))

# compute quantiles
n.quantile <- 1e5
q <- lapply(SOWs, quantile, prob=seq(0, 1, by=1/n.quantile), names=FALSE)
q.melt <- melt(q)
q.melt <- ddply(q.melt, c("L1"), cbind, q=1/(1-seq(0, 1, by=1/n.quantile)))
#q.melt <- ddply(q.melt, c("L1"), rbind, c(0, q=1))
#q.melt[q.melt$value == 0,'L1'] <- seq(1,4)
#q.melt[q.melt$value == 0,'q'] <- 1

rl.melt <- q.melt
rl.melt[,'q'] <- 1/rl.melt[,'q']

# melt for plotting
SOW.melt <- melt(SOWs)

rl.lin <- t(ldply(rl.100, function(y) unlist(lapply(q, function(x) 1/(1-which.min(abs(x-y)*3.28)/n.quantile)))))
print(rl.lin)
rl.lin2 <- data.frame(model=1:4, rl=rl.lin[,1])
rl.melt <- melt(rl.lin2, id.var='model')

# plot pdfs
p1 <- ggplot(SOW.melt) + geom_density(aes(x=value*3.28, color=as.factor(L1)), show.legend=FALSE) +
  stat_density(aes(x=value*3.28, color=as.factor(L1)), geom='line', position='identity') +
  theme_pubr(base_size=7)  + guides(colour=guide_legend(ncol=2)) +
  scale_y_continuous('Probability density', expand=c(0.001, 0.001)) +
  scale_x_continuous('Water height anomaly relative to 2015 (ft)', limits=c(0,12), expand=c(0.01, 0.01), breaks=seq(0, 11, 1)) +
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(),
        legend.key.size=unit(1.5, 'lines'), legend.position='top', legend.direction='horizontal',
        axis.line=element_line(arrow=arrow(length=unit(0.05, 'inches'), type='closed')), plot.title=element_text(size=9),
        text=element_text(family="ArialMT")) +
  coord_flip() +
  scale_color_brewer('', palette='Dark2', labels=c('Linear extrapolation\nwith historical surge record',
                                                   'No accelerated ice sheet melting\nand stationary surge',
                                                   'Possibly accelerated ice sheet melting\nand stationary surge',
                                                   'Possibly accelerated ice sheet melting\nand possible changes to surge')) +
  geom_text(aes(x=11.5, y=0.1, label='A)', family='ArialMT', fontface='bold'), size=3.214)
#  ggtitle(expression(paste(bold('A)'), 'Water height distribution')))

p2 <- ggplot() + geom_line(data=q.melt, aes(y=value*3.28, x=q, color=as.factor(L1))) +
  theme_pubr(base_size=7) + guides(colour=guide_legend(ncol=2)) +
  geom_segment(aes(y=8.06, yend=8.06, x=min(rl.lin2[,'rl']), xend=100), linetype=3) +
  geom_point(data=rl.lin2, aes(y=8.06, x=rl, color=as.factor(model)), size=1.5) +
  geom_segment(data=rl.lin2[2:4,], aes(y=0, yend=8.06, x=rl, xend=rl), linetype=3) +
  scale_x_continuous('Return period (yr)', expand=c(0.01, 0.01), limits=c(1, 125), breaks=c(seq(0, 200, by=100), round(rl.lin2[2:4,'rl']))) +
  geom_vline(aes(xintercept=100), linetype=2) +
  scale_y_continuous('Water Level (ft)', limits=c(0, 12), expand=c(0.01, 0.01), breaks=seq(0, 11, 1)) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        legend.key.size=unit(1.5, 'lines'), legend.position='top', legend.direction='horizontal',
        axis.line=element_line(arrow=arrow(length=unit(0.05, 'inches'), type='closed')), plot.title=element_text(size=9),
        text=element_text(family="ArialMT")) +
  scale_color_brewer('Model Assumptions', palette='Dark2', labels=c('Linear extrapolation\nwith historical surge record',
                                                   'No accelerated ice sheet melting\nand stationary surge',
                                                   'Accelerated ice sheet melting\nand stationary surge',
                                                   'Accelerated ice sheet melting\nand changes to surge')) +
  geom_text(aes(y=11.5, x=5, label='B)', family='ArialMT', fontface='bold'), size=3.214)
#  ggtitle(expression(paste(bold('B)'), 'Water height return levels')))

# plot lengths of 100-year return
p3 <- ggplot(melt(rl.100)) + geom_segment(aes(y=0, x=L1/2, xend=L1/2, yend=value*3.28, color=as.factor(L1)), lineend='round') +
  theme_pubr(base_size=7) + guides(colour=guide_legend(ncol=2)) +
  geom_text(aes(y=value*3.28+0.75, x=L1/2, color=as.factor(L1), label=paste(round(value*3.28,1), 'ft', sep=' '), family='ArialMT'), size=2, angle=-90) +
  scale_x_continuous('', limits=c(0.25, 2.25)) +
  scale_y_continuous('Maximum Water Level Relative to 2015 Mean Sea Level (ft)', limits=c(0,12), expand=c(0.01, 0.01), breaks=seq(0, 11, 1)) +
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.line.x=element_blank(), axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        legend.key.size=unit(1.5, 'lines'), legend.position='top', legend.direction='horizontal',
        axis.line=element_line(arrow=arrow(length=unit(0.05, 'inches'), type='closed')),
        legend.text=element_text(size=8), plot.title=element_text(size=9),
        text=element_text(family="ArialMT")) +
  scale_color_brewer('Model assumptions', palette='Dark2', labels=c('Linear extrapolation\nwith historical surge record',
                                                   'No accelerated ice sheet melting\nand stationary surge',
                                                   'Possibly accelerated ice sheet melting\nand stationary surge',
                                                   'Possibly accelerated ice sheet melting\nand possible changes to surge')) +
  geom_text(aes(y=11.5, x=0.4, label='C)', family='ArialMT', fontface='bold'), size=3.214)
#  ggtitle(expression(paste(bold('C)'), '100-year return level')))


plots <- list(p1, p2, p3)
grobs <- list()
heights <- list()
g <- ggplotGrob(plots[[2]] + theme(legend.position="top"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]] + theme(legend.position='none'))
  heights[[i]] <- grobs[[i]]$heights[2:10]
}
maxheight <- do.call(grid::unit.pmax, heights)
for (i in 1:length(grobs)){
  grobs[[i]]$heights[2:10] <- as.list(maxheight)
}

p <- do.call("grid.arrange", list(grobs=grobs, nrow = 1, widths=c(1, 1.5, 0.5)))

pdf(file.path(plot.dir, 'Norfolk_prelim_plot.pdf'), width=6, height=3, fonts=c('sans', 'serif', 'ArialMT'))
grid.arrange(
  legend,
  p,
  nrow = 2,
  heights = grid::unit.c(lheight, unit(1, "npc") - lheight))
dev.off()

png(file=file.path(plot.dir, 'Norfolk_prelim_plot.png'), type='cairo-png', width=6, height=3, units='in', res=600, family='ArialMT')
grid.arrange(
  legend,
  p,
  nrow = 2,
  heights = grid::unit.c(lheight, unit(1, "npc") - lheight))
dev.off()

()

