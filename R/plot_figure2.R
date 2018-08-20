library(ggplot2)
library(reshape2)
library(ggpubr)

# load data
plot.dir <- 'figures'
out.dir <- 'output'
rl <- readRDS(file.path(out.dir, 'rl100.rds'))
rl.melt <- melt(rl)
rl.melt$value <- rl.melt$value*3.28

# construct levee heights
# signpost years
year.WAIS <- 2045
year.ns <- 2060

year.h1 <- 2030

# starting and ending heights
levee.start <- 8
levee.end <- 9.5
levee.WAIS <- 9
levee.ns <- 9.5
levee.h1 <- 8.5

levee <- data.frame(height=c(rep(levee.start, (year.h1-2020+5)/5), rep(levee.h1, (year.WAIS-year.h1)/5), rep(levee.WAIS, (year.ns-year.WAIS)/5), rep(levee.ns, (2070-year.ns)/5)), year=seq(2020, 2070, 5))

# plot return levels
p <- ggplot(rl.melt) + geom_line(aes(x=as.numeric(L1), y=value, color=as.factor(L2)), linetype=2) + geom_line(data=levee, aes(x=year, y=height), color='black', linetype=1) + theme_bw() + scale_x_continuous('Year', limits=c(2020, 2065)) + scale_y_continuous('Water/Levee Height Relative to 2015 MSL (ft)', limits=c(5, 10), breaks=seq(0, 10, 1)) + guides(colour=guide_legend(ncol=2)) + theme(legend.key.size=unit(1.5, 'lines'), legend.position='top', legend.direction='horizontal') + geom_vline(aes(xintercept=year.WAIS), linetype=3, color='black') + geom_vline(aes(xintercept=year.ns), linetype=3, color='black') + geom_text(aes(x=year.WAIS+0.5, y=6, angle=-90), size=2, label='WAIS Collapse Observed') + geom_text(aes(x=year.ns+0.5, y=6, angle=-90), size=2, label='Surge Nonstationary Detected') + theme_pubr(base_size=7) + scale_color_brewer('Levee Design Assumptions', palette='Dark2', labels=c('linear extrapolation\nwith historical surge record',
                                                   'no accelerated ice sheet melting\nand stationary surge',
                                                   'possibly accelerated ice sheet melting\nand stationary surge',
                                                   'possibly accelerated ice sheet melting\nand possible changes to surge'))
                                                   
pdf(file.path(plot.dir, 'Norfolk_adapt_plan.pdf'), width=4, height=4)
p
dev.off()

bitmap(file=file.path(plot.dir, 'Norfolk_adapt_plan.png'), type='png256', width=4, height=4, units='in', res=600)
p
dev.off()