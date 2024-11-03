#!/usr/bin/env Rscript
library(ggplot2)
library(cowplot)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)


plot_theme <- theme(panel.background = element_rect(fill = NaN, colour = NaN),
                    plot.background = element_rect(fill=NaN, colour = NaN),
                    axis.line = element_line(),
                    strip.background = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(), 
                    strip.text = element_text(face = "bold", size = 11, angle = 0, hjust = 0.5),
                    strip.background.x = element_blank(),
                    text = element_text(size = 12, color = "black"),
                    axis.text.x = element_text(size = 12, angle = 0, color = "black"),
                    axis.title = element_text(size=12),
                    axis.text.y = element_text(size=12, color = "black"),
                    plot.margin=unit(c(0,0,0,0), "cm"),
                    legend.position = "None",
                    legend.key = element_blank(),
                    legend.text = element_text(size=12),
                    plot.title = element_text(hjust = 0.5))

table <- fread(args[1])
table[,read_status:="NA"]
table[bait_count!=0, read_status:="PASS"]
table[bait_count==0, read_status:="FAIL"]
summary_table_density_top <- table[,.(count=length(bait_count)), by=c('read_status', 'read_length')]
summary_table_density_right <- table[,.(count=length(bait_count)), by=c('read_status', 'bait_alignment')]

read2align_length <- ggplot()+plot_theme+
  geom_point(data=table, aes(read_length, bait_alignment, col=read_status))+
  xlab("Simulated read length")+
  ylab("Alignment lenght")+
  scale_color_manual(values = c('red', 'black'))

density <- ggplot()+plot_theme+
  geom_point(data=summary_table_density_top, aes(read_length, count, col=read_status))+
  theme(legend.position = 'top',
        axis.text.x = element_blank())+
  scale_color_manual(name='Read Status', values = c('red', 'black'))+
  xlab("")

density2 <- ggplot()+plot_theme+
  geom_point(data=summary_table_density_right, aes(count, bait_alignment, col=read_status))+
  scale_color_manual(values = c('red', 'black'))+
  theme(axis.text.y = element_blank())+
  ylab("")+
  annotate(x = 10, y = 30, geom = 'text', 
         label=paste(table[,1-(sum(read_status=='FAIL')/nrow(table))], ' of the simulated reads match a bait with:\n\t >91% identity \n\t >50 nt. alignment', sep=''), 
         hjust = 0, vjust = 0.5, size=4)

histogram <- ggplot()+plot_theme+
  geom_histogram(data=table, aes(bait_count), col='black', binwidth = 1)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(plot.margin=unit(c(1,0,0,0), "cm"))+
  xlab("Number of baits")  


suppressWarnings({
final_plot <- plot_grid(density, histogram, read2align_length, density2, ncol = 2, align = 'v', axis = 'lr')
})

suppressWarnings({
ggsave("validation_plot.pdf", plot = final_plot, width = 8, height = 5, device = "pdf")
})
