#Open required libraries
library(reshape2)
library(ggplot2)
library(ggpubr)

#Make sure that there is no open graphs
graphics.off()

size = 24

#Plots of time-gap effect
results<-read.table("../results/t_int.dat", header=T)

p1<-ggplot(results, aes(x=t_int, y=survival, color=factor(T)))+geom_line(linewidth = 1.5)+scale_color_manual(values=c("blue", "orange", "red"), labels=c("39 °C", "41 °C", "43 °C"))+xlab("Time interval (h)")+ylab("Survival Fraction")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.title=element_blank(),legend.box.background=element_rect(linewidth = 1.5),legend.position=c(0.1,0.12), legend.direction="vertical",legend.margin=margin(1,1,1,0, unit='mm'),legend.spacing.x=unit(0, "mm"),legend.spacing.y=unit(0, "mm"),strip.text.x = element_text(size = size))+facet_grid(cols=vars("mu==0.027~h^-1"),labeller = label_parsed)+xlim(-4,4)+ylim(-.02,0.45)+geom_hline(yintercept=0.42390030409889934, linewidth=1.5)+annotate("text", x=0, y=0.45, label="Only RT", size=8)
p2<-ggplot(results, aes(x=t_int, y=survival_high, color=factor(T)))+geom_line(linewidth = 1.5)+scale_color_manual(values=c("blue", "orange", "red"), labels=c("39 °C", "41 °C", "43 °C"))+xlab("Time interval (h)")+ylab("Survival Fraction")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.title=element_blank(),legend.box.background=element_rect(linewidth = 1.5),legend.position=c(0.1,0.12), legend.direction="vertical",legend.margin=margin(1,1,1,0, unit='mm'),legend.spacing.x=unit(0, "mm"),legend.spacing.y=unit(0, "mm"),strip.text.x = element_text(size = size))+facet_grid(cols=vars("mu==0.5~h^-1"),labeller = label_parsed)+xlim(-4,4)+ylim(-.02,0.45)+geom_hline(yintercept=0.42390030409889934, linewidth=1.5)+annotate("text", x=0, y=0.45, label="Only RT", size=8)


pdf("Figure_1.pdf",width=16,height=6)
ggarrange(p1, p2, font.label = list(size = size),ncol = 2, nrow = 1)
graphics.off()