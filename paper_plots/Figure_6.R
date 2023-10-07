#Open required libraries
library(reshape2)
library(ggplot2)
library(ggpubr)

#Make sure that there is no open graphs
graphics.off()

size = 24

#Plots of time interval
results<-read.table("../results/TCP_t_int.dat", header=T)

p1<-ggplot(results, aes(x=interval, y=TCP, color=factor(T)))+geom_point()+scale_color_manual(values=c("blue", "orange", "red"), labels=c("39 °C", "41 °C", "43 °C"))+xlab("Mean Time Interval (min)")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.title=element_blank(),legend.box.background=element_rect(linewidth = 1.5),legend.position=c(0.12,0.14), legend.direction="vertical",legend.margin=margin(1,1,1,0, unit='mm'),legend.spacing.x=unit(0, "mm"),legend.spacing.y=unit(0, "mm"),strip.text.x = element_text(size = size))+facet_grid(cols=vars("mu==0.027~h^-1"),labeller = label_parsed)+scale_x_continuous(limits=c(10,240), breaks=c(10,50,100,150,200,240))+ylim(.55,1.)
p2<-ggplot(results, aes(x=interval, y=TCP_high, color=factor(T)))+geom_point()+scale_color_manual(values=c("blue", "orange", "red"), labels=c("39 °C", "41 °C", "43 °C"))+xlab("Mean Time Interval (min)")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.title=element_blank(),legend.box.background=element_rect(linewidth = 1.5),legend.position=c(0.125,0.14), legend.direction="vertical",legend.margin=margin(1,1,1,0, unit='mm'),legend.spacing.x=unit(0, "mm"),legend.spacing.y=unit(0, "mm"),strip.text.x = element_text(size = size))+facet_grid(cols=vars("mu==0.5~h^-1"),rows=vars("Without~direct~killing"),labeller = label_parsed)+scale_x_continuous(limits=c(10,240),breaks=c(10,50,100,150,200,240))+ylim(.55,1.)
p3<-ggplot(results, aes(x=interval, y=TCP_kill, color=factor(T)), show.legend = FALSE)+geom_point()+scale_color_manual(values=c("blue", "orange", "red"), labels=c("39 °C", "41 °C", "43 °C"))+xlab("Mean Time Interval (min)")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.title=element_blank(),legend.box.background=element_rect(linewidth = 1.5),legend.position=c(0.12,0.13), legend.direction="vertical",legend.margin=margin(1,1,1,0, unit='mm'),legend.spacing.x=unit(0, "mm"),legend.spacing.y=unit(0, "mm"),strip.text.x = element_text(size = size))+facet_grid(labeller = label_parsed)+scale_x_continuous(limits=c(10,240),breaks=c(10,50,100,150,200,240))+ylim(.55,1.)
p4<-ggplot(results, aes(x=interval, y=TCP_high_kill, color=factor(T)), show.legend = FALSE)+geom_point()+scale_color_manual(values=c("blue", "orange", "red"), labels=c("39 °C", "41 °C", "43 °C"))+xlab("Mean Time Interval (min)")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.title=element_blank(),legend.box.background=element_rect(linewidth = 1.5),legend.position=c(0.125,0.13), legend.direction="vertical",legend.margin=margin(1,1,1,0, unit='mm'),legend.spacing.x=unit(0, "mm"),legend.spacing.y=unit(0, "mm"),strip.text.x = element_text(size = size))+facet_grid(rows=vars("With~direct~killing"),labeller = label_parsed)+scale_x_continuous(limits=c(10,240),breaks=c(10,50,100,150,200,240))+ylim(.55,1.)


pdf("Figure_6.pdf",width=16,height=12)
ggarrange(p1+ xlab(""), p2 + xlab(""), p3, p4, font.label = list(size = size),ncol = 2, nrow = 2)
graphics.off()