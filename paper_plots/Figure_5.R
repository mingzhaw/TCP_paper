#Open required libraries
library(reshape2)
library(ggplot2)
library(ggpubr)

#Make sure that there is no open graphs
graphics.off()

size = 24

#Plots of temperature
results<-read.table("../results/TCP_T.dat", header=T)

p1<-ggplot(results, aes(x=Temp_min, y=TCP))+geom_point(color="blue")+xlab("Temperature / °C")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank(),strip.text.x = element_text(size = size))+facet_grid(cols=vars("Minimum~Temperature"),labeller = label_parsed)+xlim(39,43)+ylim(.7,1.)
p2<-ggplot(results, aes(x=Temp, y=TCP))+geom_point(color="blue")+xlab("Temperature / °C")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank(),strip.text.x = element_text(size = size))+facet_grid(cols=vars("Mean~Temperature"),labeller = label_parsed)+xlim(39,43)+ylim(.7,1.)
p3<-ggplot(results, aes(x=Temp_max, y=TCP))+geom_point(color="blue")+xlab("Temperature / °C")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank(),strip.text.x = element_text(size = size))+facet_grid(cols=vars("Maximum~Temperature"), rows=vars("Without~direct~killing"),labeller = label_parsed)+xlim(39,43)+ylim(.7,1.)

p4<-ggplot(results, aes(x=Temp_min, y=TCP_kill))+geom_point(color="blue")+xlab("Temperature (°C)")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank(),strip.text.x = element_text(size = size))+facet_grid(labeller = label_parsed)+xlim(39,43)+ylim(.7,1.)
p5<-ggplot(results, aes(x=Temp, y=TCP_kill))+geom_point(color="blue")+xlab("Temperature (°C)")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank(),strip.text.x = element_text(size = size))+facet_grid(labeller = label_parsed)+xlim(39,43)+ylim(.7,1.)
p6<-ggplot(results, aes(x=Temp_max, y=TCP_kill))+geom_point(color="blue")+xlab("Temperature (°C)")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank(),strip.text.x = element_text(size = size))+facet_grid(rows=vars("With~direct~killing"),labeller = label_parsed)+xlim(39,43)+ylim(.7,1.)


pdf("Figure_5.pdf",width=24,height=12)
ggarrange(p1+ xlab(NULL), p2 + xlab(NULL), p3+ xlab(NULL), p4, p5, p6, font.label = list(size = size),ncol = 3, nrow = 2)
graphics.off()