#Open required libraries
library(reshape2)
library(ggplot2)
library(ggpubr)

#Make sure that there is no open graphs
graphics.off()

#Plots of dose-control curves
size = 24
results<-read.table("../results/dose_control_het.dat", header=T)

p1<-ggplot(results, aes(x=D, y=TCP))+geom_line(color="blue",linewidth=1.5)+xlab("Dose (Gy)")+ylab("Mean TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),strip.text.x = element_text(size = size))+facet_grid(cols=vars("No~heterogeneity"),labeller = label_parsed)+ylim(.0,1.0)+scale_x_continuous(limits=c(0,120),breaks=c(0,20,40,60,80,100,120))
p2<-ggplot(results, aes(x=D))+geom_line(aes(y=TCP_het_alpha, color="Radiosensitivy"),linewidth=1.5)+geom_line(aes(y=TCP_het_vol, color="Volume"),linewidth=1.5)+scale_color_manual(values=c("Radiosensitivy" = "red", "Volume" = "darkgreen"))+xlab("Dose (Gy)")+ylab("Mean TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.title=element_blank(),legend.box.background=element_rect(linewidth = 1.5),legend.position=c(0.75,0.12), legend.direction="vertical",legend.margin=margin(1,1,1,0, unit='mm'),legend.spacing.x=unit(0, "mm"),legend.spacing.y=unit(0, "mm"),strip.text.x = element_text(size = size))+facet_grid(cols=vars("Heterogeneity"),labeller = label_parsed)+ylim(.0,1.)+scale_x_continuous(limits=c(0,120),breaks=c(0,20,40,60,80,100,120))


pdf("Figure_3.pdf",width=16,height=6)
ggarrange(p1, p2, font.label = list(size = size),ncol = 2, nrow = 1)
graphics.off()

#Plots of temperature (heterogeneity)
results<-read.table("../results/TCP_T_het.dat", header=T)

p1<-ggplot(results, aes(x=Temp, y=TCP))+geom_point(color="blue")+xlab("Mean Temperature (°C)")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),strip.text.x = element_text(size = size))+facet_grid(cols=vars("No~Heterogeneity"),labeller = label_parsed)+xlim(39,43)+ylim(.0,1.)
p2<-ggplot(results, aes(x=Temp, y=TCP_5))+geom_point(color="blue")+xlab("Mean Temperature (°C)")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),strip.text.x = element_text(size = size))+facet_grid(cols=vars('Heterogeneity~(c.v.~"="~"5%")'),labeller = label_parsed)+xlim(39,43)+ylim(.0,1.)
p3<-ggplot(results, aes(x=Temp, y=TCP_50))+geom_point(color="blue")+xlab("Mean Temperature (°C)")+ylab("TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),strip.text.x = element_text(size = size))+facet_grid(cols=vars('Heterogeneity~(c.v.~"="~"50%")'),labeller = label_parsed)+xlim(39,43)+ylim(.0,1.)


pdf("Figure_10.pdf",width=24,height=6)
ggarrange(p1, p2, p3, labels = c("A", "B","C"), font.label = list(size = size),ncol = 3, nrow = 1)
graphics.off()

#Plots of temperature (heterogeneity, TCP mean)
results<-read.table("../results/TCP_mean_T_het.dat", header=T)

p1<-ggplot(results, aes(x=Temp, y=TCP))+geom_line(color="blue", linewidth =1.5)+xlab("Mean Temperature (°C)")+ylab("Mean TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),strip.text.x = element_text(size = size))+facet_grid(cols=vars("No~Heterogeneity"),labeller = label_parsed)+xlim(39.5,42.5)+ylim(.5,1.)+geom_hline(yintercept=0.5533459312634031, linewidth=1.5)+annotate("text", x=42., y=0.58, label="Only RT", size=8)
p2<-ggplot(results, aes(x=Temp, y=TCP_5))+geom_line(color="blue", linewidth =1.5)+xlab("Mean Temperature (°C)")+ylab("Mean TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),strip.text.x = element_text(size = size))+facet_grid(cols=vars('Heterogeneity~(c.v.~"="~"5%")'),labeller = label_parsed)+xlim(39.5,42.5)+ylim(.5,1.)+geom_hline(yintercept=0.5110489601188253, linewidth=1.5)+annotate("text", x=42., y=0.54, label="Only RT", size=8)
p3<-ggplot(results, aes(x=Temp, y=TCP_50))+geom_line(color="blue", linewidth =1.5)+xlab("Mean Temperature (°C)")+ylab("Mean TCP")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),strip.text.x = element_text(size = size))+facet_grid(cols=vars('Heterogeneity~(c.v.~"="~"50%")'),labeller = label_parsed)+xlim(39.5,42.5)+ylim(.5,1.)+geom_hline(yintercept=0.5819018273528764, linewidth=1.5)+annotate("text", x=42., y=0.61, label="Only RT", size=8)


pdf("Figure_11.pdf",width=24,height=6)
ggarrange(p1, p2, p3, labels = c("A", "B","C"), font.label = list(size = size),ncol = 3, nrow = 1)
graphics.off()