library(ggplot2)
library(gridExtra)

## Manuscript Figure 3
plot1 <- ggplot(NewTools, aes(x=Sensitivity, y=Precision)) +
  xlim(0,1) +
  ylim(0,1 ) +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=bquote(Recall^est), y=bquote(Precision^est), title = "Individual Tools")+
  geom_text(label= NewTools$Label, fontface=ifelse(NewTools$Cohort == "HTMCP","bold", "plain"), size=6)

plot2 <- ggplot() +
  xlim(0, 1)+
  ylim(0,1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x=bquote(Recall^est), y=bquote(Precision^est), title = "Combinations and FFPolish")+
  geom_point(data = NewCombi,aes(x=Sensitivity, y=Precision), size=5, shape=ifelse(NewCombi$Cohort == "HTMCP",16, 1)) +
  geom_point(data = NewCombi1, aes(x=Sensitivity, y=Precision), size=5, shape=ifelse(NewCombi1$Cohort == "HTMCP",15, 0))+
  geom_point(data = Newffpe, aes(x=Sensitivity, y=Precision), size=5, shape=ifelse(Newffpe$Cohort == "HTMCP",17, 2))
  
complete <- grid.arrange(plot1, plot2, ncol=2)


## Supplementary Section 5 - Correlations 
library("ggpubr")
a <- ggscatter(CorrResultsBLGSP, x = "Precision", y = "Duplicate_reads" , 
               color = "Tool",
               add = "reg.line", conf.int = FALSE, 
               cor.coef = FALSE, cor.method = "pearson",
               xlab = "Precision", ylab = "Duplicate Reads", 
               title = "Correlation of Precision and Duplicate Reads BLGSP") + 
  stat_cor(aes(color = Tool), label.x = 0.25)   
plot(a)
a <- ggscatter(CorrResultsBLGSP, x = "Sensitivity", y = "mean_mapping_quality" , 
               color = "Tool",
               add = "reg.line", conf.int = FALSE, 
               cor.coef = FALSE, cor.method = "pearson",
               xlab = "Sensitivity", ylab = "Mean Mapping Quality", 
               title = "Correlation of Sensitivity and Mean Mapping Quality BLGSP")
a + stat_cor(aes(color = Tool), label.x = 0.25)

a <- ggplot(PrecHTMCP, aes(x=Metric, y=Correlation.Coefficient..r., fill = Tool)) + 
  geom_boxplot(fill="white", outlier.color=NA)+ geom_dotplot(binaxis='y', stackdir='center', 
                                                             dotsize=1.2, alpha=0.7,
                                                             position=position_dodge(width=0.5),
                                                             colour=ifelse(PrecHTMCP$p.value < 0.05, "red", "black"),
                                                             stroke=ifelse(PrecHTMCP$p.value < 0.05, "2", "1")) +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text.x = element_blank(),
    axis.title.x=element_blank(), axis.title.y=element_blank()) +
  labs(x="", y="Correlation Coefficient (r)") + ylim(-1,1) +
  ggtitle("Precision HTMCP")
a
b <- ggplot(SensHTMCP, aes(x=Metric, y=Correlation.Coefficient..r., fill = Tool)) + 
  geom_boxplot(fill="white", outlier.color=NA)+ geom_dotplot(binaxis='y', stackdir='center', 
                                                             dotsize=1.2, alpha=0.7,
                                                             position=position_dodge(width=0.5),
                                                             colour=ifelse(SensHTMCP$p.value < 0.05, "red", "black"),
                                                             stroke=ifelse(SensHTMCP$p.value < 0.05, "2", "1")) +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    axis.text.x = element_blank(),
    axis.title.x=element_blank(), axis.title.y=element_blank()
  )  + 
  labs(x="", y="Correlation Coefficient (r)") + ylim(-1,1) +
  ggtitle("Sensitivity HTMCP")
b
c <- ggplot(PrecBLGSP, aes(x=Metric, y=Correlation.Coefficient..r., fill = Tool)) + 
  geom_boxplot(fill="white", outlier.color=NA)+ geom_dotplot(binaxis='y', stackdir='center', 
                                                             dotsize=1.2, alpha=0.7,
                                                             position=position_dodge(width=0.5),
                                                             colour=ifelse(PrecBLGSP$p.value < 0.05, "red", "black"),
                                                             stroke=ifelse(PrecBLGSP$p.value < 0.05, "2", "1")) +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
    axis.text.x = element_blank(),
    axis.title.x=element_blank())+
  labs(x="", y="Correlation Coefficient (r)") + ylim(-1,1) +
  ggtitle("Precision BLGSP")
c
d <- ggplot(SensBLGSP, aes(x=Metric, y=Correlation.Coefficient..r., fill = Tool)) + 
  geom_boxplot(fill="white",outlier.color=NA)+ geom_dotplot(binaxis='y', stackdir='center', 
                                                            dotsize=1.2, alpha=0.7,
                                                            position=position_dodge(width=0.5),
                                                            colour=ifelse(SensBLGSP$p.value < 0.05, "red", "black"),
                                                            stroke=ifelse(SensBLGSP$p.value < 0.05, "2", "1")) +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.text.x = element_blank(),
    axis.title.x=element_blank(), axis.title.y=element_blank())   + 
  labs(x="", y="Correlation Coefficient (r)") + ylim(-1,1) +
  ggtitle("Sensitivity BLGSP")
d
figure <- ggarrange(a, b, c, d,
                    labels = c("A ", "B", "C", "D"),
                    ncol = 2, nrow = 2, common.legend = TRUE )
figure


## Supplementary Section 6 - Individual Tools 

BLGSP_MatchedPlot <- ggscatterhist(BLGSP_Matched, x="Specificity", y="Sensitivity", size = 3, color="Tool", 
                                   alpha=0.9, margin.plot = "histogram", 
                                   margin.params = list(fill="Tool", color="black", size=0.2), 
                                   bins =10, main.plot.size = 1, margin.plot.size = 1, xlim=c(0,1),
                                   title = "BLGSP FFPE with matched normal")
plot(BLGSP_MatchedPlot)


HTMCP_MatchedPlot <- ggscatterhist(HTMCP_Matched, x="Specificity", y="Sensitivity", size = 3, color="Tool", 
                                   alpha=0.9, margin.plot = "histogram", 
                                   margin.params = list(fill="Tool", color="black", size=0.2), 
                                   bins =10, main.plot.size = 1, margin.plot.size = 1, xlim=c(0,1),
                                   title = "HTMCP FFPE with matched normal")
plot(HTMCP_MatchedPlot)

## Supplementary Section 7 - Results in comparison with Brienen et. al

plot(FinalCombi_HTMCP$Sensitivity, FinalCombi_HTMCP$Precision, xlim=c(0,1), ylim=c(0,1), xlab="Median sensitivity",
     ylab="Median precision", main="FFPE SNV Caller - Combined Results HTMCP",
     cex=ifelse(FinalCombi_HTMCP$Type=='Indi', 3, 1),
     pch=ifelse(FinalCombi_HTMCP$Type=='Indi', 19, 4),
     col=c('#E64B35FF', '#4DBBD5FF', '#00A087FF',
           '#3C5488FF', '#F39B7FFF', '#7E6148FF',
           '#A9A9A9', '#A9A9A9', '#A9A9A9',
           '#A9A9A9', '#A9A9A9', '#A9A9A9',
           '#696969', '#696969', '#696969'))


Groups <- c("LoFreq", "Strelka2",
            "Shimmer", "Virmid", "Mutect2", "Best F1 scores", "Best sensitivity or precision")
legend(x=-0.03, y=0.9, bty="n", legend=Groups, 
       fill=c('#E64B35FF', '#00A087FF',
              '#3C5488FF', '#F39B7FFF', '#7E6148FF',
              '#696969', '#DCDCDC'), cex=1.2)

plot(Finalcombi$Sensitivity, Finalcombi$Precision, axes=TRUE, xlim=c(0,1), ylim=c(0,1), xlab="Median sensitivity",
     ylab="Median precision",
     asp=0, cex=ifelse(Finalcombi$Type=='Indi', 3, 1),
     pch=ifelse(Finalcombi$Type=='Indi', 19, 4),
     col=c('#E64B35FF', '#4DBBD5FF', '#00A087FF',
           '#3C5488FF', '#F39B7FFF', '#7E6148FF',
           '#A9A9A9', '#A9A9A9', '#A9A9A9',
           '#A9A9A9', '#A9A9A9', '#A9A9A9',
           '#696969', '#696969', '#696969'))

ggscatter(data=Finalcombi, x="Sensitivity", y="Precision")

scatterDF <- na.omit(Scatter)
FinalcombiDF <- na.omit(Finalcombi)

a <- ggplot(data=FinalcombiDF, aes(x=Sensitivity, y=Precision) )+ 
  geom_point(colour = c('#E64B35FF', '#00A087FF', '#3C5488FF', 
                        '#F39B7FFF', '#7E6148FF',
                        '#A9A9A9', '#A9A9A9', '#A9A9A9',
                        '#A9A9A9', '#A9A9A9', '#A9A9A9',
                        '#696969', '#696969', '#696969',
                        '#ae22ff'), 
             pch=c(19, 19, 19,
                   19, 19,
                   15, 15, 15,
                   15, 15, 15,
                   15, 15, 15,
                   19),
             size=c(6,6,6,
                    6,6,
                    2,2,2,
                    2,2,2,
                    2,2,2,
                    6), 
             alpha=c(1, 1,1,
                     1,1,
                     1,1,1,
                     1,1,1,
                     1,1,1,
                     0.8)) +
  geom_scatterpie(aes(x=Sensitivity, y=Precision, group=Region, r=0.015), alpha=0.8,data=scatterDF, 
                  cols=c("LoFreq", "Strelka2", "Shimmer", "Virmid", "Mutect2"),color=NA) +
  scale_fill_manual(values = c('#E64B35FF', '#7E6148FF', '#3C5488FF', '#00A087FF', 
                               '#F39B7FFF')) + coord_equal(xlim = c(0,1), ylim = c(0,1))+ theme_classic()+theme(panel.grid.major = element_blank(), 
                                                                                                                panel.grid.minor = element_blank(),
                                                                                                                panel.background = element_rect(colour = "black", size=0.7)) +
  ggtitle("BLGSP")
a

scatterDF <- na.omit(scatter_HTMCP)
FinalcombiDF <- na.omit(FinalCombi_HTMCP)

b <- ggplot(data=FinalcombiDF, aes(x=Sensitivity, y=Precision) )+ 
  geom_point(colour = c('#E64B35FF', '#00A087FF', '#3C5488FF', 
                        '#F39B7FFF', '#7E6148FF',
                        '#A9A9A9', '#A9A9A9', '#A9A9A9',
                        '#A9A9A9', '#A9A9A9', '#A9A9A9',
                        '#696969', '#696969', '#696969',
                        '#ae22ff'), 
             pch=c(19, 19, 19,
                   19, 19,
                   15, 15, 15,
                   15, 15, 15,
                   15, 15, 15,
                   19),
             size=c(6,6,6,
                    6,6,
                    2,2,2,
                    2,2,2,
                    2,2,2,
                    6), 
             alpha=c(1, 1, 1,
                     1,1,
                     1,1,1,
                     1,1,1,
                     1,1,1,
                     0.8)) +
  geom_scatterpie(aes(x=Sensitivity, y=Precision, group=Region, r=0.015), alpha=0.8,data=scatterDF, 
                  cols=c("LoFreq", "Strelka2", "Shimmer", "Virmid", "Mutect2"),color=NA) +
  scale_fill_manual(values = c('#E64B35FF', '#7E6148FF', '#3C5488FF', '#00A087FF', 
                               '#F39B7FFF')) + coord_equal(xlim = c(0,1), ylim = c(0,1))+ theme_classic()+theme(panel.grid.major = element_blank(), 
                                                                                                                panel.grid.minor = element_blank(),
                                                                                                                panel.background = element_rect(colour = "black", size=0.7))+
  ggtitle("HTMCP")
b
figure1 <- ggarrange(a, b,
                     legend = "none" )
figure1