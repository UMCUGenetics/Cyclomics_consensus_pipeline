library(ggplot2)
library(gridExtra)
library(grid)
library(stringr)

#### Functions ####
args = commandArgs(trailingOnly=TRUE)
input=args[1]

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)

  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                 ncol = cols, nrow = ceiling(numPlots/cols))
}

if (numPlots == 1) {
print(plots[[1]])

} else {
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

for (i in 1:numPlots) {
  matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

  print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                  layout.pos.col = matchidx$col))
 }
}
}

#### START  ####
args = commandArgs(trailingOnly=TRUE)
input=args[1]
name<-(substr(args[1], 0 ,str_length(args[1])-4 ))
df = read.table(input,sep="\t", header=TRUE) 
max_reads=args[2]

df_lines <- head(subset(df,B>=0 & I>=0), max_reads)
max_x=max(df_lines$B)
max_y=max(df_lines$I)
if (max_x>max_y){ 
max<-max_x
} else {
max<-max_y
}

p1 <- ggplot(df_lines, aes(x=I, y=B)) +
  geom_point(alpha = 1/20) +
  labs(title = paste("Backbone vs Insert (max_values=",max_reads,")")) + 
  xlab("#Insert")+
  ylab("#Backbone")+
  scale_x_continuous(limits=c(0, max))+
  scale_y_continuous(limits=c(0, max))


sub_ratio<-subset(subset(df,BI>=0.8), BI<=1.2)
sub_df1<-subset(df,BI>0)
perc_BI <-(nrow(sub_ratio))/(nrow(sub_df1))*100
max_x=10
bin_size=0.1
max_y =(layer_scales(ggplot(data=sub_df1, aes(sub_df1$BI)) + geom_histogram(breaks=seq(0, max_x , bin_size)))$y$range$range)[2]
p2 <- ggplot(data=sub_df1, aes(sub_df1$BI)) + 
  geom_histogram(breaks=seq(0, max_x, bin_size), col="black", fill="gray", alpha = .8) + 
  labs(title="Histogram B/I ratio") +
  labs(x="#B/I ratio", y="Count") +
  scale_x_continuous(breaks=seq(0, max_x, 1),limits=c(0, max_x))+
  #scale_x_log10(breaks=c(0.01,0.1,1,10,100,1000),limits=c(0.01,1000))+
  geom_vline(aes(xintercept=median(sub_df1$BI)),color="red", linetype="dashed", size=0.5)+
  annotate("label", x = max_x/2, y = (max_y*0.9), label = paste("Percentage >=0.8 and <=1.2","\n",round(perc_BI, digits = 2),"%")) 

max_x=40
sub_df2<-subset(df,MeanBaseQualityMappedRead>0)
sub_ratio<-subset(subset(df,MeanBaseQualityMappedRead>=13))
bin_size=1
max_y = (layer_scales(ggplot(data=sub_df2, aes(sub_df2$MeanBaseQualityMappedRead)) + geom_histogram(breaks=seq(0, max_x , bin_size)))$y$range$range)[2]
perc_Q13<-(nrow(sub_ratio))/(nrow(sub_df2))*100
p3 <- ggplot(data=sub_df2, aes(sub_df2$MeanBaseQualityMappedRead)) + 
  geom_histogram(breaks=seq(0, max_x, bin_size), col="black", fill="gray", alpha = .8) + 
  labs(title="Histogram Mean BaseQuality") +
  labs(x="Mean BaseQuality", y="Count") +
  scale_x_continuous(breaks=seq(0, max_x, 2),limits=c(1, max_x))  + 
  geom_vline(aes(xintercept=median(sub_df2$MeanBaseQualityMappedRead)),color="red", linetype="dashed", size=0.5)+
  annotate("label", x = max_x/1.5 , y = (max_y*0.9), label = paste("Percentage >=Q13","\n",round(perc_Q13, digits = 2),"%")) 


sub_df_ins1<-subset(df,DifInserts==1)
sub_df_insall<-subset(subset(df,DifInserts>0))
perc_perfectIns= (nrow(sub_df_ins1))/(nrow(sub_df_insall))*100
max_y=nrow(sub_df_ins1)
p4 <- ggplot(data=df, aes(df$DifInserts)) + 
  geom_histogram(breaks=seq(0, 10, 1), col="black", fill="gray", alpha = .8) + 
  labs(title="Histogram #Different Insert Fragments") +
  labs(x="#Insert Fragments", y="Count") +
  scale_x_continuous(breaks=seq(0, 10, 1),limits=c(0, 10))+
  geom_vline(aes(xintercept=median(subset(df,DifInserts>0)$DifInserts)),color="red", linetype="dashed", size=0.5)+
  annotate("label", x = 5 , y = (max_y*0.9), label = paste("Percentage Insert = 1","\n",round(perc_perfectIns, digits = 2),"%"))

sub_df_bb1<-subset(df,DifBackbone==1)
sub_df_bbmore1<-subset(df,DifBackbone>1)
sub_df_bball<-subset(subset(df,DifBackbone>0))

perc_perfectBb= (nrow(sub_df_bb1))/(nrow(sub_df_bball))*100
perc_more1Bb= (nrow(sub_df_bbmore1))/(nrow(sub_df_bball))*100
max_y=nrow(sub_df_bb1)
p5 <- ggplot(data=df, aes(df$DifBackbone)) +
  geom_histogram(breaks=seq(0, 10, 1), col="black", fill="gray", alpha = .8) +
  labs(title="Histogram #Different Backbone Fragments") +
  labs(x="#Backbone Fragments", y="Count") +
  scale_x_continuous(breaks=seq(0, 10, 1),limits=c(0, 10))+
  geom_vline(aes(xintercept=median(subset(df,DifBackbone>0)$DifBackbone)),color="red", linetype="dashed", size=0.5)+
  annotate("label", x = 5 , y = (max_y*0.9), label = paste("Percentage Backbone = 1","\n",round(perc_perfectBb, digits = 2),"%"))

p6 <- ggplot(data=df, aes(df$MADinsert)) + 
  geom_histogram(breaks=seq(0, 100, 1), col="black", fill="gray", alpha = .8) + 
  labs(title="Histogram MAD insert startsite (in-target only)") +
  labs(x="MAD insert", y="Count")+
  scale_x_continuous(breaks=seq(0, 100, 10),limits=c(0, 100))+
  geom_vline(aes(xintercept=median(subset(df,MADinsert>0)$MADinsert)),color="red", linetype="dashed", size=0.5)

pdf(paste(name,".pdf",sep=""),width=10,height=10)
multiplot(p1,p3,p2,p4,p5,p6,cols=2)
dev.off()




