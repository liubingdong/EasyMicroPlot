#' @import ggplot2
#' @import vegan
#' @import ape
#' @import grid
#' @import dplyr
#' @import multcomp
#' @import patchwork
#' @import fs
#' @import stringr
#' @import ggiraph
#' @import plotly
#' @import agricolae

options(dplyr.summarise.inform = FALSE)

alpha_caculate <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  #Chao1 <- est[2, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  result <- data.frame(Richness, Shannon, Simpson)
  return(result)
}

alpha_plot <- function(data,design){
  data=data
  try(mapping<-read.table(paste0(design),header = T),silent = T)
  mapping$Group=as.factor(mapping$Group)
  real_sample=Reduce(intersect,list(mapping$SampleID,rownames(data)))
  mapping=mapping[mapping$SampleID%in%real_sample,]
  data$SampleID=rownames(data)
  data=full_join(design,data,by='SampleID')
  try(data<-subset(data,select = -c(SampleID,barcode,primer,Description)),silent=T)
  group_name=as.character(unique(data$Group))
  group_combn=combn(group_name,2)
  compare=alply(group_combn,2)
  

  Richness=ggplot(data=data,aes(x=Group,y=data[,'Richness'] ))+geom_boxplot(aes(fill=Group),colour="black",notch=F,outlier.colour=NA)+geom_point(position = "jitter",color="black",alpha=1)+
      ylab('Richness')+stat_compare_means(comparisons = compare,method="t.test")+
      theme_bw()
  Shannon=ggplot(data=data,aes(x=Group,y=data[,'Shannon'] ))+geom_boxplot(aes(fill=Group),colour="black",notch=F,outlier.colour=NA)+geom_point(position = "jitter",color="black",alpha=1)+
    ylab('Shannon')+stat_compare_means(comparisons = compare,method="t.test")+
    theme_bw()
  Simpson=ggplot(data=data,aes(x=Group,y=data[,'Simpson'] ))+geom_boxplot(aes(fill=Group),colour="black",notch=F,outlier.colour=NA)+geom_point(position = "jitter",color="black",alpha=1)+
    ylab('Simpson')+stat_compare_means(comparisons = compare,method="t.test")+
    theme_bw()
  
  deposit=list()
  deposit$Richness=Richness
  deposit$Shannon=Shannon
  deposit$Simpson=Simpson
  return(deposit)
}