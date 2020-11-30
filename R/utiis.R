
pca_boxplot=function(data,design,group_level=c('default'),method=c('Tukey'),distance=c('bray'),palette=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
                       "#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666")){
  library(vegan)
  library(ape)
  library(ggplot2)
  library(grid)
  library(dplyr)
  options(dplyr.summarise.inform = FALSE)
  library(multcomp)
  library(patchwork)
  library(fs)
  library(stringr)

  data=data
  mapping=read.table(paste0(design),header = T)
  mapping$Group=as.factor(mapping$Group)
  # 检查输入的group 新顺序是否符合要求，来决定是否使用默认字符串顺序
  name_group=unique(mapping$Group)
  group_level_check=all(name_group%in%group_level)&all(group_level%in%name_group)
  if (group_level_check==T) {
    mapping$Group<-factor(mapping$Group,levels = group_level)
  }else{
      if (group_level != 'default') {
        warning('group level can not match, and pipe will follow the default level !')
      }
    }
  groups <- data.frame(mapping$SampleID,mapping$Group)
  colnames(groups) <- c("V1","V2")
  length=length(unique(as.character(groups$V1)))
  times1=length%/%8
  res1=length%%8
  times2=length%/%5
  res2=length%%5
  col1=rep(1:8,times1)
  col=c(col1,1:res1)
  pich1=rep(c(21:24),times2)
  pich=c(pich1,15:(15+res2))
  ###########bray
  ###PCoA分析--轴12
  data <- vegdist(data,method=distance)
  pcoa<- pcoa(data, correction = "none", rn = NULL)
  PC1 = pcoa$vectors[,1]
  PC2 = pcoa$vectors[,2]
  PC3 = pcoa$vectors[,3]
  plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2,PC3,groups$V2)
  colnames(plotdata) <-c("sample","PC1","PC2","PC3","Group")
  pc1 <-round(pcoa$values$Relative_eig[1]*100,digits = 2)
  pc2 <-round(pcoa$values$Relative_eig[2]*100,digits = 2)
  pc3 <-round(pcoa$values$Relative_eig[3]*100,digits = 2)
  #plotdata$Group <- factor(plotdata$Group,levels = name_group)

  #PC1和PC2的显著性检验(PC1,PC2,PC3进行组间差异检验)
  yf <- plotdata
  yd1 <- yf %>% group_by(Group) %>% summarise(Max = max(PC1))
  yd2 <- yf %>% group_by(Group) %>% summarise(Max = max(PC2))
  yd3 <- yf %>% group_by(Group) %>% summarise(Max = max(PC3))
  yd1$Max <- yd1$Max + max(yd1$Max)*0.1
  yd2$Max <- yd2$Max + max(yd2$Max)*0.1
  yd3$Max <- yd3$Max + max(yd2$Max)*0.1
  fit1 <- aov(PC1~Group,data = plotdata)    #ANOVA检验≥3组样本
  tuk1<-glht(fit1,linfct=mcp(Group=method))#glht多重比较-Tukey(模型)
  res1 <- cld(tuk1,alpah=0.05)
  fit2 <- aov(PC2~Group,data = plotdata)
  tuk2<-glht(fit2,linfct=mcp(Group=method))
  res2 <- cld(tuk2,alpah=0.05)
  fit3 <- aov(PC3~Group,data = plotdata)
  tuk3<-glht(fit3,linfct=mcp(Group=method))
  res3 <- cld(tuk3,alpah=0.05)

  test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,PC3 = res3$mcletters$Letters,
                     yd1 = yd1$Max,yd2 = yd2$Max,yd3 = yd3$Max,Group = yd1$Group)
  test$Group <- factor(test$Group,levels = name_group)

  #相须图绘制
  p1 <- ggplot(plotdata,aes(Group,PC1)) +
    geom_boxplot(aes(fill = Group),outlier.colour = NA) +
    geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
              size = 7,color = "black",fontface = "bold") +
    coord_flip() +geom_point(position = "jitter",color="black",alpha=1)+
    scale_fill_manual(values=palette) +
    theme_bw()+
    theme(axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_text(colour='black',size=20,face = "bold"),
          axis.text.x=element_blank(),
          legend.position = "none")

  p2 <- ggplot(plotdata,aes(Group,PC2),outlier.colour = NA) +
    geom_boxplot(aes(fill = Group)) +
    geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
              size = 7,color = "black",fontface = "bold") +geom_point(position = "jitter",color="black",alpha=1)+
    scale_fill_manual(values=palette) +
    theme_bw()+
    theme(axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(colour='black',size=20,angle = 45,
                                   vjust = 1,hjust = 1,face = "bold"),
          axis.text.y=element_blank(),
          legend.position = "none")

  p3 <- ggplot(plotdata,aes(Group,PC3),outlier.colour = NA) +
    geom_boxplot(aes(fill = Group)) +
    geom_text(data = test,aes(x = Group,y = yd3,label = PC3),
              size = 7,color = "black",fontface = "bold") +geom_point(position = "jitter",color="black",alpha=1)+
    scale_fill_manual(values=palette) +
    theme_bw()+
    theme(axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(colour='black',size=20,angle = 45,
                                   vjust = 1,hjust = 1,face = "bold"),
          axis.text.y=element_blank(),
          legend.position = "none")

  #PCoA结果图绘制
  p12<-ggplot(plotdata, aes(PC1, PC2)) +
    geom_point(aes(fill=Group),size=8,pch = 21)+
    scale_fill_manual(values=palette,name = "Group")+
    xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) +
    ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
    xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
    ylim(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range) +
    theme(text=element_text(size=30))+
    geom_vline(aes(xintercept = 0),linetype="dotted")+
    geom_hline(aes(yintercept = 0),linetype="dotted")+
    theme(panel.background = element_rect(fill='white', colour='black'),
          panel.grid=element_blank(),
          axis.title = element_text(color='black',size=34),
          axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=34),
          axis.title.y=element_text(colour='black', size=34),
          axis.text=element_text(colour='black',size=28),
          legend.title=element_text(size = 24,face = "bold"),
          legend.text=element_text(size=20),
          legend.key=element_blank(),legend.position = c('left'),
          legend.background = element_rect(colour = "black"),
          legend.key.height=unit(1,"cm")) +
    guides(fill = guide_legend(ncol = 1))


  p13<-ggplot(plotdata, aes(PC1, PC3)) +
    geom_point(aes(fill=Group),size=8,pch = 21)+
    scale_fill_manual(values=palette,name = "Group")+
    xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) +
    ylab(paste("PC3 ( ",pc3,"%"," )",sep=""))+
    xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
    ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
    theme(text=element_text(size=30))+
    geom_vline(aes(xintercept = 0),linetype="dotted")+
    geom_hline(aes(yintercept = 0),linetype="dotted")+
    theme(panel.background = element_rect(fill='white', colour='black'),
          panel.grid=element_blank(),
          axis.title = element_text(color='black',size=34),
          axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=34),
          axis.title.y=element_text(colour='black', size=34),
          axis.text=element_text(colour='black',size=28),
          legend.title=element_text(size = 24,face = "bold"),
          legend.text=element_text(size=20),
          legend.key=element_blank(),legend.position = c('left'),
          legend.background = element_rect(colour = "black"),
          legend.key.height=unit(1,"cm")) +
    guides(fill = guide_legend(ncol = 1))

  p23<-ggplot(plotdata, aes(PC2, PC3)) +
    geom_point(aes(fill=Group),size=8,pch = 21)+
    scale_fill_manual(values=palette,name = "Group")+
    xlab(paste("PC2 ( ",pc1,"%"," )",sep="")) +
    ylab(paste("PC3 ( ",pc3,"%"," )",sep=""))+
    xlim(ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range) +
    ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
    theme(text=element_text(size=30))+
    geom_vline(aes(xintercept = 0),linetype="dotted")+
    geom_hline(aes(yintercept = 0),linetype="dotted")+
    theme(panel.background = element_rect(fill='white', colour='black'),
          panel.grid=element_blank(),
          axis.title = element_text(color='black',size=34),
          axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=34),
          axis.title.y=element_text(colour='black', size=34),
          axis.text=element_text(colour='black',size=28),
          legend.title=element_text(size = 24,face = "bold"),
          legend.text=element_text(size=20),
          legend.key=element_blank(),legend.position = c('left'),
          legend.background = element_rect(colour = "black"),
          legend.key.height=unit(1,"cm")) +
    guides(fill = guide_legend(ncol = 1))




  #PERMANOVA分析
  otu.adonis=adonis(data~V2,data = groups,distance = "bray")

  p5 <- ggplot() +
    geom_text(aes(x = -0.5,y = 0.6,
                  label = paste("PERMANOVA:\ndf = ",
                                otu.adonis$aov.tab$Df[1],
                                "\nR2 = ",round(otu.adonis$aov.tab$R2[1],4),
                                "\np-value = ",otu.adonis$aov.tab$`Pr(>F)`[1],sep = "")),
              size = 7) +
    theme_bw() +
    xlab("") + ylab("") +
    theme(panel.grid=element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())

  #图像拼接-使用patchwork包将4幅图拼在一起
  p12 <- p1 + p5 + p12 + p2 +
    plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
  p23 <- p2 + p5 + p23 + p3 +
    plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
  p13 <- p1 + p5 + p13 + p3 +
    plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
  deposit=list()
  deposit$p12=p12
  deposit$p23=p23
  deposit$p13=p13
  return(deposit)
}


modify_data=function(data,design, min_relative,min_odd) {
  library(dplyr)
  otu_origin=data
  mapping=design
  # 将第一列转为行名
  rownames(otu_origin)<-otu_origin[,1]
  otu_origin<-otu_origin[,-1]

  # 过滤掉小于千分之一的数据
  otu_origin[otu_origin<min_relative]=0

  #先矩阵化，再倒置
  otu_t=as.data.frame(t(as.matrix(otu_origin)))

  # 修改第一列
  otu_t=data.frame( SampleID=rownames(otu_t),otu_t)
  rownames(otu_t)=c(1:nrow(otu_t))

  #与mapping文件合并，获得数据切分依据
  otu_merge=full_join(otu_t,mapping,by="SampleID")
  otu_group_split=split(otu_merge,otu_merge$Group)

  name_group=unique(mapping$Group)

  # 创建ID空集
  data_for_filter=list()

  for (j in 1:length(name_group)) {
    # 循环读取各组数据
    otu=otu_group_split[[j]]
    # 删除前面合并的不必要的数据
    otu=subset(otu,select = -c(barcode,primer,Group,Description))

    # 更改数据名称为V，便于下游制图
    colnames(otu)=c("SampleID",paste("V", 1:c(ncol(otu)-1), sep = ""))

    # 将第一列转换为rownames
    otu=data.frame(otu,row.names = 1)

    # 利用mapping文件获得此次分组名称
    mapping=read.table("mapping.txt",header = T,sep = "\t")

    id_group=mapping$SampleID%in%rownames(otu)

    group_name=mapping$Group[id_group]

    name=as.character(group_name[1])

    # 统计为0的菌种比例
    no_zero_odd=apply(otu, 2, function(m){sum(m!=0)/nrow(otu)})

    # 过滤掉低占比的菌种，可根据情况设置20%，30%
    otu=otu[,no_zero_odd>=min_odd]
    data_for_filter[[j]]=colnames(otu)
  }

  ### 根据每个单独过滤组的变量取并集，得到所需的序列
  idx=Reduce(union,data_for_filter)

  # 获得根据过滤条件的适合随机森林的数据RF
  RF_raw=subset(otu_merge,select = -c(SampleID,barcode,primer,Description))
  ## 获得变量info
  var_names=colnames(RF_raw[,-ncol(RF_raw)])
  ## 重命名
  colnames(RF_raw)<-paste("V",1:ncol(RF_raw),sep="")

  RF_filter=RF_raw[,idx]
  # 得到所需数据
  RF=data.frame(Group=otu_merge$Group,RF_filter)

  # 变量存储
  deposit=list()
  deposit$mapping_file=mapping
  deposit$otu_group_split=otu_group_split
  deposit$union_id=idx
  deposit$var_names=var_names
  deposit$filtered_data=RF
  return(deposit)
}


data_filter=function(dir,min_relative,min_ratio,design,adjust=F,output=F,pattern=''){
  file_name<-list.files(path =dir ,pattern = pattern)
  file_name=file_name[which(file_name!=design)]
  file_num=length(file_name)
  sep_num=c()
  tax_total=c('phylum','class','order','family','genus','species')
  mapping=read.table(paste0(design),header = T)
  deposit=list()
  for (i in c(1:file_num)){
    file_data=read.table(file_name[i],sep='\t',header=T,check.names = F,row.names= 1)
    sep_num[i]=str_count(rownames(file_data)[1],pattern = ';')
  }
  if (length(unique(sep_num))<file_num) {
    warning('Data level can not match well, plz check files in different microbial level!')
    warning('Duplicated level detcted!')
    warning('input this to acquire more error information: View(info_data)')
  }else{
    for (i in c(1:file_num)){
      tax_level=tax_total[sep_num[i]]
      data=read.table(file_name[i], sep="\t", header=T,check.names = F)
      SampleID=colnames(data)[-1]
      filter_result=modify_data(data=data,design = mapping,min_relative = min_relative,min_odd =min_ratio )
      sub_data=filter_result$filtered_data
      sub_data=data.frame(SampleID=SampleID,sub_data)

      # check empty data
      sub_data_check=subset(sub_data,select=-c(SampleID,Group))
      if (ncol(sub_data_check) != 0) {
        id_na<-apply(sub_data_check,1,max)==0
        
        empty_data_name=paste0(tax_level,'_empty_data')
        empty_data=data.frame()
        try(empty_data<-sub_data[id_na,],silent=T)
        if (length(rownames(empty_data))>0) {
          warning('Empty data detected!')
          warning(paste0('input this to acquire more error information: View(<output>$empty_data$',empty_data_name,')'))
          deposit$empty_data[[paste0(empty_data_name)]]=empty_data
        }
        
        # 设置极小值校正特征 使空行可以继续
        if ((adjust==T)&(length(rownames(empty_data))>0)) {
          empty_ID=sub_data$SampleID[id_na]
          empty_ID_adjust=paste0(empty_ID,'_adjust')
          for (k in c(1:length(empty_ID_adjust))) {
            sub_data[empty_ID_adjust[k]]=rep(0,nrow(sub_data))
            idx<-sub_data$SampleID==empty_ID[k]
            sub_data[empty_ID_adjust[k]][idx,]=0.00001
          }
        }
        
        if (output==T){
          dir.create(tax_total[sep_num[i]])
          write.csv(sub_data, file =paste0(tax_level,"_",min_ratio*100,"%.csv"),row.names = F)
          write.csv(filter_result$var_names, file =paste0(tax_level,"_",min_ratio*100,"%_info.csv"))
          file_move(paste0(tax_level,"_",min_ratio*100,"%.csv"),tax_level)
          file_move(paste0(tax_level,"_",min_ratio*100,"%_info.csv"),tax_level)
        }
        deposit$filter_data[[paste0(tax_level)]]=sub_data
      }else{
        print(paste0(file_name[i],' has been filtered into empty data, please reset the filter parameter !'))
        deposit$filter_data[[paste0(tax_level)]]='empty data'
      }
    }
  }
  return(deposit)
}

beta_plot=function(dir,group_level=c('default'),min_relative = 0,min_ratio = 0,design ,adjust = F,pattern = '',output = F,
                   method='Tukey',distance = 'bray',palette=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
                                                              "#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666") ){
  library(vegan)
  library(ape)
  library(ggplot2)
  library(grid)
  library(dplyr)
  options(dplyr.summarise.inform = FALSE)
  library(multcomp)
  library(patchwork)
  library(fs)
  library(stringr)
  deposit=list()
  deposit$result=data_filter(dir = dir,min_relative = min_relative,min_ratio = min_ratio,design = design,adjust = adjust,pattern = pattern,output = output)
  for (i in names(deposit$result$filter_data)){
    data=data.frame()
    try(data<-subset(deposit$result$filter_data[[i]],select=-c(SampleID,Group)),silent=T)
    if (ncol(data) != 0) {
      deposit$plot[[i]]<-pca_boxplot(data =data ,design = design,group_level=group_level)
    }
  }
  return(deposit)
}






