#' @import ggplot2
#' @import randomForest
#' @import plyr
#' @import dplyr




RFCVSEED= function(RF,seed_start=123,ntree=1000,core=1,kfold=5,rep=10,RF_importance=1,step=1,each_ouput=F){
  
  value=c("MeanDecreaseAccuracy","MeanDecreaseGini")
  value=value[RF_importance]
  
  result=list()
  RF_var=list()
  RF_seed_plot=list()
  # R升级4.0后需要指明因子
  RF$Group=as.factor(RF$Group)
  RFCV= function(i){
    # 删除id为i的行，创建训练集
    # 选id为i的行，创建训练集
    trainingset <- subset(data, id !=i)
    testset <- subset(data, id==i)
    #运行一个随机森林模型
    set.seed(seed)
    mymodel <- randomForest(trainingset$Group ~ ., data = trainingset[,-ncol(trainingset)], ntree = ntree,keep.forest=T,importance=TRUE)
    #去掉回应列1
    testdata=subset(testset,select =-c(Group,id) )
    temp <- as.data.frame(predict(mymodel, testdata))
    idx=(temp[,1]==testset$Group)
    num=rep(1,length(idx))
    sum(num[idx])/length(idx)
  }
  
  seed_confirm=c()
  seed_check=seed_start
  while (length(seed_confirm)<rep) {
    set.seed(seed_check)
    num=length(unique(sample(1:kfold, nrow(RF), replace = TRUE)))
    if (num==kfold) {
      seed_confirm=append(seed_confirm,seed_check)
    }
    seed_check=seed_check+1
  }
  
  for (seed in seed_confirm) {
    set.seed(seed)
    mymodel <- randomForest(RF$Group ~ ., data = RF, ntree = ntree,keep.forest=T,importance=TRUE)
    importance=as.data.frame(importance(mymodel))
    importance_order=importance[order(importance[,value]),]
    if (each_ouput == T){
    pdf(paste0("seed",seed,"_varimplot.pdf"))
    varImpPlot(mymodel)
    dev.off()
    }
    # 获取倒序剔除特征值顺序
    del_order=row.names(importance_order)
    RF_var[[seed-seed_start+1]]=rev(del_order)
    
    
    accuracy_step=list()
    
    # 决定step_order的，注意不要削减至空数据集 将无法训练
    data_filter_decide<-c( (ncol(RF)-1) /step)==round((ncol(RF)-1)/step)
    if (data_filter_decide&&step!=1) {
      step_order=seq(0,ncol(RF)-3,by=step)
    }else{
      step_order=seq(0,ncol(RF)-2,by=step)
    }
    
    
    progress.bar <- create_progress_bar("text")
    progress.bar$init(length(step_order))
    
    for (j in 1:length(step_order)) {
      # 过滤数据
      if (step_order[j]==0) { data=RF
      }else{
        id_retain=!(colnames(RF)%in%del_order[1:step_order[j]])
        data=RF[id_retain]}
      # 并行开始
      set.seed(seed)
      k=kfold
      data$id <- sample(1:k, nrow(data), replace = TRUE)
      #####
      #cl <- makeCluster(core)
      #registerDoParallel(cl)
      #accuracy=foreach(i=1:k,.packages="randomForest") %dopar% RFCV(i)
      #print(mean(unlist(accuracy)))
      #accuracy_step[[j]]=mean(unlist(accuracy))
      #stopCluster(cl)
      #####
      accuracy=c()
      for (i in 1:2) {
        accuracy=append(accuracy,RFCV(i))
      }
      accuracy_step[[j]]= mean(accuracy)
      progress.bar$step()
    }
    result[[seed-seed_start+1]]=accuracy_step
    
    
    print(paste0("RF_CV_",seed," is done!"))
    
  }
  names(result)=paste( "Seed_", seed_confirm, sep="")
  names(RF_var)=paste( "Seed_", seed_confirm, sep="")
  deposit=list()
  deposit$CV_accuracy=result
  deposit$RF_var=RF_var
  deposit$Seed_num=seed_confirm
  deposit$Trails_num=length(seed_confirm)
  deposit$step_order=step_order
  deposit$Feature_num=ncol(RF)-1
  return(deposit)
}





RFCV_plot= function(data,y_break=1,cutoff_colour=c("red"),palette=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
                                                                    "#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666")){
  # 确定色板颜色
  if (length(palette)==1) {
    palette=rep(palette,data$Trails_num)
  }
  deposit=list()
  Test=list()
  for (i in 1:data$Trails_num) {
    Test[[i]]=rep(paste0("Trail_",i),length(data$step_order))
  }
  
  Test=unlist(Test)
  
  xlab=rep(data$Feature_num-data$step_order,data$Trails_num)
  y=data$CV_accuracy
  
  for (i in 1:data$Trails_num) {y[[i]]=unlist(y[[i]])
  }
  
  ylab=unlist(y)
  f=data.frame(xlab,ylab,Test)
  
  y2=as.data.frame(y)
  total_mean=apply(y2, 1, mean)
  total_sd=apply(y2, 1, sd)
  total_points=1-total_mean+total_sd
  
  max_mean=max(unlist(total_mean))
  cut_off=data$Feature_num - data$step_order[which(unlist(total_points)%in%min(unlist(total_points)))]
  
  opt_mean=mean(as.numeric(y2[cut_off,]))
  
  p=ggplot(f,aes(x=xlab,y=(1-ylab),fill=Test,color=Test))+geom_line()+xlab( "Number of var")+ylab("CV Error")
  
  curve_plot=p+geom_vline(xintercept = cut_off,colour=cutoff_colour)+scale_x_reverse(breaks=seq(0, data$Feature_num, y_break))+
    labs(title=paste0("Max mean accuracy is ",max_mean," and ",
                      "Optional mean accuracy is ",opt_mean,"  ; ","Optional number is ",cut_off))+ theme_bw()+scale_colour_manual(values=palette)
  
  opt_random_num=colnames(y2)[as.logical(y2[cut_off,]==max(y2[cut_off,]))]
  
  # 根据cut_off 确认变量交并集
  var_select=list()
  for (i in names(data$RF_var)) {
    var_select[[i]]=data$RF_var[[i]][1:cut_off]
  }
  
  intersect_num=Reduce(intersect,var_select)
  union_num=Reduce(union,var_select)
  
  deposit$curve_plot=curve_plot
  deposit$opt_num=cut_off
  deposit$opt_random_num=opt_random_num
  deposit$intersect_num=intersect_num
  deposit$union_num=union_num
  return(deposit)
}




