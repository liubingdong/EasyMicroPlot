# EasyMicroPlot-- A easy R script plot  for Microbiome  analysis

EasyMicroPlot aims to be an **easy-to-use Microbiome  analysis pipeline** that accomplishes the core tasks of metagenomic analysis from start to finish: Data filter, Alpha diversity, Beta diversity, Co-occurence, Venn plot, Heat map, Key taxa screen,etc. Additionally, EasyMicroPlot could generated beautiful and qualified picture and tables for dissertation and paper. EasyMicroPlot is meant to be a fast and simple approach before you delve deeper into parameterization of your analysis. Each individual module of EasyMicroPlot is a standalone program, which means you can use only the modules you are interested in for your data.

## PACKAGES DEPENDENCIES 
* vegan (>= 2.5-6)
* ape (>= 5.3) 
* grid (>= 3.5.1)
* dplyr (>= 1.0.2)
* multcomp (>= 1.4-14)
* patchwork (>= 1.0.1)
* fs (>= 1.5.0)
* stringr (>= 1.4.0)
* plotly (>= 4.9.2.1)
* ggiraph (>= 0.7.0)
* agricolae（>= 1.3-3)
* randomForest (>= 4.6-14)
* doParallel (>= 1.0.15)
* parallel (>= 3.5.1)

## INSTALLATION

	if(! require("devtools")) install.packages("devtools")
	library(devtools)
	install_github("liubingdong/EasyMicroPlot",subdir='Version_0.4.2',
					upgrade = 'never')


				
## Main FUNCTION

* data_filter
* pca_boxplot
* beta_plot


## USAGE

**Easy to run (Recomened):**

	beta_result=beta_plot(group_level=c('Ctrl','CtrlD', 'HFD', 'HFDDH'),
	                      dir = '.',min_relative = 0.001,min_ratio = 0.85,
	                      design = 'mapping.txt',adjust = T,pattern = 'L',
	                      html_out = F,output=F,seed = 124)

**Step by step :**

```
## filtering the data
data=data_filter(dir = '.',design = 'mapping.txt',min_relative = 0.005,
				min_ratio = 0.7,adjust = T,pattern = 'L4',output=F)
## modifying the data,make SampleID as rowname
rownames(data)<-data[,1]
data<-data[,-1]
data<-subset(data,select=-c(Group))
## plot all
beta=pca_boxplot(data =data ,design = 'mapping.txt',seed=12,
					method='HSD',distance='bray')
## gereate html
beta=pca_boxplot(data =data ,design = 'mapping.txt',seed=12,
                    method='HSD',distance='bray')

## generate pdf
## better set seed as the fomer to make pdf as same as html output
pdf('sp_1_2.pdf',height = 15,width = 15)
beta_result$plot$species$p12
dev.off()
```

## UPDATE
* Version_0.4.3

  1 . **Add randomForest togther with N-fold cross validation** 
	 
	```
	generate filter data
	re=data_filter(dir = '.',min_relative = 0.001,min_ratio = 0.8,
	                 design = 'mapping.txt',adjust = F,pattern = 'L7')              
	modify data               
	rf=re$filter_data$species
	rf=subset(rf,select = -c(SampleID))
	RFCV
	result<-RFCVSEED(rep = 6,RF = rf,seed_start = 123,ntree = 10,core = 1,
	                              kfold = 5,RF_importance = 1,step = 1)
	plot=RFCV_plot(data = result,y_break = 1,
	                        cutoff_colour = 'red',palette = 'black')	
	```
 

 
   2 . **Add alpha caculation function**  
	```
	re=data_filter(dir = '.',min_relative = 0.001,min_ratio = 0.3,
				design = 'mapping.txt',adjust = F,pattern = 'L7')
	data=re$filter_data$species
	data=subset(data,select = -c(Group))
	rownames(data)<-data[,1]
	data<-data[,-1]
	data=round(data*1000000,0))
	alpha_result=alpha_caculate(data)
	alpha_pic=alpha_plot(data =alpha_result,design = 'mapping.txt' )
	```



* Version_0.4.2


  1 . **Rebuilde multible comparison** :
  	
  	```
	HSD: Multiple comparisons, Tukey
	
	LSD: Least significant difference
	
	duncan: Duncan's new multiple range test
	
	scheffe: Multiple comparisons, scheffe
	
	REGW: Ryan, Einot and Gabriel and Welsch multiple range test
	
	SNK: Student-Newman-Keuls
	```
  2 .  **Add more waring information and debug some tepo**
  
  ```
  When unequal replication in Tukey comparisons , now it could show correct warning.
  
  Add method infomation to the name of html files.
  
  ```







