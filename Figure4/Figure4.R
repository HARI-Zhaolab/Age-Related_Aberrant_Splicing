setwd('/home/cmx/senes_splicing/')
rm(list = ls())
options(stringsAsFactors = F)

library(dplyr)
library(TCGAbiolinks)

load("./outdata/step1_cluster_os.Rdata")
# query <- GDCquery(project = "TCGA-BRCA",
#                   data.category = "Copy Number Variation",
#                   data.type = "Masked Copy Number Segment")
# 
# BRCA_CNV_download <- GDCprepare(query = query, save = TRUE, save.filename = "./outdata/BRCA_CNV_download.rda")

getwd()
A=load("./outdata/BRCA_CNV_download.rda")
tumorCNV <- eval(parse(text = A))

#改名
tumorCNV=tumorCNV[,2:7]
tumorCNV=tumorCNV[,c('Sample','Chromosome','Start','End','Num_Probes','Segment_Mean')]
#write.table(tumorCNV,file = 'BRCA_CNV.txt',sep = '\t',quote = F,row.names = F)
##只要带01A的
tumorCNV_01A <- tumorCNV[grepl("01A", tumorCNV$Sample), ]
#write.table(tumorCNV_01A,file = 'segment_file.txt',sep = '\t',quote = F,row.names = F)

colnames(tumorCNV_01A)

#分三组
##第一组：
# tumorCNV_01A$cluster <- substr(tumorCNV_01A$Sample,1,12)
# tumorCNV_01A <- tumorCNV_01A[which(tumorCNV_01A$cluster %in% clinical_OS[which(clinical_OS$group1 == "Cluster1"),]$sample),]
# tumorCNV_01A <- tumorCNV_01A[,-7]
# write.table(tumorCNV_01A,file = './outdata/MaskedCopyNumberSegment1.txt',sep = '\t',quote = F,row.names = F)

# tumorCNV_01A$cluster <- substr(tumorCNV_01A$Sample,1,12)
# tumorCNV_01A <- tumorCNV_01A[which(tumorCNV_01A$cluster %in% clinical_OS[which(clinical_OS$group1 == "Cluster2"),]$sample),]
# tumorCNV_01A <- tumorCNV_01A[,-7]
# write.table(tumorCNV_01A,file = './outdata/MaskedCopyNumberSegment2.txt',sep = '\t',quote = F,row.names = F)

tumorCNV_01A$cluster <- substr(tumorCNV_01A$Sample,1,12)
tumorCNV_01A <- tumorCNV_01A[which(tumorCNV_01A$cluster %in% clinical_OS[which(clinical_OS$group1 == "Cluster3"),]$sample),]
tumorCNV_01A <- tumorCNV_01A[,-7]
write.table(tumorCNV_01A,file = './outdata/MaskedCopyNumberSegment3.txt',sep = '\t',quote = F,row.names = F)
###必须过滤freqcnv=FALSE
filemarker<-read.table("./indata/snp6.na35.remap.hg38.subset.txt")
colnames(filemarker) <- filemarker[1,]       ## 列名定义为第一行
filemarker <- filemarker[-1,]    ## 删除第一行实现将第一行转换为列名
filemarker<- filemarker[grepl("FALSE",filemarker$freqcnv), ]
filemarker<-filemarker[,c(1,2,3)]
colnames(filemarker)<-c("Marker_name","Chromosome","Marker_position")
# write.table(filemarker,file = './outdata/marker_file.txt',sep = '\t',quote = F,row.names = F)

####cluster1学习网上的绘图方式——G-score####
library(maftools)

options(stringsAsFactors = F)
#读入GISTIC文件
pathway <- "./outdata/cluster1_cnv/"
all.lesions <- paste0(pathway,"all_lesions.conf_90.txt")
amp.genes <- paste0(pathway,"amp_genes.conf_90.txt")
del.genes <- paste0(pathway,"del_genes.conf_90.txt")
scores.gis <- paste0(pathway,"scores.gistic")

laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, 
                         gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)
##绘图
#genome plot
gisticChromPlot(gistic = laml.gistic, markBands="all", ref.build = 'hg38')

#人类的全基因组序列包BSgenome.Hsapiens.UCSC.hg38，这里面就有各个基因组的位置和长度信息。
#然后再继续学习下就知道BSgenome也是一个对象，可以通过特定函数提取信息
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), # 染色体名字
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)# 每条染色体长度
)

df$chromNum <- 1:length(df$chromName) # 染色体名字纯数字版，为了和scores.gistic里面的对应

# 我们用的原始CNV文件是没有性染色体信息的
df <- df[1:22,] # 只要前22条染色体信息

df

str(df)


#计算从0开始的每条染色体位置坐标，就是简单的线段长度加减法，不过对于我这种好久不搞数学的人来说也是很费脑子的！
#在scores.gistic这个文件里，第一条染色体位置是从0开始的，所以不用怎么改，但是第2条染色体的Start的坐标，
#应该是再加上第一条染色体的长度才是我们需要的，以此类推，不断相加！
#所以我们先计算下每条染色体从0开始的起始坐标是多少！第一条染色体起始位置就是0，第二条起始位置是第一条长度的位置，第3条是前两条长度的位置，以此类推！
# 小发现，在用cumsum前要把int变成numeric
df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) # 染色体累加长度

# 得到每条染色体从0开始的起始坐标
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])

# 计算每条染色体中间位置坐标，用来最后加文字
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle

df

# 如果你不知道用哪个函数读取，多试几次就知道了！
scores <- read.table("./outdata/cluster1_cnv/scores.gistic",sep="\t",header=T,stringsAsFactors = F)

head(scores)
##   Type Chromosome   Start     End X.log10.q.value.  G.score average.amplitude
## 1  Amp          1 3302046 3371973                0 0.021837          0.339205
## 2  Amp          1 3375059 3380822                0 0.021099          0.354653
## 3  Amp          1 3381074 3449929                0 0.020520          0.368893
## 4  Amp          1 3451503 3503571                0 0.022561          0.392242
## 5  Amp          1 3505022 4071958                0 0.022036          0.376773
## 6  Amp          1 4072066 4090117                0 0.022599          0.374114
##   frequency
## 1  0.036810
## 2  0.036810
## 3  0.034765
## 4  0.034765
## 5  0.032720
## 6  0.032720
#每一个G.score都对应一个坐标，这样才能画图，
#但其实每一个Amp或者Del是一个区间，为了方便，我们就用起始坐标代替了，当然你也可以用中点、结束点表示
chromID <- scores$Chromosome

scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
#得到的G.score都是正数，需要把Del的G.score变成负数。
range(scores$G.score)
## [1] 0.000000 0.564793

scores[scores$Type == "Del", "G.score"] <- scores[scores$Type == "Del", "G.score"] * -1

range(scores$G.score)
## [1] -0.564793  0.280607
library(ggsci)
# ggplot(scores, aes(StartPos, G.score))+
#   geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
#   scale_fill_lancet(guide=guide_legend(reverse = T))+
#   geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
#   geom_text(data = df,aes(x=chromMidelePosFrom0,y=0.2,label=chromName))+
#   scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9))+
#   ylim(-0.3,0.3)+
#   theme_minimal()

###增加错落感
df$ypos <- rep(c(-0.3,0.3),11)
ggplot(scores, aes(StartPos, G.score))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName),size=6,size = 6)+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-0.3,0.3)+
  theme_bw()+
  #ggplot2_theme_bw去掉网格线和背景色
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16)
  )
####cluster1绘图——Frequency####
#总的
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), # 染色体名字
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)# 每条染色体长度
)

df$chromNum <- 1:length(df$chromName) 

# 我们用的原始CNV文件是没有性染色体信息的
df <- df[1:22,] # 只要前22条染色体信息
df
str(df)

df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) # 染色体累加长度

# 得到每条染色体从0开始的起始坐标
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])

# 计算每条染色体中间位置坐标，用来最后加文字
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle

df

# 如果你不知道用哪个函数读取，多试几次就知道了！
scores <- read.table("./outdata/cluster1_cnv/scores.gistic",sep="\t",header=T,stringsAsFactors = F)
head(scores)

chromID <- scores$Chromosome

scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
#得到的G.score都是正数，需要把Del的G.score变成负数。
range(scores$frequency)
## [1] 0.000000 0.564793

scores[scores$Type == "Del", "frequency"] <- scores[scores$Type == "Del", "frequency"] * -1

range(scores$frequency)
scores$frequency<-scores$frequency*70
range(scores$frequency)
## [1] -0.564793  0.280607
###增加错落感
df$ypos <- rep(c(-52,51),11)
ggplot(scores, aes(StartPos, frequency))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName),size=6,size = 6)+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-52,51)+
  theme_bw()+
  #ggplot2_theme_bw去掉网格线和背景色
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16))

####cluster2学习网上的绘图方式——G-score####
library(maftools)

options(stringsAsFactors = F)
#读入GISTIC文件
pathway <- "./outdata/cluster2_cnv/"
all.lesions <- paste0(pathway,"all_lesions.conf_90.txt")
amp.genes <- paste0(pathway,"amp_genes.conf_90.txt")
del.genes <- paste0(pathway,"del_genes.conf_90.txt")
scores.gis <- paste0(pathway,"scores.gistic")

laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, 
                         gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)
##绘图
#genome plot
gisticChromPlot(gistic = laml.gistic, markBands="all", ref.build = 'hg38')

#人类的全基因组序列包BSgenome.Hsapiens.UCSC.hg38，这里面就有各个基因组的位置和长度信息。
#然后再继续学习下就知道BSgenome也是一个对象，可以通过特定函数提取信息
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), # 染色体名字
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)# 每条染色体长度
)

df$chromNum <- 1:length(df$chromName) # 染色体名字纯数字版，为了和scores.gistic里面的对应

# 我们用的原始CNV文件是没有性染色体信息的
df <- df[1:22,] # 只要前22条染色体信息

df

str(df)


#计算从0开始的每条染色体位置坐标，就是简单的线段长度加减法，不过对于我这种好久不搞数学的人来说也是很费脑子的！
#在scores.gistic这个文件里，第一条染色体位置是从0开始的，所以不用怎么改，但是第2条染色体的Start的坐标，
#应该是再加上第一条染色体的长度才是我们需要的，以此类推，不断相加！
#所以我们先计算下每条染色体从0开始的起始坐标是多少！第一条染色体起始位置就是0，第二条起始位置是第一条长度的位置，第3条是前两条长度的位置，以此类推！
# 小发现，在用cumsum前要把int变成numeric
df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) # 染色体累加长度

# 得到每条染色体从0开始的起始坐标
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])

# 计算每条染色体中间位置坐标，用来最后加文字
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle

df

# 如果你不知道用哪个函数读取，多试几次就知道了！
scores <- read.table("./outdata/cluster2_cnv/scores.gistic",sep="\t",header=T,stringsAsFactors = F)

head(scores)
##   Type Chromosome   Start     End X.log10.q.value.  G.score average.amplitude
## 1  Amp          1 3302046 3371973                0 0.021837          0.339205
## 2  Amp          1 3375059 3380822                0 0.021099          0.354653
## 3  Amp          1 3381074 3449929                0 0.020520          0.368893
## 4  Amp          1 3451503 3503571                0 0.022561          0.392242
## 5  Amp          1 3505022 4071958                0 0.022036          0.376773
## 6  Amp          1 4072066 4090117                0 0.022599          0.374114
##   frequency
## 1  0.036810
## 2  0.036810
## 3  0.034765
## 4  0.034765
## 5  0.032720
## 6  0.032720
#每一个G.score都对应一个坐标，这样才能画图，
#但其实每一个Amp或者Del是一个区间，为了方便，我们就用起始坐标代替了，当然你也可以用中点、结束点表示
chromID <- scores$Chromosome

scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
#得到的G.score都是正数，需要把Del的G.score变成负数。
range(scores$G.score)
## [1] 0.000000 0.564793

scores[scores$Type == "Del", "G.score"] <- scores[scores$Type == "Del", "G.score"] * -1

range(scores$G.score)
## [1] -0.564793  0.280607
library(ggsci)
# ggplot(scores, aes(StartPos, G.score))+
#   geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
#   scale_fill_lancet(guide=guide_legend(reverse = T))+
#   geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
#   geom_text(data = df,aes(x=chromMidelePosFrom0,y=0.2,label=chromName))+
#   scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9))+
#   ylim(-0.3,0.3)+
#   theme_minimal()

###增加错落感
df$ypos <- rep(c(-0.3,0.3),11)
ggplot(scores, aes(StartPos, G.score))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName),size=6,size=6)+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-0.3,0.3)+
  theme_bw()+
  #ggplot2_theme_bw去掉网格线和背景色
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16)
  )
####cluster2绘图——Frequency####
#总的
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), # 染色体名字
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)# 每条染色体长度
)

df$chromNum <- 1:length(df$chromName) 

# 我们用的原始CNV文件是没有性染色体信息的
df <- df[1:22,] # 只要前22条染色体信息
df
str(df)

df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) # 染色体累加长度

# 得到每条染色体从0开始的起始坐标
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])

# 计算每条染色体中间位置坐标，用来最后加文字
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle

df

# 如果你不知道用哪个函数读取，多试几次就知道了！
scores <- read.table("./outdata/cluster2_cnv/scores.gistic",sep="\t",header=T,stringsAsFactors = F)
head(scores)

chromID <- scores$Chromosome

scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
#得到的G.score都是正数，需要把Del的G.score变成负数。
range(scores$frequency)
## [1] 0.000000 0.564793

scores[scores$Type == "Del", "frequency"] <- scores[scores$Type == "Del", "frequency"] * -1

range(scores$frequency)
scores$frequency<-scores$frequency*70
range(scores$frequency)
## [1] -0.564793  0.280607
###增加错落感
df$ypos <- rep(c(-54.5,56.5),11)
ggplot(scores, aes(StartPos, frequency))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName),size=6)+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-54.5,56.5)+
  theme_bw()+
  #ggplot2_theme_bw去掉网格线和背景色
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16))



####cluster3学习网上的绘图方式——G-score####
library(maftools)

options(stringsAsFactors = F)
#读入GISTIC文件
pathway <- "./outdata/cluster3_cnv/"
all.lesions <- paste0(pathway,"all_lesions.conf_90.txt")
amp.genes <- paste0(pathway,"amp_genes.conf_90.txt")
del.genes <- paste0(pathway,"del_genes.conf_90.txt")
scores.gis <- paste0(pathway,"scores.gistic")

laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, 
                         gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)
##绘图
#genome plot
gisticChromPlot(gistic = laml.gistic, markBands="all", ref.build = 'hg38')

#人类的全基因组序列包BSgenome.Hsapiens.UCSC.hg38，这里面就有各个基因组的位置和长度信息。
#然后再继续学习下就知道BSgenome也是一个对象，可以通过特定函数提取信息
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), # 染色体名字
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)# 每条染色体长度
)

df$chromNum <- 1:length(df$chromName) # 染色体名字纯数字版，为了和scores.gistic里面的对应

# 我们用的原始CNV文件是没有性染色体信息的
df <- df[1:22,] # 只要前22条染色体信息

df

str(df)


#计算从0开始的每条染色体位置坐标，就是简单的线段长度加减法，不过对于我这种好久不搞数学的人来说也是很费脑子的！
#在scores.gistic这个文件里，第一条染色体位置是从0开始的，所以不用怎么改，但是第2条染色体的Start的坐标，
#应该是再加上第一条染色体的长度才是我们需要的，以此类推，不断相加！
#所以我们先计算下每条染色体从0开始的起始坐标是多少！第一条染色体起始位置就是0，第二条起始位置是第一条长度的位置，第3条是前两条长度的位置，以此类推！
# 小发现，在用cumsum前要把int变成numeric
df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) # 染色体累加长度

# 得到每条染色体从0开始的起始坐标
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])

# 计算每条染色体中间位置坐标，用来最后加文字
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle

df

# 如果你不知道用哪个函数读取，多试几次就知道了！
scores <- read.table("./outdata/cluster3_cnv/scores.gistic",sep="\t",header=T,stringsAsFactors = F)

head(scores)
##   Type Chromosome   Start     End X.log10.q.value.  G.score average.amplitude
## 1  Amp          1 3302046 3371973                0 0.021837          0.339205
## 2  Amp          1 3375059 3380822                0 0.021099          0.354653
## 3  Amp          1 3381074 3449929                0 0.020520          0.368893
## 4  Amp          1 3451503 3503571                0 0.022561          0.392242
## 5  Amp          1 3505022 4071958                0 0.022036          0.376773
## 6  Amp          1 4072066 4090117                0 0.022599          0.374114
##   frequency
## 1  0.036810
## 2  0.036810
## 3  0.034765
## 4  0.034765
## 5  0.032720
## 6  0.032720
#每一个G.score都对应一个坐标，这样才能画图，
#但其实每一个Amp或者Del是一个区间，为了方便，我们就用起始坐标代替了，当然你也可以用中点、结束点表示
chromID <- scores$Chromosome

scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
#得到的G.score都是正数，需要把Del的G.score变成负数。
range(scores$G.score)
## [1] 0.000000 0.564793

scores[scores$Type == "Del", "G.score"] <- scores[scores$Type == "Del", "G.score"] * -1

range(scores$G.score)
## [1] -0.564793  0.280607
library(ggsci)
# ggplot(scores, aes(StartPos, G.score))+
#   geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
#   scale_fill_lancet(guide=guide_legend(reverse = T))+
#   geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
#   geom_text(data = df,aes(x=chromMidelePosFrom0,y=0.2,label=chromName))+
#   scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9))+
#   ylim(-0.3,0.3)+
#   theme_minimal()

###增加错落感
df$ypos <- rep(c(-0.3,0.3),11)
ggplot(scores, aes(StartPos, G.score))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName),size=6)+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-0.3,0.3)+
  theme_bw()+
  #ggplot2_theme_bw去掉网格线和背景色
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16)
  )
####cluster3绘图——Frequency####
#总的
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), # 染色体名字
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)# 每条染色体长度
)

df$chromNum <- 1:length(df$chromName) 

# 我们用的原始CNV文件是没有性染色体信息的
df <- df[1:22,] # 只要前22条染色体信息
df
str(df)

df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) # 染色体累加长度

# 得到每条染色体从0开始的起始坐标
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])

# 计算每条染色体中间位置坐标，用来最后加文字
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle

df

# 如果你不知道用哪个函数读取，多试几次就知道了！
scores <- read.table("./outdata/cluster3_cnv/scores.gistic",sep="\t",header=T,stringsAsFactors = F)
head(scores)

chromID <- scores$Chromosome

scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
#得到的G.score都是正数，需要把Del的G.score变成负数。
range(scores$frequency)
## [1] 0.000000 0.564793

scores[scores$Type == "Del", "frequency"] <- scores[scores$Type == "Del", "frequency"] * -1

range(scores$frequency)
scores$frequency<-scores$frequency*70
range(scores$frequency)
## [1] -0.564793  0.280607
###增加错落感
df$ypos <- rep(c(-49,61),11)
ggplot(scores, aes(StartPos, frequency))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName),size=6)+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-49,61)+
  theme_bw()+
  #ggplot2_theme_bw去掉网格线和背景色
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16))


