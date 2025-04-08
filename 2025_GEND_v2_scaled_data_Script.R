

### R CMD BATCH /u/project/kp1/jlangerm/Projects/Gomperts_SMG/Analysis/2024-10-10_BG15_organoid_VoxOutput/2024.10.10_BG15_Vox_V8_Part1_Script.R
###Part 1 script

library(Seurat)
library(ggplot2)


Dir.Name <- "BG15_organoid_VoxOutput_part2"
dir.create(paste0("/u/project/kp1/jlangerm/Projects/Gomperts_SMG/Analysis/", Sys.Date(),"_",Dir.Name,"/"))
outdir <- c(paste0("/u/project/kp1/jlangerm/Projects/Gomperts_SMG/Analysis/", Sys.Date(),"_",Dir.Name,"/"))


Exp.Name <- "BG15_SMGOrgo"


BG15_Merge <- readRDS("/u/project/kp1/jlangerm/Projects/Gomperts_SMG/Analysis/2024-07-24_BG15_First_Look2/BG15_SMGOrgo.v2.object.obj")

markers<-read.table("/u/project/kp1/jlangerm/Projects/Gomperts_SMG/Analysis/2024-10-10_BG15_organoid_VoxOutput/BG15_SMGOrgo.Aggrecluster.metadata.txt")


normalized.data <- as.data.frame(BG15_Merge@assays$RNA$data)




####
## Start VOX workflow

library(dplyr)
###Have to rerun tag table now with new cluster names
tag.table <- NULL
tag.table$rn <- row.names(normalized.data)
for( Sample1 in sort(unique(markers$Cluster))  ){
  print(Sample1)
  
  if(length(row.names(markers[markers$Cluster %in% Sample1,]) )>1){
    
    tag.temp <-cbind.data.frame(rn=row.names(normalized.data),
                                Sample1=rowMeans(normalized.data[ ,colnames(normalized.data) %in% row.names(markers[markers$Cluster %in% Sample1,]) ]) )
    tag.table <- cbind.data.frame(tag.table, tag.temp$Sample1)
    names(tag.table)[ncol(tag.table)] <- paste(Sample1)} else{
      
      tag.temp <- normalized.data[ ,colnames(normalized.data) %in% row.names(markers[markers$Cluster %in% Sample1,]) ]
      tag.table$temp <- tag.temp
      names(tag.table)[ncol(tag.table)] <- paste(Sample1)
    }
}
row.names(tag.table) <- tag.table$rn
tag.table$rn <- NULL


### Allow for any merging and then finish here

new.sample.df<-as.data.frame(t(t(tag.table)/colSums(tag.table)))
new.sample.df$rn <- NULL









##########
##Builds variance tag


gene.info <- cbind.data.frame(Gene=row.names(new.sample.df), Mean=rowMeans(new.sample.df))
gene.info$Stdev <- apply(new.sample.df,1, FUN=sd)
gene.info$Variance <- gene.info$Stdev/gene.info$Mean
gene.info$Bin.Mean <- cut(gene.info$Mean,5, labels=F)

gene.info2 <- gene.info
tissue.names <-names(new.sample.df)
for(Tissue in tissue.names){
  gene.info2$input <- new.sample.df[,Tissue]
  names(gene.info2)[ncol(gene.info2)] <- Tissue
}

scaled.tag<-(gene.info2[,c(6:ncol(gene.info2))]-gene.info2$Mean)/gene.info2$Stdev
scaled.tag <- scaled.tag[!row.names(scaled.tag) %in% (row.names(scaled.tag[is.na(scaled.tag),])) ,]

# ggplot(scaled.tag)+geom_density(aes(C_23))

#########quin Binning
scaled.bin.tag <- scaled.tag
scaled.bin.tag[scaled.bin.tag>=2]<- 5
scaled.bin.tag[scaled.bin.tag>=1 & scaled.bin.tag<(2)] <- 4
scaled.bin.tag[scaled.bin.tag<1 & scaled.bin.tag>=(-1)] <- 3
scaled.bin.tag[scaled.bin.tag<(-1) & scaled.bin.tag>=(-2)] <- 2
scaled.bin.tag[scaled.bin.tag<(-2)] <- 1
#########







###Try to get above 10% of genes here or you will have a very large N1-noise Grouping
variance.tag <- gene.info[gene.info$Gene %in% row.names(scaled.bin.tag),]
variance.tag$Scaled_transcode <-apply( scaled.bin.tag , 1 , paste , collapse = "" )
transcode.table <-as.data.frame(table(variance.tag$Scaled_transcode))
nrow(transcode.table)
mean(transcode.table$Freq)
transcode.table<- transcode.table[with(transcode.table, order(Freq, decreasing=T)),]
head(transcode.table)


transcode.table$Trans_num<- paste0("T", c(1:nrow(transcode.table)))
variance.tag$Trans_num <- transcode.table[match(variance.tag$Scaled_transcode , transcode.table$Var1),]$Trans_num

###How many genes are included in freq>4
included.percent<-100*nrow(variance.tag[variance.tag$Trans_num %in% transcode.table[transcode.table$Freq>3,]$Trans_num,])/nrow(variance.tag)
print( paste0("includes ",  round(included.percent,2),"% of genes"))



variance.tag2 <- variance.tag


####use a softened dictionary and avoid the correlation issue
##

# variance.tag2 <-variance.tag
variance.tag2$Trans_num2 <- variance.tag2$Trans_num
n<-1


length(unique(variance.tag$Scaled_transcode))
length(unique(variance.tag$Trans_num))










############
### Corr Round 1
transcode.table2<-as.data.frame(table(variance.tag2$Trans_num2))

nrow(transcode.table2)-nrow(transcode.table)

trans.means <- cbind.data.frame(rn= names(new.sample.df))
##Recalculate transcode means 
for(Transcode in transcode.table2$Var1){
  
  print(Transcode)
  top.genes<- variance.tag2[variance.tag2$Trans_num2==Transcode,]$Gene
  # freq.amount <- transcode.table[transcode.table$Var1==Transcode,]$Freq
  trans.means$Gene <- colMeans(new.sample.df[top.genes,])
  # plot.data <- cbind.data.frame(Sample= names(scaled.Stam), Gene= colMeans(scaled.Stam[top.genes,]), group=1)
  names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
  # rank_n=rank_n+1
}
row.names(trans.means) <- trans.means$rn
trans.means$rn <- NULL
exp.level.tester <- as.data.frame(t(trans.means))


reasonably.expressed.transcodes<-row.names(exp.level.tester[ rowSums(exp.level.tester)>0.00003,])
nrow(variance.tag2[variance.tag2$Trans_num2 %in% reasonably.expressed.transcodes,])/nrow(variance.tag2)
nrow(variance.tag2[variance.tag2$Trans_num2 %in% transcode.table2[transcode.table2$Freq>9,]$Var1,])
##How many members are in low freq timecodes
as.data.frame(table(transcode.table2[!transcode.table2$Var1 %in% reasonably.expressed.transcodes,]$Freq))
table(transcode.table2$Freq)


library(reshape2)
###Precalc
tester.cor<-cor(trans.means)
diag(tester.cor)<-0
# 
tester.melt <- melt(tester.cor)
tester.melt<- tester.melt[tester.melt$value>0.5,]
tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]
# head(tester.melt)

nrow(tester.melt[tester.melt$value>0.90,])
nrow(tester.melt)
##Collapse based on gentle some cutoff

variance.tag2$Trans_num3 <- variance.tag2$Trans_num2
transcode.table2$Trans_num3 <- transcode.table2$Var1 



library(reshape2)
for(Transcode in unique(transcode.table2[transcode.table2$Freq<9,]$Var1)){
  print(paste0("low freq transcode ",Transcode))
  tester.melt2<-tester.melt[tester.melt$Var1==Transcode,]
  if(max(tester.melt2$value)>0.9){
    
    change.target  <- tester.melt2$Var2[1]
    print(paste0("changing to " , change.target))
    variance.tag2[variance.tag2$Trans_num2==Transcode,]$Trans_num3 <- as.character(change.target)
    transcode.table2[transcode.table2$Var1==Transcode,]$Trans_num3 <- as.character(change.target)
    
  } else {
    print("Nothing yet")
  }
}

variance.tag2$Transcode <- variance.tag2$Trans_num3
dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag2$Trans_num3)){
  print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    print("next")}
  dictionary.spread<-Transcode
  tester.df <- variance.tag2[,c("Trans_num3","Trans_num2")]
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num3 %in% dictionary.spread,]$Trans_num3)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num2 %in% dictionary.spread,]$Trans_num3)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num3 %in% dictionary.spread,]$Trans_num2)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num2 %in% dictionary.spread,]$Trans_num2)
    dictionary.spread<-unique(dictionary.spread)
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag2[variance.tag2$Trans_num3 %in% dictionary.spread,]$Transcode <- paste0("TB", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}
variance.tag2$Trans_num3<- variance.tag2$Transcode
print(paste0("Now have ", length(unique(variance.tag2$Trans_num3)), " transcodes compared to ", length(unique(variance.tag2$Trans_num2))) )









#### Corr round 2


transcode.table2<- as.data.frame(table(variance.tag2$Trans_num3))

trans.means <- cbind.data.frame(rn= names(new.sample.df))
##### Round 3 - correlation
Transcode <- "TA 10053"
##Recalculate transcode means 
for(Transcode in transcode.table2$Var1){
  
  print(Transcode)
  top.genes<- variance.tag2[variance.tag2$Trans_num3==Transcode,]$Gene
  trans.means$Gene <- colMeans(new.sample.df[top.genes,])
  names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
}
row.names(trans.means) <- trans.means$rn
trans.means$rn <- NULL
exp.level.tester <- as.data.frame(t(trans.means))

reasonably.expressed.transcodes<-row.names(exp.level.tester[ rowSums(exp.level.tester)>0.00003,])
nrow(variance.tag2[variance.tag2$Trans_num3 %in% reasonably.expressed.transcodes,])/nrow(variance.tag2)
##How many members are in low freq timecodes
as.data.frame(table(transcode.table2[!transcode.table2$Var1 %in% reasonably.expressed.transcodes,]$Freq))
table(transcode.table2$Freq)
nrow(variance.tag2[variance.tag2$Trans_num3 %in% transcode.table2[transcode.table2$Freq>9,]$Var1,])




###Precalc
tester.cor<-cor(trans.means)
diag(tester.cor)<-0
tester.melt <- melt(tester.cor)
tester.melt<- tester.melt[tester.melt$value>0.5,]
tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]



##Collapse based on gentle some cutoff
variance.tag2$Trans_num4 <- variance.tag2$Trans_num3
transcode.table2$Trans_num4 <- transcode.table2$Var1 

library(reshape2)
Transcode<- "TA 1105"
for(Transcode in unique(transcode.table2[transcode.table2$Freq<10,]$Var1)){
  print(paste0("low freq transcode ",Transcode))
  tester.melt2<-tester.melt[tester.melt$Var1==Transcode,]
  if(max(tester.melt2$value)>0.8){
    change.target  <- tester.melt2$Var2[1]
    print(paste0("changing to " , change.target))
    variance.tag2[variance.tag2$Trans_num3==Transcode,]$Trans_num4 <- as.character(change.target)
    transcode.table2[transcode.table2$Var1==Transcode,]$Trans_num4 <- as.character(change.target)
  } else {
    print("Nothing yet")
  }
}

variance.tag2$Transcode <- variance.tag2$Trans_num4
dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag2$Trans_num4)){
  print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    print("next")}
  dictionary.spread<-Transcode
  tester.df <- variance.tag2[,c("Trans_num4","Trans_num3")]
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num4 %in% dictionary.spread,]$Trans_num4)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num3 %in% dictionary.spread,]$Trans_num4)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num4 %in% dictionary.spread,]$Trans_num3)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num3 %in% dictionary.spread,]$Trans_num3)
    
    dictionary.spread<-unique(dictionary.spread)
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag2[variance.tag2$Trans_num4 %in% dictionary.spread,]$Transcode <- paste0("TC", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}
variance.tag2$Trans_num4<- variance.tag2$Transcode
print(paste0("Now have ", length(unique(variance.tag2$Trans_num4)), " transcodes"))







###############
#### Corr round 3
transcode.table2<- as.data.frame(table(variance.tag2$Trans_num4))
trans.means <- cbind.data.frame(rn= names(new.sample.df))
##Recalculate transcode means 
for(Transcode in transcode.table2$Var1){
  print(Transcode)
  top.genes<- variance.tag2[variance.tag2$Trans_num4==Transcode,]$Gene
  trans.means$Gene <- colMeans(new.sample.df[top.genes,])
  names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
}
row.names(trans.means) <- trans.means$rn
trans.means$rn <- NULL
exp.level.tester <- as.data.frame(t(trans.means))

reasonably.expressed.transcodes<-row.names(exp.level.tester[ rowSums(exp.level.tester)>0.00003,])
nrow(variance.tag2[variance.tag2$Trans_num4 %in% reasonably.expressed.transcodes,])/nrow(variance.tag2)
##How many members are in low freq timecodes
as.data.frame(table(transcode.table2[!transcode.table2$Var1 %in% reasonably.expressed.transcodes,]$Freq))
table(transcode.table2$Freq)
nrow(variance.tag2[variance.tag2$Trans_num4 %in% transcode.table2[transcode.table2$Freq>9,]$Var1,])

###Precalc
tester.cor<-cor(trans.means)
diag(tester.cor)<-0
tester.melt <- melt(tester.cor)
tester.melt<- tester.melt[tester.melt$value>0.5,]
tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]

##Collapse based on gentle some cutoff

variance.tag2$Trans_num5 <- variance.tag2$Trans_num4
transcode.table2$Trans_num5 <- transcode.table2$Var1 

library(reshape2)
Transcode<- "TA 1105"
for(Transcode in unique(transcode.table2[transcode.table2$Freq<9,]$Var1)){
  print(paste0("low freq transcode ",Transcode))
  tester.melt2<-tester.melt[tester.melt$Var1==Transcode,]
  if(max(tester.melt2$value)>0.7){
    change.target  <- tester.melt2$Var2[1]
    print(paste0("changing to " , change.target))
    variance.tag2[variance.tag2$Trans_num4==Transcode,]$Trans_num5 <- as.character(change.target)
    transcode.table2[transcode.table2$Var1==Transcode,]$Trans_num5 <- as.character(change.target)
  } else {
    print("Nothing yet")
  }
}

##Clean up
rm(tester.cor)
rm(tester.melt)
variance.tag2$Transcode <- variance.tag2$Trans_num5

dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag2$Trans_num5)){
  print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    print("next")}
  dictionary.spread<-Transcode
  tester.df <- variance.tag2[,c("Trans_num5","Trans_num4")]
  
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num5 %in% dictionary.spread,]$Trans_num5)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num4 %in% dictionary.spread,]$Trans_num5)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num5 %in% dictionary.spread,]$Trans_num4)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num4 %in% dictionary.spread,]$Trans_num4)
    dictionary.spread<-unique(dictionary.spread)
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag2[variance.tag2$Trans_num5 %in% dictionary.spread,]$Transcode <- paste0("TD", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}
variance.tag2$Trans_num5<- variance.tag2$Transcode
print(paste0("Now have ", length(unique(variance.tag2$Trans_num5)), " transcodes"))








##Do we do round 5 is we're at 2k??

###### Round 5  collapse


#### Corr round 5
transcode.table2<- as.data.frame(table(variance.tag2$Trans_num5))
trans.means <- cbind.data.frame(rn= names(new.sample.df))
##Recalculate transcode means 
for(Transcode in transcode.table2$Var1){
  print(Transcode)
  top.genes<- variance.tag2[variance.tag2$Trans_num5==Transcode,]$Gene
  trans.means$Gene <- colMeans(new.sample.df[top.genes,])
  names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
}
row.names(trans.means) <- trans.means$rn
trans.means$rn <- NULL
exp.level.tester <- as.data.frame(t(trans.means))

reasonably.expressed.transcodes<-row.names(exp.level.tester[ rowSums(exp.level.tester)>0.00003,])
nrow(variance.tag2[variance.tag2$Trans_num5 %in% reasonably.expressed.transcodes,])/nrow(variance.tag2)
##How many members are in low freq timecodes
as.data.frame(table(transcode.table2[!transcode.table2$Var1 %in% reasonably.expressed.transcodes,]$Freq))
table(transcode.table2$Freq)
nrow(variance.tag2[variance.tag2$Trans_num5 %in% transcode.table2[transcode.table2$Freq>9,]$Var1,])

###Precalc
tester.cor<-cor(trans.means)
diag(tester.cor)<-0
tester.melt <- melt(tester.cor)
tester.melt<- tester.melt[tester.melt$value>0.5,]
tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]

##Collapse based on gentle some cutoff
variance.tag2$Trans_num6 <- variance.tag2$Trans_num5
transcode.table2$Trans_num6 <- transcode.table2$Var1 

library(reshape2)
Transcode<- "TA 1105"
for(Transcode in unique(transcode.table2[transcode.table2$Freq<10,]$Var1)){
  print(paste0("low freq transcode ",Transcode))
  tester.melt2<-tester.melt[tester.melt$Var1==Transcode,]
  if(max(tester.melt2$value)>0.6){
    change.target  <- tester.melt2$Var2[1]
    print(paste0("changing to " , change.target))
    variance.tag2[variance.tag2$Trans_num5==Transcode,]$Trans_num6 <- as.character(change.target)
    transcode.table2[transcode.table2$Var1==Transcode,]$Trans_num6 <- as.character(change.target)
  } else {
    print("Nothing yet")
  }
}
rm(tester.cor)
rm(tester.melt)
variance.tag2$Transcode <- variance.tag2$Trans_num6

dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag2$Trans_num6)){
  print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    print("next")}
  
  dictionary.spread<-Transcode
  tester.df <- variance.tag2[,c("Trans_num6","Trans_num5")]
  
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num6 %in% dictionary.spread,]$Trans_num6)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num5 %in% dictionary.spread,]$Trans_num6)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num6 %in% dictionary.spread,]$Trans_num5)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num5 %in% dictionary.spread,]$Trans_num5)
    
    dictionary.spread<-unique(dictionary.spread)
    
    n=n+1
    
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag2[variance.tag2$Trans_num6 %in% dictionary.spread,]$Transcode <- paste0("TE", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}
variance.tag2$Trans_num6<- variance.tag2$Transcode
print(paste0("Now have ", length(unique(variance.tag2$Trans_num6)), " transcodes"))







# variance.tag2 <- read.table("/u/project/kp1/jlangerm/Projects/Gomperts_SMG/Analysis/2024-10-11_BG15_organoid_VoxOutput/BG15_SMGOrgo_VOX_v8.Gene.Table.to.CorrMergeFinal_HalfWay.txt",
#   sep="\t")


#####################
###### Noise collapse
###This has to be tweaked to be more useful 

variance.tag2$Transcode <- variance.tag2$Trans_num6
variance.tag2$Trans_num_final <- variance.tag2$Trans_num6

## Now recalculate transcode means
trans.means <- cbind.data.frame(rn= names(new.sample.df))
# trans.means$rn <- variance.test$Gene
# rank_n=1
for(Transcode in unique(variance.tag2$Transcode)){
  print(Transcode)
  top.genes<- variance.tag2[variance.tag2$Transcode==Transcode,]$Gene
  # freq.amount <- transcode.table[transcode.table$Var1==Transcode,]$Freq
  trans.means$Gene <- colMeans(new.sample.df[top.genes,])
  # plot.data <- cbind.data.frame(Sample= names(scaled.Stam), Gene= colMeans(scaled.Stam[top.genes,]), group=1)
  names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
  # rank_n=rank_n+1
}
row.names(trans.means) <- trans.means$rn

transcode.table2 <- as.data.frame(table(variance.tag2$Transcode))

library(reshape2)
trans.means$rn <- NULL
transcode.cor <- cor(trans.means)

transcode.cor[1:10,1:10]
diag(transcode.cor)<-0
tester.melt <- melt(transcode.cor)
tester.melt<-tester.melt[tester.melt$Var2 %in% transcode.table2[transcode.table2$Freq>=10,]$Var1,]
tester.melt<- tester.melt[tester.melt$value>0.0,]
tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]

tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]

nrow(transcode.table2[transcode.table2$Freq<10,])
# unique(tester.melt$Var2)
# unique(tester.melt$Var1)


# ggplot(tester.melt[tester.melt$Var1 %in% transcode.table2[transcode.table2$Freq<10,]$Var1,])+geom_bar(aes(Var1, value), stat="identity",fill="black")+
#   theme_classic()+
#   geom_hline(yintercept = mean(tester.melt$value), color="red")+
#   theme(axis.text.x=element_blank())
# ggsave(paste0(outdir, Exp.Name, ".KFinal.CorrTest.data.bargraph.png"),
#        height = 3, width=4)


# nrow(transcode.table2[transcode.table2$Freq>10,])

# Group <- "HG137"
variance.tag2$Transcode2 <- variance.tag2$Transcode

# variance.tag2[variance.tag2$Gene=="CTNNB1",]
# transcode.table2[transcode.table2$Var1=="TG32",]

final.cor.cutoff<- median(tester.melt[tester.melt$Var1 %in% transcode.table2[transcode.table2$Freq<10,]$Var1,]$value)
nrow(tester.melt[tester.melt$value>final.cor.cutoff,])
nrow(tester.melt[tester.melt$value>0.4,])
transcode.table2[transcode.table2$Var1 %in% tester.melt[tester.melt$value<0.4,]$Var1,]
final.cor.cutoff <- 0.40

# tester.melt[tester.melt$Var1=="TE3853",]
# transcode.table2[transcode.table2$Var1=="TE3853",]


###WAIT What happens to the noise group here where did it go 
noise.group <- paste0("TH",length(unique(variance.tag2$Transcode))+1  )
# variance.tag2$Transcode <- variance.tag2$HighGroup
for(Group in transcode.table2[transcode.table2$Freq<6,]$Var1){
  print(Group)
  if( tester.melt[tester.melt$Var1==Group,]$value> final.cor.cutoff ){
    
    change.target<- tester.melt[tester.melt$Var1==Group,]$Var2
    variance.tag2[variance.tag2$Transcode==Group,]$Transcode <- as.character(change.target)
    # transcode.table2[transcode.table2$Var1==Transcode,]$Trans_num6 <- as.character(change.target)
    
  } else {
    change.target<- noise.group
    variance.tag2[variance.tag2$Transcode==Group,]$Transcode <- as.character(change.target)
  }
}
noisy.genes <-variance.tag2[variance.tag2$Transcode== noise.group,]$Gene
# 
# unique(variance.tag2[variance.tag2$Transcode==noise.group,]$Trans_num6)

dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag2$Transcode)){
  print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    print("next")}
  ##Maybe I can build a growing dictionary of referenced states
  
  dictionary.spread<-Transcode
  tester.df <- variance.tag2[,c("Transcode","Trans_num_final")]
  
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Transcode %in% dictionary.spread,]$Transcode)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num_final %in% dictionary.spread,]$Transcode)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Transcode %in% dictionary.spread,]$Trans_num_final)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num_final %in% dictionary.spread,]$Trans_num_final)
    
    dictionary.spread<-unique(dictionary.spread)
    
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  
  
  variance.tag2[variance.tag2$Transcode %in% dictionary.spread,]$Transcode2 <- paste0("TI", tb.counter)
  
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}

tail(dictionary.summary)

table(variance.tag2$Transcode2)
table(variance.tag2$Transcode)


names(variance.tag2)
variance.tag2$TMerge_Final<- variance.tag2$Transcode2

variance.tag2$Trans_num_final <-NULL
# variance.tag2$KGroupFinal<-gsub("G","KF",variance.tag2$KGroupFinal)
# 
print(paste0("Now have ", length(unique(variance.tag2$TMerge_Final)), " transcodes"))


table(variance.tag2$TMerge_Final)












## Drop out low enriched vox groups to noise
new.sample.group.data <- cbind.data.frame(Sample=names(new.sample.df))
Group<- "TI77"
for(Group in sort(unique(variance.tag2$TMerge_Final))){
  print(Group)
  gene.list <- (variance.tag2[variance.tag2$TMerge_Final==Group,]$Gene)
  new.sample.group.data <- cbind.data.frame( new.sample.group.data ,Score=colSums(new.sample.df[gene.list,]) ) 
  names(new.sample.group.data)[ncol(new.sample.group.data)]<- paste0(Group)
}
new.sample.group.data$Sample <- NULL
new.sample.group.data <-((new.sample.group.data)/rowSums(new.sample.group.data))
rowSums(new.sample.group.data)

enrich.data <- NULL
enrich.data <- cbind.data.frame(Group=names(new.sample.group.data))
for(Sample in names(new.sample.df)){
  group.mean.tester <- new.sample.group.data
  group.mean.tester$Sample <- NULL
  diff.groups.tester<-(group.mean.tester[Sample,]+0.01)/(colMeans(group.mean.tester[row.names(group.mean.tester)!=Sample,])+0.01)
  diff.groups.tester<-unlist(diff.groups.tester)
  enrich.data<-cbind(enrich.data, Data=diff.groups.tester)
  names(enrich.data)[ncol(enrich.data)] <- Sample
  
}

row.names(enrich.data) <- enrich.data$Group
enrich.data$Group<-NULL
head(enrich.data)


ml3.table<-cbind.data.frame(Group=names(apply(enrich.data,1, max)),
                            Max=unlist(apply(enrich.data,1, max)),
                            Min=unlist(apply(enrich.data,1, min)))
ml3.table$Max_Up <- (ml3.table$Max-1)
ml3.table$Min_Up <- (1-ml3.table$Min)
ml3.table$AvgEnrich <- (ml3.table$Max_Up+ ml3.table$Min_Up)/2
ml3.table$Enrich_Bin <- cut(ml3.table$AvgEnrich, labels=F, breaks=20)


ml3.table <- ml3.table[with(ml3.table, order(-Enrich_Bin)),]


##Precalc for specificity
tmp.norm <- normalized.data[,colnames(normalized.data) %in% row.names(markers)]
tmp.norm[tmp.norm>0]<-1
###THIS CAN TAKE AWHILE FOR BIG DFs
tmp.norm[1:10,1:10]

variance.tag2$TMerge_Specificity <- 0
variance.tag2$TMerge_Correlation <- 0
# Group<- "TI129"

for(Group in unique(variance.tag2$TMerge_Final)){
  print(Group)
  genes.of.interest<- variance.tag2[variance.tag2$TMerge_Final==Group,]$Gene
  norm.data.chunk <- normalized.data[row.names(normalized.data) %in% genes.of.interest,]
  
  score.piece<-colMeans(norm.data.chunk)
  
  score.piece[score.piece>0]
  score.cor <- cor(t(norm.data.chunk), score.piece)
  
  ##Specificity calc
  # particular.enriched.groups <- names(enrich.data)[enrich.data[Group,] >1.1 ]
  
  # cells.to.take<- row.names(markers[markers$Cluster %in% particular.enriched.groups,])
  # subtmp.norm <- tmp.norm[,colnames(tmp.norm) %in% cells.to.take]
  # subtmp.norm<- subtmp.norm[row.names(subtmp.norm) %in% genes.of.interest,]
  # gene.tester <- rowMeans(subtmp.norm)
  # gene.tester<- cbind.data.frame(Gene=names(gene.tester), Score_in=gene.tester)
  
  # subtmp.norm <- tmp.norm[,!colnames(tmp.norm) %in% cells.to.take]
  # subtmp.norm<- subtmp.norm[row.names(subtmp.norm) %in% genes.of.interest,]
  # gene.tester2 <- rowMeans(subtmp.norm)
  # gene.tester2<- cbind.data.frame(Gene=names(gene.tester2), Score_out=gene.tester2)
  
  # final.gene.tester<-cbind(gene.tester, gene.tester2)
  # final.gene.tester$Specificity <- (final.gene.tester$Score_in +0.01)/(final.gene.tester$Score_out+0.01)
  
  final.gene.tester<- cbind.data.frame(Gene=row.names(score.cor),Specificity=0, Score_Corr =as.numeric(score.cor))
  
  # final.gene.tester$Score_Corr <- as.numeric(score.cor)
  # names(final.gene.tester)[ncol(final.gene.tester)]<-"Score_Corr"
  
  variance.tag2[variance.tag2$TMerge_Final==Group,]$TMerge_Specificity <- final.gene.tester[match(variance.tag2[variance.tag2$TMerge_Final==Group,]$Gene,final.gene.tester$Gene),]$Specificity
  variance.tag2[variance.tag2$TMerge_Final==Group,]$TMerge_Correlation <- final.gene.tester[match(variance.tag2[variance.tag2$TMerge_Final==Group,]$Gene,final.gene.tester$Gene),]$Score_Corr
  
}

# variance.tag2[is.na(variance.tag2)]<-0

##Plot to reflect these metalayer data

library(dplyr)

variance.tag2 %>% group_by(TMerge_Final) %>% summarize(Corr=mean(TMerge_Correlation), Specificity=mean(TMerge_Specificity)) -> metalayer.table


metalayer.table$Enrichment <- ml3.table[match(metalayer.table$TMerge_Final,ml3.table$Group),]$AvgEnrich

metalayer.table$CorrEnr_Score <- metalayer.table$Corr * metalayer.table$Enrichment

variance.tag2$TMerge_Enrichment <- metalayer.table[match(variance.tag2$TMerge_Final,metalayer.table$TMerge_Final),]$Enrichment
variance.tag2$TMerge_CorrEnrScore <- metalayer.table[match(variance.tag2$TMerge_Final,metalayer.table$TMerge_Final),]$CorrEnr_Score





###### Noise collapse
tmp.table <-metalayer.table[with(metalayer.table, order(CorrEnr_Score)),]
groups.to.break<-tmp.table[tmp.table$CorrEnr_Score<(mean(tmp.table$CorrEnr_Score)/5),]$TMerge_Final

as.data.frame(tmp.table[tmp.table$CorrEnr_Score<(mean(tmp.table$CorrEnr_Score)/5),])


######BREAK POINT
##Just letting these groups go out into the next point 
# genes.to.test <-c(noisy.genes)
##########The noise is too crappy and causes problems

###Method where I went ahead and broke them
genes.to.test <-c(noisy.genes, variance.tag2[variance.tag2$TMerge_Final %in% groups.to.break,]$Gene)

variance.tag2$Transcode <- variance.tag2$TMerge_Final
filt.variance.tag <- variance.tag2[!variance.tag2$Gene %in% genes.to.test,]

## Now recalculate transcode means
trans.means <- cbind.data.frame(rn= names(new.sample.df))
# trans.means$rn <- variance.test$Gene
# rank_n=1
for(Transcode in unique(filt.variance.tag$Transcode)){
  print(Transcode)
  top.genes<- filt.variance.tag[filt.variance.tag$Transcode==Transcode,]$Gene
  # freq.amount <- transcode.table[transcode.table$Var1==Transcode,]$Freq
  trans.means$Gene <- colMeans(new.sample.df[top.genes,])
  # plot.data <- cbind.data.frame(Sample= names(scaled.Stam), Gene= colMeans(scaled.Stam[top.genes,]), group=1)
  names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
  # rank_n=rank_n+1
}
row.names(trans.means) <- trans.means$rn
trans.means$rn <-NULL

transcode.table2 <- as.data.frame(table(variance.tag2$Transcode))

trans.means[1:5,1:5]
identical(colnames(tag.table), row.names(trans.means))
tester.join <- cbind.data.frame(trans.means,t(tag.table[row.names(tag.table) %in% genes.to.test,]))

gene.cor <- cor(tester.join)
diag(gene.cor)<-0

mean(apply(gene.cor[rownames(gene.cor) %in% genes.to.test,!colnames(gene.cor) %in% genes.to.test],1,max))

#############
n<-1; m<-1
for(Gene in genes.to.test){
  tmp.gene.cor <- gene.cor[Gene,]
  tmp.gene.cor<- tmp.gene.cor[!colnames(gene.cor) %in% genes.to.test]
  if(tmp.gene.cor[which.max(tmp.gene.cor)]>0.5){
    change.target <- names(tmp.gene.cor)[which.max(tmp.gene.cor)]
    
    variance.tag2[variance.tag2$Gene==Gene,]$Transcode <- change.target
    # print(paste0("changing group for ", Gene))
    n<-n+1
  } else {
    variance.tag2[variance.tag2$Gene==Gene,]$Transcode <- "Noise"
    # print(paste0("noise for group ", Gene))
    m <-m+1
  } 
} 

print(paste0(n , " merged and ", m, " noise calls"))

new.noisy.genes<-variance.tag2[variance.tag2$Transcode=="Noise",]$Gene

table(variance.tag2$Transcode)



###Have lost macrophage group here to noise




library(pheatmap)
ranged.means <- t(t(trans.means)/apply(trans.means, 2, max))
ordering.groups.df <- pheatmap(ranged.means, silent=T)

pheatmap(ranged.means,
         filename=paste0(outdir, Exp.Name, ".PostCorr.MeanExpGroups.heatmap.png"),
         height=8, width=10)

# ordering.groups.df$tree_row$labels[ordering.groups.df$tree_row$order]

sg.table<- cbind.data.frame(Group=ordering.groups.df$tree_col$labels[ordering.groups.df$tree_col$order], Order=c(1:ncol(ranged.means)))

# sg.table<- cbind.data.frame(Group=ordering.groups.df$tree_col$labels, Order=ordering.groups.df$tree_col$order)
sg.table<- sg.table[with(sg.table, order(Order)),]
sg.table$Group_New <- paste0("F" ,c(1:nrow(sg.table)) )

##Separate low noise



low.noise.genes<-new.noisy.genes[variance.tag2[variance.tag2$Gene %in% new.noisy.genes,]$Mean<(mean(variance.tag2$Mean)/10)]
variance.tag2[variance.tag2$Gene %in% low.noise.genes,]$Transcode <- "NK1"


# 
# ## Trajecory Score Shenanigans
# plot.data <- markers
# pc.data <- as.data.frame(BG12_Merge@reductions$pca@cell.embeddings)
# # row.names(pc.data) <- row.names(BG12_Merge)
# # $PC_1
# plot.data$PC1 <- pc.data[colnames(BG12_Merge) %in% row.names(plot.data),]$PC_1
# plot.data$PC2 <- pc.data[colnames(BG12_Merge) %in% row.names(plot.data),]$PC_2
# plot.data$PC3 <- pc.data[colnames(BG12_Merge) %in% row.names(plot.data),]$PC_3
# plot.data$PC4 <- pc.data[colnames(BG12_Merge) %in% row.names(plot.data),]$PC_4
# 
# pca.traj.df <- plot.data[,c("PC1","PC2","PC3","PC4")]
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# pca.traj.df <- data.frame(apply(pca.traj.df, 2, range01))
# pca.traj.df$Traj1 <-pca.traj.df$PC1 - pca.traj.df$PC2
# pca.traj.df$Traj2 <-pca.traj.df$PC3 - pca.traj.df$PC4
# pca.traj.df$TrajectoryScore <- pca.traj.df$Traj1+pca.traj.df$Traj2/2
# plot.data$TrajectoryScore <- pca.traj.df$TrajectoryScore
# 
# plot1 <- ggplot(plot.data)+
#   theme_classic()+
#   geom_point(aes(X,Y,color=TrajectoryScore),size=0.5)+
#   scale_color_gradientn(colors=c("blue","lightblue","grey55","orange","red"))
# 
# ggsave(plot=plot1, paste0(outdir, Exp.Name, ".PC_Trajectory.UMAP.png"),
#        height=4, width=5)
# 
# markers$TrajectoryScore <- pca.traj.df$TrajectoryScore
# 
# 
# library(dplyr)
# markers %>% group_by(Cluster) %>% summarize(Traj=mean(TrajectoryScore)) -> cluster.trajectory.levels
# 
# traj.gene.score<-colSums(t(tag.table[row.names(tag.table) %in% new.noisy.genes[!new.noisy.genes %in% low.noise.genes],]) *cluster.trajectory.levels$Traj)
# traj.gene.df<-cbind.data.frame(Gene=names(traj.gene.score), Score=round(traj.gene.score,4))
# ##DYNLL2 NPTX2
# 
# traj.gene.df$TrajBin <- ntile(traj.gene.df$Score,3)
# table(traj.gene.df$TrajBin)
# ##aside to test if these look okay below
# traj.gene.df$TrajBin <- paste0("NK",traj.gene.df$TrajBin+1)
# 
# length(traj.gene.df$TrajBin)
# identical(variance.tag2[variance.tag2$Gene %in% new.noisy.genes,]$Gene, traj.gene.df$Gene)

variance.tag2[variance.tag2$Gene %in% new.noisy.genes[!new.noisy.genes %in% low.noise.genes],]$Transcode <-"NK2"

sg.noise.table<-cbind.data.frame(Group=c("NK1","NK2"), Order=c(0,0),
                                 Group_New=c("NF1", "NF2") )

sg.table<-rbind(sg.table, sg.noise.table)

# sg.table$Group_New<- factor(sg.table$Group_New, levels=sg.table$Group_New)

variance.tag2$Vox_TMerge_Final <- sg.table[match(variance.tag2$Transcode,sg.table$Group),]$Group_New

print(paste0("Now have ", length(unique(variance.tag2$Vox_TMerge_Final)), " transcodes"))


table(variance.tag2$Vox_TMerge_Final)



options(scipen = 999)

write.table(variance.tag2, paste0(outdir, Exp.Name,"_VOX_v8.Gene.Table.to.CorrMergeFinal_HalfWay.txt"),
            sep="\t", quote=F)




# variance.tag2<- read.table("/u/project/kp1/jlangerm/Projects/Gomperts_IPF/Chandani/Analysis/2024-10-02_BGC4_Vox.v8/BGC4_VOX_v8.Gene.Table.to.CorrMergeFinal_HalfWay.txt",sep="\t")







################################
###Plot group expression on UMAP


dir.create(paste0(outdir, "Vox_TMerge_Final.Expression/"))
outdir2<- paste0(outdir, "Vox_TMerge_Final.Expression/")

# variance.tag$ML1_Index <- paste0(variance.tag$MetaLayer2,"-",variance.tag$MetaLayer1)

transcode.table <- as.data.frame(table(variance.tag2$Vox_TMerge_Final))
# n=1
# Exp.Name<- "YS3.Filt"
library(ggplot2)
for(Group in unique(variance.tag2$Vox_TMerge_Final)){
  print(Group)
  plot.data <- markers
  genes.of.interest <- variance.tag[variance.tag2$Vox_TMerge_Final==Group,]$Gene
  plot.data$Score <- colMeans(normalized.data[row.names(normalized.data) %in% genes.of.interest,colnames(normalized.data) %in% row.names(plot.data)])
  
  amount.of.genes<-transcode.table[transcode.table$Var1==Group,]$Freq
  
  plot.data <- plot.data[with(plot.data, order(Score)),]
  ggplot(plot.data)+
    # facet_wrap(~Type)+
    geom_point(data=markers[,c("X","Y")], aes(X,Y), color="grey55")+
    geom_point( aes(X,Y, color=Score) )+
    theme_classic()+
    labs(color=Group)+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ggtitle(paste0("VOX ",Group, " with ", amount.of.genes," genes"))+
    scale_colour_gradientn(colours=c("lightcyan3", "gold", "orange","orangered","purple"))
  ggsave(filename=paste0(outdir2,  Exp.Name, ".Vox.avgExp.TMerge_Final_",Group,".UMAP.png"),
         height=4, width=6)
}









#############Revisit Transcodes once for merge

variance.tag <- variance.tag2[!variance.tag2$Vox_TMerge_Final %in% c("NF1", "NF2"),]


###Make this even more general 
variance.tag$Transcode <- variance.tag$Vox_TMerge_Final


## Now recalculate transcode means
trans.means <- cbind.data.frame(rn= names(new.sample.df))
rank_n=1
for(Transcode in unique(variance.tag$Transcode)){
  print(Transcode)
  top.genes<- variance.tag[variance.tag$Transcode %in% Transcode,]$Gene
  trans.means$Gene <- colMeans(new.sample.df[top.genes,])
  names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
}
row.names(trans.means) <- trans.means$rn
trans.means$rn <- NULL

trans.means<-t(trans.means)



group.gene.info <- cbind.data.frame(Group=row.names(trans.means), Mean=rowMeans(trans.means))
group.gene.info$Stdev <- apply(trans.means,1, FUN=sd)


group.scaled.tag<-(trans.means-group.gene.info$Mean)/group.gene.info$Stdev
#########quin Binning
scaled.bin.tag <- group.scaled.tag
scaled.bin.tag[scaled.bin.tag>=2]<- 5
scaled.bin.tag[scaled.bin.tag>=1 & scaled.bin.tag<(2)] <- 4
scaled.bin.tag[scaled.bin.tag<1 & scaled.bin.tag>=(-1)] <- 3
scaled.bin.tag[scaled.bin.tag<(-1) & scaled.bin.tag>=(-2)] <- 2
scaled.bin.tag[scaled.bin.tag<(-2)] <- 1
#########


transcode.table<-as.data.frame(apply( scaled.bin.tag , 1 , paste , collapse = "" ) ) 
transcode.table$Former_Transcode <- row.names(transcode.table)
colnames(transcode.table)[1] <- "NewCode"
transcode.table2 <- as.data.frame(table(transcode.table$NewCode))
transcode.table2$NewName <- paste0("GA", c(1:nrow(transcode.table2)))
transcode.table$NewName <- transcode.table2[match(transcode.table$NewCode,transcode.table2$Var1),]$NewName 

variance.tag$TGroup_1 <-transcode.table[match(variance.tag$Transcode,transcode.table$Former_Transcode),]$NewName

print(paste0("Now have ", length(unique(variance.tag$TGroup_1))," vs ",length(unique(variance.tag2$Vox_TMerge_Final)), " transcodes"))





variance.tag[!complete.cases(variance.tag),]














########### tsne attack r1
### Clustering vox groups in multidimensional reduction space Round 1


variance.tag$Transcode <- variance.tag$TGroup_1
variance.tag$Trans_num_final <- variance.tag$TGroup_1

##Memory flag for tsne
if(ncol(normalized.data)>100000){
  ###memory issue catch route
  markers2<-markers
  markers2$Cluster <- as.factor(kmeans(markers2[,c("X","Y")], iter.max=30, nstart=5, centers=round(nrow(markers)/5))$cluster) 
  
  
  mem.tag.table<-NULL
  library(dplyr)
  ###Have to rerun tag table now with new cluster names
  mem.tag.table <- NULL
  mem.tag.table$rn <- row.names(normalized.data)
  for( Sample1 in sort(unique(markers2$Cluster))  ){
    print(Sample1)
    
    if(length(row.names(markers2[markers2$Cluster %in% Sample1,]) )>1){
      
      tag.temp <-cbind.data.frame(rn=row.names(normalized.data),
                                  Sample1=rowMeans(normalized.data[ ,colnames(normalized.data) %in% row.names(markers2[markers2$Cluster %in% Sample1,]) ]) )
      mem.tag.table <- cbind.data.frame(mem.tag.table, tag.temp$Sample1)
      names(mem.tag.table)[ncol(mem.tag.table)] <- paste(Sample1)} else{
        
        tag.temp <- normalized.data[ ,colnames(normalized.data) %in% row.names(markers2[markers2$Cluster %in% Sample1,]) ]
        mem.tag.table$temp <- tag.temp
        names(mem.tag.table)[ncol(mem.tag.table)] <- paste(Sample1)
      }
  }
  row.names(mem.tag.table) <- mem.tag.table$rn
  mem.tag.table$rn <- NULL
  
  direct.means <- cbind.data.frame(rn= colnames(mem.tag.table))
  rank_n=1
  for(Transcode in unique(variance.tag$Transcode)){
    print(Transcode)
    top.genes<- variance.tag[variance.tag$Transcode==Transcode,]$Gene
    direct.means$Gene <- colMeans(mem.tag.table[row.names(mem.tag.table) %in% top.genes,])
    names(direct.means)[ncol(direct.means)] <- paste0(Transcode)
  }
  row.names(direct.means) <- direct.means$rn
  direct.means$rn <- NULL
  
  
}else{

direct.means <- cbind.data.frame(rn= row.names(markers))
rank_n=1
for(Transcode in unique(variance.tag$Transcode)){
  print(Transcode)
  top.genes<- variance.tag[variance.tag$Transcode==Transcode,]$Gene
  direct.means$Gene <- colMeans(normalized.data[row.names(normalized.data) %in% top.genes,])
  names(direct.means)[ncol(direct.means)] <- paste0(Transcode)
}
row.names(direct.means) <- direct.means$rn
direct.means$rn <- NULL
}

direct.means.norm <- as.data.frame(t(t(direct.means)/colSums(direct.means)))


##lowering perplexity did not affect RAM requirement, must be irbla size
# output.tnse <- tsne(t(direct.means.norm), perplexity=20)

##rounding this df did not affect size of object, seems the core problem is 100k cell vs cell  
# direct.means.norm[1:10,1:10]
# min(direct.means.norm[direct.means.norm>0])
# 1.706533e-07
# 0.0000001706553
# direct.means.norm2<-round(direct.means.norm, digits = 9)
# direct.means.norm2<-signif(direct.means.norm, digits = 3)
# min(direct.means.norm2[direct.means.norm2>0])
# object.size(direct.means.norm);object.size(direct.means.norm2)




library(tsne)
output.tnse <- tsne(t(direct.means.norm), perplexity=50)
output.tsne <- as.data.frame(output.tnse)
output.tsne$Group <- names(direct.means.norm)


diff.vector<-abs( abs(output.tsne[1,]$V1 - output.tsne$V1)+abs(output.tsne[1,]$V2 - output.tsne$V2) )
diff.vector[diff.vector==0] <-  1

collected.mindist <- NULL
for(n in c(1:nrow(output.tsne))){
  diff.vector<-abs( abs(output.tsne[n,]$V1 - output.tsne$V1)+abs(output.tsne[n,]$V2 - output.tsne$V2) )
  
  diff.vector[diff.vector==0] <-  1
  temp.diff <- cbind.data.frame(Group=output.tsne[n,]$Group,MinDiff=min(diff.vector), CloseGroup=output.tsne[which.min(diff.vector),]$Group[1])
  # which.min(abs( (output.tsne[1,]$V1 - output.tsne$V1)+(output.tsne[1,]$V2 - output.tsne$V2) ))
  collected.mindist<- rbind(collected.mindist, temp.diff)
}

cutoff.for.tsne <- (mean(collected.mindist$MinDiff)-sd(collected.mindist$MinDiff)/2)
pos.close.groups<-collected.mindist[collected.mindist$MinDiff<cutoff.for.tsne,]$Group

plot1 <- ggplot(output.tsne)+
  geom_point(aes(V1,V2))+
  geom_point(data=output.tsne[output.tsne$Group %in% pos.close.groups,], aes(V1,V2), color="red")+
  geom_text(aes(V1,V2,label=Group), hjust=0.25)+
  theme_classic()
ggsave(plot=plot1, paste0(outdir, Exp.Name, ".VoxTsne.CollapseRound1_withCloseGroups.png"),
       height=8, width=10)

for(Group in collected.mindist[collected.mindist$MinDiff<cutoff.for.tsne,]$Group){
  variance.tag[variance.tag$Transcode==Group,]$Transcode<- collected.mindist[collected.mindist$Group==Group,]$CloseGroup
}

###Now just need dictionary attack
dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag$Transcode)){
  print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    print("next")}
  ##Maybe I can build a growing dictionary of referenced states
  dictionary.spread<-Transcode
  tester.df <- variance.tag[,c("Transcode","Trans_num_final")]
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Transcode %in% dictionary.spread,]$Transcode)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num_final %in% dictionary.spread,]$Transcode)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Transcode %in% dictionary.spread,]$Trans_num_final)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num_final %in% dictionary.spread,]$Trans_num_final)
    
    dictionary.spread<-unique(dictionary.spread)
    
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag[variance.tag$Transcode %in% dictionary.spread,]$Transcode2 <- paste0("DA", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}
length(unique(variance.tag$Transcode2));length(unique(variance.tag$Trans_num_final))
variance.tag$DRGroup_1 <- variance.tag$Transcode2











################
### TSNE grouping Round 2
length(unique(variance.tag$Transcode2));length(unique(variance.tag$Trans_num_final))
##Needs catch
if( (length(unique(variance.tag$Trans_num_final))-length(unique(variance.tag$Transcode2)))<5 ){
  print("finishing")
} else{

variance.tag$Transcode <- variance.tag$DRGroup_1
variance.tag$Trans_num_final <- variance.tag$DRGroup_1


##Memory flag for tsne
if(ncol(normalized.data)>100000){
  
  direct.means <- cbind.data.frame(rn= colnames(mem.tag.table))
  rank_n=1
  for(Transcode in unique(variance.tag$Transcode)){
    print(Transcode)
    top.genes<- variance.tag[variance.tag$Transcode==Transcode,]$Gene
    direct.means$Gene <- colMeans(mem.tag.table[row.names(mem.tag.table) %in% top.genes,])
    names(direct.means)[ncol(direct.means)] <- paste0(Transcode)
  }
  row.names(direct.means) <- direct.means$rn
  direct.means$rn <- NULL
  
  
}else{
  
direct.means <- cbind.data.frame(rn= row.names(markers))
rank_n=1
for(Transcode in unique(variance.tag$Transcode)){
  print(Transcode)
  top.genes<- variance.tag[variance.tag$Transcode==Transcode,]$Gene
  direct.means$Gene <- colMeans(normalized.data[row.names(normalized.data) %in% top.genes,])
  names(direct.means)[ncol(direct.means)] <- paste0(Transcode)
}
row.names(direct.means) <- direct.means$rn
direct.means$rn <- NULL
}

direct.means.norm <- as.data.frame(t(t(direct.means)/colSums(direct.means)))


library(tsne)
output.tnse <- tsne(t(direct.means.norm), perplexity=50)
##Does work just takes awhile, maybe 20 min for 20k+ cells

output.tsne <- as.data.frame(output.tnse)
output.tsne$Group <- names(direct.means.norm)

###Calculate min distances

diff.vector<-abs( abs(output.tsne[1,]$V1 - output.tsne$V1)+abs(output.tsne[1,]$V2 - output.tsne$V2) )
diff.vector[diff.vector==0] <-  1

collected.mindist <- NULL
for(n in c(1:nrow(output.tsne))){
  diff.vector<-abs( abs(output.tsne[n,]$V1 - output.tsne$V1)+abs(output.tsne[n,]$V2 - output.tsne$V2) )
  
  diff.vector[diff.vector==0] <-  1
  temp.diff <- cbind.data.frame(Group=output.tsne[n,]$Group,MinDiff=min(diff.vector), CloseGroup=output.tsne[which.min(diff.vector),]$Group[1])
  collected.mindist<- rbind(collected.mindist, temp.diff)
}

cutoff.for.tsne <- (mean(collected.mindist$MinDiff)-sd(collected.mindist$MinDiff)/2)
pos.close.groups<-collected.mindist[collected.mindist$MinDiff<cutoff.for.tsne,]$Group


plot1 <- ggplot(output.tsne)+
  geom_point(aes(V1,V2))+
  geom_point(data=output.tsne[output.tsne$Group %in% pos.close.groups,], aes(V1,V2), color="red")+
  geom_text(aes(V1,V2,label=Group), hjust=0.25)+
  theme_classic()
ggsave(plot=plot1, paste0(outdir, Exp.Name, ".VoxTsne.Round2.CloseGroups.png"),
       height=8, width=10)

##Sanity check
# unique(variance.tag[variance.tag$Transcode %in% c("DA2", "DA26","DA50","DA57"),]$CGroup_1)

for(Group in collected.mindist[collected.mindist$MinDiff<cutoff.for.tsne,]$Group){
  variance.tag[variance.tag$Transcode==Group,]$Transcode<- collected.mindist[collected.mindist$Group==Group,]$CloseGroup
}

###Now just need dictionary attack

dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag$Transcode)){
  print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    print("next")}
  
  dictionary.spread<-Transcode
  tester.df <- variance.tag[,c("Transcode","Trans_num_final")]
  
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Transcode %in% dictionary.spread,]$Transcode)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num_final %in% dictionary.spread,]$Transcode)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Transcode %in% dictionary.spread,]$Trans_num_final)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num_final %in% dictionary.spread,]$Trans_num_final)
    
    dictionary.spread<-unique(dictionary.spread)
    
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag[variance.tag$Transcode %in% dictionary.spread,]$Transcode2 <- paste0("DB", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}

length(unique(variance.tag$Transcode2));length(unique(variance.tag$Trans_num_final))
## 120 to 94 wo
variance.tag$DRGroup_2 <- variance.tag$Transcode2

}
###or iterate








################
### TSNE grouping Round 3

##Needs catch
if( (length(unique(variance.tag$Trans_num_final))-length(unique(variance.tag$Transcode2)))<5 ){
  print("finishing")
} else{
  
  
variance.tag$Transcode <- variance.tag$DRGroup_2
variance.tag$Trans_num_final <- variance.tag$DRGroup_2
##Memory flag for tsne
if(ncol(normalized.data)>100000){
  
  direct.means <- cbind.data.frame(rn= colnames(mem.tag.table))
  rank_n=1
  for(Transcode in unique(variance.tag$Transcode)){
    print(Transcode)
    top.genes<- variance.tag[variance.tag$Transcode==Transcode,]$Gene
    direct.means$Gene <- colMeans(mem.tag.table[row.names(mem.tag.table) %in% top.genes,])
    names(direct.means)[ncol(direct.means)] <- paste0(Transcode)
  }
  row.names(direct.means) <- direct.means$rn
  direct.means$rn <- NULL
  
  
}else{
direct.means <- cbind.data.frame(rn= row.names(markers))
rank_n=1
for(Transcode in unique(variance.tag$Transcode)){
  print(Transcode)
  top.genes<- variance.tag[variance.tag$Transcode==Transcode,]$Gene
  direct.means$Gene <- colMeans(normalized.data[row.names(normalized.data) %in% top.genes,])
  names(direct.means)[ncol(direct.means)] <- paste0(Transcode)
}
row.names(direct.means) <- direct.means$rn
direct.means$rn <- NULL
}
direct.means.norm <- as.data.frame(t(t(direct.means)/colSums(direct.means)))


library(tsne)
output.tnse <- tsne(t(direct.means.norm), perplexity=50)
##Does work just takes awhile, maybe 20 min for 20k+ cells

output.tsne <- as.data.frame(output.tnse)
output.tsne$Group <- names(direct.means.norm)


# plot1 <- ggplot(output.tsne)+
#   geom_point(aes(V1,V2))+
#   geom_text(aes(V1,V2,label=Group), hjust=0.25)+
#   theme_classic()
# ggsave(plot=plot1, paste0(outdir, Exp.Name, ".VoxTsne.Test.png"),
#        height=8, width=10)


#### Okay this actually looks pretty good


###Calculate min distances


diff.vector<-abs( abs(output.tsne[1,]$V1 - output.tsne$V1)+abs(output.tsne[1,]$V2 - output.tsne$V2) )

diff.vector[diff.vector==0] <-  1


# plot1 <- ggplot(output.tsne[output.tsne$V1<0 & output.tsne$V2>0 & output.tsne$V1>(-1.5) & output.tsne$V2<1.5,])+
#   geom_point(aes(V1,V2))+
#   geom_text(aes(V1,V2,label=Group), hjust=0.25)+
#   theme_classic()
# ggsave(plot=plot1, paste0(outdir, Exp.Name, ".VoxTsne.Test_subset.tsne.png"),
#        height=8, width=10)

collected.mindist <- NULL
for(n in c(1:nrow(output.tsne))){
  diff.vector<-abs( abs(output.tsne[n,]$V1 - output.tsne$V1)+abs(output.tsne[n,]$V2 - output.tsne$V2) )
  
  diff.vector[diff.vector==0] <-  1
  temp.diff <- cbind.data.frame(Group=output.tsne[n,]$Group,MinDiff=min(diff.vector), CloseGroup=output.tsne[which.min(diff.vector),]$Group[1])
  collected.mindist<- rbind(collected.mindist, temp.diff)
}

cutoff.for.tsne <- (mean(collected.mindist$MinDiff)-sd(collected.mindist$MinDiff)/2)
pos.close.groups<-collected.mindist[collected.mindist$MinDiff<cutoff.for.tsne,]$Group


plot1 <- ggplot(output.tsne)+
  geom_point(aes(V1,V2))+
  geom_point(data=output.tsne[output.tsne$Group %in% pos.close.groups,], aes(V1,V2), color="red")+
  geom_text(aes(V1,V2,label=Group), hjust=0.25)+
  theme_classic()
ggsave(plot=plot1, paste0(outdir, Exp.Name, ".VoxTsne.Round3.CloseGroups.png"),
       height=8, width=10)


##Sanity check
# unique(variance.tag[variance.tag$Transcode %in% c("DA2", "DA26","DA50","DA57"),]$CGroup_1)
# unique(variance.tag[variance.tag$Transcode %in% c("DA49", "DA68"),]$CGroup_1)
# unique(variance.tag[variance.tag$Transcode %in% c("DA15", "DA12"),]$CGroup_1)


for(Group in collected.mindist[collected.mindist$MinDiff<cutoff.for.tsne,]$Group){
  
  variance.tag[variance.tag$Transcode==Group,]$Transcode<- collected.mindist[collected.mindist$Group==Group,]$CloseGroup
  
}

###Now just need dictionary attack

dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag$Transcode)){
  print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    print("next")}
  ##Maybe I can build a growing dictionary of referenced states
  
  dictionary.spread<-Transcode
  tester.df <- variance.tag[,c("Transcode","Trans_num_final")]
  
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Transcode %in% dictionary.spread,]$Transcode)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num_final %in% dictionary.spread,]$Transcode)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Transcode %in% dictionary.spread,]$Trans_num_final)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num_final %in% dictionary.spread,]$Trans_num_final)
    
    dictionary.spread<-unique(dictionary.spread)
    
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag[variance.tag$Transcode %in% dictionary.spread,]$Transcode2 <- paste0("DC", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}

length(unique(variance.tag$Transcode2));length(unique(variance.tag$Trans_num_final))
## 120 to 94 wo
variance.tag$DRGroup_3 <- variance.tag$Transcode2

}



###We stabilized at 103




################
### TSNE grouping Round 4
##Needs catch
if( (length(unique(variance.tag$Trans_num_final))-length(unique(variance.tag$Transcode2)))<5 ){
  print("finishing")
} else{

variance.tag$Transcode <- variance.tag$DRGroup_3
variance.tag$Trans_num_final <- variance.tag$DRGroup_3
##Memory flag for tsne
if(ncol(normalized.data)>100000){
  
  direct.means <- cbind.data.frame(rn= colnames(mem.tag.table))
  rank_n=1
  for(Transcode in unique(variance.tag$Transcode)){
    print(Transcode)
    top.genes<- variance.tag[variance.tag$Transcode==Transcode,]$Gene
    direct.means$Gene <- colMeans(mem.tag.table[row.names(mem.tag.table) %in% top.genes,])
    names(direct.means)[ncol(direct.means)] <- paste0(Transcode)
  }
  row.names(direct.means) <- direct.means$rn
  direct.means$rn <- NULL
  
  
}else{
direct.means <- cbind.data.frame(rn= row.names(markers))
rank_n=1
for(Transcode in unique(variance.tag$Transcode)){
  print(Transcode)
  top.genes<- variance.tag[variance.tag$Transcode==Transcode,]$Gene
  direct.means$Gene <- colMeans(normalized.data[row.names(normalized.data) %in% top.genes,])
  names(direct.means)[ncol(direct.means)] <- paste0(Transcode)
}
row.names(direct.means) <- direct.means$rn
direct.means$rn <- NULL
}

direct.means.norm <- as.data.frame(t(t(direct.means)/colSums(direct.means)))


library(tsne)
output.tnse <- tsne(t(direct.means.norm), perplexity=50)
##Does work just takes awhile, maybe 20 min for 20k+ cells

output.tsne <- as.data.frame(output.tnse)
output.tsne$Group <- names(direct.means.norm)

###Calculate min distances


diff.vector<-abs( abs(output.tsne[1,]$V1 - output.tsne$V1)+abs(output.tsne[1,]$V2 - output.tsne$V2) )
diff.vector[diff.vector==0] <-  1


collected.mindist <- NULL
for(n in c(1:nrow(output.tsne))){
  diff.vector<-abs( abs(output.tsne[n,]$V1 - output.tsne$V1)+abs(output.tsne[n,]$V2 - output.tsne$V2) )
  
  diff.vector[diff.vector==0] <-  1
  temp.diff <- cbind.data.frame(Group=output.tsne[n,]$Group,MinDiff=min(diff.vector), CloseGroup=output.tsne[which.min(diff.vector),]$Group[1])
  collected.mindist<- rbind(collected.mindist, temp.diff)
}

cutoff.for.tsne <- (mean(collected.mindist$MinDiff)-sd(collected.mindist$MinDiff)/2)
pos.close.groups<-collected.mindist[collected.mindist$MinDiff<cutoff.for.tsne,]$Group


plot1 <- ggplot(output.tsne)+
  geom_point(aes(V1,V2))+
  geom_point(data=output.tsne[output.tsne$Group %in% pos.close.groups,], aes(V1,V2), color="red")+
  geom_text(aes(V1,V2,label=Group), hjust=0.25)+
  theme_classic()
ggsave(plot=plot1, paste0(outdir, Exp.Name, ".VoxTsne.Round4.CloseGroups.png"),
       height=8, width=10)


##Sanity check

for(Group in collected.mindist[collected.mindist$MinDiff<cutoff.for.tsne,]$Group){
  variance.tag[variance.tag$Transcode==Group,]$Transcode<- collected.mindist[collected.mindist$Group==Group,]$CloseGroup
}

###Now just need dictionary attack

dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag$Transcode)){
  print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    print("next")}
  ##Maybe I can build a growing dictionary of referenced states
  
  dictionary.spread<-Transcode
  tester.df <- variance.tag[,c("Transcode","Trans_num_final")]
  
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Transcode %in% dictionary.spread,]$Transcode)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num_final %in% dictionary.spread,]$Transcode)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Transcode %in% dictionary.spread,]$Trans_num_final)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num_final %in% dictionary.spread,]$Trans_num_final)
    
    dictionary.spread<-unique(dictionary.spread)
    
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag[variance.tag$Transcode %in% dictionary.spread,]$Transcode2 <- paste0("DD", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}

length(unique(variance.tag$Transcode2));length(unique(variance.tag$Trans_num_final))
variance.tag$DRGroup_4 <- variance.tag$Transcode2
}




################
### TSNE grouping Round 5

##Needs catch
if( (length(unique(variance.tag$Trans_num_final))-length(unique(variance.tag$Transcode2)))<5 ){
  print("finishing")
} else{

variance.tag$Transcode <- variance.tag$DRGroup_4
variance.tag$Trans_num_final <- variance.tag$DRGroup_4
##Memory flag for tsne
if(ncol(normalized.data)>100000){
  
  direct.means <- cbind.data.frame(rn= colnames(mem.tag.table))
  rank_n=1
  for(Transcode in unique(variance.tag$Transcode)){
    print(Transcode)
    top.genes<- variance.tag[variance.tag$Transcode==Transcode,]$Gene
    direct.means$Gene <- colMeans(mem.tag.table[row.names(mem.tag.table) %in% top.genes,])
    names(direct.means)[ncol(direct.means)] <- paste0(Transcode)
  }
  row.names(direct.means) <- direct.means$rn
  direct.means$rn <- NULL
  
  
}else{
  
direct.means <- cbind.data.frame(rn= row.names(markers))
rank_n=1
for(Transcode in unique(variance.tag$Transcode)){
  print(Transcode)
  top.genes<- variance.tag[variance.tag$Transcode==Transcode,]$Gene
  direct.means$Gene <- colMeans(normalized.data[row.names(normalized.data) %in% top.genes,])
  names(direct.means)[ncol(direct.means)] <- paste0(Transcode)
}
row.names(direct.means) <- direct.means$rn
direct.means$rn <- NULL
}

direct.means.norm <- as.data.frame(t(t(direct.means)/colSums(direct.means)))


library(tsne)
output.tnse <- tsne(t(direct.means.norm), perplexity=50)
##Does work just takes awhile, maybe 20 min for 20k+ cells

output.tsne <- as.data.frame(output.tnse)
output.tsne$Group <- names(direct.means.norm)

###Calculate min distances


diff.vector<-abs( abs(output.tsne[1,]$V1 - output.tsne$V1)+abs(output.tsne[1,]$V2 - output.tsne$V2) )
diff.vector[diff.vector==0] <-  1


collected.mindist <- NULL
for(n in c(1:nrow(output.tsne))){
  diff.vector<-abs( abs(output.tsne[n,]$V1 - output.tsne$V1)+abs(output.tsne[n,]$V2 - output.tsne$V2) )
  
  diff.vector[diff.vector==0] <-  1
  temp.diff <- cbind.data.frame(Group=output.tsne[n,]$Group,MinDiff=min(diff.vector), CloseGroup=output.tsne[which.min(diff.vector),]$Group[1])
  collected.mindist<- rbind(collected.mindist, temp.diff)
}

cutoff.for.tsne <- (mean(collected.mindist$MinDiff)-sd(collected.mindist$MinDiff)/2)
pos.close.groups<-collected.mindist[collected.mindist$MinDiff<cutoff.for.tsne,]$Group


plot1 <- ggplot(output.tsne)+
  geom_point(aes(V1,V2))+
  geom_point(data=output.tsne[output.tsne$Group %in% pos.close.groups,], aes(V1,V2), color="red")+
  geom_text(aes(V1,V2,label=Group), hjust=0.25)+
  theme_classic()
ggsave(plot=plot1, paste0(outdir, Exp.Name, ".VoxTsne.Round5.CloseGroups.png"),
       height=8, width=10)


##Sanity check

for(Group in collected.mindist[collected.mindist$MinDiff<cutoff.for.tsne,]$Group){
  variance.tag[variance.tag$Transcode==Group,]$Transcode<- collected.mindist[collected.mindist$Group==Group,]$CloseGroup
}

###Now just need dictionary attack

dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag$Transcode)){
  print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    print("next")}
  ##Maybe I can build a growing dictionary of referenced states
  
  dictionary.spread<-Transcode
  tester.df <- variance.tag[,c("Transcode","Trans_num_final")]
  
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Transcode %in% dictionary.spread,]$Transcode)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num_final %in% dictionary.spread,]$Transcode)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Transcode %in% dictionary.spread,]$Trans_num_final)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Trans_num_final %in% dictionary.spread,]$Trans_num_final)
    
    dictionary.spread<-unique(dictionary.spread)
    
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag[variance.tag$Transcode %in% dictionary.spread,]$Transcode2 <- paste0("DE", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}

length(unique(variance.tag$Transcode2));length(unique(variance.tag$Trans_num_final))
variance.tag$DRGroup_5 <- variance.tag$Transcode2

}









##bring back in noise
noise.variance.tag2 <- variance.tag2[variance.tag2$Vox_TMerge_Final %in% c("NF1","NF2","NF3","NF4"),]
###lets automate this a bit
length(names(variance.tag)[!names(variance.tag) %in% names(variance.tag2)])
noise.variance.tag <- as.data.frame(matrix(noise.variance.tag2$Vox_TMerge_Final, ncol = length(names(variance.tag)[!names(variance.tag) %in% names(variance.tag2)]), nrow = nrow(noise.variance.tag2)))
names(noise.variance.tag) <- names(variance.tag)[!names(variance.tag) %in% names(variance.tag2)]
noise.variance.tag2<- cbind(noise.variance.tag2, noise.variance.tag)
vox.tag <- rbind(variance.tag, noise.variance.tag2)






vox.tag$Transcode <- vox.tag$Transcode2




write.table(vox.tag, paste0(outdir, Exp.Name,"_VOX_v8.Gene.Table.to.lastRD_nearEnd.txt"),
            sep="\t", quote=F)





# vox.tag <- read.table("/u/project/kp1/jlangerm/Projects/Gomperts_SMG/Analysis/2024-10-11_BG15_organoid_VoxOutput/BG15_SMGOrgo_VOX_v8.VoxGroups.txt",
#                       sep="\t")
# vox.tag$Transcode2 <- vox.tag$DRGroup_5
# vox.tag$Transcode <- vox.tag$Transcode2







###Final Pass

##Depending on what was picked have to rerun above fingerprint? Just paste here....
new.sample.group.data <- cbind.data.frame(Sample=names(new.sample.df))
Group<- "G1"
for(Group in sort(unique(vox.tag$Transcode))){
  print(Group)
  gene.list <- (vox.tag[vox.tag$Transcode==Group,]$Gene)
  new.sample.group.data <- cbind.data.frame( new.sample.group.data ,Score=colSums(new.sample.df[gene.list,]) ) 
  names(new.sample.group.data)[ncol(new.sample.group.data)]<- paste0(Group)
}
new.sample.group.data$Sample <- NULL
new.sample.group.data <-((new.sample.group.data)/rowSums(new.sample.group.data))
rowSums(new.sample.group.data)

enrich.data <- NULL
enrich.data <- cbind.data.frame(Group=names(new.sample.group.data))
for(Sample in names(new.sample.df)){
  group.mean.tester <- new.sample.group.data
  group.mean.tester$Sample <- NULL
  diff.groups.tester<-(group.mean.tester[Sample,]+0.01)/(colMeans(group.mean.tester[row.names(group.mean.tester)!=Sample,])+0.01)
  diff.groups.tester<-unlist(diff.groups.tester)
  enrich.data<-cbind(enrich.data, Data=diff.groups.tester)
  names(enrich.data)[ncol(enrich.data)] <- Sample
  
}

row.names(enrich.data) <- enrich.data$Group
enrich.data$Group<-NULL

ml3.table<-cbind.data.frame(Group=names(apply(enrich.data,1, max)),
                            Max=unlist(apply(enrich.data,1, max)),
                            Min=unlist(apply(enrich.data,1, min)))
ml3.table$Max_Up <- (ml3.table$Max-1)
ml3.table$Min_Up <- (1-ml3.table$Min)
ml3.table$AvgEnrich <- (ml3.table$Max_Up+ ml3.table$Min_Up)/2
ml3.table$Enrich_Bin <- cut(ml3.table$AvgEnrich, labels=F, breaks=20)


ml3.table <- ml3.table[with(ml3.table, order(-Enrich_Bin)),]

# height=4, width=5)
plot1<-ggplot(ml3.table)+
  geom_bar(aes(Group, AvgEnrich), stat="identity",fill="black")+
  theme_classic()+
  geom_hline(yintercept = median(ml3.table$AvgEnrich), color="red")+
  theme(axis.text.x=element_blank())
ggsave(plot=plot1,paste0(outdir, Exp.Name, ".VoxGroup.Enrichment.data.bargraph.png"),
       height = 3, width=4)






library(pheatmap)
ranged.means <- t(t(new.sample.group.data)/apply(new.sample.group.data, 2, max))
ordering.groups.df <- pheatmap(ranged.means, silent=T)

pheatmap(ranged.means,
         # clustering_method = "mcquitty",
         filename=paste0(outdir, Exp.Name, ".VoxNearFinal.MeanExpGroups.heatmap.png"),
         height=8, width=20)

# unique(vox.tag$VoxGroup)
# vox.tag[vox.tag$Transcode="DD35",]

# ordering.groups.df$tree_row$labels[ordering.groups.df$tree_row$order]
# ordering.groups.df$tree_col$labels[ordering.groups.df$tree_col$order]

###I dont think this is written right
# sg.table<- cbind.data.frame(Group=ordering.groups.df$tree_col$labels, Order=ordering.groups.df$tree_col$order)

sg.table<- cbind.data.frame(Group=ordering.groups.df$tree_col$labels[ordering.groups.df$tree_col$order], Order=c(1:length(unique(vox.tag$Transcode))))

sg.table[grepl("NF", sg.table$Group),]$Order<- sg.table[grepl("NF", sg.table$Group),]$Order +1000

# sg.table[grepl("TI", sg.table$Group),]$Order<- sg.table[grepl("TI", sg.table$Group),]$Order +1000


sg.table<- sg.table[with(sg.table, order(Order)),]
sg.table$Group_New <- paste0("V" ,c(1:nrow(sg.table)) )

##Tag noice
sg.table[grepl("NF", sg.table$Group),]$Group_New <- paste0("Z.",sg.table[grepl("NF", sg.table$Group),]$Group_New)
# sg.table[grepl("TI", sg.table$Group),]$Group_New <- "Noise1"
vox.tag$VoxGroup <- sg.table[match(vox.tag$Transcode,sg.table$Group),]$Group_New

length(unique(vox.tag$VoxGroup))
table(vox.tag$VoxGroup)


###Preserve for Ontorun
# vox.tag$OldVox <- vox.tag$VoxGroup


vox.tag$VoxGroup<- factor(vox.tag$VoxGroup, levels=unique(sg.table$Group_New))
vox.tag<-vox.tag[with(vox.tag, order(vox.tag$VoxGroup)),]

# sg.table[sg.table$SG %in% c("NF1", "NF2","NF3","NF4"),]$GG_New <- paste0("ZNoise",sg.table[sg.table$SG %in% c("NF1", "NF2","NF3","NF4"),]$GG_New)



vox.tag$MetaLayerEnrichment <- ml3.table[match(vox.tag$Transcode, ml3.table$Group),]$AvgEnrich



# vox.tag$MetaLayerSpecificity <- 0
vox.tag$MetaLayerCorrelation <- 0

for(Group in unique(vox.tag$VoxGroup)){
  print(Group)
  genes.of.interest<- vox.tag[vox.tag$VoxGroup==Group,]$Gene
  norm.data.chunk <- normalized.data[row.names(normalized.data) %in% genes.of.interest,]
  
  score.piece<-colMeans(normalized.data[row.names(normalized.data) %in% genes.of.interest, ])
  score.cor <- cor(t(norm.data.chunk), score.piece)
  
  
  final.gene.tester<- cbind.data.frame(Gene=row.names(score.cor),Specificity=0, Score_Corr =as.numeric(score.cor))
  
  vox.tag[vox.tag$VoxGroup==Group,]$MetaLayerCorrelation <- final.gene.tester[match(vox.tag[vox.tag$VoxGroup==Group,]$Gene,final.gene.tester$Gene),]$Score_Corr
  
}

vox.tag[is.na(vox.tag)]<-0


library(dplyr)

vox.tag %>% group_by(VoxGroup) %>% summarize(Corr=mean(MetaLayerCorrelation)) -> metalayer.table

metalayer.table$VoxGroup <- factor(metalayer.table$VoxGroup, levels=unique(vox.tag$VoxGroup))

metalayer.table$Enrichment <- vox.tag[match(metalayer.table$VoxGroup,vox.tag$VoxGroup),]$MetaLayerEnrichment
metalayer.table$Trajectory <- vox.tag[match(metalayer.table$VoxGroup,vox.tag$VoxGroup),]$MetaLayerTrajectory





options(scipen = 999)

write.table(vox.tag, paste0(outdir, Exp.Name,"_VOX_v8.VoxGroups.txt"),
            sep="\t", quote=F)





















################################
###Plot group expression on UMAP


dir.create(paste0(outdir, "VoxGroup.Expression/"))
outdir2<- paste0(outdir, "VoxGroup.Expression/")

# variance.tag$ML1_Index <- paste0(variance.tag$MetaLayer2,"-",variance.tag$MetaLayer1)

transcode.table <- as.data.frame(table(vox.tag$VoxGroup))
# n=1
# Exp.Name<- "YS3.Filt"
library(ggplot2)
for(Group in unique(vox.tag$VoxGroup)){
  print(Group)
  plot.data <- markers
  genes.of.interest <- vox.tag[vox.tag$VoxGroup==Group,]$Gene
  plot.data$Score <- colMeans(normalized.data[row.names(normalized.data) %in% genes.of.interest,colnames(normalized.data) %in% row.names(plot.data)])
  
  amount.of.genes<-transcode.table[transcode.table$Var1==Group,]$Freq
  
  plot.data <- plot.data[with(plot.data, order(Score)),]
  ggplot(plot.data)+
    # facet_wrap(~Type)+
    geom_point(data=markers[,c("X","Y")], aes(X,Y), color="grey55")+
    geom_point( aes(X,Y, color=Score) )+
    theme_classic()+
    labs(color=Group)+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ggtitle(paste0("VOX ",Group, " with ", amount.of.genes," genes"))+
    scale_colour_gradientn(colours=c("lightcyan3", "gold", "orange","orangered","purple"))
  ggsave(filename=paste0(outdir2,  Exp.Name, ".Vox.avgExp.VoxGroup_",Group,".UMAP.png"),
         height=4, width=6)
}





########################################
######## Create Outputs for metascape - Keeps order


##All Genes
Program.Table3 <- vox.tag
Program.Table3$FinalSet <- Program.Table3$VoxGroup
Program.Table3$Genes <- Program.Table3$Gene
Network.Source<- "VOX_Groups_All"


if(length(unique(Program.Table3$FinalSet))>30){
  print("Too many networks, dividing into parts")
  Cut.Amount <- c(1:ceiling( length(unique(Program.Table3$FinalSet))/30))
  Take.Amount <- 1
  Network.Cut.Table <- cbind.data.frame(Network=unique(Program.Table3$FinalSet),Cut.Identity=0)
  for(Cutn in c(1: max(Cut.Amount))){
    top.amount <- min( c(Take.Amount+round(nrow(Network.Cut.Table)/Cut.Amount), nrow(Network.Cut.Table) ) )
    Network.Cut.Table[c(Take.Amount:top.amount),]$Cut.Identity<-  Cutn
    Take.Amount <- Take.Amount+ round(nrow(Network.Cut.Table)/ max(Cut.Amount))+1
  }
  Program.Table3$Cut.Identity <- Network.Cut.Table[match(Program.Table3$FinalSet, Network.Cut.Table$Network),]$Cut.Identity
  for(Cut.Identity in unique(Program.Table3$Cut.Identity)){
    print(paste0("part ", Cut.Identity))
    Metascape.Multi <- NULL
    for(Set in unique(Program.Table3[Program.Table3$Cut.Identity==Cut.Identity,]$FinalSet)){
      tmp.Prog.Table <- Program.Table3[Program.Table3$FinalSet==Set,]
      Metascape.tmp<-cbind.data.frame( Name=paste0(Set),Genes=paste(as.character(tmp.Prog.Table$Genes), collapse=",")  )
      Metascape.Multi <- rbind(Metascape.Multi, Metascape.tmp)
    }
    print("done parsing")
    write.table(Metascape.Multi, paste0(outdir, Exp.Name,".",Network.Source,".Part",Cut.Identity,".Metanetworks.MultiGO.List.txt"), sep="\t", col.names =F, row.names = F,quote=F)
    
  }
}else{
  ##One loop version
  print("Network Numbers Okay")
  Metascape.Multi <- NULL
  for(Set in unique(Program.Table3$Group)){
    Program.Table.temp <- Program.Table3[Program.Table3$Group==Set,]
    Metascape.tmp<-cbind.data.frame( Name=paste0(Set),Genes=paste(as.character(Program.Table.temp$Genes), collapse=",")  )
    Metascape.Multi <- rbind(Metascape.Multi, Metascape.tmp)
  }
  print("done full")
  write.table(Metascape.Multi, paste0(outdir, Exp.Name,".VoxGroups.MultiGO.List.txt"), sep="\t", col.names =F, row.names = F,quote=F)
}







