
####################
##### Dependencies 

# install.packages("RTools") ###this is a propagate dependency for installation
library(propagate)
library(ggplot2)
library(dplyr)



################
##### Settings

########PATHS
dir.create("/outdir/")
outdir <- ("/outdir/")

####Name your sample
drop <- "Experiment Title"


######### Data


####Cell metadata
cell.design <- (Your Metadata File)
###rownames must match colnames of normalized.data
###Please provide reduced dimensionality coordinates as columns named "X" and "Y"
###Must be a dataframe





### digital gene expression matrix - must be a dataframe
normalized.data <- (Your GEX Dataframe)
###rownames must be Gene names

normalized.data2 <- normalized.data


#####OPTIONAL
###Adjust data to eliminate maximum for each gene - this removes "burst" signals from creating spurious networks, comment if unwanted
#############May want to comment if cell type of interest is extremely few cells
normalized.data2 <- normalized.data
max.store<-apply(normalized.data2, 1, which.max)
for(Gene in names(max.store)){
  normalized.data2[Gene,max.store[Gene]] <-0
}
normalized.data <- normalized.data2



###Must pass checks to proceed
identical(row.names(cell.design), colnames(normalized.data)) ##Must be TRUE
nrow(normalized.data[duplicated(row.names(normalized.data)),]) ##Must be 0





#############################################
############ Calculate GEND

###This is windows specific command, does nothing in *nix
memory.limit(10000000) 





#######Part 1 - Correlations
###Estimated time 1-4 hours

###If you cannot get propagate to install, can replace this with just "cor" function. Extremely memory intensive though
tf.cor <- bigcor(t(normalized.data2))


tf.cor <- tf.cor[1:nrow(tf.cor), 1:ncol(tf.cor)]
diag(tf.cor) = 0
colnames(tf.cor)<-row.names(normalized.data)
row.names(tf.cor)<-row.names(normalized.data)


library(reshape2)
corr <- melt(tf.cor)
names(corr) <- c('Source', 'Target', 'Corr')

##### CUTOFF for all base networks - No need to change, can make more restricted below without affecting anything
###Based on general scRNAseq Corr distributions, adjustment here unlikely to change outputs
cutoff = 0.08
corr2 <- subset(corr, Corr>cutoff)
names(corr2)[3] = 'Edge'

corr2 <- corr2[with(corr2, order(-Edge)),]
corr2 <- corr2[!duplicated(corr2$Source),]
corr2 <- corr2[with(corr2, order(Source)),]

network.cor.master2 <- cbind.data.frame(Genes=corr2$Source, Corr=corr2$Edge, Set=corr2$Target)
write.table(network.cor.master2, paste0(outdir, drop, ".GEND.Networks.Raw.Cor.",cutoff,".txt"),row.names=FALSE, quote=FALSE, sep="\t")
### you can skip the corr step once you have this Raw output, if you want to tweak downstream


#########






#####Cluster Data -
#### This is MUCH FASTER if you cluster on a precomputed reduced dimensionality, makes sure those vectors are in cell.design and columns X and Y

####
print("Phase 2 - Collapsing")
print("Phase 2 - Part 1 - Kmeans of expression space")



#####Set the number of cells to use for kmeans clustering. 25 for small datasets, 50 for 1000-10000, 100 for 10k+
num.of.cells <- 25


####K clustering on reduced dimensionality
names(cell.design)

cell.design$Cluster <- as.factor(kmeans(cell.design[,c("X","Y")], round(nrow(cell.design)/num.of.cells))$cluster)
ggplot(cell.design)+
  geom_text(aes(X,Y,label=Cluster, color=as.factor(Cluster)))+
  # ThemeOpts()+
  theme(legend.position = "")
  ggsave(paste0(outdir, drop, ".GEND.Kmeans.Tsne.png"), height=15, width=18)




 
  
  #########Use k on expression space, optimized
  # ##Makes 75mer cell groups, approximately
  # approx.amount <- round(nrow(cell.design)/(25*50))
  # cell.design$RandomCluster <- rep(1:approx.amount, nrow(cell.design))[1:nrow(cell.design)]
  # cell.design$KCluster<-0
  # for(Cluster in unique(cell.design$RandomCluster)){
  #   print(Cluster)
  #   cell.design[cell.design$RandomCluster==Cluster,]$KCluster <- as.factor(kmeans(t(normalized.data[, names(normalized.data) %in% row.names(cell.design[cell.design$RandomCluster==Cluster,])]),  25 )$cluster)
  # }
  # cell.design$Cluster<-  paste0(cell.design$RandomCluster,"00", cell.design$KCluster)
  ########################
  
  
  
  ########################
  ##If looking corr raw up from previous run
  # network.cor.master2 <- read.table(paste0("Your.Network.Raw.Cor.0.08.txt"), header=TRUE)
  ###########################
  
  
  
  #################Set Cutoff manually OR 
  ##Strict Cutoff - Higher = Small distinct networks, Lower = Big Networks More Genes
  # cutoff <-0.10
  
  
  ##########Determine correlation empirically
  
  target.threshold <- nrow(normalized.data)/2
  cutoff <-0.4
  nrow(network.cor.master2[network.cor.master2$Corr>=cutoff,])
  
  repeat {
    # print(network.regress[reg.number,3])
    print(cutoff)
    cutoff <- cutoff - 0.01
    # exit if the condition is met
    if ( nrow(network.cor.master2[network.cor.master2$Corr>=cutoff,]) > target.threshold) break
  }
  cutoff

  
  
###############################
###### Network building
  ##Time estimate: 20min-40min
  
  
  network.master.list <-network.cor.master2[network.cor.master2$Corr>=cutoff,]
  network.master.list <-network.master.list[with(network.master.list, order(-Corr)),]
  network.master.list2 <- network.master.list[!duplicated(network.master.list$Gene),]
  network.master.list2$RegSet <- 0
  # length(unique(network.master.list2$Set))
  
  nrow((network.master.list2))
    ###Is this about 5000-10000 genes
  
  n=1
  # gene="X"
  for(gene in (unique(network.master.list2$Set))){
    print(gene)
    if(unique(network.master.list2[network.master.list2$Set==gene,]$RegSet)==0){
      
      ##inital call
      network.master.list2[network.master.list2$Set==gene,]$RegSet <- n
      ## this part spreads the plague  1 edge
      for(gene2 in network.master.list2[network.master.list2$Set==gene,]$Gene){
        if(nrow(network.master.list2[network.master.list2$Set==gene2,])>0){
          network.master.list2[network.master.list2$Set==gene2,]$RegSet <- n
          ## second round 2 edges
          for(gene3 in network.master.list2[network.master.list2$Set==gene2,]$Gene){
            if(nrow(network.master.list2[network.master.list2$Set==gene3,])>0){
              network.master.list2[network.master.list2$Set==gene3,]$RegSet <- n
              ## third round 3 edges
              for(gene4 in network.master.list2[network.master.list2$Set==gene3,]$Gene){
                if(nrow(network.master.list2[network.master.list2$Set==gene4,])>0){
                  network.master.list2[network.master.list2$Set==gene4,]$RegSet <- n
                  for(gene5 in network.master.list2[network.master.list2$Set==gene4,]$Gene){
                    if(nrow(network.master.list2[network.master.list2$Set==gene5,])>0){
                      network.master.list2[network.master.list2$Set==gene5,]$RegSet <- n
                      ###5 edges
                      for(gene6 in network.master.list2[network.master.list2$Set==gene5,]$Gene){
                        if(nrow(network.master.list2[network.master.list2$Set==gene6,])>0){
                          network.master.list2[network.master.list2$Set==gene6,]$RegSet <- n
                        }
                      }# gene 6
                    }# if gene5
                  }#gene 5
                }#if gene4
              }#gene4
            }#if gene3
          }#gene3
          
        } #if gene2
        
      } #gene2 forloop 
      n=n+1
    } else {
      print("decided")
    }
  } #first loop
  
  length(unique(network.master.list2$RegSet))
  
  network.cor.master2 <- network.master.list2
  
  #####
  ###Collapse similar networks based on expression
  meganetwork.data <-NULL
  meganetwork.data <- cbind.data.frame(Cluster=unique(cell.design$Cluster))
  for (j in unique(network.cor.master2$RegSet)){
    cluster.mean <-NULL
    print(paste0(j, " Calc Cluster Exp"))
    
    
    network.master.list <- network.cor.master2[network.cor.master2$RegSet==j,]$Genes
    network.master.list <- unique(network.master.list)
    
    for (k in unique(cell.design$Cluster)){
      
      cluster.data <- normalized.data[,colnames(normalized.data) %in% row.names(cell.design[cell.design$Cluster==k,]), drop=F]
      
      if(ncol(cluster.data)>1){ ###Need this modifyer because kmeans is giving me 1 cell clusters
        cluster.mean.temp  <- cbind.data.frame(Set=j, Cluster=k, Score=mean(colMeans( cluster.data[row.names(cluster.data) %in% network.master.list,])))
      }else{##mod
        cluster.mean.temp  <- cbind.data.frame(Set=j, Cluster=k, Score=mean(cluster.data[[1]]) ) ##mod
      }##mod
      
      
      cluster.mean <- rbind(cluster.mean , cluster.mean.temp)
    }
    cluster.mean <- as.data.frame(cluster.mean[,c(3)])
    colnames(cluster.mean) <- j
    meganetwork.data <- cbind.data.frame(meganetwork.data, cluster.mean)
  }
  row.names(meganetwork.data)<- meganetwork.data$Cluster
  meganetwork.data$Cluster <-NULL
  
  ######Meganetwork collapse
  meganetwork.cor <- as.data.frame(cor(meganetwork.data ))
  diag(meganetwork.cor) <- 0
  sort(meganetwork.cor[[1]])
  ####Collapse Cutoff
  net.cor.cutoff <- 0.80
  
  
  network.cor.master<- NULL
  for (j in names(meganetwork.cor)){
    print(j)
    if (length(row.names(meganetwork.cor[meganetwork.cor[[j]]>=net.cor.cutoff, ,drop=F]))>0){
      network.correlation2<-cbind.data.frame(Genes=row.names(meganetwork.cor[meganetwork.cor[[j]]>=net.cor.cutoff, ,drop=F]),
                                             Corr=meganetwork.cor[meganetwork.cor[[j]]>=net.cor.cutoff,,drop=FALSE][[j]], Set=j)
      
      network.cor.master <- rbind(network.cor.master, network.correlation2)
    } else {
      print ("None")
    }
  }
  network.cor.master<-as.data.frame(network.cor.master)
  nrow(network.cor.master)
  length(unique(network.cor.master$Set))
  
  
  ####Network again
  network.master.list <-network.cor.master
  network.master.list2 <-network.master.list[with(network.master.list, order(-Corr)),]
  network.master.list2 <- network.master.list2[!duplicated(network.master.list2$Gene),]
  network.master.list2$RegSet <- 0
  length(unique(network.master.list2$Set))
  n=1
  # gene="X"
  for(gene in (unique(network.master.list2$Set))){
    print(gene)
    if(unique(network.master.list2[network.master.list2$Set==gene,]$RegSet)==0){
      
      ##inital call
      network.master.list2[network.master.list2$Set==gene,]$RegSet <- n
      ## this part spreads the plague  1 edge
      for(gene2 in network.master.list2[network.master.list2$Set==gene,]$Gene){
        if(nrow(network.master.list2[network.master.list2$Set==gene2,])>0){
          network.master.list2[network.master.list2$Set==gene2,]$RegSet <- n
          ## second round 2 edges
          for(gene3 in network.master.list2[network.master.list2$Set==gene2,]$Gene){
            if(nrow(network.master.list2[network.master.list2$Set==gene3,])>0){
              network.master.list2[network.master.list2$Set==gene3,]$RegSet <- n
              ## third round 3 edges
              for(gene4 in network.master.list2[network.master.list2$Set==gene3,]$Gene){
                if(nrow(network.master.list2[network.master.list2$Set==gene4,])>0){
                  network.master.list2[network.master.list2$Set==gene4,]$RegSet <- n
                  for(gene5 in network.master.list2[network.master.list2$Set==gene4,]$Gene){
                    if(nrow(network.master.list2[network.master.list2$Set==gene5,])>0){
                      network.master.list2[network.master.list2$Set==gene5,]$RegSet <- n
                      ###5 edges
                      for(gene6 in network.master.list2[network.master.list2$Set==gene5,]$Gene){
                        if(nrow(network.master.list2[network.master.list2$Set==gene6,])>0){
                          network.master.list2[network.master.list2$Set==gene6,]$RegSet <- n
                        }
                      }# gene 6
                    }# if gene5
                  }#gene 5
                }#if gene4
              }#gene4
            }#if gene3
          }#gene3
          
        } #if gene2
        
      } #gene2 forloop 
      n=n+1
    } else {
      print("decided")
    }
  } #first loop
  length(unique(network.master.list2$RegSet))
  length(unique(network.master.list2$Genes))
  
  
  #####

  network.mega.regress <- network.cor.master2
  
  network.master.list2$RegSet2<-0
  for(Group in unique(network.master.list2$RegSet)){
    
    network.master.list2[network.master.list2$RegSet==Group,]$RegSet2 <- network.master.list2[network.master.list2$RegSet==Group,]$Set[1]
  }
  
  network.mega.regress$FinalSet <- network.master.list2$RegSet2[match(network.mega.regress$RegSet,network.master.list2$Genes)]
  network.mega.regress[is.na(network.mega.regress$FinalSet),]$FinalSet <-network.mega.regress[is.na(network.mega.regress$FinalSet),]$RegSet
  
 
  network.size.table <- as.data.frame(table(network.mega.regress$FinalSet))
  network.mega.regress$Network.Size <- network.size.table$Freq[match(network.mega.regress$FinalSet, network.size.table$Var1)]
  
  ##
  
  
  
  # ######################No Reason to do stragglers if doing 2nd pass
  # ###############
  # ### Straggler Networks Merge
  # ##Small genes networks that may be formed by recursive attraction or noise. Comment out if you do not want these forcibly pooled into other networks
  # ##
  # network.mega.regress.small <- network.mega.regress[network.mega.regress$Network.Size<6,]
  # sort(apply(meganetwork.cor[,names(meganetwork.cor) %in% unique(network.mega.regress.small$RegSet)],2,max) )
  # allowable.floor.cor = 0.8
  # for (set.number in unique(network.mega.regress.small$FinalSet)){
  #   # if(nrow(network.mega.regress[network.mega.regress$FinalSet==set.number,])<5){
  #   print(paste0("Calculating new max correlation for small network ", set.number))
  #   # n=n+1
  #   set.number2 <- network.mega.regress.small[network.mega.regress.small$FinalSet==set.number,]$RegSet[1]
  #   
  #   reg.number <- row.names(meganetwork.cor)[apply(meganetwork.cor[,set.number2,drop=F] , 2,which.max)]
  #   final.number <- network.mega.regress[network.mega.regress$RegSet==reg.number,]$FinalSet[1]
  #   print(paste0("Found Corr ", max(meganetwork.cor[,set.number2]), " to Prime Network ", final.number))
  #   
  #   
  #   if(max(meganetwork.cor[,set.number]) > allowable.floor.cor){ 
  #     print("Merging")
  #     # print(network.mega.regress[network.mega.regress$RegSet==row.names(meganetwork.cor)[[1]], ]$FinalSet[[1]])
  #     network.mega.regress[network.mega.regress$FinalSet==set.number,]$FinalSet <- final.number
  #   } else {
  #     print("Below Acceptable Threshold, Not Merging")
  #   }
  # }
  #################################### Straggler merge
  ####################################
  
  
  

  
##################  
####Grade/Filter Networks
  network.mega.regress$Pct.Cells.Exp <- 0
  pct.test.df <- normalized.data[network.mega.regress$Genes, ]
  pct.test.df[pct.test.df>0]<- 1
  network.mega.regress$Pct.Cells.Exp <-rowSums(pct.test.df)/ncol(normalized.data)
  

  
  
  #### Final output prep
  network.mega.regress$FinalSet <- as.numeric(network.mega.regress$FinalSet)
  
  ###Checking network size and cutting off small networks
  network.size.table <- as.data.frame(table(network.mega.regress$FinalSet))
  network.mega.regress$Network.Size <- network.size.table$Freq[match(network.mega.regress$FinalSet, network.size.table$Var1)]
  network.mega.regress <- network.mega.regress[network.mega.regress$Network.Size>10,]
  

  ###Create network intercorrelation table to determine approximate gene-gene agreement
  network.size.table <- as.data.frame(table(network.mega.regress$FinalSet))
  network.size.table$InsideCorr=0
  network.size.table$OutsideCorr=0
  network.size.table$Avg.Pct.Exp <- 0
  network.size.table$CellAdjust <- 100
  network.size.table$SD.Out <- 0
  network.size.table$SD.In <- 0
  network.size.table$SD.All <- 0
  
  for(Program in unique(network.mega.regress$FinalSet)){
    markerlist4u <- network.mega.regress[network.mega.regress$FinalSet==Program,]$Genes
    markers<-cell.design
    markers$Score <- colMeans(normalized.data[row.names(normalized.data) %in% markerlist4u,])
    
    thresh.cutoff <- (max(markers$Score)*0.50)
    
    network.size.table[network.size.table$Var1==Program,]$Avg.Pct.Exp <- mean(network.mega.regress[network.mega.regress$FinalSet==Program,]$Pct.Cells.Exp)
    
    ##Sd
    cell.net <- normalized.data [row.names(normalized.data) %in% markerlist4u,]
    network.size.table[network.size.table$Var1==Program,]$SD.All <- mean(unlist(sd.obj))
    
    cell.net <- cell.net[,names(cell.net) %in% row.names(markers[markers$Score > thresh.cutoff,]),drop=FALSE ]
    
    if(ncol(cell.net)>1){
      cell.net <- cell.net[, colSums(cell.net)>0]
      cell.cor <-cor(cell.net)
      network.size.table[network.size.table$Var1==Program,]$InsideCorr<- mean(cell.cor)
      
      cell.net<-cell.net+0.0001
      sd.obj<-lapply(as.data.frame(t(cell.net/rowSums(cell.net))), sd)
      network.size.table[network.size.table$Var1==Program,]$SD.In <- mean(unlist(sd.obj))
      
      
    }else{ ###if only 1 cell qualifies, will just take top 10 cells?
      
      markers<- markers[with(markers, order(-Score)),]
      cell.net <- normalized.data [row.names(normalized.data) %in% markerlist4u,]
      cell.net <- cell.net[,names(cell.net) %in% row.names(markers[1:10,]),drop=FALSE ]
      cell.net <- cell.net[, colSums(cell.net)>0]
      
      ###Add a random pseudocount to break 0 distributions
      cell.net <- cell.net +(runif(length(cell.net),0.0001, 0.0002))
      
      cell.cor <-cor(cell.net)
      diag(cell.cor) <- 0
      network.size.table[network.size.table$Var1==Program,]$InsideCorr<- mean(cell.cor)
      
      test.df<-as.data.frame(table(colnames(cell.net)[apply(cell.net ,1,which.max)]) )
      test.df <- test.df[with(test.df, order(-Freq)),]
      network.size.table[network.size.table$Var1==Program,]$CellAdjust<- 100*test.df$Freq[1]/sum(test.df$Freq)
      
      sd.obj<-lapply(as.data.frame(t(cell.net/rowSums(cell.net))), sd)
      network.size.table[network.size.table$Var1==Program,]$SD.In<- mean(unlist(sd.obj))
    }
    
    table(names(cell.net) %in% row.names(markers[markers$Score > thresh.cutoff,]))
    
    
    cell.net <- normalized.data [row.names(normalized.data) %in% markerlist4u,]
    cell.net <- cell.net[,!names(cell.net) %in% row.names(markers[markers$Score > thresh.cutoff,]) ]
    cell.net <- cell.net[, colSums(cell.net)>0]
    cell.cor <-cor(cell.net)
    network.size.table[network.size.table$Var1==Program,]$OutsideCorr<- mean(cell.cor)
   
    cell.net<-cell.net+0.0001
    sd.obj<-lapply(as.data.frame(t(cell.net/rowSums(cell.net))), sd)
    network.size.table[network.size.table$Var1==Program,]$SD.Out <- mean(unlist(sd.obj))
    
    
  }
  
  network.size.table$DiffCorr <- network.size.table$InsideCorr - network.size.table$OutsideCorr
 
  network.size.table$Adjust.Corr <- network.size.table$DiffCorr * (network.size.table$Freq) * (network.size.table$CellAdjust)
  
  ##added log term to adjust.corr2 to soften pct exp scaling
  ########This doesn't work - reverting + abs(min(log(network.size.table$Avg.Pct.Exp)) )
  network.size.table$Adjust.Corr2 <- network.size.table$Adjust.Corr * (network.size.table$Avg.Pct.Exp)

  
  
  network.size.table<- network.size.table[with(network.size.table, order(-Adjust.Corr2)),]
  network.size.table$LikelyName <- c(1:nrow(network.size.table))
  
  
  network.mega.regress$DiffCorr <- network.size.table$DiffCorr[match(network.mega.regress$FinalSet, network.size.table$Var1)]
  network.mega.regress$Adjust.Corr <- network.size.table$Adjust.Corr2[match(network.mega.regress$FinalSet, network.size.table$Var1)]
  
  network.mega.regress <- network.mega.regress[with(network.mega.regress, order(-Adjust.Corr)),]
  
  
  network.mega.regress$TempSet <- 0
  temp.set=1
  for(Set.Number in unique(network.mega.regress$FinalSet) ){
    network.mega.regress[network.mega.regress$FinalSet==Set.Number,]$TempSet <- temp.set
    temp.set=temp.set+1
  }
  
  network.mega.regress$FinalSet <-  network.mega.regress$TempSet
  network.mega.regress$TempSet <- NULL

  network.mega.regress <- network.mega.regress[with(network.mega.regress, order(FinalSet)),]
  
  
  
  
  ###Quality judgement (may be skewed)
  network.mega.regress[is.na(network.mega.regress)] <- 0
  network.mega.regress$Net.Quality <- c("Noisy")
  network.mega.regress[network.mega.regress$Adjust.Corr> quantile(network.mega.regress$Adjust.Corr)[3] ,]$Net.Quality <- "Questionable"
  network.mega.regress[network.mega.regress$Adjust.Corr> quantile(network.mega.regress$Adjust.Corr)[4],]$Net.Quality <- "Good"
  
  
  ##Finish and write
  write.table(network.mega.regress, paste0(outdir, drop, ".GEND.Networks.Output_Cor.",cutoff,".txt"),row.names=FALSE, quote=FALSE, sep="\t")
  print("GEND Finished!")
  
  
  Program.Table <- network.mega.regress
  # Program.Table$Gene <- row.names(Program.Table)
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  ###########################
  ############
  #### Part 2b - 2nd pass to pool lower quality networks

  network.cor.master2 <- cbind.data.frame(Genes=corr2$Source, Corr=corr2$Edge, Set=corr2$Target)
  
  
  ##remove first pass genes
  network.cor.master2 <- network.cor.master2[!network.cor.master2$Genes %in% Program.Table$Gene,]
  network.cor.master2 <- network.cor.master2[!network.cor.master2$Set %in% Program.Table$Gene,]
  
  # ##########Determine correlation empirically
  # target.threshold <- nrow(normalized.data)/3
  # cutoff <-0.4
  # nrow(network.cor.master2[network.cor.master2$Corr>=cutoff,])
  # 
  # repeat {
  #   # print(network.regress[reg.number,3])
  #   print(cutoff)
  #   cutoff <- cutoff - 0.01
  #   # exit if the condition is met
  #   if ( nrow(network.cor.master2[network.cor.master2$Corr>=cutoff,]) > target.threshold) break
  # }
  
  ###### Network building
  
  ##Take everything left basically
  network.master.list <-network.cor.master2
  network.master.list <-network.master.list[with(network.master.list, order(-Corr)),]
  network.master.list2 <- network.master.list[!duplicated(network.master.list$Gene),]
  network.master.list2$RegSet <- 0
  length(unique(network.master.list2$Set))
  
  n=1
  # gene="X"
  for(gene in (unique(network.master.list2$Set))){
    print(gene)
    if(unique(network.master.list2[network.master.list2$Set==gene,]$RegSet)==0){
      
      ##inital call
      network.master.list2[network.master.list2$Set==gene,]$RegSet <- n
      ## this part spreads the plague  1 edge
      for(gene2 in network.master.list2[network.master.list2$Set==gene,]$Gene){
        if(nrow(network.master.list2[network.master.list2$Set==gene2,])>0){
          network.master.list2[network.master.list2$Set==gene2,]$RegSet <- n
          ## second round 2 edges
          for(gene3 in network.master.list2[network.master.list2$Set==gene2,]$Gene){
            if(nrow(network.master.list2[network.master.list2$Set==gene3,])>0){
              network.master.list2[network.master.list2$Set==gene3,]$RegSet <- n
              ## third round 3 edges
              for(gene4 in network.master.list2[network.master.list2$Set==gene3,]$Gene){
                if(nrow(network.master.list2[network.master.list2$Set==gene4,])>0){
                  network.master.list2[network.master.list2$Set==gene4,]$RegSet <- n
                  for(gene5 in network.master.list2[network.master.list2$Set==gene4,]$Gene){
                    if(nrow(network.master.list2[network.master.list2$Set==gene5,])>0){
                      network.master.list2[network.master.list2$Set==gene5,]$RegSet <- n
                      ###5 edges
                      for(gene6 in network.master.list2[network.master.list2$Set==gene5,]$Gene){
                        if(nrow(network.master.list2[network.master.list2$Set==gene6,])>0){
                          network.master.list2[network.master.list2$Set==gene6,]$RegSet <- n
                        }
                      }# gene 6
                    }# if gene5
                  }#gene 5
                }#if gene4
              }#gene4
            }#if gene3
          }#gene3
          
        } #if gene2
        
      } #gene2 forloop 
      n=n+1
    } else {
      print("decided")
    }
  } #first loop
  
  length(unique(network.master.list2$RegSet))
  
  network.cor.master2 <- network.master.list2
  
  #####
  ###Collapse similar networks based on expression
  meganetwork.data <-NULL
  meganetwork.data <- cbind.data.frame(Cluster=unique(cell.design$Cluster))
  for (j in unique(network.cor.master2$RegSet)){
    cluster.mean <-NULL
    print(paste0(j, " Calc Cluster Exp"))
    
    
    network.master.list <- network.cor.master2[network.cor.master2$RegSet==j,]$Genes
    network.master.list <- unique(network.master.list)
    
    for (k in unique(cell.design$Cluster)){
      
      cluster.data <- normalized.data[,colnames(normalized.data) %in% row.names(cell.design[cell.design$Cluster==k,]), drop=F]
      
      if(ncol(cluster.data)>1){ ###Need this modifyer because kmeans is giving me 1 cell clusters
        cluster.mean.temp  <- cbind.data.frame(Set=j, Cluster=k, Score=mean(colMeans( cluster.data[row.names(cluster.data) %in% network.master.list,])))
      }else{##mod
        cluster.mean.temp  <- cbind.data.frame(Set=j, Cluster=k, Score=mean(cluster.data[[1]]) ) ##mod
      }##mod
      
      
      cluster.mean <- rbind(cluster.mean , cluster.mean.temp)
    }
    cluster.mean <- as.data.frame(cluster.mean[,c(3)])
    colnames(cluster.mean) <- j
    meganetwork.data <- cbind.data.frame(meganetwork.data, cluster.mean)
  }
  row.names(meganetwork.data)<- meganetwork.data$Cluster
  meganetwork.data$Cluster <-NULL
  
  ######Meganetwork collapse
  meganetwork.cor <- as.data.frame(cor(meganetwork.data ))
  diag(meganetwork.cor) <- 0
  sort(meganetwork.cor[[1]])
  ####Collapse Cutoff
  net.cor.cutoff <- 0.80
  
  
  network.cor.master<- NULL
  for (j in names(meganetwork.cor)){
    print(j)
    # j="KRT8"
    # row.names(network.correlation[network.correlation[[j]]>=0.80, ,drop=F])
    if (length(row.names(meganetwork.cor[meganetwork.cor[[j]]>=net.cor.cutoff, ,drop=F]))>0){
      network.correlation2<-cbind.data.frame(Genes=row.names(meganetwork.cor[meganetwork.cor[[j]]>=net.cor.cutoff, ,drop=F]),
                                             Corr=meganetwork.cor[meganetwork.cor[[j]]>=net.cor.cutoff,,drop=FALSE][[j]], Set=j)
      
      network.cor.master <- rbind(network.cor.master, network.correlation2)
    } else {
      print ("None")
    }
  }
  network.cor.master<-as.data.frame(network.cor.master)
  nrow(network.cor.master)
  length(unique(network.cor.master$Set))
  
  
  ####Networking again
  network.master.list <-network.cor.master
  network.master.list2 <-network.master.list[with(network.master.list, order(-Corr)),]
  network.master.list2 <- network.master.list2[!duplicated(network.master.list2$Gene),]
  network.master.list2$RegSet <- 0
  length(unique(network.master.list2$Set))
  n=1
  # gene="X"
  for(gene in (unique(network.master.list2$Set))){
    print(gene)
    if(unique(network.master.list2[network.master.list2$Set==gene,]$RegSet)==0){
      
      ##inital call
      network.master.list2[network.master.list2$Set==gene,]$RegSet <- n
      ## this part spreads the plague  1 edge
      for(gene2 in network.master.list2[network.master.list2$Set==gene,]$Gene){
        if(nrow(network.master.list2[network.master.list2$Set==gene2,])>0){
          network.master.list2[network.master.list2$Set==gene2,]$RegSet <- n
          ## second round 2 edges
          for(gene3 in network.master.list2[network.master.list2$Set==gene2,]$Gene){
            if(nrow(network.master.list2[network.master.list2$Set==gene3,])>0){
              network.master.list2[network.master.list2$Set==gene3,]$RegSet <- n
              ## third round 3 edges
              for(gene4 in network.master.list2[network.master.list2$Set==gene3,]$Gene){
                if(nrow(network.master.list2[network.master.list2$Set==gene4,])>0){
                  network.master.list2[network.master.list2$Set==gene4,]$RegSet <- n
                  for(gene5 in network.master.list2[network.master.list2$Set==gene4,]$Gene){
                    if(nrow(network.master.list2[network.master.list2$Set==gene5,])>0){
                      network.master.list2[network.master.list2$Set==gene5,]$RegSet <- n
                      ###5 edges
                      for(gene6 in network.master.list2[network.master.list2$Set==gene5,]$Gene){
                        if(nrow(network.master.list2[network.master.list2$Set==gene6,])>0){
                          network.master.list2[network.master.list2$Set==gene6,]$RegSet <- n
                        }
                      }# gene 6
                    }# if gene5
                  }#gene 5
                }#if gene4
              }#gene4
            }#if gene3
          }#gene3
          
        } #if gene2
        
      } #gene2 forloop 
      n=n+1
    } else {
      print("decided")
    }
  } #first loop
  length(unique(network.master.list2$RegSet))
  length(unique(network.master.list2$Genes))
  
  
  ######
  
  network.mega.regress <- network.cor.master2
  
  network.master.list2$RegSet2<-0
  # network.master.list2<- network.master.list2[with(network.master.list2, order(RegSet, Set)),]
  for(Group in unique(network.master.list2$RegSet)){
    
    network.master.list2[network.master.list2$RegSet==Group,]$RegSet2 <- network.master.list2[network.master.list2$RegSet==Group,]$Set[1]
  }
  
  network.mega.regress$FinalSet <- network.master.list2$RegSet2[match(network.mega.regress$RegSet,network.master.list2$Genes)]
  network.mega.regress[is.na(network.mega.regress$FinalSet),]$FinalSet <-network.mega.regress[is.na(network.mega.regress$FinalSet),]$RegSet
  
  network.size.table <- as.data.frame(table(network.mega.regress$FinalSet))
  network.mega.regress$Network.Size <- network.size.table$Freq[match(network.mega.regress$FinalSet, network.size.table$Var1)]
  
  
  
  ####Grade/Filter Networks
  network.mega.regress$Pct.Cells.Exp <- 0
  pct.test.df <- normalized.data[network.mega.regress$Genes, ]
  pct.test.df[pct.test.df>0]<- 1
  network.mega.regress$Pct.Cells.Exp <-rowSums(pct.test.df)/ncol(normalized.data)
  
  
  
  
  #### Final output prep
  network.mega.regress$FinalSet <- as.numeric(network.mega.regress$FinalSet)
  
  ###Checking network size and cutting off small networks
  network.size.table <- as.data.frame(table(network.mega.regress$FinalSet))
  network.mega.regress$Network.Size <- network.size.table$Freq[match(network.mega.regress$FinalSet, network.size.table$Var1)]
  network.mega.regress <- network.mega.regress[network.mega.regress$Network.Size>9,]
  
  
  ###Create network intercorrelation table to determine approximate gene-gene agreement
  network.size.table <- as.data.frame(table(network.mega.regress$FinalSet))
  network.size.table$InsideCorr=0
  network.size.table$OutsideCorr=0
  network.size.table$Avg.Pct.Exp <- 0
  network.size.table$CellAdjust <- 100
  network.size.table$SD.Out <- 0
  network.size.table$SD.In <- 0
  network.size.table$SD.All <- 0
  
  for(Program in unique(network.mega.regress$FinalSet)){
    markerlist4u <- network.mega.regress[network.mega.regress$FinalSet==Program,]$Genes
    markers<-cell.design
    markers$Score <- colMeans(normalized.data[row.names(normalized.data) %in% markerlist4u,])
    
    thresh.cutoff <- (max(markers$Score)*0.50)
    
    network.size.table[network.size.table$Var1==Program,]$Avg.Pct.Exp <- mean(network.mega.regress[network.mega.regress$FinalSet==Program,]$Pct.Cells.Exp)
    
    ##Sd
    cell.net <- normalized.data [row.names(normalized.data) %in% markerlist4u,]
    network.size.table[network.size.table$Var1==Program,]$SD.All <- mean(unlist(sd.obj))
    
    cell.net <- cell.net[,names(cell.net) %in% row.names(markers[markers$Score > thresh.cutoff,]),drop=FALSE ]
    
    if(ncol(cell.net)>1){
      cell.net <- cell.net[, colSums(cell.net)>0]
      cell.cor <-cor(cell.net)
      network.size.table[network.size.table$Var1==Program,]$InsideCorr<- mean(cell.cor)
      
      cell.net<-cell.net+0.0001
      sd.obj<-lapply(as.data.frame(t(cell.net/rowSums(cell.net))), sd)
      network.size.table[network.size.table$Var1==Program,]$SD.In <- mean(unlist(sd.obj))
      
      
    }else{ ###if only 1 cell qualifies, will just take top 10 cells?
      
      markers<- markers[with(markers, order(-Score)),]
      cell.net <- normalized.data [row.names(normalized.data) %in% markerlist4u,]
      cell.net <- cell.net[,names(cell.net) %in% row.names(markers[1:10,]),drop=FALSE ]
      cell.net <- cell.net[, colSums(cell.net)>0]
      
      ###Add a random pseudocount to break 0 distributions
      cell.net <- cell.net +(runif(length(cell.net),0.0001, 0.0002))
      
      cell.cor <-cor(cell.net)
      diag(cell.cor) <- 0
      network.size.table[network.size.table$Var1==Program,]$InsideCorr<- mean(cell.cor)
      
      test.df<-as.data.frame(table(colnames(cell.net)[apply(cell.net ,1,which.max)]) )
      test.df <- test.df[with(test.df, order(-Freq)),]
      network.size.table[network.size.table$Var1==Program,]$CellAdjust<- 100*test.df$Freq[1]/sum(test.df$Freq)
      
      sd.obj<-lapply(as.data.frame(t(cell.net/rowSums(cell.net))), sd)
      network.size.table[network.size.table$Var1==Program,]$SD.In<- mean(unlist(sd.obj))
    }
    
    table(names(cell.net) %in% row.names(markers[markers$Score > thresh.cutoff,]))
    
    
    cell.net <- normalized.data [row.names(normalized.data) %in% markerlist4u,]
    cell.net <- cell.net[,!names(cell.net) %in% row.names(markers[markers$Score > thresh.cutoff,]) ]
    cell.net <- cell.net[, colSums(cell.net)>0]
    cell.cor <-cor(cell.net)
    network.size.table[network.size.table$Var1==Program,]$OutsideCorr<- mean(cell.cor)
    
    cell.net<-cell.net+0.0001
    sd.obj<-lapply(as.data.frame(t(cell.net/rowSums(cell.net))), sd)
    network.size.table[network.size.table$Var1==Program,]$SD.Out <- mean(unlist(sd.obj))
    
    
  }
  network.size.table$DiffCorr <- network.size.table$InsideCorr - network.size.table$OutsideCorr
  
  network.size.table$Adjust.Corr <- network.size.table$DiffCorr * (network.size.table$Freq) * (network.size.table$CellAdjust)
  
  ##added log term to adjust.corr2 to soften pct exp scaling
  network.size.table$Adjust.Corr2 <- network.size.table$Adjust.Corr * (network.size.table$Avg.Pct.Exp)
  
  network.size.table<- network.size.table[with(network.size.table, order(-Adjust.Corr2)),]
  network.size.table$LikelyName <- c(1:nrow(network.size.table))
  
  
  network.mega.regress$DiffCorr <- network.size.table$DiffCorr[match(network.mega.regress$FinalSet, network.size.table$Var1)]
  network.mega.regress$Adjust.Corr <- network.size.table$Adjust.Corr2[match(network.mega.regress$FinalSet, network.size.table$Var1)]
  
  network.mega.regress <- network.mega.regress[with(network.mega.regress, order(-Adjust.Corr)),]
  
  
  network.mega.regress$TempSet <- 0
  temp.set=1
  for(Set.Number in unique(network.mega.regress$FinalSet) ){
    network.mega.regress[network.mega.regress$FinalSet==Set.Number,]$TempSet <- temp.set
    temp.set=temp.set+1
  }
  
  network.mega.regress$FinalSet <-  network.mega.regress$TempSet
  network.mega.regress$TempSet <- NULL
  
  network.mega.regress <- network.mega.regress[with(network.mega.regress, order(FinalSet)),]
  
  
  
  
  ###Quality judgement (may be skewed)
  network.mega.regress[is.na(network.mega.regress)] <- 0
  network.mega.regress$Net.Quality <- c("Noisy")
  network.mega.regress[network.mega.regress$Adjust.Corr> quantile(network.mega.regress$Adjust.Corr)[3] ,]$Net.Quality <- "Questionable"
  network.mega.regress[network.mega.regress$Adjust.Corr> quantile(network.mega.regress$Adjust.Corr)[4],]$Net.Quality <- "Good"
  
  
  
  ##Finish and write
  # write.table(network.mega.regress, paste0(outdir, drop, ".GEND.Networks.Output_Cor.",cutoff,".txt"),row.names=FALSE, quote=FALSE, sep="\t")

  
  
  # Program.Table <- network.mega.regress
  # Program.Table$Gene <- row.names(Program.Table)
  
  
  
 
  names(network.mega.regress)
  names(Program.Table)
  
  network.mega.regress$FinalSet <- paste0("B", network.mega.regress$FinalSet)
  
  Program.Table$FinalSet <- paste0("A", Program.Table$FinalSet)
  
  Program.Table<- rbind(Program.Table, network.mega.regress)
  
  
  ###This is the target output file for downstream analysis
  write.table(Program.Table, paste0(outdir, drop,".GEND.with.Secondary.Pass.down.to.",cutoff,".table.txt"),
              sep="\t", quote=F)
  
  
  print("GEND Finished!")
  
  
  
  #################
  ######## Phase 3 - Display
  
  markers<-cell.design
  
 
  ###Can also just read in
  # Program.Table <- read.table(paste0("C:/Terminal/R TSNE/Drop10 Stca plus Virus/2020-01-17/Networks with 0.25 corr 08 collapse final kmeans cluster scaled data/Drop7.MetaNetworks.w.2ndaryCollapse.Cor.0.25.txt"), header=1)
  
  
  Program.Table$Set<-Program.Table$FinalSet
  

  
  
  marker.names <- names(markers)
  for(Set in unique(Program.Table$Set)){
    genes.of.interest <- Program.Table[Program.Table$Set==Set,]$Genes
    if(nrow(normalized.data[row.names(normalized.data) %in% genes.of.interest,])>0){
      markers <- cbind(markers, Score=colMeans(normalized.data[row.names(normalized.data) %in% genes.of.interest,]))
    } else{
      markers <- cbind(markers, Score=0)
    }
    marker.names <-c(marker.names, paste0("Network.",Set))
  }
  names(markers)<-marker.names
  Program.Sets<-na.omit(marker.names[(ncol(cell.design)+1):length(marker.names)])
  
  for(Program in unique(Program.Sets)){
    Program.number <- sapply(strsplit(as.character(Program), "\\."), "[[",2)
    Program.quality <- unique(Program.Table[Program.Table$Set==Program.number,]$Net.Quality)
    
    print(Program)
    markers2<-markers[,c("X","Y",Program)]
    names(markers2)<-c("X","Y","Score")
    markers2 <- markers2[with(markers2, order(Score)),]
    ggplot(markers2)+
      # facet_wrap(~Type)+
      geom_point(data=markers[,c("X","Y")], aes(X,Y), color="grey93")+
      geom_point( aes(X,Y, color=Score) )+
      # ThemeOpts()+
      labs(color=Program)+
      ggtitle(paste0(drop, " ", Program, " with ", nrow(Program.Table[Program.Table$Set==Program.number,])," genes (",Program.quality,")" ))+
      scale_colour_gradientn(colours=c("lightcyan3", "gold", "orange","orangered"))
      ggsave(filename=paste0(outdir, drop, ".Programs.marker.TSNE.Prog.",Program,".png"),
             height=4, width=5)
  }
  
  ####Filters small networks out for display. 
  Program.Table <- Program.Table[Program.Table$Adjust.Corr>0.05,]
  write.table(Program.Table, paste0(outdir, drop, ".GEND.Networks.Output.Filtered_Cor.",cutoff,".txt"),row.names=FALSE, quote=FALSE, sep="\t")
  
  
  ###Use GEND output for heatmap
  cell.design$Cluster2 <- as.factor(kmeans(cell.design[,c("X","Y")], round(nrow(cell.design)/(num.of.cells*5)))$cluster)
  unique(cell.design$Cluster2)
  ggplot(cell.design)+
    geom_text(aes(X,Y,label=Cluster2, color=as.factor(Cluster2)), size=2)+
    geom_text(data=cell.design[!duplicated(cell.design$Cluster2),],aes(X,Y,label=paste0("C", Cluster2)), color="black")+
    # ThemeOpts()+
    theme(legend.position = "")+
    ggsave(paste0(outdir, drop, ".KClusters.For.GEND.Heatmap_UMAP.png"), height=10, width=12)
  markers<-cell.design
  tag.table <- NULL
  tag.table$rn <- row.names(normalized.data)
  for( Sample in sort(unique(cell.design$Cluster2))){
    print(Sample)
    tag.temp <-cbind.data.frame(Sample=rowMeans(normalized.data[,colnames(normalized.data) %in% row.names(markers[markers$Cluster2==Sample,]) ]) )
    tag.table <- cbind.data.frame(tag.table, tag.temp$Sample)
    names(tag.table)[ncol(tag.table)] <- paste0("C", Sample)
  }
  row.names(tag.table) <- tag.table$rn
  tag.table$rn <-NULL
  ###Expression df
  tag.table2<-NULL
  for(Signature in unique(Program.Table$FinalSet)){
    print(Signature)
    genes.of.interest<- Program.Table[Program.Table$FinalSet==Signature,]$Genes
    
    tag.temp<- colMeans(tag.table[row.names(tag.table) %in% genes.of.interest,])
    tag.temp <-data.frame(t(tag.temp))
    row.names(tag.temp)<-paste0("Net_",Signature)
    tag.table2 <- rbind(tag.table2,tag.temp) 
  }
  tag.table2
  library(pheatmap)
  # rowSums(tag.table3)
  tag.table3 <- tag.table2/rowSums(tag.table2)
  pheatmap(tag.table2, 
           cluster_cols =T, cluster_rows = F,
           main = "Network Expression on KClusters",
           # display_numbers = TRUE, number_format = "%.1f",
           filename=paste0(outdir, drop, ".KClusters.GEND.Networks.Heatmap.png"),
           height=10, width=10)
  library(RColorBrewer)
  pheatmap(tag.table3, 
           cluster_cols =T, cluster_rows = F,
           color= colorRampPalette(c("dodgerblue", "yellow", "orange", "red", "purple", "purple"))(50),
           # display_numbers = TRUE, number_format = "%.1f",
           main = "Row Normalized Network Expression on KClusters",
           filename=paste0(outdir, drop, ".KClusters.GEND.Networks_RowNorm.Heatmap.png"),
           height=10, width=10)
  
  print("Plotting finished!")
  
  
  ########################################
  ######## Create Outputs for metascape - Keeps order
  Program.Table <- Program.Table[Program.Table$Adjust.Corr>0.05,]
  
  Network.Source<- drop
  if(length(unique(Program.Table$FinalSet))>30){print("Too many networks, dividing into parts")
    Cut.Amount <- c(1:ceiling( length(unique(Program.Table$FinalSet))/30))
    Take.Amount <- 1
    Network.Cut.Table <- cbind.data.frame(Network=unique(Program.Table$FinalSet),Cut.Identity=0) 
    for(Cutn in c(1: max(Cut.Amount))){
      top.amount <- min( c(Take.Amount+round(nrow(Network.Cut.Table)/Cut.Amount), nrow(Network.Cut.Table) ) )
      Network.Cut.Table[c(Take.Amount:top.amount),]$Cut.Identity<-  Cutn
      Take.Amount <- Take.Amount+ round(nrow(Network.Cut.Table)/ max(Cut.Amount))+1
    }
    Program.Table$Cut.Identity <- Network.Cut.Table[match(Program.Table$FinalSet, Network.Cut.Table$Network),]$Cut.Identity
    for(Cut.Identity in unique(Program.Table$Cut.Identity)){
      print(paste0("part ", Cut.Identity))
      Metascape.Multi <- NULL
      for(Set in unique(Program.Table[Program.Table$Cut.Identity==Cut.Identity,]$FinalSet)){
        Program.Table2 <- Program.Table[Program.Table$FinalSet==Set,]
        Metascape.tmp<-cbind.data.frame( Name=paste0("Network.",Set),Genes=paste(as.character(Program.Table2$Genes), collapse=",")  )
        Metascape.Multi <- rbind(Metascape.Multi, Metascape.tmp)
      }
      print("done parsing")
      write.table(Metascape.Multi, paste0(outdir, drop,".",Network.Source,".Part",Cut.Identity,".Metanetworks.MultiGO.List.txt"), sep="\t", col.names =F, row.names = F,quote=F)
    }
  }else{
    ##One loop version
    print("Network Numbers Okay")
    Metascape.Multi <- NULL
    for(Set in unique(Program.Table$FinalSet)){
      Program.Table2 <- Program.Table[Program.Table$FinalSet==Set,]
      Metascape.tmp<-cbind.data.frame( Name=paste0("Network.",Set),Genes=paste(as.character(Program.Table2$Genes), collapse=",")  )
      Metascape.Multi <- rbind(Metascape.Multi, Metascape.tmp)
    }
    print("done full")
    write.table(Metascape.Multi, paste0(outdir, drop,".",Network.Source,".Metanetworks.MultiGO.List.txt"), sep="\t", col.names =F, row.names = F,quote=F)
  }
  
  print("Done!")
