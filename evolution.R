setwd("/media/_EXTend2018/RStudio/evolution/")
library(ggplot2)
library(rWikiPathways)
library(org.Hs.eg.db)
library(dplyr)
library(Rtsne)
library(wordcloud2)
library(org.Hs.eg.db)

pathway_plot <- function(pathway){
  colnames(pathway) <- c("Pathway_ID","FDR")
  pathway$logp <- -log2(pathway$FDR)
  pathway$size=pathway$logp/10
  len=length(pathway$Pathway_ID)
  require(ggplot2)
  p <- ggplot(pathway,aes(x=logp,y=seq(len,1,by=-1)))+
    geom_point(aes(size = size,colour="red"))+
    scale_y_continuous(labels = pathway$Pathway_ID,breaks = seq(len,1,by = -1))+
    labs(x="-logP",y="pathways",title="logP plot")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y=element_text(size=20),
          legend.text = element_text(size = 12))
  p
}

Read.GeneSets.db2 <- function (gs.db, thres.min = 2, thres.max = 2000) {
  ## read gmt files
  temp <- readLines(gs.db)
  temp <- strsplit(temp, '\t')
  temp.size.G <- sapply(temp, function(x) length(x)-2)
  
  ## filter gene sets according to size
  rm.idx <- which(temp.size.G < thres.min | temp.size.G > thres.max)
  if(length(rm.idx) > 0){
    temp <- temp[-rm.idx]
    temp.size.G <- temp.size.G[-rm.idx]
  }
  
  max.Ng <- length(temp)         ## number of signature sets
  temp.size.G <- sapply(temp, function(x) length(x)-2)
  max.size.G <- max(temp.size.G) ## maximal size
  
  gs <- lapply(temp, function(x)x[3:length(x)])
  gs.names <- sapply(temp, function(x)x[1])
  gs.desc <- sapply(temp, function(x)x[2])
  size.G <- temp.size.G
  names(gs) <- names(gs.names) <- names(gs.desc) <- names(size.G) <- gs.names
  
  
  return(list(N.gs = max.Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc,
              size.G = size.G, max.N.gs = max.Ng))
}
Entrez2Symbol <- function(x){
  my.entrez <- x
  hs <- org.Hs.eg.db
  map.table <- select(hs, 
                      keys = my.entrez,
                      columns = c("ENTREZID", "SYMBOL"),
                      keytype = "ENTREZID")
  return(map.table$SYMBOL)
}

African.raw <- read.table("African.txt",header = T,sep = "\t")
African.gene.freq <- data.frame(gene = row.names(African.raw),freq = apply(African.raw, 1, sum))
wordcloud2(African.gene.freq,shape = "circle")
African <- African.raw[apply(African.raw, 1, sum)>(0.05*ncol(African.raw)),]
African.loss <- apply(African, 2,sum)
African.genes <- row.names(African)
# write.table(as.matrix(row.names(African)),row.names = F,col.names = F,quote = F,file = "African.genes.txt")
summary(African.loss)

African.indi.genes.list <- apply(African.raw, 2, function(x){row.names(African.raw)[x>0]})
all.pathway.geneset <- Read.GeneSets.db2("wikipathways-20181010-gmt-Homo_sapiens.gmt.txt")$gs
names(all.pathway.geneset) <- sapply(names(all.pathway.geneset), function(x){unlist(strsplit(x,"%"))[1]})
all.pathway.geneset.symbol <- sapply(all.pathway.geneset, function(x){Entrez2Symbol(x)})


pathway <- read.table("pathway.txt",header = F,sep = "\t")
pathway_plot(pathway)

Ashkenazi_Jewish.raw <- read.table("Ashkenazi Jewish.txt",header = TRUE,sep = "\t")
Ashkenazi_Jewish <- Ashkenazi_Jewish.raw[apply(Ashkenazi_Jewish.raw, 1, sum)>(0.05*ncol(Ashkenazi_Jewish.raw)),]
Ashkenazi_Jewish.gene.freq <- data.frame(gene = row.names(Ashkenazi_Jewish.raw),freq = apply(Ashkenazi_Jewish.raw, 1, sum))
wordcloud2(Ashkenazi_Jewish.gene.freq,shape = "circle")
Ashkenazi_Jewish.loss <- apply(Ashkenazi_Jewish, 2,sum)
Ashkenazi_Jewish.genes <- row.names(Ashkenazi_Jewish)
# write.table(as.matrix(row.names(Ashkenazi_Jewish)),row.names = F,col.names = F,quote = F,file = "Ashkenazi_Jewish.genes.txt")
summary(Ashkenazi_Jewish.loss)


East_Asian.raw <- read.table("East Asian.txt",header = T,sep = "\t")
East_Asian <- East_Asian.raw[apply(East_Asian.raw, 1, sum)>(0.05*ncol(East_Asian.raw)),]
East_Asian.gene.freq <- data.frame(gene = row.names(East_Asian.raw),freq = apply(East_Asian.raw, 1, sum))
wordcloud2(East_Asian.gene.freq,shape = "circle")
East_Asian.loss <- apply(East_Asian, 2,sum)
East_Asian.genes <- row.names(East_Asian)
# write.table(as.matrix(row.names(East_Asian)),row.names = F,col.names = F,quote = F,file = "East_Asian.genes.txt")
summary(East_Asian.loss)


European.Finnish.raw <- read.table("European (Finnish).txt",header = T,sep = "\t")
European.Finnish <- European.Finnish.raw[apply(European.Finnish.raw, 1, sum)>(0.05*ncol(European.Finnish.raw)),]
European.Finnish.gene.freq <- data.frame(gene = row.names(European.Finnish.raw),freq = apply(European.Finnish.raw, 1, sum))
wordcloud2(European.Finnish.gene.freq,shape = "circle")
European.Finnish.loss <- apply(European.Finnish, 2,sum)
European.Finnish.genes <- row.names(European.Finnish)
# write.table(as.matrix(row.names(European.Finnish)),row.names = F,col.names = F,quote = F,file = "European.Finnish.genes.txt")
summary(European.Finnish.loss)


European.Non.Finnish.raw <- read.table("European (Non-Finnish).txt",header = T,sep = "\t")
European.Non.Finnish <- European.Non.Finnish.raw[apply(European.Non.Finnish.raw, 1, sum)>(0.05*ncol(European.Non.Finnish.raw)),]
European.Non.Finnish.gene.freq <- data.frame(gene = row.names(European.Non.Finnish.raw),freq = apply(European.Non.Finnish.raw, 1, sum))
wordcloud2(European.Non.Finnish.gene.freq,shape = "circle")
European.Non.Finnish.loss <- apply(European.Non.Finnish, 2,sum)
European.Non.Finnish.genes <- row.names(European.Non.Finnish)
# write.table(as.matrix(row.names(European.Non.Finnish)),row.names = F,col.names = F,quote = F,file = "European.Non.Finnish.genes.txt")
summary(European.Non.Finnish.loss)


Latino.raw <- read.table("Latino.txt",header = T,sep = "\t")
Latino <- Latino.raw[apply(Latino.raw, 1, sum)>(0.05*ncol(Latino.raw)),]
Latino.gene.freq <- data.frame(gene = row.names(Latino.raw),freq = apply(Latino.raw, 1, sum))
wordcloud2(Latino.gene.freq,shape = "circle")
Latino.loss <- apply(Latino, 2,sum)
Latino.genes <- row.names(Latino)
# write.table(as.matrix(row.names(Latino)),row.names = F,col.names = F,quote = F,file = "Latino.genes.txt")
summary(Latino.loss)


South_Asian.raw <- read.table("South Asian.txt",header = T,sep = "\t")
South_Asian <- South_Asian.raw[apply(South_Asian.raw, 1, sum)>(0.05*ncol(South_Asian.raw)),]
South_Asian.gene.freq <- data.frame(gene = row.names(South_Asian.raw),freq = apply(South_Asian.raw, 1, sum))
wordcloud2(South_Asian.gene.freq,shape = "circle")
South_Asian.loss <- apply(South_Asian, 2,sum)
South_Asian.genes <- row.names(South_Asian)
# write.table(as.matrix(row.names(South_Asian)),row.names = F,col.names = F,quote = F,file = "South_Asian.genes.txt")
summary(South_Asian.loss)


all.mean <- c(mean(African.loss),mean(Ashkenazi_Jewish.loss),mean(East_Asian.loss),mean(European.Finnish.loss),
              mean(European.Non.Finnish.loss),mean(Latino.loss),mean(South_Asian.loss))
names(all.mean) <- c("African","Ashkenazi_Jewish","East_Asian","European.Finnish","European.Non.Finnish",
                     "Latino","South_Asian")
all.mean <- round(all.mean,2)
data.mean <- data.frame(name = names(all.mean),mean = all.mean)
ggplot(data.mean, aes(name,mean,fill=name)) + 
  geom_bar(stat="identity", width = 0.5,position ="dodge")+
  geom_text(aes(x=name,y=mean+mean/3,label=mean),position=position_dodge(.7),size=8)+
  scale_y_continuous(trans="log2")+
  labs(x="Population",y="Mean loss genes")+
  theme(strip.text.x = element_text(face="bold",size = 12),
        strip.text.y = element_text(face="bold",size = 12),
        axis.text.x=element_text(face="bold",size=12,angle = 0),
        axis.text.y=element_text(face="bold",size=12),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15),
        legend.text = element_text(size = 12))


####unique genes
African.unique <- setdiff(African.genes,c(Ashkenazi_Jewish.genes,East_Asian.genes,European.Finnish.genes,
                                          European.Non.Finnish.genes,Latino.genes,South_Asian.genes))
Ashkenazi_Jewish.unique <- setdiff(Ashkenazi_Jewish.genes,c(African.genes,East_Asian.genes,European.Finnish.genes,
                                                            European.Non.Finnish.genes,Latino.genes,South_Asian.genes))
East_Asian.unique <- setdiff(East_Asian.genes,c(African.genes,Ashkenazi_Jewish.genes,European.Finnish.genes,
                                                            European.Non.Finnish.genes,Latino.genes,South_Asian.genes))
European.Finnish.unique <- setdiff(European.Finnish.genes,c(African.genes,Ashkenazi_Jewish.genes,East_Asian.genes,
                                                European.Non.Finnish.genes,Latino.genes,South_Asian.genes))
European.Non.Finnish.unique <- setdiff(European.Non.Finnish.genes,c(African.genes,Ashkenazi_Jewish.genes,East_Asian.genes,
                                                                    European.Finnish.genes,Latino.genes,South_Asian.genes))
Latino.unique <- setdiff(Latino.genes,c(African.genes,Ashkenazi_Jewish.genes,East_Asian.genes,
                                        European.Finnish.genes,European.Non.Finnish.genes,South_Asian.genes))
South_Asian.unique <- setdiff(South_Asian.genes,c(African.genes,Ashkenazi_Jewish.genes,East_Asian.genes,
                                                  European.Finnish.genes,European.Non.Finnish.genes,Latino.genes))

# write.table(as.matrix(African.unique),row.names = F,col.names = F,quote = F,file = "African.unique.txt")
# write.table(as.matrix(Ashkenazi_Jewish.unique),row.names = F,col.names = F,quote = F,file = "Ashkenazi_Jewish.unique.txt")
# write.table(as.matrix(East_Asian.unique),row.names = F,col.names = F,quote = F,file = "East_Asian.unique.txt")
# write.table(as.matrix(European.Finnish.unique),row.names = F,col.names = F,quote = F,file = "European.Finnish.unique.txt")
# write.table(as.matrix(European.Non.Finnish.unique),row.names = F,col.names = F,quote = F,file = "European.Non.Finnish.unique.txt")
# write.table(as.matrix(Latino.unique),row.names = F,col.names = F,quote = F,file = "Latino.unique.txt")
# write.table(as.matrix(South_Asian.unique),row.names = F,col.names = F,quote = F,file = "South_Asian.unique.txt")


all.uniq <- c(length(African.unique),length(Ashkenazi_Jewish.unique),length(East_Asian.unique),length(European.Finnish.unique),
              length(European.Non.Finnish.unique),length(Latino.unique),length(South_Asian.unique))
all.uniq.mean <- c(length(African.unique)/ncol(African.raw),length(Ashkenazi_Jewish.unique)/ncol(Ashkenazi_Jewish.raw),length(East_Asian.unique)/ncol(East_Asian.raw),
                   length(European.Finnish.unique)/ncol(European.Finnish.raw),length(European.Non.Finnish.unique)/ncol(European.Non.Finnish.raw),length(Latino.unique)/ncol(Latino.raw),
                   length(South_Asian.unique)/ncol(South_Asian.raw))
all.uniq.mean <- round(all.uniq.mean,2)
data.uniq <- data.frame(name = names(all.mean),uniq = all.uniq,uniqmean = all.uniq.mean)
data.uniq <- dplyr::mutate(data.uniq,loguniq= log(uniq+1),loguniqmean = log(uniqmean+1))

ggplot(data.uniq, aes(name,loguniqmean,fill=name)) + 
  geom_bar(stat="identity", width = 0.5,position ="dodge")+
  geom_text(aes(x=name,y=loguniqmean+loguniqmean/15,label=uniqmean),position=position_dodge(.7),size=8)+
  labs(x="Population",y="Unique loss genes")+
  theme(strip.text.x = element_text(face="bold",size = 12),
        strip.text.y = element_text(face="bold",size = 12),
        axis.text.x=element_text(face="bold",size=12,angle = 0),
        axis.text.y=element_text(face="bold",size=12),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15),
        legend.text = element_text(size = 12))

African$gene <- row.names(African)
Ashkenazi_Jewish$gene <- row.names(Ashkenazi_Jewish)
East_Asian$gene <- row.names(East_Asian)
European.Finnish$gene <- row.names(European.Finnish)
European.Non.Finnish$gene <- row.names(European.Non.Finnish)
Latino$gene <- row.names(Latino)
South_Asian$gene <- row.names(South_Asian)

all.African.genes <- dplyr::left_join(African,Ashkenazi_Jewish,by="gene") %>% left_join(East_Asian,by="gene") %>% left_join(European.Finnish,by="gene") %>% 
  left_join(European.Non.Finnish,by="gene") %>% left_join(Latino,by="gene") %>% left_join(South_Asian,by="gene")
all.African.genes[is.na(all.African.genes)] <- 0
all.African.chisq <- t(rbind(apply(as.matrix(African[,-ncol(African)]), 1, table),
                             apply(as.matrix(all.African.genes[,-ncol(African)]), 1, table)))
all.African.chisq.pval <- apply(all.African.chisq, 1, function(x){tab=matrix(c(x[2:1],x[4:3]-x[2:1]),2);chisq.test(tab)$p.value})
all.African.chisq.padj <- p.adjust(all.African.chisq.pval)
table(all.African.chisq.padj<0.05)

all.Ashkenazi_Jewish.genes <- dplyr::left_join(Ashkenazi_Jewish,African,by="gene") %>% left_join(East_Asian,by="gene") %>% left_join(European.Finnish,by="gene") %>% 
  left_join(European.Non.Finnish,by="gene") %>% left_join(Latino,by="gene") %>% left_join(South_Asian,by="gene")
all.Ashkenazi_Jewish.genes[is.na(all.Ashkenazi_Jewish.genes)] <- 0
all.Ashkenazi_Jewish.chisq <- t(rbind(apply(as.matrix(Ashkenazi_Jewish[,-ncol(Ashkenazi_Jewish)]), 1, table),
                                      apply(as.matrix(all.Ashkenazi_Jewish.genes[,-ncol(Ashkenazi_Jewish)]), 1, table)))
all.Ashkenazi_Jewish.chisq.pval <- apply(all.Ashkenazi_Jewish.chisq, 1, function(x){tab=matrix(c(x[2:1],x[4:3]-x[2:1]),2);chisq.test(tab)$p.value})
all.Ashkenazi_Jewish.chisq.padj <- p.adjust(all.Ashkenazi_Jewish.chisq.pval)
table(all.Ashkenazi_Jewish.chisq.padj<0.05)

all.East_Asian.genes <- dplyr::left_join(East_Asian,African,by="gene") %>% left_join(Ashkenazi_Jewish,by="gene") %>% left_join(European.Finnish,by="gene") %>% 
  left_join(European.Non.Finnish,by="gene") %>% left_join(Latino,by="gene") %>% left_join(South_Asian,by="gene")
all.East_Asian.genes[is.na(all.East_Asian.genes)] <- 0
all.East_Asian.chisq <- t(rbind(apply(as.matrix(East_Asian[,-ncol(East_Asian)]), 1, table),
                                      apply(as.matrix(all.East_Asian.genes[,-ncol(East_Asian)]), 1, table)))
all.East_Asian.chisq.pval <- apply(all.East_Asian.chisq, 1, function(x){tab=matrix(c(x[2:1],x[4:3]-x[2:1]),2);chisq.test(tab)$p.value})
all.East_Asian.chisq.padj <- p.adjust(all.East_Asian.chisq.pval)
table(all.East_Asian.chisq.padj<0.05)

all.European.Finnish.genes <- dplyr::left_join(European.Finnish,African,by="gene") %>% left_join(Ashkenazi_Jewish,by="gene") %>% left_join(East_Asian,by="gene") %>% 
  left_join(European.Non.Finnish,by="gene") %>% left_join(Latino,by="gene") %>% left_join(South_Asian,by="gene")
all.European.Finnish.genes[is.na(all.European.Finnish.genes)] <- 0
all.European.Finnish.chisq <- t(rbind(apply(as.matrix(European.Finnish[,-ncol(European.Finnish)]), 1, table),
                                apply(as.matrix(all.European.Finnish.genes[,-ncol(European.Finnish)]), 1, table)))
all.European.Finnish.chisq.pval <- apply(all.European.Finnish.chisq, 1, function(x){tab=matrix(c(x[2:1],x[4:3]-x[2:1]),2);chisq.test(tab)$p.value})
all.European.Finnish.chisq.padj <- p.adjust(all.European.Finnish.chisq.pval)
table(all.European.Finnish.chisq.padj<0.05)

all.European.Non.Finnish.genes <- dplyr::left_join(European.Non.Finnish,African,by="gene") %>% left_join(Ashkenazi_Jewish,by="gene") %>% left_join(East_Asian,by="gene") %>% 
  left_join(European.Finnish,by="gene") %>% left_join(Latino,by="gene") %>% left_join(South_Asian,by="gene")
all.European.Non.Finnish.genes[is.na(all.European.Non.Finnish.genes)] <- 0
all.European.Non.Finnish.chisq <- t(rbind(apply(as.matrix(European.Non.Finnish[,-ncol(European.Non.Finnish)]), 1, table),
                                      apply(as.matrix(all.European.Non.Finnish.genes[,-ncol(European.Non.Finnish)]), 1, table)))
all.European.Non.Finnish.chisq.pval <- apply(all.European.Non.Finnish.chisq, 1, function(x){tab=matrix(c(x[2:1],x[4:3]-x[2:1]),2);chisq.test(tab)$p.value})
all.European.Non.Finnish.chisq.padj <- p.adjust(all.European.Non.Finnish.chisq.pval)
table(all.European.Non.Finnish.chisq.padj<0.05)

all.Latino.genes <- dplyr::left_join(Latino,African,by="gene") %>% left_join(Ashkenazi_Jewish,by="gene") %>% left_join(East_Asian,by="gene") %>% 
  left_join(European.Finnish,by="gene") %>% left_join(European.Non.Finnish,by="gene") %>% left_join(South_Asian,by="gene")
all.Latino.genes[is.na(all.Latino.genes)] <- 0
all.Latino.chisq <- t(rbind(apply(as.matrix(Latino[,-ncol(Latino)]), 1, table),
                            apply(as.matrix(all.Latino.genes[,-ncol(Latino)]), 1, table)))
all.Latino.chisq.pval <- apply(all.Latino.chisq, 1, function(x){tab=matrix(c(x[2:1],x[4:3]-x[2:1]),2);chisq.test(tab)$p.value})
all.Latino.chisq.padj <- p.adjust(all.Latino.chisq.pval)
table(all.Latino.chisq.padj<0.05)

all.South_Asian.genes <- dplyr::left_join(South_Asian,African,by="gene") %>% left_join(Ashkenazi_Jewish,by="gene") %>% left_join(East_Asian,by="gene") %>% 
  left_join(European.Finnish,by="gene") %>% left_join(European.Non.Finnish,by="gene") %>% left_join(Latino,by="gene")
all.South_Asian.genes[is.na(all.South_Asian.genes)] <- 0
all.South_Asian.chisq <- t(rbind(apply(as.matrix(South_Asian[,-ncol(South_Asian)]), 1, table),
                            apply(as.matrix(all.South_Asian.genes[,-ncol(South_Asian)]), 1, table)))
all.South_Asian.chisq.pval <- apply(all.South_Asian.chisq, 1, function(x){tab=matrix(c(x[2:1],x[4:3]-x[2:1]),2);chisq.test(tab)$p.value})
all.South_Asian.chisq.padj <- p.adjust(all.South_Asian.chisq.pval)
table(all.South_Asian.chisq.padj<0.05)

all.sign.num <- c(table(all.African.chisq.padj<0.05)[2],table(all.Ashkenazi_Jewish.chisq.padj<0.05)[2],table(all.East_Asian.chisq.padj<0.05)[2],
                  table(all.European.Finnish.chisq.padj<0.05)[2],table(all.European.Non.Finnish.chisq.padj<0.05)[2],table(all.Latino.chisq.padj<0.05)[2],
                  table(all.South_Asian.chisq.padj<0.05)[2])
data.sign <- data.frame(name = names(all.mean),sign = all.sign.num)
data.sign <- dplyr::mutate(data.sign,logsign = log(sign+1))
ggplot(data.sign, aes(name,logsign,fill=name)) + 
  geom_bar(stat="identity", width = 0.5,position ="dodge")+
  geom_text(aes(x=name,y=logsign+logsign/15,label=sign),position=position_dodge(.7),size=8)+
  scale_y_continuous(trans="log2")+
  labs(x="Population",y="Significant loss genes")+
  theme(strip.text.x = element_text(face="bold",size = 12),
        strip.text.y = element_text(face="bold",size = 12),
        axis.text.x=element_text(face="bold",size=12,angle = 0),
        axis.text.y=element_text(face="bold",size=12),
        axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15),
        legend.text = element_text(size = 12))

MutualExclusivity <- function(Part.mutual.exclusity.mat){
  E <- t(Part.mutual.exclusity.mat)
  g <- colnames(E)
  #s <- intersect(colnames(E),as.character(specific.gene.list[[2]]))
  resu <- NULL
  for(i in 1:length(g)){
    for(j in i:length(g)){
      if(i!=j){
        f <- fisher.test(rbind(E[,i],E[,j]))
        resu <- rbind(resu,cbind(geneA=g[i],geneB=g[j],oddsRatio=f$estimate,pvalue=f$p.value))
      }
    }
  }
  # return(resu)
  # some formatting
  res.new <- as.data.frame(resu)
  res.new$oddsRatio <- as.numeric(as.character(res.new$oddsRatio))
  res.new$pvalue <- as.numeric(as.character(res.new$pvalue))
  row.names(res.new) <- NULL
  # use p.adjust to correct for multi testing using a FDR
  res2 <- cbind(res.new,fdr=p.adjust(res.new$pvalue,"fdr"))
  # change the FDR in labels for plotting
  res2$stars <- cut(res2$fdr,
                    breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
                    label=c("***", "**", "*", ""))
  res2$events <- cut(res2$oddsRatio,
                     breaks=c(0, 0.1, 0.5, 2, 10, Inf),
                     label=c("SE", "TE", "NO", "TC","SC"),
                     include.lowest = T)
  return(res2)
}

MutualExclusivity_plot <- function(Part.mutual.exclusity.mat,cc){
  E <- t(Part.mutual.exclusity.mat)
  g <- colnames(E)
  #s <- intersect(colnames(E),as.character(specific.gene.list[[2]]))
  resu <- NULL
  for(i in 1:length(g)){
    for(j in i:length(g)){
      if(i!=j){
        f <- fisher.test(rbind(E[,i],E[,j]))
        resu <- rbind(resu,cbind(geneA=g[i],geneB=g[j],oddsRatio=f$estimate,pvalue=f$p.value))
      }
    }
  }
  # return(resu)
  # some formatting
  res.new <- as.data.frame(resu)
  res.new$oddsRatio <- as.numeric(as.character(res.new$oddsRatio))
  res.new$pvalue <- as.numeric(as.character(res.new$pvalue))
  row.names(res.new) <- NULL
  # use p.adjust to correct for multi testing using a FDR
  res2 <- cbind(res.new,fdr=p.adjust(res.new$pvalue,"fdr"))
  # change the FDR in labels for plotting
  res2$stars <- cut(res2$fdr,
                    breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
                    label=c("***", "**", "*", ""))
  res2$events <- cut(res2$oddsRatio,
                     breaks=c(0, 0.1, 0.5, 2, 10, Inf),
                     label=c("SE", "TE", "NO", "TC","SC"),
                     include.lowest = T)
  # return(res2)
  # plot with ggplot 2
  require(ggplot2)
  require(cowplot) # not necessary but the plot is nicer
  p <- ggplot(res2, aes(geneA, geneB)) +
    geom_tile(aes(fill = events),colour = "white") +
    geom_text(aes(label=stars), color="black", size=5) +
    xlab("Gene A") +
    ylab("Gene B") +
    scale_fill_brewer(palette = "Set3") +
    ggtitle(paste(cc,"Mutual exclusivity",sep = "_")) +
    theme(legend.key.size = unit(0.4, "cm"),
          axis.text.y=element_text(size=9,face = "bold"),
          axis.text.x=element_text(angle=70,size=9,hjust=1,face = "bold",
                                   margin = margin(0,0,0,20)))+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5,
                                    margin = margin(l=100,r=50,t=10,b=10),
                                    face = "bold", colour = "black",size=10))
  p
}
MutualExclusivity_plot(all.Ashkenazi_Jewish.chisq[,1:2],"Ashkenazi Jewish")
African.mutual.res <- MutualExclusivity(all.African.chisq[,1:2])
Ashkenazi_Jewish.mutual.res <- MutualExclusivity(all.Ashkenazi_Jewish.chisq[,1:2])
East_Asian.mutual.res <- MutualExclusivity(all.East_Asian.chisq[,1:2])
European.Finnish.mutual.res <- MutualExclusivity(all.European.Finnish.chisq[,1:2])
European.Non.Finnish.mutual.res <- MutualExclusivity(all.European.Non.Finnish.chisq[,1:2])
Latino.mutual.res <- MutualExclusivity(all.Latino.chisq[,1:2])
South_Asian.mutual.res <- MutualExclusivity(all.South_Asian.chisq[,1:2])


all.loss.mat <- full_join(African,Ashkenazi_Jewish,by="gene") %>% full_join(East_Asian,by="gene") %>% full_join(European.Finnish,by="gene") %>% 
  full_join(European.Non.Finnish,by="gene") %>% full_join(Latino,by="gene") %>% full_join(South_Asian,by="gene")
all.loss.mat[is.na(all.loss.mat)] <- 0
all.loss.mat <- all.loss.mat[,-ncol(African)]
tsne_model <- Rtsne(as.matrix(t(all.loss.mat)), check_duplicates=FALSE, pca=TRUE, perplexity=20,
                    theta=0.5, dims=2)

## getting the two dimension matrix
d_tsne <- as.data.frame(tsne_model$Y) %>% mutate(Individual = colnames(all.loss.mat),
                                                 Population = rep(c("African","Ashkenazi Jewish","East Asian","European.Finnish",
                                                                "European.Non.Finnish","Latino","South Asian"),
                                                              c(ncol(African.raw),ncol(Ashkenazi_Jewish.raw),ncol(East_Asian.raw),ncol(European.Finnish.raw),
                                                                ncol(European.Non.Finnish.raw),ncol(Latino.raw),ncol(South_Asian.raw))))

# options(repr.plot.width=5, repr.plot.height=3.5)
ggplot(d_tsne, aes(x=V1, y=V2)) +  
  geom_point(size=2.5,aes(colour = Population)) +
  xlab("") + ylab("") +
  labs(title = "t-SNE") +
  theme_light(base_size=10) +
  theme(plot.title = element_text( face="bold",hjust = 0.5,size = 20),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 20))


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rWikiPathways", version = "3.8")
                    
                    
                    
