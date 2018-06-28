##################################################################################################################
##################################################################################################################
# if you lost your .RData got to line 253 #### maybe bug here!!!!!!
##################################################################################################################
##################################################################################################################
source ("http://www.bioconductor.org/biocLite.R")
setwd("/Users/Giulia/Desktop/WP1/")
library(openxlsx)
library(matrixStats)
library(RColorBrewer)
library(ggplot2)
#GSK3B P49841
#GSK3A P49840
#LDHA P00338
#PSMB8 P28062 (from pandey article The draft human proteome)
# note: If you do not have one of these packages installed, install by source("https://bioconductor.org/biocLite.R") and then biocLite("packageName")

###### import all Expr data
ATCCcellExpr <- read.xlsx(xlsxFile ="ExprProt/ATCC_celllines?original PrOEF-150320-OPL1022-TLL-Expression.Pancreatic.Cell.Lines(exp 46).xlsx",sheet=5)
PDX10Expr <- read.xlsx(xlsxFile = "ExprProt/10PDX_PrOEF-150717-OPL1022-TLL-Protein.Profiling.PDAC.PDX(AMC),human.vs.mouse.contribution(exp 53)-tp160411.xlsx",sheet = 7)
PDX5CellLinesExpr <- read.xlsx(xlsxFile = "ExprProt/5PDXandPRIMARYCELLLINESProtein.Profiling.PDX.models&Cell.Lines.E_Giovannetti(exp49)-tp171026.xlsx",sheet = 7)
ClinicalExpr <- read.xlsx("../LCM_Nov/Final completeLCM cox regression data version 3-tp171130TL11122017.xlsx",sheet = 4)
ClinicalExpr5gel <- read.xlsx("/Users/Giulia/Desktop/INKA/input_files/Prot.Expr./ClinicalSamples/PrOEF-180327-OPL1022-TRAS_TLL-Protein.Profiling-Human.PDAC.Tumour.Tissues.xlsx",sheet = 6)

###### import all PHOSPHO data
ATCCcellPHO <- read.xlsx(xlsxFile = "TiOX/ATCCcelllineTLL PrOEF-150501-sOPL1022-TLL-TiOx.Pancreatic.Cell.Lines(exp 46.part2.repeat).xlsx",sheet = 8)
PDX10PHO <- read.xlsx(xlsxFile = "TiOX/SUBPDXtllPrOEF-151007-1022-TLL-TiOx-Pancreatic.Cancer.PDX.Models.AMC(exp65)-tp151110.xlsx",sheet = 8)
PDX5CellLinesPHO <- read.xlsx(xlsxFile = "TiOX/5PDX5PRIMARYPrOEF-150507-OPL1022-TLL-TiOx PDX models E_Giovanetti(exp50) TLL22-6-2015-tp150914.xlsx",sheet = 8)

##################################################################################################################
                                          #PROTEIN EXPRESSION#
##################################################################################################################
###################################### ATCC cell lines ########################################
colnames(ATCCcellExpr)
ATCCexprLFQ <- ATCCcellExpr[,26:39]
ATCCexprLFQ$`HPDE(Epithelial)`=NULL
ATCCexprLFQ$`LFQ.intensity.(J).HPAF-II.Repl2`=NULL
ATCCexprLFQ$`LFQ.intensity.(K).HPAF-II.Repl3`=NULL
ATCCexprLFQ$`LFQ.intensity.(N).HPDE(Epithelial)`=NULL
cnames <- sapply(strsplit(colnames(ATCCexprLFQ),"[.]"), `[`, 4)
colnames(ATCCexprLFQ) <- cnames
ATCCexprLFQ$"Protein" <- sapply(strsplit(ATCCcellExpr$Protein.IDs.short,";"), `[`, 1)
ATCCexprLFQ$"Gene" <- sapply(strsplit(ATCCcellExpr$Gene.name,";"), `[`, 1)



################################################# # PDX 10 ##################################################################
                                           
colnames(PDX10Expr)
PDX10ExprLFQ <- PDX10Expr[,23:33]
cnames <- sapply(strsplit(colnames(PDX10ExprLFQ),"[.]"), `[`, 6)
colnames(PDX10ExprLFQ) <- cnames
PDX10ExprLFQ$"Protein" <- PDX10Expr$First.accession
PDX10ExprLFQ$"Gene" <- sapply(strsplit(PDX10Expr$Gene.names,";"), `[`, 1)



####################################### PDX 5 and Primary Cell Lines #############################################################################
                                    
colnames(PDX5CellLinesExpr)
PDX5CellLinesExprLFQ <- PDX5CellLinesExpr[,30:47]


colnames(PDX5CellLinesExprLFQ) <- c("PDX_PDAC-1","PDX_PDAC-2","PDX_PDAC-3","PDX_PDAC-4","PDX_PDAC-5","Tx=short","Tx=long","Tx=long","Tx=long","HPAF-II,","PDAC-1","PDAC-2","PDAC-3","PDAC-4","PDAC-5","pp161","pp391","HPAF-II,")
PDX5CellLinesExprLFQ$"Protein" <- PDX5CellLinesExpr$Protein.IDs.short
PDX5CellLinesExprLFQ$"Gene" <- sapply(strsplit(PDX5CellLinesExpr$Gene.name,";"), `[`, 1)

####################################### Clinical Samples 5 gel band #################################################################
colnames(ClinicalExpr5gel)
ClinicalExpr5gelLFQ <- ClinicalExpr5gel[,55:74]
cnames <- sapply(strsplit(colnames(ClinicalExpr5gelLFQ),"[.]"), `[`, 3)
colnames(ClinicalExpr5gelLFQ) <- cnames
ClinicalExpr5gelLFQ$"Protein" <- ClinicalExpr5gel$First.accession
ClinicalExpr5gelLFQ$"Gene" <- sapply(strsplit(ClinicalExpr5gel$Gene.names,";"), `[`, 1)

##################################################################################################################
                                                #PHOSPHOPROTEIN#
##################################################################################################################

############################################### ATCC cell lines ########################################
colnames(ATCCcellPHO)

ATCCcellPHOLFQ <- ATCCcellPHO[,27:40]
cnames <- sapply(strsplit(colnames(ATCCcellPHOLFQ),"[.]"), `[`, 5)
colnames(ATCCcellPHOLFQ) <- cnames
ATCCcellPHOLFQ$"Protein" <- ATCCcellPHO$Proteins
ATCCcellPHOLFQ$"Gene" <- ATCCcellPHO$Gene.Names
###################################### 10 PDX ########################################

colnames(PDX10PHO)
PDX10PHOLFQ <- PDX10PHO[,33:53]
cnames <- sapply(strsplit(colnames(PDX10PHOLFQ),"[.]"), `[`, 6)
colnames(PDX10PHOLFQ) <- cnames
PDX10PHOLFQ$"Protein" <- PDX10PHO$Proteins
PDX10PHOLFQ$"Gene" <- PDX10PHO$Gene.Names


###################################### 5 PDX AND PRIMARY CELL LINES ########################################
colnames(PDX5CellLinesPHO)
PDX5CellLinesPHOLFQ <- PDX5CellLinesPHO[,36:59]
cnames <- sapply(strsplit(colnames(PDX5CellLinesPHOLFQ),"[.]"), `[`, 6)
PDX5CellLinesPHOLFQ$"Protein" <- PDX5CellLinesPHO$Proteins
PDX5CellLinesPHOLFQ$"Gene" <- PDX5CellLinesPHO$Gene.Names


##################################################################################################################

##########################################        PREPARATION TO ComBat      ####################################################

##################################################################################################################

# include Human Tumours

ClinicalExprLFQ <-ClinicalExpr[grep("LFQ.intensity",colnames(ClinicalExpr))] 
ClinicalExprBulk <- ClinicalExprLFQ[grep("PDAC",colnames(ClinicalExprLFQ))]
ClinicalExprBulk$`LFQ.intensity.(NN).ExpII.33_PDAC_21_2014.BioRepl2`=NULL
ClinicalExprBulk$`LFQ.intensity.(PP).ExpII.33_PDAC_26_2014.TechRepl2`=NULL
ClinicalExprBulk$"Protein" <- ClinicalExpr$First.accession
ClinicalExprBulk$"Gene" <- ClinicalExpr$`HGNC+.Gene.Names`
ClinicalExprBulk$Protein <- sapply(strsplit(ClinicalExprBulk$Protein,";"), `[`, 1)
ClinicalExprBulk$Gene <- sapply(strsplit(ClinicalExprBulk$Gene,";"), `[`, 1)


######################################## primary analysis on TiOX ###############################

## CHECK TOT. NUMBER OF PROTEIN IN TIOX
ATCCcellPHOLFQ$Protein <- sapply(strsplit(ATCCcellPHOLFQ$Protein,";"), `[`, 1)
TotalATCCphosProtein <- length(unique(ATCCcellPHOLFQ$Protein))
PDX10PHOLFQ$Protein <- sapply(strsplit(PDX10PHOLFQ$Protein,";"), `[`, 1)
TotalPDX10phosProtein <- length(unique(PDX10PHOLFQ$Protein))
PDX5CellLinesPHOLFQ$Protein <- sapply(strsplit(PDX5CellLinesPHOLFQ$Protein,";"), `[`, 1)
TotalPDX5phosProtein <- length(unique(PDX5CellLinesPHOLFQ$Protein))

##CHECK TOTAL NUMBER OF GENES IN TIOX
ATCCcellPHO$Gene.Names <- sapply(strsplit(ATCCcellPHO$Gene.Names,";"), `[`, 1)
TotalATCCphosGenes <- length(unique(ATCCcellPHO$Gene.Names))
PDX10PHO$Gene.Names <- sapply(strsplit(PDX10PHO$Gene.Names,";"), `[`, 1)
TotalPDX10phosGenes <- length(unique(PDX10PHO$Gene.Names))
PDX5CellLinesPHO$Gene.Names <- sapply(strsplit(PDX5CellLinesPHO$Gene.Names,";"), `[`, 1)
TotalPDX5phosGenes <- length(unique(PDX5CellLinesPHO$Gene.Names))

#phospho intersection
peptideIntersect <- intersect(intersect(ATCCcellPHO$Sequence,PDX10PHO$Sequence),PDX5CellLinesPHO$Sequence)
geneIntersection <- intersect(intersect(ATCCcellPHO$Gene.Names,PDX10PHO$Gene.Names),PDX5CellLinesPHO$Gene.Names)
proteinIntersection <- intersect(intersect(ATCCcellPHOLFQ$Protein,PDX10PHOLFQ$Protein),PDX5CellLinesPHOLFQ$Protein)

phosphoInteraction <- data.frame(length(peptideIntersect),length(geneIntersection),length(proteinIntersection))
colnames(phosphoInteraction) <- c("Peptide_level","Gene_level","Protein_level")
pdf("Integration/PhophoInteraction.pdf")
b1 <- barplot(as.matrix(phosphoInteraction),ylim = c(0,5000),col=c("orange"),axes=FALSE,main = "Datasets interaction on different levels")
text(b1, y = as.numeric(phosphoInteraction), label =as.character(phosphoInteraction), pos = 3, cex = 1.2, col = "black")
dev.off()


ATCCphosphopeptides <- ATCCcellPHO$Sequence[!is.na(ATCCcellPHO$Sequence)]
PDX10phosphopeptides <- PDX10PHO$Sequence[!is.na(PDX10PHO$Sequence)]
PDX5CellLinesphosphopeptides <- PDX5CellLinesPHO$Sequence[!is.na(PDX5CellLinesPHO$Sequence)]


TotalforTiOX1 <- data.frame(ATCC=TotalATCCphosGenes,PDX10=TotalPDX10phosGenes,PDX5andPRIMARY=TotalPDX5phosGenes)
TotalforTiOX2 <- data.frame(ATCC=TotalATCCphosProtein,PDX10=TotalPDX10phosProtein,PDX5andPRIMARY=TotalPDX5phosProtein)
TotalforTiOX3 <- data.frame(ATCC=length(unique(ATCCcellPHO$Sequence)),PDX10=length(unique(PDX10PHO$Sequence)),PDX5andPRIMARY=length(unique(PDX5CellLinesPHO$Sequence)))

TotalforTiOX <- rbind.data.frame(TotalforTiOX1,TotalforTiOX2,TotalforTiOX3)
rownames(TotalforTiOX) <- c("Genes","Proteins","Peptides")
coul = brewer.pal(3, "Pastel2") 
# Grouped barplot
textforbarplot <- c(3439,3545,11089,3462,3567,10917,2312,2366,5933)

pdf("Integration/SummaryTiOX.pdf")
b2 <- barplot(as.matrix(TotalforTiOX), col=coul,border="white", font.axis=1, beside=T, legend=rownames(TotalforTiOX), xlab="", font.lab=2, ylim = c(0,12000),main = "Summary for TiOX data")
text(b2, textforbarplot+0.3, round(textforbarplot, 1),cex=1,pos=3)

dev.off()





#venn diagram

library(VennDiagram)

venn.diagram(
  x = list(ATCCphosphopeptides,PDX10phosphopeptides,PDX5CellLinesphosphopeptides),
  category.names = c("ATCC" , "PDX10" , "5PDX+PrimaryCell"),
  filename = '#14_venn_diagramm.png',
  output = TRUE ,
  imagetype="png" ,
  height = 780, 
  width = 780, 
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = c('yellow', 'purple', 'green'),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

######################################## primary analysis on EXPRESSION ###############################


######## REMOVE USELESS SAMPLES #####

PDX10ExprLFQ$`HPAF-II`=NULL
PDX10ExprLFQ$p53M=NULL

PDX5CellLinesExprCorrect <- PDX5CellLinesExprLFQ[,-c(6,7,8,9,10,18)]
PDX5CellLinesExprCorrect$"Protein" <- sapply(strsplit(PDX5CellLinesExprLFQ$Protein,"[|]"), `[`, 2)


###### check how many protiens you have per sample in clinical samples.

par(mar=c(4,4,4,1))

pdf("samplesize_dataset_LFQ_plusSS.pdf",width = 15)
par(mar=c(6,1,1,1))
vectorSSClinic <- apply(ClinicalExpr5gelLFQ[,1:20],2,function(x) sum(x > 0))
vectorSSPDX5CellLines <- apply(PDX5CellLinesExprCorrect[,1:12],2,function(x) sum(x > 0))
vectorSSATCC <- apply(ATCCexprLFQ[,1:11],2,function(x) sum(x > 0))
vectorSSPDX10 <- apply(PDX10ExprLFQ[,1:9],2,function(x) sum(x > 0))
vectorSSClinicSHOT <- apply(ClinicalExprBulk[,1:16],2,function(x) sum(x > 0))
entire_ss <- c(vectorSSATCC,vectorSSPDX5CellLines,vectorSSPDX10,vectorSSClinic,vectorSSClinicSHOT)
xx <- barplot(entire_ss,ylim = c(0,6000),col = c(rep("red",11),rep("blue",5),rep("black",7),rep("darkgreen",9),rep("yellow",20),rep("coral",16)),xaxt = 'n', xlab = '', width = 0.85,ylab="LFQ",
              main = "Sample size datasets")
## Add text at top of bars
text(x = xx, y = entire_ss, label = entire_ss, pos = 3, cex = 0.8, col = "black")
## Add x-axis labels 
axis(1, at=xx, labels=names(entire_ss), tick=FALSE, las=2, line=-0.5, cex.axis=1)

legend("topright", 
       legend = c("ATCC_cellLines" , "SubcutanousPDX" , "OrthotopicPDX","PrimaryCellLines","BulkTumors5gel","BulkTumorsSingleShot"), 
       col = c("red" , "darkgreen" , "blue" , "black", "yellow","coral"), 
       pch = c(20,20,20,20,20), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0.000001, 0.001))
dev.off()


boxplot(log10(ClinicalExpr5gelLFQ[,1:20]+1))
colnames(mergedALLExpr)
pdf("totalBoxplotwithzero.pdf",width = 20)
boxplot(log10(mergedALLExpr+1),col=c(rep("red",11),rep("darkgreen",9),rep("blue",5),rep("black",7),rep("yellow",20)),xaxt = "n",  xlab = "")
#axis(-1, labels = FALSE)
text(x =  seq_along(colnames(mergedALLExpr)), y = par("usr")[3] - 0.2, srt = 45, adj = 1,
     labels = colnames(mergedALLExpr), xpd = TRUE)
dev.off()

barplot(apply(ClinicalExpr5gelLFQ[,1:20],2,function(x) sum(x > 0)),ylim = c(0,4000),col = "coral")
#expression intersection

geneIntersection <- intersect(intersect(intersect(ATCCexprLFQ$Gene,PDX10ExprLFQ$Gene),PDX5CellLinesExprCorrect$Gene),ClinicalExpr5gelLFQ$Gene)
proteinIntersection <- intersect(intersect(intersect(ATCCexprLFQ$Protein,PDX10ExprLFQ$Protein),PDX5CellLinesExprCorrect$Protein),ClinicalExpr5gelLFQ$Protein)

pdf("vennDiagram5gel.pdf",width = 10)

vennOverLeaf <- venn.diagram(list(ATCCexprLFQ$Protein,PDX10ExprLFQ$Protein,PDX5CellLinesExprCorrect$Protein,ClinicalExpr5gelLFQ$Protein),NULL,
                         col = c("yellow", "green", "cyan", "coral"),
                         fill=c("white", "white", "white", "white"), 
                         alpha=c(0.5,0.5,0.5,0.5), 
                         cex = 1.5, 
                         cat.fontface=2,
                         #cat.fontfamily = rep("arial", 4),
                         category.names=c("ATCC cell lines = 7806", "SubPDX = 6114", "OrthoPDXandPrimary = 7035", "Clinical = 5763"))

grid.draw(vennOverLeaf)




dev.off()




fATCC <- ATCCexprLFQ[ATCCexprLFQ$Protein %in% proteinIntersection,]

f10PDX <- PDX10ExprLFQ[PDX10ExprLFQ$Protein %in% proteinIntersection,]
fPrimary <- PDX5CellLinesExprCorrect[PDX5CellLinesExprCorrect$Protein %in% proteinIntersection,]
fClinical <- ClinicalExpr5gelLFQ[ClinicalExpr5gelLFQ$Protein %in% proteinIntersection,]


fClinical$Gene=NULL
fPrimary$Gene=NULL
f10PDX$Gene=NULL
fATCC$Gene=NULL
merge1Expr <- merge(fATCC,f10PDX, by="Protein")
merge2Expr <- merge(merge1Expr,fPrimary,by="Protein")
mergedALLExpr <- as.data.frame(merge(merge2Expr,fClinical,by="Protein"))


rownames(mergedALLExpr) <- mergedALLExpr$Protein
mergedALLExpr$Protein=NULL


colnames(mergedALLExpr)[33] <- c("Hum.Tum.II")
colnames(mergedALLExpr)[34] <- c("Hum.Tum.JJ")
colnames(mergedALLExpr)[35] <- c("Hum.Tum.KK")
colnames(mergedALLExpr)[36] <- c("Hum.Tum.LL")
colnames(mergedALLExpr)[37] <- c("Hum.Tum.MM")
colnames(mergedALLExpr)[38] <- c("Hum.Tum.OO")
colnames(mergedALLExpr)[39] <- c("Hum.Tum.QQ")
colnames(mergedALLExpr)[40] <- c("Hum.Tum.RR")
colnames(mergedALLExpr)[41] <- c("Hum.Tum.SS")
colnames(mergedALLExpr)[42] <- c("Hum.Tum.TT")
colnames(mergedALLExpr)[43] <- c("Hum.Tum.UU")
colnames(mergedALLExpr)[44] <- c("Hum.Tum.VV")
colnames(mergedALLExpr)[45] <- c("Hum.Tum.WW")
colnames(mergedALLExpr)[46] <- c("Hum.Tum.XX")
colnames(mergedALLExpr)[47] <- c("Hum.Tum.YY")
colnames(mergedALLExpr)[48] <- c("Hum.Tum.ZZ")


write.table(file = "Expressionmerged1.txt",mergedALLExpr,quote = F)
# new data
write.table(file = "5gel_Expressionmerged1.txt",mergedALLExpr,quote = F)


# in case you have lost your .RData
mergedALLExpr_check <- read.table("Expressionmerged1.txt")
#new data with 5 gel band clinical
mergedALLExpr_check <- read.table("5gel_Expressionmerged1.txt")

# preparing for ComBat in Genepattern
gct_file <- mergedALLExpr
gct_file$"Description" <- ""
combat_gct <- cbind.data.frame(rownames(gct_file),gct_file$Description,gct_file[,1:48])
colnames(combat_gct)[1] <- "NAME"
colnames(combat_gct)[2] <- "Description"
write.table(file = "Expressionmerged1Combat.txt",combat_gct,quote = F,sep = "\t")


TmergedALLExpr <- as.data.frame(t(mergedALLExpr))



TmergedALLExpr$"Datatype" <- c(rep("ATCC_cellLines",11),rep("SubcutanousPDX",9),rep("OrthotopicPDX",5),rep("PrimaryCellLines",7),rep("BulkTumors",20))
TmergedALLExpr$"Batch" <- c(rep("Batch1",11),rep("Batch2",9),rep("Batch3",12),rep("Batch4",20))

sample_info_combat <- TmergedALLExpr
sample_info_combat$"Array" <- ""
sample_info <- cbind.data.frame(sample_info_combat$Array,rownames(TmergedALLExpr),TmergedALLExpr$Batch)
colnames(sample_info) <- c("Array","Sample","Batch")
write.table(file = "sample.info.combat.txt",sample_info,sep = "\t",row.names = F,quote = F)



checkpoint <- log10(TmergedALLExpr[,c(1:2697)])
checknew<- do.call(data.frame,lapply(checkpoint, function(x) replace(x, is.infinite(x),NA)))
checknew$Sample <- TmergedALLExpr$Datatype
checknew$Batch <- TmergedALLExpr$Batch

###### Apparently is not possible to avoid NA in hierarchical clustering, so I'm try the new imputation method: for each row, 
# create a normal distribution with mean the minimum number divide by 2 and SD=SD of that row. From this distribution select random number.

######

normalize_df <- function(df){
  
  pep_log <- log10(df)
  pep_log[pep_log == -Inf] <- NA
  
  return(pep_log)
}

set.seed("1203")
imputation_per_col <- function(Col){
  
  Col_mean <- min(Col, na.rm=TRUE)
  Col_sd <- sd(Col,na.rm = T)
  
  Col[is.na(Col)]<-rnorm(length(Col[is.na(Col)]), mean = Col_mean, sd = Col_sd)
  
  return(Col)
}

trial <- checknew[,1:2697]
result_trial <- as.data.frame(apply(trial, 2, function(x) imputation_per_col(x)))
dist=dist(result_trial , diag=TRUE)
hc=hclust(dist)
dhc=as.dendrogram(hc)
specific_leaf=dhc[[1]][[1]][[1]]
specific_leaf
attributes(specific_leaf)
i=0
colLab<<-function(n){
  if(is.leaf(n)){
    
    #I take the current attributes
    a=attributes(n)
    
    #I deduce the line in the original data, and so the treatment and the specie.
    ligne=match(attributes(n)$label,rownames(TmergedALLExpr))
    batch=TmergedALLExpr[ligne,2699];
    if(batch=="Batch1"){col_batch="blue"};if(batch=="Batch2"){col_batch="red"};if(batch=="Batch3"){col_batch="black"};if(batch=="Batch4"){col_batch="yellow"}
    datatype=TmergedALLExpr[ligne,2698];
    if(datatype=="ATCC_cellLines"){col_datatype="red"};if(datatype=="SubcutanousPDX"){col_datatype="Darkgreen"};if(datatype=="OrthotopicPDX"){col_datatype="blue"};if(datatype=="PrimaryCellLines"){col_datatype="black"};if(datatype=="BulkTumors"){col_datatype="chocolate"}
    
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=1,pch=20,col=col_batch,lab.col=col_datatype,lab.font=1,lab.cex=1))
  }
  return(n)
}
dL <- dendrapply(dhc, colLab) #### ---> there is a bug here and I don't know what does it means!!!!!!

transformed <- normalize_df(base_file)
imputated2 <- imputation(transformed)


####### Try GSEA to understand if there is a different biology between the 2 groups of clinical samples that are clusterizing apart
### expected to find more stroma genes in the aside clustering
############################# LOADING R PACKAGES AND PRE-COOKED R OBJECTS ######
library(stringr)
library(ggplot2)
library(ggthemes)
library(data.table)
library(gplots)
library(openxlsx)
library(gage)
#library(gageData)
library(qusage)
library(stringr)
library(parallel)
library(fmsb)
library(RColorBrewer)
library(openxlsx)

MSigDbGeneSets <- list(
  hallMarkGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/h.all.v6.1.symbols.gmt"),
  positionalGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c1.all.v6.1.symbols.gmt"),
  allCuratedGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c2.all.v6.1.symbols.gmt"),
  chemicalGenticPertubationGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c2.cgp.v6.1.symbols.gmt"),
  allCanonicalPathwaysGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c2.cp.v6.1.symbols.gmt"),
  bioCartaGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c2.cp.biocarta.v6.1.symbols.gmt"),
  keggGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c2.cp.kegg.v6.1.symbols.gmt"),
  reactomeGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c2.cp.reactome.v6.1.symbols.gmt"),
  allMotifGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c3.all.v6.1.symbols.gmt"),
  microRNAGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c3.mir.v6.1.symbols.gmt"),
  transcriptionFactorGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c3.tft.v6.1.symbols.gmt"),
  allComputationalGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c4.all.v6.1.symbols.gmt"),
  cancerGeneNeighborhoodsGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c4.cgn.v6.1.symbols.gmt"),
  cancerModules = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c4.cm.v6.1.symbols.gmt"),
  allGOGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c5.all.v6.1.symbols.gmt"),
  goBPGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c5.bp.v6.1.symbols.gmt"),
  goCCGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c5.cc.v6.1.symbols.gmt"),
  goMFGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c5.mf.v6.1.symbols.gmt"),
  allOncogenicSignaturesGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c6.all.v6.1.symbols.gmt"),
  allImmunologicSignaturesGeneSet = read.gmt("/Users/Giulia/Desktop/TCPA_GSEA/Input files/MSigDb_v6.1/c7.all.v6.1.symbols.gmt"),
  BindeaSignature=read.gmt("/Users/Giulia/Dropbox/LCM data and paper/GSEA/GSEA poor prognosis vs good prognosis/custom.bindea.gmt"),
  SenBabaogluSignature=read.gmt("/Users/Giulia/Dropbox/LCM data and paper/GSEA/GSEA poor prognosis vs good prognosis/custom.senbabaoglu.gmt")
)

gseaAnalyser <- function(dataSet,comps,geneSetCollection){
  dataSet = dataSet
  
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, c(library(data.table),library(gage)))
  
  l1 <- list()
  
  for(i in 1: length(comps)){
    targetGroup <- grep(comps[i], colnames(dataSet))
    others <- grep(comps[i], colnames(dataSet), invert = T)
    
    clusterExport(cl, varlist =  c("dataSet","geneSetCollection","targetGroup","others"),envir=environment())
    
    l2 <- parLapply(cl,geneSetCollection, function(x) { 
      gage(dataSet, gsets = x, ref = others, samp = targetGroup, compare = "unpaired",set.size = c(8,1200))
    })
    
    l1[[comps[i]]] <- l2
    
  }
  
  stopCluster(cl)
  return(l1)
}


gseaRadarPlotter <- function(dataObject,geneSet,plotTitle, direction = "greater", top = 5){
  grps_colors <- c("#377EB8","#4DAF4A","#E41A1C")
  grps <- names(dataObject)
  MyMerge <- function(x, y){
    dt <- merge(x, y, by= "Title", all= T)
    return(dt)
  }
  add.alpha <- function(col, alpha=1){
    
    apply(sapply(col, col2rgb)/255, 2,
          function(x) {
            rgb(x[1], x[2], x[3], alpha=alpha)
          }
    )
  }
  
  switch(direction,
         greater = {
           dt1 <- Reduce(MyMerge,lapply(dataObject, function(x){
             
             dt0 <- as.data.table(x[[geneSet]]$greater[,1:5])
             dt0[,Title := unlist(lapply(str_split(rownames(x[[geneSet]]$greater),pattern = "_"), function(x){
               paste(x[-1],collapse = " ")
             })) ]
             dt0[,negLogQval := round(log10(1/q.val),2)]
             dt0 <- dt0[,.(Title,negLogQval)]
             
           }))
           names(dt1) <- c("Title",grps)
           dt1[is.na(dt1)] <- 0
           
           dt1 <- dt1[ Title %in% c(sapply(grps, function(x) {
             dt1[head(order(dt1[,.SD,.SDcol = x],decreasing = T),top),Title]})),]
           
           dt1[,max := round(ceiling(max(dt1[,-1],na.rm = T)+10),-1)]
           dt1[,min := 0]
           dt1 <- dt1[,.SD,.SDcols = unlist(lapply(c("Title","max","min",grps),function(x) match(x,names(dt1)))) ]
           
           dt2 <- as.data.frame(t(dt1[,-1]))
           colnames(dt2) <- dt1[,Title]
           
           radarchart(dt2, axistype=1 , seg = 5, title = paste0("GSEA: ",plotTitle,"$",direction),
                      pcol= add.alpha(grps_colors[1:length(grps)],alpha = 0.8),plty = 1,
                      pfcol=add.alpha(grps_colors[1:length(grps)],alpha = 0.4), plwd=1 ,
                      cglcol="gray90", cglty=1, axislabcol="grey", cglwd=0.8,caxislabels=seq(0,100,20),
                      vlcex=0.6 )
           
           legend(1.5,1, legend= grps , 
                  seg.len=0.5, pch=1, bty="n" ,
                  lwd=3, y.intersp=1, horiz=FALSE, 
                  col=add.alpha(grps_colors[1:length(grps)],alpha = 0.4))
           
           
           
         },
         less = {
           dt1 <- Reduce(MyMerge,lapply(dataObject, function(x){
             
             dt0 <- as.data.table(x[[geneSet]]$less[,1:5])
             dt0[,Title := unlist(lapply(str_split(rownames(x[[geneSet]]$less),pattern = "_"), function(x){
               paste(x[-1],collapse = " ")
             })) ]
             dt0[,negLogQval := round(log10(1/q.val),2)]
             dt0 <- dt0[,.(Title,negLogQval)]
             
           }))
           names(dt1) <- c("Title",grps)
           dt1[is.na(dt1)] <- 0
           
           dt1 <- dt1[ Title %in% c(sapply(grps, function(x) {
             dt1[head(order(dt1[,.SD,.SDcol = x],decreasing = T),top),Title]})),]
           
           dt1[,max := round(ceiling(max(dt1[,-1],na.rm = T)+10),-1)]
           dt1[,min := 0]
           dt1 <- dt1[,.SD,.SDcols = unlist(lapply(c("Title","max","min",grps),function(x) match(x,names(dt1)))) ]
           
           dt2 <- as.data.frame(t(dt1[,-1]))
           colnames(dt2) <- dt1[,Title]
           
           radarchart(dt2, axistype=1 , seg = 5, title = paste0("GSEA: ",plotTitle,"$",direction),
                      pcol= add.alpha(grps_colors[1:length(grps)],alpha = 0.8),plty = 1,
                      pfcol=add.alpha(grps_colors[1:length(grps)],alpha = 0.4), plwd=1 ,
                      cglcol="gray90", cglty=1, axislabcol="grey", cglwd=0.8,caxislabels=seq(0,100,20),
                      vlcex=0.6 )
           
           legend(1.5,1, legend= grps , 
                  seg.len=0.5, pch=1, bty="n" ,
                  lwd=3, y.intersp=1, horiz=FALSE, 
                  col=add.alpha(grps_colors[1:length(grps)],alpha = 0.4))
           
         })
  
  
  
  
  
  
}

exp.moreStroma <- ClinicalExpr5gelLFQ[,c("(O)","(P)","(I)","(N)","(H)","(J)","(C)","(A)","(D)","(Q)","(F)","(M)")]
exp.moreTumor <- ClinicalExpr5gelLFQ[,c("(E)","(R)","(B)","(G)","(S)","(T)","(K)","(L)")]
colnames(exp.moreStroma) <-c("moreStroma(O)","moreStroma(P)","moreStroma(I)","moreStroma(N)","moreStroma(H)","moreStroma(J)","moreStroma(C)","moreStroma(A)","moreStroma(D)","moreStroma(Q)","moreStroma(F)","moreStroma(M)") 
colnames(exp.moreTumor) <- c("moreTumor(E)","moreTumor(R)","moreTumor(B)","moreTumor(G)","moreTumor(S)","moreTumor(T)","moreTumor(K)","moreTumor(L)")
combinedforGSEAonClinical <- cbind.data.frame(exp.moreStroma,exp.moreTumor)
rownames(combinedforGSEAonClinical) <- make.names(ClinicalExpr5gelLFQ$Gene,unique = T)
system.time(
  gseaResult <- gseaAnalyser(dataSet = as.matrix(combinedforGSEAonClinical),
                             comps = c("moreStroma","moreTumor"),
                             geneSetCollection = MSigDbGeneSets))

pdf("/Users/Giulia/Desktop/meglio2.pdf",width = 11.69, height = 8.27)
for (i in c(1,5,6,7,10,13,14,15,16,17,18,20)){
  print(names(MSigDbGeneSets)[i])
  gseaRadarPlotter(dataObject = gseaResult,geneSet = names(MSigDbGeneSets)[i],plotTitle = names(MSigDbGeneSets)[i])
  gseaRadarPlotter(dataObject = gseaResult,geneSet = names(MSigDbGeneSets)[i],plotTitle = names(MSigDbGeneSets)[i],direction = "less")
  
}
dev.off()



####### HIERARCHICAL CLUSTERING ########

dist=dist(log10(TmergedALLExpr[ , c(1:2697)]+1) , diag=TRUE)
hc=hclust(dist)
dhc=as.dendrogram(hc)
specific_leaf=dhc[[1]][[1]][[1]]
specific_leaf
attributes(specific_leaf)

i=0
colLab<<-function(n){
  if(is.leaf(n)){
    
    #I take the current attributes
    a=attributes(n)
    
    #I deduce the line in the original data, and so the treatment and the specie.
    ligne=match(attributes(n)$label,rownames(TmergedALLExpr))
    batch=TmergedALLExpr[ligne,2699];
    if(batch=="Batch1"){col_batch="blue"};if(batch=="Batch2"){col_batch="red"};if(batch=="Batch3"){col_batch="black"};if(batch=="Batch4"){col_batch="yellow"}
    datatype=TmergedALLExpr[ligne,2698];
    if(datatype=="ATCC_cellLines"){col_datatype="red"};if(datatype=="SubcutanousPDX"){col_datatype="Darkgreen"};if(datatype=="OrthotopicPDX"){col_datatype="blue"};if(datatype=="PrimaryCellLines"){col_datatype="black"};if(datatype=="BulkTumors"){col_datatype="chocolate"}
    
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=1,pch=20,col=col_batch,lab.col=col_datatype,lab.font=1,lab.cex=1))
  }
  return(n)
}
dL <- dendrapply(dhc, colLab)
pdf("Integration/dendrogram3.pdf",width = 12)
plot(dL , main="Integration in Expression proteomics")
legend("topright", 
       legend = c("Exp1" , "Exp2" ,  "Exp3" , "Exp4" ,"ATCC_cellLines" , "SubcutanousPDX" , "OrthotopicPDX","PrimaryCellLines","BulkTumors"), 
       col = c("blue", "red" , "black" , "yellow" , "red", "Darkgreen","blue","black","chocolate"), 
       pch = c(20,20,20,20,4,4,4,4,4), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0.000001, 0.001))

dev.off()


##### simplistic visualization of the dendrogram #####

simpleDendr <- TmergedALLExpr
simpleDendr$Batch=NULL

dist=dist(TmergedALLExpr[ , c(1:2697)], diag=TRUE)

hc=hclust(dist)
dhc=as.dendrogram(hc)
specific_leaf=dhc[[1]][[1]][[1]]
specific_leaf
attributes(specific_leaf)

i=0
colLab<<-function(n){
  if(is.leaf(n)){
    
    #I take the current attributes
    a=attributes(n)
    
    #I deduce the line in the original data, and so the treatment and the specie.
    ligne=match(attributes(n)$label,rownames(TmergedALLExpr))
    datatype=TmergedALLExpr[ligne,2698];
    if(datatype=="ATCC_cellLines"){col_datatype="red"};if(datatype=="SubcutanousPDX"){col_datatype="Darkgreen"};if(datatype=="OrthotopicPDX"){col_datatype="blue"};if(datatype=="PrimaryCellLines"){col_datatype="black"};if(datatype=="BulkTumors"){col_datatype="yellow"}
    
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=1,pch=20,col=col_datatype,lab.font=1,lab.cex=1))
  }
  return(n)
}
dL <- dendrapply(dhc, colLab)
pdf("Integration/Simple_dendrogram_noLOG.pdf",width = 12)
plot(dL , main="Integration of ATCC, Primary, Sub_PDX, Ortho_PDX, Human_Tumours \n in PDAC proteomics")
legend("topright", 
       legend = c("ATCC_cellLines" , "SubcutanousPDX" , "OrthotopicPDX","PrimaryCellLines","BulkTumors"), 
       col = c("red", "Darkgreen","blue","black","yellow"), 
       pch = c(20,20,20,20,20), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0.000001, 0.001))

dev.off()


######## PCA is not looking nice!!!
##### PCA ######
PCA_all <- TmergedALLExpr[,1:2697]
PCA_all$SampleType <- c(rep("ATCC",11),rep("SubPDX",9),rep("OrthoPDX",5),rep("PrimaryPDX",7),rep("HumanTumor",20))
PCA_all$col <- c(rep("red",11),rep("green",9),rep("pink",5),rep("blue",7),rep("orange",20))

autoplot(prcomp(PCA_all[,1:2697]), data = PCA_all, colour = "SampleType", main="",x=1, y=2,frame=F, label = TRUE, label.size =5)
dev.off()
set.seed(1)
KM <- PCA_all[,1:2697]
cnames <- sapply(strsplit(colnames(KM),"[-]"), `[`, 1)
colnames(KM) <- cnames
cnames <-  sapply(strsplit(colnames(KM),"[.]"), `[`, 1)
colnames(KM) <- cnames

#### remove duplicated rows
dedup <- KM[ , -which(table(names(KM)) >1)]
autoplot(kmeans(dedup, 3), data = dedup)


###### try other types of clusterings #####
# Ward Hierarchical Clustering
d <- dist(KM, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")
####### don't run too much memory ######
# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit <- pvclust(KM, method.hclust="ward",
               method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)
#################### end of don't run ########


# K-Means Clustering with 5 clusters
fit <- kmeans(KM, 5)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster) 
clusplot(KM, fit$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(KM, fit$cluster)



################### skip combat? compute correlation in all the datasets ########

colnames(mergedALLExpr)
Clinical_cor <- mergedALLExpr[,grep("Hum.Tum",colnames(mergedALLExpr))]
SubPDX_cor <- mergedALLExpr[,c("X93","X72","X84","X90","X96","X67","X76","p28","p53tb")]
OrthoPDX_cor<- mergedALLExpr[,grep("PDX",colnames(mergedALLExpr))]
Primary_cor <- mergedALLExpr[,c("pp391","PDAC.1","PDAC.3","pp161","PDAC.4","PDAC.2","PDAC.5")]
ATCC_cor <- mergedALLExpr[,c("AsPC.1","Capan.2","CFPAC.1","HPAF.II","BxPC.3","HPAC","PL45","Miapaca.2","Suit.2","PANC.1","Hs766t")]

A <- cor(t(Primary_cor))
B <- cor(t(SubPDX_cor))
C <- cor(t(OrthoPDX_cor))
D <- cor(t(ATCC_cor))


library(mixOmics)
pdf("ATCCandPrimarycor.pdf")
imgCor(ATCC_cor, Primary_cor, X.var.names = F, Y.var.names = FALSE, interactive.dev=FALSE)
dev.off()

pdf("ATCCandOrthoPDXcor.pdf")
imgCor(ATCC_cor, OrthoPDX_cor, X.var.names = F, Y.var.names = FALSE, interactive.dev=FALSE)
dev.off()


pdf("ATCCandSubPDXcor.pdf")
imgCor(ATCC_cor, SubPDX_cor, X.var.names = F, Y.var.names = FALSE, interactive.dev=FALSE)
dev.off()

pdf("ATCCandCliniccor.pdf")
imgCor(ATCC_cor, Clinical_cor, X.var.names = F, Y.var.names = FALSE, interactive.dev=FALSE)
dev.off()

pdf("PrimaryndSubcor.pdf")
imgCor(Primary_cor, SubPDX_cor, X.var.names = F, Y.var.names = FALSE, interactive.dev=FALSE)
dev.off()

pdf("PrimaryndOrthoPDXcor.pdf")
imgCor(Primary_cor, OrthoPDX_cor, X.var.names = F, Y.var.names = FALSE, interactive.dev=FALSE)
dev.off()

pdf("PrimaryndClinicalcor.pdf")
imgCor(Primary_cor, Clinical_cor, X.var.names = F, Y.var.names = FALSE, interactive.dev=FALSE)
dev.off()

pdf("OrthondSubcor.pdf")
imgCor(OrthoPDX_cor, SubPDX_cor, X.var.names = F, Y.var.names = FALSE, interactive.dev=FALSE)
dev.off()


GOODcorr <- cbind.data.frame(ATCC_cor,Primary_cor,OrthoPDX_cor,SubPDX_cor)
POORcorr <- cbind.data.frame(Clinical_cor)

pdf("GoodVsPoor.pdf")
imgCor(GOODcorr, POORcorr, X.var.names = F, Y.var.names = FALSE, interactive.dev=FALSE)
dev.off()

################## data imputation ##################
require(DMwR)
require(imputation)
require(VIM)
ATCCforImp <- ATCC_cor
ATCCforImp[ATCCforImp==0] <- NA
ATCC_imputated <- kNN(ATCCforImp,k=7)
boxplot(log2(ATCC_imputated[,1:11]),col = "yellow")


PrimaryforImp <- Primary_cor
PrimaryforImp[PrimaryforImp==0] <- NA
Primary_imputated <- kNN(PrimaryforImp,k=7)
boxplot(log2(Primary_imputated[,1:7]),col = "yellow")

SubPDXforImp <- SubPDX_cor
SubPDXforImp[SubPDXforImp==0] <- NA
SubPDXforImputated <- kNN(SubPDXforImp,k=7)

OrthoforImp <- OrthoPDX_cor
OrthoforImp[OrthoforImp==0] <- NA
OrthoforImputated<- kNN(OrthoforImp,k=5)

Clinical_forImp <- ClinicalExpr5gelLFQ[,1:20]
Clinical_forImp[Clinical_forImp==0] <- NA
Clinical_Imputated <- kNN(Clinical_forImp,k=20)

forImputatedClustering <- cbind.data.frame(ATCC_imputated[,1:11])


###### BOXPLOT after imputation #######
colorBOX <- c(rep("darkseagreen2",11),rep("bisque",7),rep("aliceblue",9),rep("darkgoldenrod1",5))
DATAforBoxplot <- cbind.data.frame(ATCC_imputated[,1:11],Primary_imputated[,1:7],SubPDXforImputated[,1:9],OrthoforImputated[,1:5])
pdf("boxplotCLandPDX.pdf")
par(mar=c(4,2,2,8))
bp <- boxplot(log2(DATAforBoxplot),col = colorBOX)
legend(35,35,c("ATCC","Primary","SubPDX","OrthoPDX"),fill = c("darkseagreen2","bisque","aliceblue","darkgoldenrod1"))

dev.off()


################## data normalization ##################
boxplot(ATCC_cor)
min(ATCC_cor)
max(ATCC_cor)
max(Primary_cor)
max(SubPDX_cor)
max(OrthoPDX_cor)

##########  density 2D plots ########## 
# Packages
library(hexbin)
library(RColorBrewer)

# Create data
ATCC <- rowMeans(log2(ATCC_imputated[,1:11]))
Primary <- rowMeans(log2(Primary_imputated[,1:7]))

# Make the plot
pdf("densityplotATCC.pdf")
bin<-hexbin(ATCC, Primary, xbins=40)
my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
plot(bin, main="" , colramp=my_colors , legend=F)
dev.off()


##### MA plot #######
library(rafalib)
M = log10(DATAforBoxplot$AsPC.1) - log10(DATAforBoxplot$pp391)
A = (log10(DATAforBoxplot$AsPC.1) + log10(DATAforBoxplot$pp391))/2
pdf("before_quantile.pdf")
splot(A, M, ylim = c(-1.5, 1.5))
abline(h = 0, col = 2, lwd = 2)
dev.off()

library(preprocessCore)
nqDATAforBoxplot <- as.data.frame(normalize.quantiles(as.matrix(DATAforBoxplot)))
colnames(nqDATAforBoxplot) <- colnames(DATAforBoxplot)
pdf("after_quantile_norm.pdf")
boxplot(log2(nqDATAforBoxplot),col = colorBOX)
dev.off()


M = log10(nqDATAforBoxplot$AsPC.1) - log10(nqDATAforBoxplot$pp391)
A = (log10(nqDATAforBoxplot$AsPC.1) + log2(nqDATAforBoxplot$pp391))/2
pdf("after_quantile.pdf")
splot(A, M, ylim = c(-1.5, 1.5))
abline(h = 0, col = 2, lwd = 2)
dev.off()


###### check for batch effect? ######
library("sva")
library(bladderbatch)
data(bladderdata)
data(bladderEset)
edata = exprs(bladderEset)
edata <- DATAforBoxplot[,1:18]
yu <- data.frame(colnames(edata),colorBOX3)
colnames(yu) <- c("Sample","cancer")

#Applying the ComBat function to adjust for known batches

batch = c(rep(1,11),rep(2,7))


modcombat = model.matrix(~1, data=yu)
combat_edata = ComBat(dat=log2(edata), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots = T)
  
pValuesComBat = f.pvalue(as.matrix(combat_edata),mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")                     
n.sv = num.sv(nqDATAforBoxplot,mod,vfilter=2000,method="leek")
svobj = sva(nqDATAforBoxplot,mod,mod0,n.sv=n.sv,vfilter=2000)


####### density plots #########
colorBOX2 <- c(rep("blue",11),rep("red",7),rep("lightgreen",9),rep("black",5))

colorBOX3 <- c(rep("blue",11),rep("red",7))

par(mfrow=c(1,1))
pdf("density_on_imputated.pdf")
par(new=F)
for (i in 1:length(colnames(DATAforBoxplot))){
  d <- density(log2(DATAforBoxplot[,i]))
  plot(d,xaxt="none",yaxt="none",xlab="",ylab="",lwd=1.5,xlim=c(18,40),main="",col=colorBOX2[i]) # plots the results
  par(new=T)
}
legend("topright",legend = c("ATCC","Primary","SubPDX","OrthoPDX"),col = c("blue","red","lightgreen","black"),lty = 1,lwd = 2,bty = "n")
# returns the density data 
dev.off()

pdf("density_on_qnormalized.pdf")
par(new=F)
for (i in 1:length(colnames(DATAforBoxplot))){
  d <- density(log2(nqDATAforBoxplot[,i]))
  plot(d,xaxt="none",yaxt="none",xlab="",ylab="",lwd=1.5,xlim=c(18,40),main="",col=colorBOX2[i]) # plots the results
  par(new=T)
}
legend("topright",legend = c("ATCC","Primary","SubPDX","OrthoPDX"),col = c("blue","red","lightgreen","black"),lty = 1,lwd = 2,bty = "n")
# returns the density data 
dev.off()


min(log2(DATAforBoxplot))






###### I cannot find batch effects  #####
###### start with WGCNA #####
library(fastcluster)
library(WGCNA)
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function on the gene expression profiles
sft <- pickSoftThreshold(log10(nqDATAforBoxplot), powerVector = powers, verbose = 5)

pdf("Network_topology.pdf")
# Plot the results:
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

sft <- pickSoftThreshold(log10(nqDATAforBoxplot[,1:27]), powerVector = powers, verbose = 5)

pdf("Network_topology_selected.pdf")
# Plot the results:
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

sft <- pickSoftThreshold(log10(nqDATAforBoxplot[,27:32]), powerVector = powers, verbose = 5)

pdf("Network_topology_PDX2807.pdf")
# Plot the results:
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#############################################   try with the entire dataset ###################################  

sft <- pickSoftThreshold(mergedALLExpr, powerVector = powers, verbose = 5)

pdf("Network_topology_PDX2807.pdf")
# Plot the results:
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

######### try normfinder #######
source("normFinder.txt")
imputNormFinder <- mergedALLExpr
rownames(imputNormFinder) <- rownames(mergedALLExpr)
imputNormFinder[2808,] <- c(rep("ATCC",11),rep("Primary",7),rep("Sub",9),rep("Ortho",5),rep("Clinical",16))
rownames(imputNormFinder)[2808] <- "Groups"

write.table(file = "inputforNormFinder.txt",imputNormFinder,quote=F)
Result=Normfinder("inputforNormFinder.txt")
Result2=Normfinder("inputforNormFinder.txt",ctVal = FALSE)


######### PCA on all samples per Protein expression #######
library(ggfortify); library(ggplot2)

PCA_all <- rbind.data.frame(mergedALLExpr,batch = c(rep("1",11),rep("2",9),rep("3",12),rep("4",20)))
pdf("pca_prebatch.pdf")
autoplot(prcomp(as.matrix(t(mergedALLExpr))), data = t(PCA_all),col="batch", label = T, label.size = 4, label.vjust=1.5, label.colour="black",label.lineheight=1, main="Principal Component Analysis")
dev.off()
GenesOnPCs = predict(prcomp(as.matrix(t(mergedALLExpr))))

###### combat data #####
human_toimp <- mergedALLExpr[,33:48]
human_toimp[human_toimp==0] <- NA
Human_imputated <- kNN(human_toimp,k=16)
newCombatdata <- cbind.data.frame(combat_edata,log2(DATAforBoxplot[19:32]),log2(Human_imputated[,1:16]))

PCA_all <- rbind.data.frame(newCombatdata,batch = c(rep("blue",11),rep("black",7),rep("red",9),rep("black",5),rep("yellow",16)))
DF <- as.data.frame(t(PCA_all))
pdf("PCA_combat_on_cellLines_only.pdf")
autoplot(prcomp(as.matrix(t(newCombatdata))), data = DF,col=DF$V2808, label = T, label.size = 3, label.vjust=1.5, label.colour="black",label.lineheight=1, main="Principal Component Analysis")
dev.off()

PCA <-prcomp(as.matrix(newCombatdata)) 
library(rgl)
plot3d(PCA$rotation[,1],PCA$rotation[,2],PCA$rotation[,3], col="red", size=3,xlab = "PCA1",ylab = "PCA2",zlab = "PCA3",label=T)
par(mfrow=c(1,1))
plot(PCA$rotation[,1],PCA$rotation[,2])
text(PCA$rotation[,1],PCA$rotation[,2], row.names(PCA$rotation), cex=0.6, pos=4, col=batch)


#################### WGCNA on Cell Lines+ PDX models    ##############
CellLinesPDX <- newCombatdata[,1:32]
sft <- pickSoftThreshold(t(CellLinesPDX), powerVector = powers, verbose = 5)
pdf("CellLinesandPDXCombat2807.pdf")
# Plot the results:
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#### adjacency ###
softPower = 3;
adjacency = adjacency(t(CellLinesPDX), power = softPower);
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath

sizeGrWindow(8,6)
pdf("merged_modulesandmax30.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


# Calculate eigengenes
MEList = moduleEigengenes(t(CellLinesPDX), colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
pdf("Clust_MO_eigengene.pdf")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()
# Call an automatic merging function
merge = mergeCloseModules(t(CellLinesPDX), dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
pdf(file = "mergedModulesdendro.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "PDXandCelllinesComBat.RData")
rownames(CellLinesPDX) <- rownames(mergedALLExpr)
turquoiseProteins <- rownames(CellLinesPDX)[moduleColors=="turquoise"]
blueProteins <- rownames(CellLinesPDX)[moduleColors=="blue"]
yellowProteins <- rownames(CellLinesPDX)[moduleColors=="yellow"]
brownProteins <- rownames(CellLinesPDX)[moduleColors=="brown"]

write.table(file="mod1.txt",turquoiseProteins,quote = F,col.names = F,row.names = F)
write.table(file="mod2.txt",blueProteins,quote = F,col.names = F,row.names = F)
write.table(file="mod3.txt",brownProteins,quote = F,col.names = F,row.names = F)
write.table(file="mod4.txt",yellowProteins,quote = F,col.names = F,row.names = F)


#### module analysis

net <- blockwiseModules(t(CellLinesPDX), power = 3, TOMType = "unsigned", minModuleSize = 5, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 0)

# Convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
pdf("module_network_CellLinesandPDX.pdf")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Frequency table of the number of genes in each module
table(mergedColors)

#You can extract, for example, the genes belonging to the turqoise module and write them to a text file as follows:
turquoiseGenes <- colnames(Y)[mergedColors=="blue"]
write(turquoiseGenes, file="blueGenes.txt")



















##############################   try WGCNA only with cell lines ############################### 

sft <- pickSoftThreshold(log10(DATAforBoxplot[,1:11]), powerVector = powers, verbose = 5)

pdf("Network_topology_cell_linesCombat2807.pdf")
# Plot the results:
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

###################### try to do ComBat the Cell lines with all genes shared between them #####

ATCCexprLFQ$Protein <- sapply(strsplit(ATCCexprLFQ$Protein,";"), `[`, 1)
PDX5CellLinesExprCorrect$Protein <- sapply(strsplit(PDX5CellLinesExprCorrect$Protein,"[|]"), `[`, 2)
proteinIntersection2 <- intersect(ATCCexprLFQ$Protein,PDX5CellLinesExprCorrect$Protein)

fATCC2 <- ATCCexprLFQ[ATCCexprLFQ$Protein %in% proteinIntersection2,]


fPrimary2 <- PDX5CellLinesExprCorrect[PDX5CellLinesExprCorrect$Protein %in% proteinIntersection2,]

mergeonlyCellLines <- merge(fATCC2,fPrimary2, by="Protein")

mergeonlyCellLines$Gene.x=NULL
mergeonlyCellLines$Gene.y=NULL
rownames(mergeonlyCellLines) <- mergeonlyCellLines$Protein
mergeonlyCellLines$Protein=NULL
new_merge_cell <- mergeonlyCellLines[,-c(12,13,14,15,16)]

##### ComBat ####

edata <- new_merge_cell
yu <- data.frame(colnames(edata),colorBOX3)
colnames(yu) <- c("Sample","cancer")

#Applying the ComBat function to adjust for known batches

batch = c(rep(1,11),rep(2,7))


modcombat = model.matrix(~1, data=yu)
combat_edata = ComBat(dat=log2(edata), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots = T)

###### pca on combat data and hclus ######
pdf("pca_onlycellLines_combat.pdf")
autoplot(prcomp(as.matrix(t(combat_edata))),col=batch, label = T, label.size = 3, label.vjust=1.5, label.colour="black",label.lineheight=1, main="Principal Component Analysis")
dev.off()
par(mfrow=c(1,1))
sampleTree <- hclust(dist(t(combat_edata)), method="average")
pdf("sample_clustering_onlycellLines_combat.pdf")
plot(sampleTree,labels=colnames(combat_edata))
dev.off()

################## imputation ####
cellLines_toimp <- combat_edata
cellLines_toimp[cellLines_toimp==0] <- NA
cellLines_imputated <- kNN(cellLines_toimp,k=18)




sft <- pickSoftThreshold(combat_edata, powerVector = powers, verbose = 5)

pdf("Network_topology_cell_linesCombat4826.pdf")
# Plot the results:
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()




library(fastcluster)
library(WGCNA)


### RISOLTO! DEVI TRASPORRE IDIOTA!!!

newPDX5 <- PDX5CellLinesExprLFQ[sapply(strsplit(PDX5CellLinesExprLFQ$Protein,"[|]"), `[`, 2) %in% proteinIntersection,]

rownames(newPDX5) <- sapply(strsplit(newPDX5$Protein,"[|]"), `[`, 2)

newnewPDX5 <- newPDX5[,1:5]

sampleTree <- hclust(dist(t(newnewPDX5)), method="average")
pdf("sample_clustering_onlycellLines_combat.pdf")
plot(sampleTree,labels=colnames(fPDX5))
dev.off()


#impute
newnewPDX5[newnewPDX5==0] <- 0.0001
#remove rows with 0
row_sub = apply(newnewPDX5, 1, function(row) all(row !=0 ))
newnewPDX5[row_sub,]

T_newnewPDX5 <- t(newnewPDX5)
sft <- pickSoftThreshold(T_newnewPDX5, powerVector = powers, verbose = 5)


pdf("fuckyouuuuuuuuu.pdf")
# Plot the results:
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()



############## Riprova con ALL #######

sft <- pickSoftThreshold(T_newnewPDX5, powerVector = powers, verbose = 5)


pdf("fuckyouuuuuuuuu.pdf")
# Plot the results:
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


######thang meeting 21st June 2018

# do the hierarchical clustering only on bulk tumors and see how do they cluster 
par(mar=c(4,4,4,4))
colnames(fClinical)
clusters <- hclust(dist(t(fClinical[, -c(21)]+1)))
pdf("/Users/Giulia/Desktop/WP1/onlyBulkHC.pdf")
plot(clusters)
dev.off()
###### why the clustering from transformed data and not transformed are so different????

####check the raw count on every sample type#### 
toselect <- ClinicalExpr[,32:66]
ClinicalSS_rawcount<- toselect[,grep("PDAC",colnames(toselect))]
vectorSSClinicSHOT_count <- apply(ClinicalSS_rawcount[,-c(6,8)],2,function(x) sum(x > 0))

colnames(PDX5CellLinesExpr)
toselect <- PDX5CellLinesExpr[,12:29]
colnames(toselect[,c(1,2,3,4,5,11,12,13,14,15,16,17)])
vectorSSPDX5CellLines_count <- apply(toselect[,c(1,2,3,4,5,11,12,13,14,15,16,17)],2,function(x) sum(x > 0))

colnames(ATCCcellExpr)
toselect <- ATCCcellExpr[,12:24]
colnames(toselect[,c(1,2,3,4,5,6,7,8,9,12,13)])
vectorSSATCC_count <- apply(toselect[,c(1,2,3,4,5,6,7,8,9,12,13)],2,function(x) sum(x > 0))

colnames(PDX10Expr)
toselect <- PDX10Expr[,12:21]
colnames(toselect)
vectorSSPDX10_count <- apply(toselect,2,function(x) sum(x > 0))

colnames(ClinicalExpr5gel)
toselect <- ClinicalExpr5gel[,15:34]
colnames(toselect)
vectorSSClinic5gel_count <- apply(toselect,2,function(x) sum(x > 0))
entire_ss_count <- c(vectorSSATCC_count,vectorSSPDX5CellLines_count,vectorSSPDX10_count,vectorSSClinic5gel_count,vectorSSClinicSHOT_count)
pdf("samplesize_dataset_Count_plusSS.pdf",width = 15)
xx <- barplot(entire_ss_count,ylim = c(0,6000),col = c(rep("red",11),rep("blue",5),rep("black",7),rep("darkgreen",10),rep("yellow",20),rep("coral",16)),xaxt = 'n', xlab = '', width = 0.85,ylab="Raw_Count",
              main = "Sample size datasets")

legend("topright", 
       legend = c("ATCC_cellLines" , "SubcutanousPDX" , "OrthotopicPDX","PrimaryCellLines","BulkTumors5gel","BulkTumorsSingleShot"), 
       col = c("red" , "darkgreen" , "blue" , "black", "yellow","coral"), 
       pch = c(20,20,20,20,20), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0.000001, 0.001))
dev.off()


##### check the zeros ####

zero_check_Clinic <- colSums(ClinicalExpr5gelLFQ[,1:20] == 0)
zero_check_PDX5Primary <- colSums(PDX5CellLinesExprCorrect[,1:12] == 0)
zero_check_ATCC <- colSums(ATCCexprLFQ[,1:11] == 0)
zero_check_PDX10 <- colSums(PDX10ExprLFQ[,1:9] == 0)
zero_check_ClinicSHOT <- colSums(ClinicalExprBulk[,1:16] == 0)
entire_ss_zero <- c(zero_check_ATCC,zero_check_PDX5Primary,zero_check_PDX10,zero_check_Clinic,zero_check_ClinicSHOT)
pdf("NumberofZEROS_plusSS.pdf",width = 15)
xx <- barplot(entire_ss_zero,ylim = c(0,6000),col = c(rep("red",11),rep("blue",5),rep("black",7),rep("darkgreen",10),rep("yellow",20),rep("coral",16)),xaxt = 'n', xlab = '', width = 0.85,ylab="Number of zeros",
              main = "Zero count datasets")

legend("topright", 
       legend = c("ATCC_cellLines" , "SubcutanousPDX" , "OrthotopicPDX","PrimaryCellLines","BulkTumors5gel","BulkTumorsSingleShot"), 
       col = c("red" , "darkgreen" , "blue" , "black", "yellow","coral"), 
       pch = c(20,20,20,20,20), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0.000001, 0.001))

dev.off()

######### check the sample sum for all samples #######

pdf("proteins_Sums.pdf",width = 15)
samples_sum <- c(rowSums(TmergedALLExpr[,1:2697]),colSums(ClinicalExprBulk[,1:16]))
xx <- barplot(log10(samples_sum),ylim = c(0,15),col = c(rep("red",11),rep("blue",5),rep("black",7),rep("darkgreen",10),rep("yellow",20),rep("coral",16)),xaxt = 'n', xlab = '', width = 0.85,ylab="Log10(Sum_of_Proteins)",
              main = "Proteins sum for each sample")
dev.off()

######## KUSTER methods #### allow only proteins that are present at least 95 % of all samples
cutoff=95
atmax <- (52/100)*cutoff
rownames(TmergedALLExpr)
ismajorthanzero <- apply(TmergedALLExpr[,1:2697],2,function(x) sum(x > 0))
selectedProteins <- ismajorthanzero[ismajorthanzero >= atmax]

SubsetTmerged <- TmergedALLExpr[,names(selectedProteins)]
SubsetTmerged$"SampleType" <- TmergedALLExpr[,2698]
dist=dist(log10(SubsetTmerged[ , -c(1)]+1) , diag=TRUE)
hc=hclust(dist)
dhc=as.dendrogram(hc)
specific_leaf=dhc[[1]][[1]][[1]]
specific_leaf
attributes(specific_leaf)

i=0
colLab<<-function(n){
  if(is.leaf(n)){
    
    #I take the current attributes
    a=attributes(n)
    
    #I deduce the line in the original data, and so the treatment and the specie.
    ligne=match(attributes(n)$label,rownames(SubsetTmerged))
    datatype=SubsetTmerged[ligne,length(names(SubsetTmerged))];
    if(datatype=="ATCC_cellLines"){col_datatype="red"};if(datatype=="SubcutanousPDX"){col_datatype="Darkgreen"};if(datatype=="OrthotopicPDX"){col_datatype="blue"};if(datatype=="PrimaryCellLines"){col_datatype="black"};if(datatype=="BulkTumors"){col_datatype="chocolate"}
    
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=1,pch=20,lab.col=col_datatype,lab.font=1,lab.cex=1))
  }
  return(n)
}
dL <- dendrapply(dhc, colLab)

plot(dL , main="Integration in Expression proteomics")
legend("topright", 
       legend = c("ATCC_cellLines" , "SubcutanousPDX" , "OrthotopicPDX","PrimaryCellLines","BulkTumors"), 
       col = c("blue", "red" , "black" , "yellow" , "red", "Darkgreen","blue","black","chocolate"), 
       pch = c(20,20,20,20,20), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0.000001, 0.001))

########## correlation between sample size and tumor % ##########
names(entire_ss)
left <- entire_ss[names(entire_ss) %in% c("(A)","(C)","(D)","(F)","(H)","(I)","(J)","(M)","(N)","(O)","(P)","(Q)")]
right <- entire_ss[names(entire_ss) %in% c("(B)","(E)","(G)","(K)","(L)","(R)","(S)","(T)")]

tumor_perc_left <- c(10,40,10,40,20,50,10,35,15,70,10,15)
tumor_perc_right <- c(50,70,40,30,40,60,10,35)
max(c(left,right))

pdf("correlation_samplesID_tumorPercentage.pdf")
plot(c(left,right),c(tumor_perc_left,tumor_perc_right), col=c(rep("blue",length(tumor_perc_left)),rep("red",length(tumor_perc_right))),pch=20,cex=2,xlab="sample IDentification",xlim=c(1800,3600),ylab="Tumor percentage",cex.lab=1.5)
text(c(left,right),c(tumor_perc_left,tumor_perc_right), c(names(left),names(right)), cex=1, pos=4, col="black")
dev.off()


###### analysis of Sunday 24th of June 2018
### what is the difference between the 2 subgroups in clinical samples?

# analysis done only on the shared proteins
exp.moreStroma_2 <- TmergedALLExpr[c("(O)","(P)","(I)","(N)","(H)","(J)","(C)","(A)","(D)","(Q)","(F)","(M)"),]
exp.moreTumor_2 <- TmergedALLExpr[c("(E)","(R)","(B)","(G)","(S)","(T)","(K)","(L)"),]
colnames(exp.moreStroma) <-c("moreStroma(O)","moreStroma(P)","moreStroma(I)","moreStroma(N)","moreStroma(H)","moreStroma(J)","moreStroma(C)","moreStroma(A)","moreStroma(D)","moreStroma(Q)","moreStroma(F)","moreStroma(M)") 
colnames(exp.moreTumor) <- c("moreTumor(E)","moreTumor(R)","moreTumor(B)","moreTumor(G)","moreTumor(S)","moreTumor(T)","moreTumor(K)","moreTumor(L)")


#### check the number of zero
zero_freq_exp_stroma <- table(colSums(exp.moreStroma_2[,1:2697] == 0))

zero_freq_exp_tumor <- table(colSums(exp.moreTumor_2[,1:2697] == 0))

zero_proteins_stroma <- which(colSums(exp.moreStroma_2[,1:2697] == 0) == 0)

zero_proteins_tumor <- which(colSums(exp.moreTumor_2[,1:2697] == 0) == 0)

interesection_zero_proteins <- intersect(names(zero_proteins_stroma),names(zero_proteins_tumor))

pdf("zero_freq_expectation.pdf")
barplot(c(zero_freq_exp_stroma,zero_freq_exp_tumor),
        col = c(rep("yellow",length(zero_freq_exp_stroma)),rep("blue",length(zero_freq_exp_tumor))))
dev.off()


require(venneuler)
pdf("venn_overlap_zero.pdf",width = 15,height = 12)
fit <- euler(c(out_of_clust=650-633, in_clust=1220-633, "out_of_clust&in_clust"=633),shape = "ellipse")
plot(fit,fills=c("yellow","blue"))
dev.off()


##### check the blue part of the venn diagram, there are some proteins in the out_of_cluster that are not expressed in the in_clust!!

blue_proteins <- setdiff(names(zero_proteins_tumor),names(zero_proteins_stroma))
write.table(file="blue_proteins_to_string.txt",blue_proteins,quote = F,row.names = F,col.names = F)


########## differential analysis ########
library(limma)
left_c <- TmergedALLExpr[c("(A)","(C)","(D)","(F)","(H)","(I)","(J)","(M)","(N)","(O)","(P)","(Q)"),1:2697]
right_c <- TmergedALLExpr[c("(B)","(E)","(G)","(K)","(L)","(R)","(S)","(T)"),1:2697]
#prepare design matrix
design <- cbind(left=c(rep(1, 12), rep(0, 8)), right=c(rep(0, 12), rep(1,8))); design

DataforLimma <- as.data.frame(t(log10(TmergedALLExpr[c("(A)","(C)","(D)","(F)","(H)","(I)","(J)","(M)","(N)","(O)","(P)","(Q)","(B)","(E)","(G)","(K)","(L)","(R)","(S)","(T)"),1:2697]+1)))
#fit linear model, according to design matrix 
half_min <- min(DataforLimma[DataforLimma > 0])/2

DataforLimma2 <- DataforLimma+half_min
fit <- lmFit(DataforLimma2, design)

#prepare contrast matrix
cont.matrix <- makeContrasts( left.vs.right =  left -  right, levels=design)


# calculate
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
fit2
#check top candidates
topTable(fit2, adjust="BH")
topTable(fit2)

tab <- topTable(fit2, n=Inf)


with(tab, plot(tab$logFC, -log10(tab$P.Value), pch=20, main="Volcano plot", xlim=c(-10,10)))

with(subset(tab, P.Value < 0.05 ), points(logFC, -log10(P.Value), pch=20, col="red"))
# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(tab, -log10(tab$P.Value)>5 & abs(tab$logFC)>5), textxy(tab$logFC, -log10(tab$P.Value), labs=rownames(tab), cex=.8))

sign_right <- subset(tab,tab$P.Value < 0.05 & tab$logFC < 0)
sign_left<- subset(tab,tab$P.Value < 0.05 & tab$logFC > 0)

# map protein to gene names
sign_left_genes <- DataforLimma[rownames(DataforLimma) %in% rownames(sign_left),]
gene_inexp.stroma <- ClinicalExpr5gelLFQ[ClinicalExpr5gelLFQ$Protein %in% rownames(sign_left_genes),]$Gene

sign_right_genes <- DataforLimma[rownames(DataforLimma) %in% rownames(sign_right),]
gene_inexp.tumor <- ClinicalExpr5gelLFQ[ClinicalExpr5gelLFQ$Protein %in% rownames(sign_right_genes),]$Gene
write.table(file = "gene_up_inexp.tumor.txt", gene_inexp.tumor,quote = F,row.names = F,col.names = F)

stroma_LCM <- read.table("/Users/Giulia/Dropbox/LCM data and paper/Data/up_in_STROMA.txt",header = T)
tumor_LCM <- read.table("/Users/Giulia/Dropbox/LCM data and paper/Data/up_in_TUMOR.txt",header = T)
# rank the gene from the highest to the lowest p.val

######## check in LCM signature ### no correlation!
for_tumor <-tumor_LCM[tumor_LCM$Gene.names %in% gene_inexp.tumor, ]
for_stroma <-stroma_LCM[stroma_LCM$Gene.names %in% gene_inexp.stroma, ]


###### check on MOFFIT since it was done on bulk #####
StromaSG_normal <- c("ZNF469","VCAN","THBS2","SULF1","SPARC","SFRP2","POSTN","MMP11","LUM","ITGA11","INHBA","GREM1","FNDC1","FN1","FAP","CTHRC1","COMP","COL5A2","COL5A1","COL3A1","COL1A2","COL1A1","COL11A1","COL10A1","CDH11")
StromaSG_activated <- c("VIT","SYNM","SCRG1","RSPO3","RERGL","RBPMS2","PTX3","PLP1","OGN","MYH11","MEOX2","LPHN3","LMOD1","IGF1","ID4","GPM6B","FABP4","DES","CDH19","ANGPTL7","ADAMTS1","ACTG2","ABCA8")
stroma_MOFFIT <- c(StromaSG_normal,StromaSG_activated)
Tumor_basallike <- c("VGLL1","UCA1","S100A2","LY6D","SPRR3","SPRR1B","LEMD1","KRT15","CTSL2","DHRS9","AREG","CST6","SERPINB3","KRT6C","KRT6A","SERPINB4","FAM83A","SCEL","FGFBP1","KRT7","KRT17","GPR87","TNS4","SLC2A1","ANXA8L2")
Tumor_classical <- c("BTNL8","FAM3D","ATAD4","AGR3","CTSE","LOC400573","LYZ","TFF2","TFF1","ANZA10","LGALS4","PLA2G10","CEACAM6","VSIG2","TSPAN8","ST6GALNAC1","AGR2","TFF3","CYP3A7","MYO1A","CLRN3","KRT20","CDH17","SPINK4","REG4")
tumor_MOFFIT <- c(Tumor_basallike,Tumor_classical)

for_tumor <- intersect(tumor_MOFFIT,gene_inexp.stroma)
for_stroma <- intersect(stroma_MOFFIT,gene_inexp.tumor)

# try top 100 genes

moffitt_top100 <- read.xlsx("/Users/Giulia/Dropbox/LCM data and paper/SUBTYPING/moffitt subtype top 100 genes.xlsx")

stroma_MOFFIT <- c(moffitt_top100[,1],moffitt_top100[,2])
tumor_MOFFIT <- c(moffitt_top100[,3],moffitt_top100[,4])

for_stroma <- intersect(stroma_MOFFIT,gene_inexp.stroma)
for_tumor <- intersect(tumor_MOFFIT,gene_inexp.tumor)


#### try Collison

COLLISSON_classical <- c("TMEM4","SDR16C5","GPRC5A","AGR2","S100P","FXYD3","ST6GALNAC1","CEACAM5","CEACAM6","TFF1","TFF3","CAPN8","FOXQ1","ELF3","ERBB3","TSPAN8","TOX3","LGALS4","PLS1","GPX2","ATP10B","MUC13")
COLLISSON_quasimesenchymal <- c("AIM2","FAM26F","GPM68","S100A2","KRT14","CAV1","LOX","SLC2A3","TWIST1","PAPPA","NT5E","CKS2","HMMA","SLC5A3","PMAIP1","PHLDA1","SLC16A1","FERMT1","HK2","AKNAK2")
COLLISSON_exocriene<- c("REG1B","REG3B","REG1A","PRLIPRP2","CEL","PNLIP","PLA2G1B","CeLA3A","CPB1","CELA3B","CTRB2","CLPS","CELA2B","PRSS2","PRSS1","GP2","SLC3A1","CTRF","SLC4A4","SPINK1")

Collison <- c(COLLISSON_classical,COLLISSON_quasimesenchymal,COLLISSON_exocriene)
for_tumor <- intersect(Collison,gene_inexp.tumor)

##### probably biological reason???? ##########

#Hypothesis: if we remove the “stromal” contamination do we have a better cluster?

library(limma)
out_c <- TmergedALLExpr[c("(A)","(C)","(D)","(F)","(H)","(I)","(J)","(M)","(N)","(O)","(P)","(Q)"),1:2697]
in_c <- TmergedALLExpr[!rownames(TmergedALLExpr) %in% c("(A)","(C)","(D)","(F)","(H)","(I)","(J)","(M)","(N)","(O)","(P)","(Q)"),1:2697]
### are there some proteins that are completely zeros? which means that they were choosen only because the presence of the out samples
SUMCOL <- colSums(in_c)
min(SUMCOL)

#prepare design matrix
design <- cbind(out_c=c(rep(1, 12), rep(0, 40)), in_c=c(rep(0, 12), rep(1,40))); design

DataforLimma <- as.data.frame(t(rbind.data.frame(TmergedALLExpr[rownames(out_c),1:2697],TmergedALLExpr[rownames(in_c),1:2697])))
#fit linear model, according to design matrix 
half_min <- min(TmergedALLExpr[,1:2697][TmergedALLExpr[,1:2697] > 0])/2

DataforLimma2 <- DataforLimma+half_min
fit <- lmFit(log10(DataforLimma2), design)

#prepare contrast matrix
cont.matrix <- makeContrasts( out_c.vs.in_c = out_c - in_c, levels=design)


# calculate
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
fit2
#check top candidates
topTable(fit2, adjust="BH")
topTable(fit2)

tab <- topTable(fit2, n=Inf)

with(tab, plot(tab$logFC, -log10(tab$P.Value), pch=20, main="Volcano plot", xlim=c(-10,10)))

with(subset(tab, P.Value < 0.05 ), points(logFC, -log10(P.Value), pch=20, col="red"))

up_out_c <- subset(tab,tab$logFC > 0 & tab$P.Value < 0.05)

up_in_c <- subset(tab,tab$logFC < 0 & tab$P.Value < 0.05)


sub_Tmerged <- TmergedALLExpr[,! colnames(TmergedALLExpr) %in% rownames(up_out_c)]
sub_Tmerged[,1245]
####### try clustering??????

sampleTree <- hclust(dist(sub_Tmerged[,2435]), method = "average")
pdf("sample_clustering.pdf")
plot(sampleTree,labels=rownames(sub_Tmerged))
dev.off()



dist=dist(log10(sub_Tmerged[,1:2435]+1) , diag=TRUE)
hc=hclust(dist)
dhc=as.dendrogram(hc)
specific_leaf=dhc[[1]][[1]][[1]]
specific_leaf
attributes(specific_leaf)

i=0
colLab<<-function(n){
  if(is.leaf(n)){
    
    #I take the current attributes
    a=attributes(n)
    
    #I deduce the line in the original data, and so the treatment and the specie.
    ligne=match(attributes(n)$label,rownames(sub_Tmerged))
    batch=sub_Tmerged[ligne,2437];
    if(batch=="Batch1"){col_batch="blue"};if(batch=="Batch2"){col_batch="red"};if(batch=="Batch3"){col_batch="black"};if(batch=="Batch4"){col_batch="yellow"}
    datatype=sub_Tmerged[ligne,2436];
    if(datatype=="ATCC_cellLines"){col_datatype="red"};if(datatype=="SubcutanousPDX"){col_datatype="Darkgreen"};if(datatype=="OrthotopicPDX"){col_datatype="blue"};if(datatype=="PrimaryCellLines"){col_datatype="black"};if(datatype=="BulkTumors"){col_datatype="chocolate"}
    
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=1,pch=20,col=col_batch,lab.col=col_datatype,lab.font=1,lab.cex=1))
  }
  return(n)
}
dL <- dendrapply(dhc, colLab)### 
pdf("Integration/dendrogram3.pdf",width = 12)
plot(dL , main="Integration in Expression proteomics")
legend("topright", 
       legend = c("Exp1" , "Exp2" ,  "Exp3" , "Exp4" ,"ATCC_cellLines" , "SubcutanousPDX" , "OrthotopicPDX","PrimaryCellLines","BulkTumors"), 
       col = c("blue", "red" , "black" , "yellow" , "red", "Darkgreen","blue","black","chocolate"), 
       pch = c(20,20,20,20,4,4,4,4,4), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0.000001, 0.001))

dev.off()


#### try heatmap ####

heatmap.2(as.matrix(log10(sub_Tmerged[,1:1245]+1)),trace = "none",scale = "col")
##### correlation

corra <- cor(as.matrix(TmergedALLExpr[,1:2697]))
heatmap.2(as.matrix(corra),trace = "none")

#[#### check protein Identification level]
proteinID<- apply(DataforLimma[,1:20],2,function(x) sum(x > 0))
barplot(proteinID,col=c(rep("yellow",12),rep("blue",8)),names.arg = "",ylab = "per sample protein groups")
### nice!
#correlate with raw spectral counts

Clinical5gel_rawcount <- ClinicalExpr5gel[,15:34]
colnames(Clinical5gel_rawcount) <- sapply(strsplit(colnames(Clinical5gel_rawcount),"[.]"), `[`, 3)
raw_spectral_count_5gel <- colSums(Clinical5gel_rawcount[,c("(A)","(C)","(D)","(F)","(H)","(I)","(J)","(M)","(N)","(O)","(P)","(Q)","(B)","(E)","(G)","(K)","(L)","(R)","(S)","(T)")])

raw_intensiti5gel <- ClinicalExpr5gel[,35:54]
colnames(raw_intensiti5gel) <- sapply(strsplit(colnames(raw_intensiti5gel),"[.]"), `[`, 3)
raw_intensiti5gel_sort <- colSums(raw_intensiti5gel[,c("(A)","(C)","(D)","(F)","(H)","(I)","(J)","(M)","(N)","(O)","(P)","(Q)","(B)","(E)","(G)","(K)","(L)","(R)","(S)","(T)")])

barplot(raw_intensiti5gel_sort)
pdf("correlation_proteinID_spectralCount.pdf")
plot(proteinID,raw_spectral_count_5gel,pch=20,cex=2,col=c(rep("yellow",12),rep("blue",8)))
text(proteinID,raw_spectral_count_5gel, names(proteinID), cex=1, pos=4, col="black")
dev.off()
wo=""
for (i in 1:length(colnames(DataforLimma[,1:12]))){
  wo <- c(wo,rownames(DataforLimma[which(DataforLimma[,i] > 0),]))
  print(wo)
}

u_wo <- unique(wo)


fo=""
for (i in 1:length(colnames(DataforLimma[,13:20]))){
  fo <- c(fo,rownames(DataforLimma[which(DataforLimma[,i] > 0),]))
  print(fo)
}

u_fo <- unique(fo)


diff <- setdiff(u_wo,u_fo)

write.table(file ="diff_prot_groups.txt",gene_diff,row.names = F,col.names = F,quote=F)
gene_diff <- ClinicalExpr5gelLFQ[ClinicalExpr5gelLFQ$Protein %in% diff,]$Gene

