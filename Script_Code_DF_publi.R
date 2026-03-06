#Script containing the code for "Heterochronic shift in gene expression leads to larval onset of spermatogenesis and miniaturized males" Rouan et al.

#load Rdata 
load("path/Data_Code_DF_publi.RData")

#contains:
#annotation transcriptome annotation
#log2.cpm.filtered.norm.df_Male
#targets sample design
#color_num corresponding cluster number to cluster colors
#module_new_num a df containing the gene_id, averaged log-cpm for a sex and satge category and the associated cluster color and numbers
#gene2GO corresponding GO terms to Gene_ID (Trinity id) obtained from annotation
#GO4 the Top20 gene ontology terms and pvalue for the 4 clusters shown in Figure5

#load packages
library(stringr)
library(stringi) 
library(tidyr)
library(WGCNA)
library(tidytext)
library(dplyr)
library(ggplot2)
library(svglite)
library(topGO)

#To perform the clustering analysis WGCNA
## Step 1
###cluster only on males larva and detect developmental heterochronies: isolate the male sample from D2 to D5 larva from the TMM transformed cpm
log2.cpm.filtered.norm.df_Male <- log2.cpm.filtered.norm.df

#Rename the column using the combination of stage and sex
colnames(log2.cpm.filtered.norm.df_Male) <- targets$group_sex

#Now select male alrval stage
log2.cpm.filtered.norm.df_MaleL <- log2.cpm.filtered.norm.df_Male[,which(colnames(log2.cpm.filtered.norm.df_Male)%in% c("D2Larva_M" ,"D3Larva_M","D4Larva_M","D5Larva_M"))]

#transpose your data for WGCNA
W_input <-t(log2.cpm.filtered.norm.df_MaleL)
W_input <- as.data.frame(W_input)

#Load the package

allowWGCNAThreads()     

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  W_input,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

#Plot the results to choose the power treshold
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)

#Plot again with the power=13 with corresponding softreshold R square of 0.768
abline(h = 0.768, col = "blue")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

#Now cluster your genes based on a power=13
picked_power = 13
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(W_input,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
#Now get the module 
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

#Save your data 

## Step 2
#Plot the gene expression for 4 selected cluster (Figure 5) from the 16 showing early offset heterochronies in males
#Cluster were initially encoded as colored
#They were changed to numbers, see color_num df for correspondence
#module_new_num contains the cluster numbers

#module_new_num contains descriptive variables from targets and average log cpm values
#obtained as followed:
#library(limma)
#log2.cpm.filtered.norm.df_Male.AVG <- avearrays(log2.cpm.filtered.norm.df_Male)
#module_new$log2_cpm_norm <- sapply(1:nrow(module_new),function(n){ 
#  g <- log2.cpm.filtered.norm.df_Male.AVG[which(rownames(log2.cpm.filtered.norm.df_Male.AVG)%in%module_new[n,"GeneID"]),which(colnames(log2.cpm.filtered.norm.df_Male.AVG)%in%module_new[n,"group_gender"])]
#  return(g) })

#List of the 4 clusters displayed in Figure 5 encoded with colors
target <- c("sienna3","yellowgreen","skyblue","blue")

#List of the remaining 14 clusters showing male heterochronies
remaining_target <- c("tan","violet","lightcyan1","skyblue3", "darkred","paleturquoise","bisque4","darkslateblue","pink","brown","navajowhite2","salmon")

#Get the subset of 4 target clusters
module_4 <- module_new_num[which(module_new_num$cluster %in% target),]

#Plot
module4 <- ggplot(aes(x =factor(stage,levels = c( "D2Larva","D3Larva","D4Larva","D5Larva","juvenile","adult") ), y = log2_cpm_norm, color= factor(gender),group=gender, fill= factor(gender)), data = module_4) + 
  stat_smooth(method = "loess") +
  #geom_line(aes(x =factor(stage, levels = c("D2Larva","D3Larva","D4Larva","D5Larva","juvenile","adult")), y = log2_cpm_norm), size=0.5,  alpha = 0.8)+
  #scale_fill_manual(values=c("#a31d20ff","#0C459A"))+
  scale_color_manual(values=c("#a31d20ff","#0C459A"))+
  # labs(title="16 clusters showing an early onset and/or early offset", subtitle = "WGCNA based on Male larvae expression ")+
  xlab("Stages")+
  facet_wrap(~factor(num, levels=c("cluster1","cluster2","cluster3","cluster4")), ncol =1)+
  theme_light()+theme(legend.position = "none", axis.text.x=element_text(size=16, angle =+45, face= "bold") ,strip.background = element_rect(fill = "white"),axis.text.y=element_text(size=16), strip.text.x.top = element_text(size= 16, colour ="grey29", face = "bold"))

#Generate a Top20 GO enrichment matrix
#Step 1

###Step TOP20
#use the generated list of correspondence between gene_id and GO terms (produced from from annotation) 

#Make a for loop to apply it to several cluster 
for (L in target){ 
  gene_universe <- rownames(log2.cpm.filtered.norm.df_MaleL)
  module_genes <- module_new_num[which(module_new_num$cluster == paste(L)), "gene_id"]
  geneList <- factor(as.integer(gene_universe %in% module_genes))
  names(geneList) <- gene_universe
  gene2GO_filtered <- gene2GO[gene_universe]
  
  my_go_data_up <- new("topGOdata", 
                       description = paste("GOtest_BP"),
                       ontology = "BP", 
                       allGenes = geneList,
                       geneSel = function(p) p == 1,
                       annot = annFUN.gene2GO,
                       gene2GO = gene2GO_filtered)
  ##Top20 table
  top20 <- GenTable(my_go_data_up, classicFisher = runTest(object= my_go_data_up, algorithm = "weight01",statistic = "fisher"), topNodes = 20)
  #save the output
  n <-paste(L,"WGCNA_module")
   write.table(top20, row.names=F, file = paste("path/GO_enr_Top20_",n[1],"_","BP.csv", sep = "_"), sep = ",",dec = ".", fileEncoding = "UTF-8")  
}

#Now import the tables and make the Top20 GO plot from Figure4
#use GO4 a fused table of the 4 cluster of interest with simplified terms

GO4_filtered <- GO4 %>%
  group_by(cluster) %>%
  mutate(
    Term_wrapped = ifelse(str_count(Full_Term, "\\S+") > 2,
                          str_wrap(Full_Term, width = 40),
                          Full_Term),
    Term_wrapped = reorder(Term_wrapped, log_p)
  ) %>%
  ungroup()

ggplot(GO4_filtered) +
  geom_col(aes(x = log_p, y = reorder_within(Term_wrapped, log_p, num), fill = num), show.legend=F) +  
  scale_y_reordered() +  # very important!
  scale_fill_manual(values=c("#a020f0ff","#a020f0ff","#a020f0ff","#a020f0ff"))+
  facet_wrap(~factor(num), scales = "free", ncol = 1, drop = TRUE) +
  labs(
    title = "Top20 GO Terms by Cluster",
    x = expression(-log[10](italic(p))),
    y = "GO Term"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8),axis.text = element_text(size = 8),
        strip.text = element_text(size = 8))
