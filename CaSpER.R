#title: Copy Number Alterations using CaSpER on bulk RNA-Seq
#authors: Nayanga Thirimanne

library(sva)
library(CaSpER)
library(GenomicRanges)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(edgeR)
install.packages("survival")
library(survival)
install.packages("ggsurvfit")
library(ggsurvfit)
install.packages("survminer")
library(survminer)

#ggplot2 specifications
plot_title_size = 20
legend_pt_size = 4
axis_text_size = 25
axis_title_size = 25
legend_text_size = 15
spacing = 1
chosen_margin = c(0.5,1,0.5,1) #top,right,bottom,left

theme_legend <- theme_blank()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "right",
        legend.justification = "left",
        legend.title = element_blank() )

theme_nolegend <- theme_blank()+
  theme(plot.title = element_text(hjust=0, vjust=0, lineheight=.8, face="bold", size=plot_title_size),
        plot.margin = unit(chosen_margin,"cm"),
        legend.text = element_text(size=legend_text_size),
        legend.key.height = unit(spacing,"cm"),
        legend.position = "none",
        legend.justification = "left",
        legend.title = element_blank() )



####
#Run CaSpER 
####

# Here we are using hippocampus and frontol cortex normal tumors as controls
# For normal samples RNA-Seq data was downloaded from GTEX and processed using the same pipelline

#normalize and batch correct meningioma datasets along with GTEX samples
df1 <- readRDS("~/Meningioma_dataset1_538/GSE139652_10/GSE139652_10samples_proteincodingGencode_rawcounts.rds")
df2 <- readRDS("~/Meningioma_dataset1_538/GSE136661_160/GSE136661_159samples_proteincodingGencode_rawcounts.rds")
df3 <- readRDS("~/Meningioma_dataset1_538/Holland_cohort_1_101/Hollandcohort1_93samples_proteincodingGencode_rawcounts.rds")
df4 <- readRDS("~/Meningioma_dataset1_538/Holland_cohort_2_108/Hollandcohort2_108samples_proteincodingGencode_rawcounts.rds")
df5 <- readRDS("~/Meningioma_dataset1_538/Felix_42/Felix_44samples_proteincodingGencode_rawcounts.rds")
df6 <- readRDS("~/Meningioma_dataset1_538/Toronto_124/Toronto_123samples_proteincodingGencode_rawcounts.rds")
df7 <- readRDS("~/Meningioma_RNASeq_dataset2/GSE101638/GSE101638_42samples_proteincodingGencode_rawcounts.rds")
df8 <- readRDS("~/Meningioma_RNASeq_dataset2/GSE183653/GSE183653_185samples_proteincodingGencode_rawcounts.rds")
df9 <- readRDS("~/Meningioma_RNASeq_dataset2/GSE151921/GSE151921_13samples_proteincodingGencode_rawcounts.rds")
df10 <- readRDS("~/Meningioma_RNASeq_dataset2/GSE85133/GSE85133_19samples_proteincodingGencode_rawcounts.rds")
df11 <- readRDS("~/Meningioma_RNASeq_dataset2/PRJNA705586/PRJNA705586_70samples_proteincodingGencode_rawcounts.rds")
df12 <- readRDS("~/Meningioma_RNASeq_dataset2/CAVATICA_CHOP_analysis/CAVATICA_25samples_proteincodingGencode_rawcounts.rds")
df13 <- readRDS("~/Meningioma_RNASeq_dataset2/HollandUW_3/HollandUW_3_78samples_wodups_rawcounts_proteincoGencode_013123.rds")
df14 <- readRDS("~/Meningioma_RNASeq_dataset2/GSE212666_HKU/GSE212666_301samples_proteincodingGencode_rawcounts.rds")
df15 <- readRDS("~/Meningioma_RNASeq_dataset2/CBTN_youngadults/CBTN_28samples_proteincodingGencode_rawcounts.rds")

# GTEX frontal cortex =30 and GTEX hippocampus =8
df16 <- readRDS("~/GTEX/GTEX_30frontalcor_proteincodingGencode_rawcounts.rds")
df17 <- readRDS("~/GTEX/GTEX_8hippocamp_proteincodingGencode_rawcounts.rds")

DF <- cbind(df1, df2, df3, df4, df5, df6, 
            df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17)
DF <- as.data.frame(DF)

batch <- c(rep(1, ncol(df1)), 
           rep(2, ncol(df2)),  
           rep(3, ncol(df3)), 
           rep(4, ncol(df4)),
           rep(5, ncol(df5)),
           rep(6, ncol(df6)),
           rep(7, ncol(df7)),
           rep(8, ncol(df8)),
           rep(9, ncol(df9)),
           rep(10, ncol(df10)),
           rep(11, ncol(df11)), 
           rep(12, ncol(df12)),
           rep(13, ncol(df13)), 
           rep(14, ncol(df14)), 
           rep(15, ncol(df15)),
           rep(16, ncol(df16)),
           rep(17, ncol(df17))
)

#batch correction
adjusted_CNV <-  ComBat_seq(as.matrix(DF), batch=batch, group=NULL)

#normalize read counts 
gtf <- rtracklayer::import("~/Reference_Genomes/hg38_release39_053122/gencode.v39.primary_assembly.annotation.gtf")
hg38gtf <- as.data.frame(gtf)
hg38gtf_sub <- subset(hg38gtf, gene_name %in% rownames(adjusted_CNV)) 
hg38gtf_sub2 <- dplyr::filter(hg38gtf_sub, hg38gtf_sub$gene_type == "protein_coding" & hg38gtf_sub$type == "gene")
hg38gtf_sub2 <- hg38gtf_sub2 %>% arrange(gene_name) %>%
  dplyr::filter(duplicated(gene_name) == FALSE)
gene_width = hg38gtf_sub2$width
gene_name = rownames(adjusted_CNV)
y <- DGEList(counts = adjusted_CNV)
rpkm = rpkm(y, gene.length = gene_width)
rownames(rpkm) = gene_name
tpm = apply(rpkm, 2, function(x){
  (x/sum(x))*10^6
})
log2_tpm_ar =log2(tpm+1)

#### 
df50 <- log2_tpm_ar
#drop samples that doesn't have CaSpER calls
drop <- c("SRR3995987", "SRR13810536", "H486", "1288", "SRR1319539", "SRR1342045", "SRR1388305", "SRR1402900", "SRR1433796")
df50.1 <- df50[, !(colnames(df50) %in% drop)]
cnvtestdata <- as.matrix(df50.1)

#Cytoband information
cytoband <- read.delim("~/cytoBand_hg38_ucsc.txt", header = FALSE)
cytoband <- data.frame(V1=gsub("chr", "", cytoband[,1]), V2=cytoband[,2], V3=cytoband[,3], V4=substring(cytoband$V4, 1, 1), stringsAsFactors=F)
start <- do.call(rbind, lapply(split(cytoband$V2, paste0(cytoband$V1, cytoband$V4)), min))
end <- do.call(rbind, lapply(split(cytoband$V3, paste0(cytoband$V1, cytoband$V4)), max))
cytoband <- data.frame(V1=gsub("p", "", gsub("q", "", rownames(start))), V2=start, V3=end, V4=rownames(start), stringsAsFactors=F)
cytoband <- cytoband [as.vector(unlist(sapply(c(1:22, "X"), function(x) which(cytoband$V1 %in% x)))), ]
cytoband$V4[grep("q", cytoband$V4)] <- "q"
cytoband$V4[grep("p", cytoband$V4)] <- "p"
rownames(cytoband) <- NULL

#Centromere information
centromere <- read.delim("~/centromere_hg38.txt", header = FALSE)

#control sample ids
controls <- read.csv("~/CaSpER_test/control.sample.id.csv", header = TRUE)
control.sample.ids <- controls$controls_GTEX

#Annotation data (positions of each gene along each chromosome)
annotation <- CaSpER::generateAnnotation(id_type="hgnc_symbol", genes=rownames(df50.1), ishg19 = F, centromere, host="useast.ensembl.org")

#subset data to only have genes as same as in annotation file
cnvtestdata <- cnvtestdata[match( annotation$Gene,rownames(cnvtestdata)), ]

#Reading BAFExtract output files
#baf files of all meningioma samples were in the compile_menin folder
loh <- CaSpER::readBAFExtractOutput(path="~/CaSpER_Meningioma/compile_menin/", sequencing.type = "bulk", suffix = "baf")

#generate dataframe of loh.name.mapping
files <- list.files(pattern = ".baf")
loh.name <- c(files)
sample.name <- c(files)
loh.name.mapping <- data.frame(loh.name, sample.name)

#creating CaSpER object
object <- CaSpER::CreateCasperObject(raw.data = cnvtestdata, loh.name.mapping = loh.name.mapping, sequencing.type = "bulk", cnv.scale =3, loh.scale = 3, matrix.type = "normalized",
                                     expr.cutoff = 0.1, annotation = annotation, method = "iterative", loh=loh, filter = "median", log.transformed = TRUE,
                                     control.sample.ids = control.sample.ids, cytoband = cytoband, genomeVersion = "hg38")

#pairwise comparison of scales from BAF and expression signals
final.objects <- CaSpER::runCaSpER(object, removeCentromere = T, cytoband = cytoband, method = "iterative")

#Large scale CNV summarization
finalChrMat <- CaSpER::extractLargeScaleEvents(final.objects, thr = 0.75)

#Visualization 
obj <- final.objects[[9]]
CaSpER::plotLargeScaleEvent2(chrMat = finalChrMat, fileName = "large.scale.events.1328samples_ChrMat_050223.pdf")
dev.off()

#Visualize only GTEX samples
finalChrMat_GTEX <- finalChrMat[rownames(finalChrMat) %in% controls$controls_GTEX, ]
CaSpER::plotLargeScaleEvent2(chrMat = finalChrMat_GTEX, fileName = "large.scale.events.1328samples_GTEX_ChrMat_050223.pdf")


