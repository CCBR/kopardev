###############################################################################################################
# Skyler Kuhn
# ccbr1003
# Calculates geneset enrichment via over-representation analysis (ORA)
# Methods: ORA (Hypergeometric test)
# USAGE: Rscript ccbr1003_ORA.R -d Acidovorax_expression_correlation_ALL_SAMPLES.txt -o '/path/to/output/' -c 'output_prefix_test-tmp'
###############################################################################################################
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library(clusterProfiler))
suppressMessages(library(argparse))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(ReactomePA))
suppressMessages(library(ggplot2))
######################################################################
# Functions
######################################################################
# Saves image as pdf (requested format)
create_pdf <- function(image, filename, x, y) {
  "Saves an image as a pdf"
  pdf(filename, width=x, height=y)
  image
  dev.off()
  return()
}
# Basic wrapper for write.table() function, add logic to take a list of df as input
write2file <- function(data, outputfilename, include_rownames=TRUE){
  "Writes data to output file"
  write.table(data, file = outputfilename, sep = "\t", row.names = include_rownames)
}
ensembl2entrez <- function(dge_results, organism, write = FALSE){
  "Converts EnsemblIDs to ENTREZID"
  # For mapping EnsemblID to ENTREZID 
  ensembl2ENTREZID <- bitr(dge_results$Gene, fromType = "ENSEMBL", toType = c("ENTREZID"), OrgDb = organism)
  # Rename column ENSEMBL to Gene before join on column name
  colnames(ensembl2ENTREZID)[colnames(ensembl2ENTREZID)=="ENSEMBL"] <- "Gene"
  # Saving mapping file to output file
  if (write){
    write2file(ensembl2ENTREZID, paste("ensembl2entrez_", gsub('\\.', '-', organism$packageName), ".txt", sep = ""), include_rownames = FALSE)
  }
  deg_res_entrez <- merge(dge_results, ensembl2ENTREZID, by="Gene")
  return(deg_res_entrez)
}
ensembl2uniprot <- function(dge_results, organism, write = FALSE){
  "Converts EnsemblIDs to UNIPROTID"
  # For mapping EnsemblID to UNIPROTID 
  ensembl2UNIPROTID <- bitr(dge_results$Gene, fromType = "ENSEMBL", toType = c("UNIPROT"), OrgDb = organism)
  # Rename column ENSEMBL to Gene before join on column name
  colnames(ensembl2UNIPROTID)[colnames(ensembl2UNIPROTID)=="ENSEMBL"] <- "Gene"
  # Saving mapping file to output file
  if (write) {
    write2file(ensembl2UNIPROTID, paste("ensembl2uniprot_", gsub('\\.', '-', organism$packageName), ".txt", sep = ""), include_rownames = FALSE)
  }
  deg_res_uniprot <- merge(dge_results, ensembl2UNIPROTID, by="Gene")
  return(deg_res_uniprot)
}
# Filters DEG results based on adjusted p-value and FC direction
directionality <- function(deg_results, pvalue_threshold = 0.05 , rho_threshold = 0){
  "Returns filtered DGE list based on pvalue and FC direction"
  # Significant up and down regulated DEG results
  sig_deg = filter(deg_results, pvalue <= pvalue_threshold & (rho >= rho_threshold | rho <= -rho_threshold))
  # Singificant up-regulated DEG results
  sig_up = filter(deg_results, pvalue <= pvalue_threshold & rho >= rho_threshold)
  # Singificant down-regulated DEG results
  sig_down = filter(deg_results, pvalue <= pvalue_threshold & rho <= -rho_threshold) 
  return(list(both=sig_deg, up=sig_up, down=sig_down))
}
# Running GO Enrichment Analysis of the significant DE genese
geneOntologyEnrichmentAnalysis <- function(genelist, organism, background_genelist, pvalue_threshold=1.0, fdr_threshold=1.0){
  "Returns list of enriched genesets from an over-representation test for the following sub-ontologies: BP, MF, CC"
  # Cellular Compartment
  goCC <- enrichGO(gene=genelist, OrgDb=organism, universe = background_genelist, keyType='ENSEMBL', ont="CC", 
                   pAdjustMethod = "BH", pvalueCutoff  = pvalue_threshold, qvalueCutoff  = fdr_threshold)
  # Biological Process
  goBP <- enrichGO(gene=genelist, OrgDb=organism, universe = background_genelist, keyType='ENSEMBL', ont="BP", 
                   pAdjustMethod = "BH", pvalueCutoff  = pvalue_threshold, qvalueCutoff  = fdr_threshold)
  # Molecular Function
  goMF <- enrichGO(gene=genelist, OrgDb=organism, universe = background_genelist, keyType='ENSEMBL', ont="MF", 
                   pAdjustMethod = "BH", pvalueCutoff  = pvalue_threshold, qvalueCutoff  = fdr_threshold)
  # Converting Geneset descriptions to upper case for plotting
  goCC@result$Description <- toupper(goCC@result$Description)
  goBP@result$Description <- toupper(goBP@result$Description)
  goMF@result$Description <- toupper(goMF@result$Description)
  return(list(CC=goCC, BP=goBP, MF=goMF))
}
######################################################################
# Main Method
######################################################################
main <- function(workdir, ora_input, mycontrast, myseed=1234) {
  # Set seed
  set.seed(myseed)
  filename = ora_input  # 'Acidovorax_expression_correlation_ALL_SAMPLES.txt'
  Contrast = mycontrast # 'Test-tmp'
  # Reading in DEG results and removing genename from Gene (where "ENSEMBLID|GENENAME")
  deg = read.table(filename, header = TRUE)
  deg$Gene = gsub('\\..*','', deg$Gene)  # regular-expression to remove everything after the '.'
  # Set working directory
  setwd(workdir)  # '/Users/kuhnsa/Desktop/Archive/ccbr1003/PCAs/Chenran_Samples'
  #####################################
  # Over-representation Analysis (ORA)
  #####################################
  cat("\nOver-representation Analysis (ORA)\n")
  # Initialize Output Directory (ORA)
  dir.create(file.path(getwd(), "ORA"), showWarnings = FALSE)
  # Background geneset for over-representation test
  background = deg$Gene
  # Filtered List of Significant Differential DGE Results for GOEA
  ## where: both = up and down regulated; up = up-regulated; down = down-regulated 
  sig_res <- directionality(deg)
  ###################
  # C5: BP, MF, CC
  ###################
  # Running GO Enrichment Analysis of the significant DE geneset for each sub-ontology
  ## Cellular Compartment (CC), Biological Process (BP), Molecular Function (MF)
  cat("Running ORA for C5: BP, MF, CC\n")
  goea_res <- geneOntologyEnrichmentAnalysis(genelist = sig_res$both$Gene, organism = org.Hs.eg.db, background_genelist = background)
  # Saving results and figures for GO term over-representation tests
  ## Cellular Compartment
  write2file(goea_res$CC@result, paste("ORA/", Contrast, "_C5CC_GO_results.txt", sep=""), include_rownames = FALSE)
  mypdf <- dotplot(goea_res$CC, showCategory=15, color="qvalue")
  ggsave(filename = paste("ORA/", Contrast, "_C5CC_GO_dotplot.pdf", sep=""), height= 8.90, width = 12.80, device = "pdf", plot = mypdf)
  ## Biological Process
  write2file(goea_res$BP@result, paste("ORA/", Contrast, "_C5BP_GO_results.txt", sep=""), include_rownames = FALSE)
  mypdf <- dotplot(goea_res$BP, showCategory=15, color="qvalue")
  ggsave(filename = paste("ORA/", Contrast, "_C5BP_GO_dotplot.pdf", sep=""), height= 8.90, width = 12.80, device = "pdf", plot = mypdf)
  ## Molecular Function 
  write2file(goea_res$MF@result, paste("ORA/", Contrast, "_C5MF_GO_results.txt", sep=""), include_rownames = FALSE)
  mypdf <- dotplot(goea_res$MF, showCategory=15, color="qvalue")
  ggsave(filename = paste("ORA/", Contrast, "_C5MF_GO_dotplot.pdf", sep=""), height= 8.90, width = 12.80, device = "pdf", plot = mypdf)  
  ###################
  # C2CP: REACTOME
  ###################
  cat("Running ORA for C2CP: Reactome\n")
  # Running ORA on Reactome genesets
  ## Significant up and down regulated genes
  sig_res_entrez <- ensembl2entrez(sig_res$both, organism = org.Hs.eg.db)
  ## Background genes set (everything)
  background_entrez <- ensembl2entrez(deg, organism = org.Hs.eg.db)
  # Using ReactomePA to find enriched genesets
  goea_c2_cp <- enrichPathway(gene=sig_res_entrez$ENTREZID, organism = 'human', universe = background_entrez$ENTREZID, 
                              pvalueCutoff=1.0, qvalueCutoff=1.0, pAdjustMethod = "fdr", readable=T, minGSSize = 20)
  # Writing results and figures to files
  write2file(goea_c2_cp@result, paste("ORA/", Contrast, "_C2CP_REACTOME_results.txt", sep=""), include_rownames = FALSE)
  mypdf <- dotplot(goea_c2_cp, showCategory=15, color="qvalue")
  ggsave(filename = paste("ORA/", Contrast, "_C2CP_REACTOME_dotplot.pdf", sep=""), height = 8.90, width = 12.80, device = "pdf", plot = mypdf)  
  ###################
  # C2CP: KEGG
  ###################
  cat("Running ORA for C2CP: KEGG\n")
  # Running ORA on KEGG genesets
  ## Significant up and down regulated genes (must convert to valid geneid ~ entrez)
  kegg_res_ora <- enrichKEGG(gene=sig_res_entrez$ENTREZID, organism = 'hsa', universe = background_entrez$ENTREZID, 
                             pvalueCutoff=1.0, qvalueCutoff=1.0, pAdjustMethod = "fdr", minGSSize = 10)
  # Saving KEGG ORA results
  write2file(kegg_res_ora@result, paste("ORA/", Contrast, "_C2CP_KEGG_results.txt", sep=""), include_rownames = FALSE)
  mypdf <- dotplot(kegg_res_ora, showCategory=15, color="qvalue")
  ggsave(filename = paste("ORA/", Contrast, "_C2CP_KEGG_dotplot.pdf", sep=""), height = 8.90, width = 12.80, device = "pdf", plot = mypdf)  
}
#####################################
# Argument Parsing
#####################################
# Pass command line args to main
# create parser object
parser <- ArgumentParser()
# Adding command-line arguments
parser$add_argument("-d", "--input_data", type="character", required=TRUE,
                    help="File containing the list of significant correlated genes and pvaules <REQUIRED>")
parser$add_argument("-o", "--output_dir", type="character", required=TRUE,
                    help="Output working diretory <REQUIRED>")
parser$add_argument("-c", "--contrast", type="character", required=TRUE,
                    help="Differential Expression Analysis Contrast. String is used as a prefix for output filenames <REQUIRED>")
args <- parser$parse_args()
main(workdir = args$output_dir, ora_input = args$input_data, mycontrast = args$contrast)
