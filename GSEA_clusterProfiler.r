# 1. Data Import
# 2. Renaming ensembl gene id versions to gene names
# 3. Loading gene sets from MSigDB databatse
# 4. Running GSEA 


# Clear environment
rm(list=ls(all.names = TRUE))

### ----------------------------- Libraries --------------------------------- ##
################################################################################
# Package to perform GSEA
library(clusterProfiler)
# Retrieving
library(msigdbr)
library(biomaRt)
# Data manuiplulation
library(tibble)
library(dplyr)





# 1. Data Import
################################################################################

# Load differential expression data
cancer_vs_normal <- readRDS("cancer_vs_normal_DElist.rds")
# Select log2FC and ensembl_id_version column
gene.list <- cancer_vs_normal$log2FoldChange



# 2. Renaming ensembl gene id versions to gene names
################################################################################

# Connect to biomart servers
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Retrieve external_gene_name of the genes using their ensembl_gene_id_version
gene_names <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id_version"),
  values = names(gene.list),
  filter = "ensembl_gene_id_version",
  mart = ensembl
)


# Merge gene.list and gene_symbols by ensembl_gene_id_version
gene.df <- gene.list %>% as.data.frame() %>% 
  {colnames(.) <- "log2FC";.} %>% 
  rownames_to_column(., var = "ensembl_gene_id_version") %>% 
  merge(., gene_names, by = "ensembl_gene_id_version")


# Filter out duplicated rows of the external_gene_name column
filt_gene.df <- gene.df %>% distinct(external_gene_name, .keep_all = T)
# Obtain log2FC with gene names
filt.list <- filt_gene.df$log2FC
names(filt.list) <- filt_gene.df$external_gene_name

# Sort gene list in a descending order
sorted.list <- filt.list %>% 
   sort(., decreasing = T)



# Ensure gene_list is a named numeric vector
if (!is.numeric(sorted.list) || is.null(names(sorted.list))){
  stop("sorted.list must be a names numeric vector.")
}


# 3. Loading gene sets from MSigDB databatse
################################################################################

# Obtain MSigDB pathways
msigdb <- msigdbr(species = "Homo sapiens")

# obtain msigdb
gda_custom <- data.frame(
  geneSymbol =  msigdb$human_gene_symbol,
  pathwayID = msigdb$gs_id,
  pathwatName = msigdb$gs_name
)

# Select pathways of interest
pathway_types <- c("KEGG_")

# filter pathways
filtered_gda_custom <- gda_custom[grepl(pathway_types, 
                                        gda_custom$pathwatName),]


# Construct the term objects as GSEA input
term2gene_custom <- filtered_gda_custom[,c("pathwayID", "geneSymbol")] 
term2name_custom <- filtered_gda_custom[,c("pathwayID", "pathwatName")]

# 4.  Running GSEA 
################################################################################

# Perform GSEA
gse <- GSEA(
  sorted.list, TERM2GENE = term2gene_custom, TERM2NAME = term2name_custom
)
# Convert gse into dataframe
gse.df <- gse %>% as.data.frame()
