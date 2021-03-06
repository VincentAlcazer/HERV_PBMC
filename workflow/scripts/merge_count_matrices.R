
#!/usr/bin/R
########################################################################
## merge_count_matrices.R
##
## hreyes Jan 2022
########################################################################
# Read in and merge two single-cell RNA-seq count matrices
# one is protein coding transcripts from starsolo
# two is transposable elements transcripts from stellarscope
#
# both matrices are sparse
########################################################################
#
#################### import libraries and set options ##################
library(optparse)
library(Matrix)
#
message("\nRequired libraries have been loaded.")
#
########################## functions ###################################
# it's dangerous to go alone! take this.
#
# function to whatever
#
########################## read in data ###################################
option_list = list(
  make_option(opt_str = c("-s", "--sample"),
              type = "character",
              help = "input sample name"),
  make_option(opt_str = c("-c", "--coding"), 
              type = "character",
              help = "directory of protein coding transcripts matrix"),
  make_option(opt_str = c("-t", "--transposable"),
              type = "character",
              help = "directory of transposable elements transcripts matrix"),
  make_option(opt_str = c("-o", "--output"), 
              type = "character", 
              help = "output directory for the merged matrix object")
)

opt <- parse_args(OptionParser(option_list=option_list))

opt$sample = snakemake@params[["sample"]]
opt$coding = snakemake@input[["protein_coding_dir"]]
opt$transposable = snakemake@input[["transposable_elements_dir"]]
# In this case I am making the output the whole file 
opt$output = snakemake@output[[1]]

if (is.null(opt$sample)){
  print_help(OptionParser(option_list=option_list))
  stop("\nThe input sample is mandatory", call.=FALSE)
}

##########################  read in protein coding matrix ########################## 
PC.mat <- readMM(file = paste0(opt$coding, "/matrix.mtx"))
message(paste0("\nProtein coding counts matrix is loaded: ", opt$coding, "/matrix.mtx"))

# read in protein coding feature names
PC.features <- read.delim(file = paste0(opt$coding, "/features.tsv"),
                          header = FALSE, stringsAsFactors = FALSE)

# read in protein coding barcodes
PC.barcodes <- read.delim(file = paste0(opt$coding, "/barcodes.tsv"),
                          header = FALSE, stringsAsFactors = FALSE)
# in the features
# V1 are ensemble gene IDs
# V2 are gene symbols
# I might only need to add a "ENSG" tag to the beginning of gene symbols, if we never use ensembl gene ID going forward
rownames(PC.mat) = paste0(PC.features$V1, "hrg", PC.features$V2)
colnames(PC.mat) = PC.barcodes$V1

rm(PC.features, PC.barcodes)
########################## read in transposable elements matrix ########################## 
TE.mat <- readMM(file = paste0(opt$transposable, "/", opt$sample, "-TE_counts.mtx"))
TE.mat <- t(TE.mat)
message(paste0("\nTransposable elements counts matrix is loaded: ", opt$transposable, "/", opt$sample, "-TE_counts.mtx"))

# read in transposable elements feature names
TE.features <- read.delim(file = paste0(opt$transposable, "/", opt$sample, "-features.tsv"),
                          header = TRUE, stringsAsFactors = FALSE)

# read in transposable elements barcodes
TE.barcodes <- read.delim(file = paste0(opt$transposable, "/", opt$sample, "-barcodes.tsv"),
                          header=TRUE, stringsAsFactors = FALSE)

# why is no feature the second one?
TE.features <- subset(x = TE.features, X0 != "__no_feature")

rownames(TE.mat) <- TE.features$X0
colnames(TE.mat) <- TE.barcodes$barcodes

rm(TE.features, TE.barcodes)
########################## merge both matrices ########################## 
# match the order of columns aka barcodes
TE.mat <- TE.mat[,colnames(PC.mat)]

PC.TE.mat <- rbind(PC.mat, TE.mat)

message("\nThe two matrices have been merged")

#rm(PC.mat, TE.mat)
########################## save object ########################## 
# save object
#saveRDS(object = PC.TE.mat, file = paste0(opt$output, opt$sample, "_matrix.Rds"))
saveRDS(object = PC.TE.mat, file = opt$output)
#
# should I save the matrix in text format (three files)? get a flag in opts
message(paste0("\nMerged matrix is saved to: ", opt$output))
