library(DESeq2)

`%notin%` <- Negate(`%in%`)

# the count matrix was generated with the following command run from within /work/5GB1c_rnaseq
# python barrelseq/shim.py extract --config-file barrelseq_config_gb.yaml --output data/extract_RAW_counts.tsv --values raw --format tsv
# read in the count matrix
countdata <- read.table("data/extract_RAW_counts.tsv", header=TRUE, row.names=1, sep='\t')
countdata <- countdata[,0:98]
countdata <- as.matrix(countdata)
#head(countdata)
#names(countdata)

# read in the sample matrix 
coldata <- read.table("data/5G_exp_metadata.tsv", header=TRUE, sep='\t', comment.char = '', row.names=3)
coldata <- coldata[,c(1, 2)]
#rownames(coldata) <- coldata[,1]
head(as.data.frame(coldata))

colnames(countdata)[colnames(countdata) %notin% rownames(coldata)]
all(colnames(countdata) %in% rownames(coldata))
all(rownames(coldata) %in% colnames(countdata))
all(rownames(coldata) == colnames(countdata))

# run the deseq2 estimators
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~exp_condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)
