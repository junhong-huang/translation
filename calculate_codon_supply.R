rm(list = ls())
library(tAI)

extract_cod <- function (trnas, anticod){
  output = data.frame(row.names = anticod)
  trnas_acod = sapply(rownames(trnas), function(x) substr(x,nchar(x)-2,nchar(x)))#Extract the last three characters.
  
  for (s in colnames(trnas)){
    output[,s] = sapply(anticod, function(x) if(any(trnas_acod==x)){mean(trnas[trnas_acod==x,s])}else{0})
  }
  
  return(output)
}

AAnormalize <- function(data,codons){
  # Compute relative data
  aa = sapply(codons,function(x) substr(x,1,nchar(x)-3))
  uniqueaa = unique(aa)
  outdata = numeric(length=length(data))
  for (n in uniqueaa){
    idx = (aa %in% n)
    total = max(data[idx],na.rm=T)
    outdata[idx] = data[idx]/total
    if (total %in% 0){
      outdata[idx] = 1.0/sum(idx)
    }
  }
  return(outdata)
}

transformdata <- function(data,transf){
  # aa_idx = regexpr("i?[A-Z][a-z]{2}[A-Z]{3}",rownames(data))==1
  # data = data[aa_idx,]
  if (transf=="log"){
    outdata = sapply(data,log)
    # Remove inf values
    outdata[outdata==-Inf] = NaN
    rownames(outdata)=rownames(data)
  }else if (transf=="arcsinh"){
    outdata = sapply(data,asinh)
    rownames(outdata)=rownames(data)
  }else if (transf=="sqrt"){
    outdata = sapply(data,sqrt)
    rownames(outdata)=rownames(data)
  }else if (transf=="rel"){
    # Compute relative data
    outdata = data.frame(matrix(ncol = ncol(data), nrow = nrow(data)),row.names = rownames(data))
    colnames(outdata)= colnames(data)
    aa = sapply(rownames(outdata),function(x) substr(x,1,nchar(x)-3))
    uniqueaa = unique(aa)
    for (n in uniqueaa){
      idx = (aa %in% n)
      idx_data = matrix(as.matrix(data[idx,]), ncol = ncol(data), nrow = sum(idx))
      total = colSums(idx_data)
      outdata[idx,] = t(apply(idx_data,1,function(x) x/total))
      iszero = (total %in% 0)
      if (any(iszero)){
        outdata[idx,iszero] = 1.0/sum(idx)
      }
    }
  }else{
    outdata=data
  }
  return(outdata)
}

## Load trna and weighted CU
# Codon table
codons = read.table("/public/home/huangjh/perlBioTools/Patho-DBiT/tRNA_viruses/data/codons_table.txt", sep="\t", header=TRUE, row.names = 1, check.names = FALSE)
trna = read.table("/public1/home/huangjh/spatial_omics/Patho-DBiT/processingData/Human_AITL_20um/translation/expression_tRNA_by_spot.tsv", sep="\t", header=TRUE, row.names = 1, check.names = FALSE)
codus = read.table("/public1/home/huangjh/spatial_omics/Patho-DBiT/processingData/Human_AITL_20um/translation/frequency_mRNA_by_codon.tsv",sep="\t", row.names = 1, header=TRUE, check.names = FALSE)
codon_demand <- read.table("/public1/home/huangjh/spatial_omics/Patho-DBiT/processingData/Human_AITL_20um/translation/demand_codon_by_spot.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)


# tRNAs
cells = colnames(trna)

# Compute mean of tissue
trna_mean = data.frame(row.names = rownames(trna))
for (t in unique(cells)) {
  trna_mean[, t] = rowMeans(trna[, cells == t, drop = FALSE], na.rm = TRUE)
}

### tRNA: a vector of length 64 with tRNA gene copy numbers
anticodon = extract_cod(transformdata(trna_mean,"sqrt"), codons$ANTICODON)


# Genomic codon usage
# Keep only columns with codon info
# codus_clean = data.frame(sapply(rownames(codus),function(x) as.numeric(codus[x,1:ncol(codus)])), row.names = colnames(codus)[1:ncol(codus)])
codus_clean=t(codus)

rownames(codus_clean) = sapply(rownames(codus_clean),function(x) paste(codons[x,"AA"],x,sep=""))


print("Calculate tAI for each codon...")
codon = extract_cod(transformdata(codus_clean,""), rownames(codons)[!(codons$AA %in% c("Stop","Met"))])

TAIs = data.frame(matrix(ncol = ncol(anticodon), nrow = ncol(codon)),row.names = colnames(codon)); colnames(TAIs) = colnames(anticodon)

### 0.5 for I-T/C/A; 0.75 for GU
# initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)
initial_s = c(0, 0, 0, 0, 0.5, 1, 1, 0.5, 1)
all_ws <- data.frame(row.names = rownames(codon))

write.table(anticodon, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, "/public1/home/huangjh/spatial_omics/Patho-DBiT/processingData/Human_AITL_20um/translation/anticodon.tsv")

# Calculate tAI foreach sample
for (sample in colnames(anticodon)){
  # Calculate relative adaptiveness values (ws)
  sample.ws = get.ws(tRNA=anticodon[,sample], s=initial_s, sking=0)
  
  sample.ws = AAnormalize(sample.ws,paste0(codons[rownames(codon),"AA"],rownames(codon)))# normalize by AA
  all_ws[,sample] = sample.ws
  
  # Calculate tAI for all CUs
  # sample.tai <- get.tai(t(codon), sample.ws)
  # TAIs[,sample] = sample.tai
}

write.table(all_ws, file = "/public1/home/huangjh/spatial_omics/Patho-DBiT/processingData/Human_AITL_20um/translation/supply_codon_by_spot.tsv", sep = "\t", quote = FALSE, col.names = NA)

# TAIs[,c("annotation","Accession","Species")] = codus[,c("annotation","Accession","Species")]
# write.table(TAIs, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, "/public1/home/huangjh/spatial_omics/Patho-DBiT/processingData/Human_AITL_20um/translation/RtAI_all_mRNAs.tsv")




print("sum tAI for each gene...")
codon_SDA <- all_ws / codon_demand
write.table(codon_SDA, file = "/public1/home/huangjh/spatial_omics/Patho-DBiT/processingData/Human_AITL_20um/translation/SDA_codon_by_spot.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

for (sample in colnames(anticodon)){
  # Calculate tAI for all CUs
  sample.tai <- get.tai(t(codon), codon_SDA[,sample] )
  TAIs[,sample] = sample.tai
}

write.table(TAIs, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, "/public1/home/huangjh/spatial_omics/Patho-DBiT/processingData/Human_AITL_20um/translation/SDA_mRNA_by_spot.tsv")



print("Group spot tAI by clsuter...")
# Read cluster_position.csv, skipping the first row if it's header info, and set column names
cluster_position <- read.table("/public1/home/huangjh/spatial_omics/Patho-DBiT/processing/public1/home/huangjh/spatial_omics/Patho-DBiT/processingData/Human_AITL_20um/translation/Human_AITL_20um/Spatial_Report/cluster_position.csv", sep = ",", skip = 1, header = FALSE, check.names = FALSE)
colnames(cluster_position) <- c("sample", "group")

# Ensure sample names are characters to avoid type mismatches
cluster_position$sample <- as.character(cluster_position$sample)
cluster_position$group <- as.character(cluster_position$group)

# Get unique groups
unique_groups <- unique(cluster_position$group)

# Create a new data frame for averages by group
TAIs_by_group <- data.frame(matrix(nrow = nrow(TAIs), ncol = length(unique_groups)), row.names = rownames(TAIs))
colnames(TAIs_by_group) <- unique_groups

# Calculate the mean for each group
for (g in unique_groups) {
  group_samples <- cluster_position$sample[cluster_position$group == g]
  # Filter to samples present in TAIs columns
  group_samples <- group_samples[group_samples %in% colnames(TAIs)]
  if (length(group_samples) > 0) {
    TAIs_by_group[, g] <- rowMeans(TAIs[, group_samples, drop = FALSE], na.rm = TRUE)
  } else {
    message("Warning: No valid samples found for group ", g, ". Assigning NA.")
    TAIs_by_group[, g] <- NA
  }
}

# Write the averaged table
write.table(TAIs_by_group, file = "/public1/home/huangjh/spatial_omics/Patho-DBiT/processing/public1/home/huangjh/spatial_omics/Patho-DBiT/processingData/Human_AITL_20um/translation/Human_AITL_20um/translation/SDA_mRNA_by_cluster.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
