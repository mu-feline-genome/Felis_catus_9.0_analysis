rm(list = ls())
options(stringsAsFactors = FALSE)

library(vcfR)
library(GenomicRanges)

# functions
convertVcfGt <- function(vcf, ann){
  gt <- extract.gt(vcf,convertNA = FALSE, return.alleles = TRUE)
  #rownames(gt) <- gsub(x = rownames(gt), pattern = "_", replacement = ":")
  ann$gt_rows <- gsub(":", "_", ann$Location)
  
  print("vcf variants not annotated:")
  print(sum(!(rownames(gt) %in% ann$gt_rows)))
  print("")
  print("annotated variants not in vcf:")
  print(sum(!(ann$gt_rows %in% rownames(gt))))
  
  ann <- ann[ann$gt_rows %in% rownames(gt),]
  
  gt.loc <- gt[ann$gt_rows,]
  gt.locA <- gsub("/.", "", gt.loc)
  gt.locB <- gsub("./", "", gt.loc)
  
  ref <- getREF(vcf)
  names(ref) <- rownames(gt)
  ref.loc <- ref[ann$gt_rows]
  
  alt.loc <- ann$Allele
  
  RR <- gt.locA == ref.loc & gt.locB == ref.loc
  R_ <- (gt.locA == ref.loc | gt.locB == ref.loc) & !RR
  
  AA <- gt.locA == alt.loc & gt.locB == alt.loc
  A_ <- (gt.locA == alt.loc | gt.locB == alt.loc) & !AA
  
  NN <- gt.locA == "." & gt.locB == "."
  N_ <- (gt.locA == "." | gt.locB == ".") & !NN
  
  DD <- gt.locA == "*" & gt.locB == "*"
  D_ <- (gt.locA == "*" | gt.locB == "*") & !DD
  
  gt.list <- list(RR = RR, R_ = R_, AA = AA, A_ = A_, NN = NN, N_ = N_, DD = DD, D_ = D_)
  
  gt.df <- data.frame(
    matrix( nrow = nrow(gt.loc), ncol = length(gt.list),
            dimnames = list(NULL, names(gt.list)))
  )
  
  for(i in names(gt.list)){
    gt.samples = col(gt.loc) * gt.list[[i]]
    gt.samples[gt.samples == 0] <- NA
    gt.samples.list <- apply(gt.samples, 1, FUN = function(x){ names(x)[!is.na(x)] })
    if(length(gt.samples.list) == 0){
      gt.samples.list = list()
      length(gt.samples.list) <- nrow(gt.df)
    }
    names(gt.samples.list) <- NULL
    gt.df[[i]] <- gt.samples.list
  }
  return(cbind(ann,gt.df))
}



colCollapse <- function(df, cols = c("Location", "Gene", "IMPACT")){
  
  link <- apply((df[cols]), FUN = paste, collapse = " ", MARGIN = 1)
  
  df0 <- data.frame(matrix(nrow = length(unique(link)),
                           ncol = ncol(df),
                           dimnames = list(NULL, colnames(df)))
  )
  
  df.long <- NULL
  for( i in colnames(df)){
    x <- df[,i]
    row <- rep(1:length(x), sapply(x, length))
    sample <- unlist(x)
    gt <- i
    if(length(row) == 0){
      row = NA
      sample = NA
    }
    
    df.long0 <- data.frame(row = row, sample = sample, gt = gt)
    
    if(length(unique(df.long0$row)) != length(x)){
      df.miss <- data.frame(row = (1:length(x))[!((1:length(x)) %in% row)],
                            sample = NA, 
                            gt = gt)
      df.long0 <- rbind(df.long0, df.miss)
    }
    
    df.long0$link <- factor(link[df.long0$row], levels = unique(link[sort(df.long0$row)]))
    
    agg <- aggregate(df.long0$sample, list(link = df.long0$link), unique, drop = FALSE)
    agg <- agg[order(agg$link),]
    
    df0[[i]] <- agg$x
    print(i)
  }
  return(df0)
}


calculateAF <- function(df, n){
  
  df0 <- df[,c("RR","R_","AA","A_","NN","N_","DD","D_")]
  
  ref <- (sapply(df0$RR, length) - is.na(df0$RR)) * 2 + (sapply(df0$R_, length) - is.na(df0$R_))
  missing <- (sapply(df0$NN, length) - is.na(df0$NN)) * 2 + (sapply(df0$N_, length) - is.na(df0$N_))
  deleted <- (sapply(df0$DD, length) - is.na(df0$DD)) * 2 + (sapply(df0$D_, length) - is.na(df0$D_))
  alt <- n * 2 - ref - missing - deleted
  
  
  df.af <- data.frame(ref, alt, missing, deleted)
  
  df.af$alleles <- df.af$ref + df.af$alt
  df.af$AF <- df.af$alt/df.af$alleles
  df.af$AF[df.af$missing + df.af$deleted == 2 * n] <- 0
  df.af$MAF <- df.af$alt/(df.af$alleles + 1)
  df.af$MAF[df.af$AF > .5 & !is.na(df.af$AF)] <- 1 - df.af$MAF[df.af$AF > .5  & !is.na(df.af$AF)]
  df.af$singleton <- NA
  df.af$singleton[df.af$alt == 1] <- sapply(df0$A_[df.af$alt == 1], "[[", 1)
  df.af$doubleton <- NA
  df.af$doubleton[df.af$alt == 2] <- sapply(df0$AA[df.af$alt == 2], "[[", 1)
  
  return(df.af)
}


indiv_variants <- function(df, search_cat){
  het = sapply(lapply(df0$A_, FUN = match, table=search_cat, nomatch = 0), sum)
  hom = sapply(lapply(df0$AA, FUN = match, table=search_cat, nomatch = 0), sum)
  
  ref = sapply(lapply(df0$R_, FUN = match, table=search_cat, nomatch = 0), sum)
  mis = sapply(lapply(df0$N_, FUN = match, table=search_cat, nomatch = 0), sum)
  del = sapply(lapply(df0$D_, FUN = match, table=search_cat, nomatch = 0), sum)
  
  indivStatus = data.frame(variant = hom + het,
                           hom, 
                           het = het * ref, 
                           het_multi = het * !ref * !(mis + del > 0), 
                           het_absent = het * (mis + del > 0))
  return(indivStatus)
}

list.contains <- function(list, x){
  sapply(lapply(x, FUN =  match, table = x), sum, na.rm = TRUE) > 0
}

list.greater <- function(list, x){
  sapply(lapply(list, FUN = '>', y=x), sum, na.rm = TRUE) > 0
}

list.less <- function(list, x){
  sapply(lapply(list, FUN = '<', y=x), sum, na.rm = TRUE) > 0
}

list.order <- function(list){
  order(sapply(list, max, na.rm = TRUE))
}

zygosity <- function(df){
  data.frame(
    hom.ref = sapply(df$RR, function(x){length(x[!is.na(x)])}),
    het = sapply(df$A_, function(x){length(x[!is.na(x)])}),
    hom.alt = sapply(df$AA, function(x){length(x[!is.na(x)])})
  )
}

ann2gr <- function(df){
  GRanges(seqnames = sapply(strsplit(df$Location, ":"), "[[", 1),
        ranges = IRanges(start = as.numeric(sapply(strsplit(df$Location, ":"), "[[", 2)),
                         width = 1))
}

varify_var <- function(df, varify_gr){
  ol <- findOverlaps(ann2gr(df), varify_gr, maxgap = 2)
  return(df[unique(queryHits(ol)),])
}

view_data <- function(df){
  df <- df[,!(colnames(df) %in% c("AA","A_", "RR", "R_", "NN", "N_", "DD", "D_"))]
  for(i in 1:ncol(df)){
      df[,i] <- sapply(df[[i]], FUN = paste, collapse = ";") 
  }
  return(df)
}


results_dir <- "/home/buckley/Documents/projects/project_files/vcf_impact_analysis_V9/Results"

# get annotation
ann <- read.table("/mnt/raid/projects/v9paper/results/annotation/v9.ensembl.slim.cds.tab", header = FALSE, 
                  col.names = c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
                                "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
                                "Codons", "Existing_variation", "IMPACT", "DISTANCE", "STRAND", "FLAGS",
                                "SYMBOL","SYMBOL_SOURCE", "HGNC_ID")
)
# get ortholog table
ortho <- read.table("~/Documents/projects/project_files/vcf_impact_analysis_V9/human_tx_ortho_pLI.tsv", sep= "\t", header = TRUE)
# merege annotation and orthologs
ann <- merge(ann, unique(ortho[,c("tx_id","human_transcript_stable_id","human_gene_stable_id","pLI")]), by.x = "Feature", by.y = "tx_id", all.x = TRUE, sort = FALSE)

# get vcf
vcf <- read.vcfR("/mnt/raid/projects/v9paper/data/vcfs/snps_only/filtered_samples/v9.ens.impact.recode.vcf")

# get refSeq genes
refSeq <- read.table("~/Documents/projects/reference_files/annotation/Felis_catus_9.0/NCBI/Felis_catus_9.0.gff3", 
                     sep = "\t", comment.char = "#", skip = 9, quote = "", fill = TRUE)
refSeq.CDS <- refSeq[refSeq$V3 == "CDS",]
refSeq.info <- strsplit(refSeq.CDS$V9, ";")
refSeq.info.len <- sapply(refSeq.info, length)
refSeq.info.df <- data.frame(row.no = rep(1:length(refSeq.info), refSeq.info.len),
                             type = sapply(strsplit(unlist(refSeq.info), "="), "[[", 1),
                             value = sapply(strsplit(unlist(refSeq.info), "="), "[[", 2)
)
refSeq.gene <- refSeq.info.df[refSeq.info.df$type == "gene",]
refSeq.protein_id <- refSeq.info.df[refSeq.info.df$type == "protein_id",]

refSeq.CDS.gr <- GRanges(seqnames = refSeq.CDS$V1,
                          ranges = IRanges(start = refSeq.CDS$V4, 
                                           end = refSeq.CDS$V5))
refSeq.CDS.gr$gene <- NA
refSeq.CDS.gr$gene[refSeq.gene$row.no] <- refSeq.gene$value
refSeq.CDS.gr$protein_id <- NA
refSeq.CDS.gr$protein_id[refSeq.protein_id$row.no] <- refSeq.protein_id$value




# analyze data
df.gt <- convertVcfGt(vcf, ann)

#df.combine <- cbind(ann, df.gt)
df.combine <- df.gt


dfCollapse <- colCollapse(df.combine)

df.all <- cbind(dfCollapse, 
                calculateAF(dfCollapse, n = 54), 
                zygosity(dfCollapse))


#high_doubleton <- df.all[!is.na(df.all$doubleton) & df.all$IMPACT == "HIGH",]
#high_doubleton_pLI <- high_doubleton[rev(list.order(high_doubleton$pLI)),]
#view_data(varify_var(high_doubleton_pLI, refSeq.CDS.gr))

#high_doubleton <- df.all[!is.na(df.all$doubleton) & df.all$IMPACT == "HIGH",]
#high_doubleton_pLI <- high_doubleton[rev(list.order(high_doubleton$pLI)),]
#view_data(varify_var(high_doubleton_pLI, refSeq.CDS.gr))



high_singleton <- df.all[!is.na(df.all$singleton) & df.all$IMPACT == "HIGH",]
high_singleton_pLI <- high_singleton[rev(list.order(high_singleton$pLI)),]
x = view_data(varify_var(high_singleton_pLI, refSeq.CDS.gr))
write.table(x = x, 
            file = paste(results_dir, "/singleton_variants_HIGH_refSeq.tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

x <- view_data((high_singleton_pLI))[1:16,]
x[order(x$Location),]
write.table(x = x, 
            file = paste(results_dir, "/singleton_variants_HIGH_ens_pLI0.90.tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




high_allAF <- df.all[df.all$IMPACT == "HIGH",]
high_allAF_pLI <- high_allAF[rev(list.order(high_allAF$pLI)),]
view_data(varify_var(high_allAF_pLI, refSeq.CDS.gr))

a <- view_data(high_allAF_pLI)
write.table(x = a, 
            file = paste(results_dir, "/allAF_variants_HIGH_refSeq.tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)





# create tables

MAF.low <- c(0,0.01,0.10,0)
MAF.high <- c(0.01,0.10,0.5,.5)
impact <- c("HIGH", "MODERATE", "LOW")
min_criteria <- df.all$MAF < 1 & !is.na(df.all$MAF) & !is.na(df.all$IMPACT) & df.all$IMPACT != "MODIFIER" & df.all$MAF > 0
consequences <- c("synonymous_variant",
                  "stop_retained_variant",
                  "missense_variant", 
                  "stop_gained",
                  "start_lost", 
                  "stop_lost")
MAF.range.low <- NULL
MAF.range.high <- NULL
impact.range <- NULL
consequences.range <- NULL
variant_count <- NULL

for(i in 1:length(MAF.low)){
  for(j in 1:length(impact)){
    for(k in 1:length(consequences)){
      variant_count <- c(variant_count,
                         nrow(df.all[df.all$MAF <= MAF.high[i] & df.all$MAF > MAF.low[i] & grepl(consequences[k], df.all$Consequence) & df.all$IMPACT == impact[j] & min_criteria,]))
      
      MAF.range.low <- c(MAF.range.low,MAF.low[i])
      MAF.range.high <- c(MAF.range.high,MAF.high[i])
      impact.range <- c(impact.range,impact[j])
      consequences.range <- c(consequences.range,consequences[k])

    }
    variant_count <- c(variant_count,
                       nrow(df.all[df.all$MAF <= MAF.high[i] & df.all$MAF > MAF.low[i] & df.all$IMPACT == impact[j] & min_criteria,]))
    
    MAF.range.low <- c(MAF.range.low,MAF.low[i])
    MAF.range.high <- c(MAF.range.high,MAF.high[i])
    impact.range <- c(impact.range,impact[j])
    consequences.range <- c(consequences.range,"All")
    
  }
  variant_count <- c(variant_count,
                     nrow(df.all[df.all$MAF <= MAF.high[i] & df.all$MAF > MAF.low[i] & min_criteria,]))
  
  MAF.range.low <- c(MAF.range.low,MAF.low[i])
  MAF.range.high <- c(MAF.range.high,MAF.high[i])
  impact.range <- c(impact.range,"All")
  consequences.range <- c(consequences.range,"All")
  
}

df.counts <- data.frame(
  MAF.range.low,
  MAF.range.high,
  impact.range,
  consequences.range,
  variant_count
)

write.table(df.counts,
            file = paste(results_dir,"/MAF_consequences_counts.tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE)


### what about permutation and snp_density
ens.gr <- GRanges(seqnames = ortho$chr,
                  ranges = IRanges(start = ortho$start, end = ortho$end),
                  tx_id = ortho$tx_id, gene_id = ortho$gene_id, ortholog = !is.na(ortho$human_gene_stable_id),
                  pLI = ortho$pLI)

# need to compress based on gene IDs!
# or is there some other way to work out density?
# is it basicaly just number of variants divided by coding sequence?

cds_len <- c(all_genes = sum(width(reduce(ens.gr))), 
             orthologs = sum(width(reduce(ens.gr[ens.gr$ortholog == TRUE]))),
             pLI = sum(width(reduce(ens.gr[!is.na(ens.gr$pLI)]))),
             tollerant = sum(width(reduce(ens.gr[ens.gr$pLI < 0.1 & !is.na(ens.gr$pLI)]))),
             medium = sum(width(reduce(ens.gr[ens.gr$pLI > 0.1 & ens.gr$pLI < 0.9 & !is.na(ens.gr$pLI)]))),
             constrained = sum(width(reduce(ens.gr[ens.gr$pLI > 0.9 & !is.na(ens.gr$pLI)]))))

gene_couunts <- c(all_genes = length(unique(ens.gr$gene_id)), 
                  orthologs = length(unique(ens.gr$gene_id[ens.gr$ortholog == TRUE])),
                  pLI = length(unique(ens.gr$gene_id[!is.na(ens.gr$pLI)])),
                  tollerant = length(unique(ens.gr$gene_id[ens.gr$pLI < 0.1 & !is.na(ens.gr$pLI)])),
                  medium = length(unique(ens.gr$gene_id[ens.gr$pLI > 0.1 & ens.gr$pLI < 0.9 & !is.na(ens.gr$pLI)])),
                  constrained = length(unique(ens.gr$gene_id[ens.gr$pLI > 0.9 & !is.na(ens.gr$pLI)])))


variant_gene_group <- data.frame(all_genes = rep(TRUE, nrow(df.all)),
                                 orthologs = !is.na(df.all$human_gene_stable_id),
                                 pLI = !is.na(df.all$pLI),
                                 tollerant = list.less(df.all$pLI, 0.1) & !list.greater(df.all$pLI, 0.1),
                                 medium = list.greater(df.all$pLI, 0.1) & list.less(df.all$pLI, 0.9),
                                 constrained = list.greater(df.all$pLI,0.9) & !list.less(df.all$pLI, 0.9))[min_criteria,]

singleton <- !is.na(df.all$singleton)[min_criteria]

impact <- df.all$IMPACT[min_criteria]

n = 10000
impact.matrix <- list(all_genes = replicate(n, 
                                            sample(impact, 
                                                   length(impact), 
                                                   replace = FALSE)),
                      orthologs = replicate(n, 
                                            sample(impact[variant_gene_group[,"orthologs"]], 
                                                   sum(variant_gene_group[,"orthologs"]), 
                                                   replace = FALSE))
)


# the loop starts here 

impacts <- c("LOW", "MODERATE", "HIGH")
gene_groups <- c("all_genes", "orthologs","pLI", "tollerant","medium", "constrained")
backgrounds <- c("all_genes", "orthologs")

df_density = NULL

for(b in backgrounds){
  for(g in gene_groups){
    for(i in impacts){
      
      if(b == "orthologs" & g == "all_genes"){
        next()
      }
      
      var_select <- variant_gene_group[variant_gene_group[,b],g]
      impact_matrix_select <- matrix(as.integer(impact.matrix[[b]] == i),ncol = n)
      impact_select <- impact[variant_gene_group[,b]] == i
      
      observed <- sum(as.integer(var_select) * as.integer(impact_select)) / cds_len[g] * 1000
      expected <- sum(var_select) * (sum(impact_select)/length(impact_select)) / cds_len[g] * 1000
      CI <- colSums(as.integer(var_select) * impact_matrix_select) / cds_len[g] * 1000
      
      hist(CI, breaks = 10, xlim = range(c(observed, CI)),
           main = paste(b,g,i))
      abline(v = expected)
      abline(v = mean(CI), col =2)
      abline(v = observed, col = 4)
      
      print(paste(b,g,i,mean(CI),expected,observed))
      
      df <- data.frame(background = b, gene_group = g, impact = i, 
                       CI99 = sort(abs(expected - CI))[.99 * length(CI)],
                       CI95 = sort(abs(expected - CI))[.95 * length(CI)],
                       expected = expected,
                       observed = observed,
                       observed_count = sum(as.integer(var_select) * as.integer(impact_select)),
                       expected_count = sum(var_select) * (sum(impact_select)/length(impact_select)),
                       cds_len = cds_len[g],
                       genes = gene_couunts[g])
      df_density = rbind(df_density, df)
    }
  }
}

write.table(df_density,
            file = paste(results_dir,"/SNP_density.tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE)


##### Singleton fraction
singleton.matrix <- list(all_genes = replicate(n,
                                            sample(singleton,
                                                   length(singleton),
                                                   replace = FALSE)),
                      orthologs = replicate(n,
                                            sample(singleton[variant_gene_group[,"orthologs"]],
                                                   sum(variant_gene_group[,"orthologs"]),
                                                   replace = FALSE))
)

# denominator keeps changing.

df_singleton = NULL

for(b in backgrounds){
  for(g in gene_groups){
    for(i in impacts){
      
      if(b == "orthologs" & g == "all_genes"){
        next()
      }
      
      var_select <- variant_gene_group[variant_gene_group[,b],g]
      impact_select <- impact[variant_gene_group[,b]] == i
      singleton_select <- singleton[variant_gene_group[,b]]
      
      observed_var <- sum(as.integer(var_select) * as.integer(impact_select)) 
      observed_singleton <- sum(as.integer(var_select) * as.integer(impact_select) * as.integer(singleton_select))
      observed <- observed_singleton/observed_var
      
      expected <- sum(singleton_select)/length(singleton_select)
      
      CI <- colSums((as.integer(var_select) * as.integer(impact_select) * singleton.matrix[[b]]))/observed_var
      
      hist(CI, breaks = 10, xlim = range(c(observed, CI)),
           main = paste(b,g,i))
      abline(v = expected)
      abline(v = mean(CI), col =2)
      abline(v = observed, col = 4)
      
      print(paste(b,g,i,mean(CI),expected,observed))
      
      df <- data.frame(background = b, gene_group = g, impact = i, 
                       CI99 = sort(abs(expected - CI))[.99 * length(CI)],
                       CI95 = sort(abs(expected - CI))[.95 * length(CI)],
                       expected = expected,
                       observed = observed,
                       observed_count = observed_singleton,
                       variant_count = observed_var,
                       total_singletons = sum(singleton_select),
                       total_variants = length(singleton_select),
                       genes = gene_couunts[g])
      df_singleton = rbind(df_singleton, df)
    }
  }
}
write.table(df_singleton,
            file = paste(results_dir,"/SNP_singleton_fraction.tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE)




# one more table intersecting AF distributions with gene groups.
MAF.low <- c(0,0.01,0.10,0)
MAF.high <- c(0.01,0.10,0.5,.5)
impact <- c("HIGH", "MODERATE", "LOW")



variant_gene_group <- data.frame(all_genes = rep(TRUE, nrow(df.all)),
                                 orthologs = !is.na(df.all$human_gene_stable_id),
                                 pLI = !is.na(df.all$pLI),
                                 tollerant = list.less(df.all$pLI, 0.1) & !list.greater(df.all$pLI, 0.1),
                                 medium = list.greater(df.all$pLI, 0.1) & list.less(df.all$pLI, 0.9),
                                 constrained = list.greater(df.all$pLI,0.9) & !list.less(df.all$pLI, 0.9))

gene_groups <- c("all_genes", "orthologs","pLI", "tollerant","medium", "constrained")

MAF.range.low <- NULL
MAF.range.high <- NULL
impact.range <- NULL
gene_group.range <- NULL
variant_count <- NULL

for(i in 1:length(MAF.low)){
  for(j in 1:length(gene_groups)){
    for(k in 1:length(impacts)){
      variant_count <- c(variant_count,
                         nrow(df.all[df.all$MAF <= MAF.high[i] & df.all$MAF > MAF.low[i] &
                                       variant_gene_group[,j] & df.all$IMPACT == impact[k] & 
                                       min_criteria,]))
      
      MAF.range.low <- c(MAF.range.low,MAF.low[i])
      MAF.range.high <- c(MAF.range.high,MAF.high[i])
      impact.range <- c(impact.range,impact[k])
      gene_group.range <- c(gene_group.range,gene_groups[j])
      
    }
    variant_count <- c(variant_count,
                       nrow(df.all[df.all$MAF <= MAF.high[i] & df.all$MAF > MAF.low[i] & variant_gene_group[,j] & min_criteria,]))
    
    MAF.range.low <- c(MAF.range.low,MAF.low[i])
    MAF.range.high <- c(MAF.range.high,MAF.high[i])
    impact.range <- c(impact.range,"All")
    gene_group.range <- c(gene_group.range,gene_groups[j])
    
  }
  variant_count <- c(variant_count,
                     nrow(df.all[df.all$MAF <= MAF.high[i] & df.all$MAF > MAF.low[i] & min_criteria,]))
  
  MAF.range.low <- c(MAF.range.low,MAF.low[i])
  MAF.range.high <- c(MAF.range.high,MAF.high[i])
  impact.range <- c(impact.range,"All")
  gene_group.range <- c(gene_group.range,"All")
  
}

df.counts <- data.frame(
  MAF.range.low,
  MAF.range.high,
  impact.range,
  gene_group.range,
  variant_count
)

write.table(df.counts,
            file = paste(results_dir,"/MAF_gene_groups.tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE)
samples <- names(table(unlist(df.all[,c("AA", "A_")])))


# SNVs in genes



variant_gene_group <- data.frame(all_genes = rep(TRUE, nrow(df.all)),
                                 orthologs = !is.na(df.all$human_gene_stable_id),
                                 pLI = !is.na(df.all$pLI),
                                 tollerant = list.less(df.all$pLI, 0.1) & !list.greater(df.all$pLI, 0.1),
                                 medium = list.greater(df.all$pLI, 0.1) & list.less(df.all$pLI, 0.9),
                                 constrained = list.greater(df.all$pLI,0.9) & !list.less(df.all$pLI, 0.9))


is <- c("LOW", "MODERATE", "HIGH")
gs <- c("all_genes", "orthologs","pLI", "tollerant","medium", "constrained")
zs <- c("AA", "A_")
MAF.low <- c(0,0.01,0.10,0)
MAF.high <- c(0.01,0.10,0.5,.5)


MAF.range.low <- NULL
MAF.range.high <- NULL
impact.range <- NULL
gene_group.range <- NULL
zygosity.range <- NULL
site_count <- NULL

#for(i in 1:length(MAF.low)){
#  for(j in 1:length(gene_groups)){
#    for(k in 1:length(impacts)){
for(af in 1:length(MAF.low)){
  for(g in gs){
    for(z in zs){
      for(i in is){
        site_count <- rbind(site_count,
                        table(factor(unlist(df.all[df.all$IMPACT == i & 
                                                     variant_gene_group[,g] & 
                                                     min_criteria &
                                                     df.all$MAF <= MAF.high[af] & 
                                                     df.all$MAF > MAF.low[af]
                                                   , z]), levels = samples)))
        MAF.range.low <- c(MAF.range.low, MAF.low[af])
        MAF.range.high <- c(MAF.range.high, MAF.high[af])
        impact.range <- c(impact.range,i)
        gene_group.range <- c(gene_group.range,g)
        zygosity.range <- c(zygosity.range,z)
      }
      site_count <- rbind(site_count,
                      table(factor(unlist(df.all[variant_gene_group[,g] & 
                                                   min_criteria &
                                                   df.all$MAF <= MAF.high[af] & 
                                                   df.all$MAF > MAF.low[af]
                                                 , z]), levels = samples)))
      MAF.range.low <- c(MAF.range.low, MAF.low[af])
      MAF.range.high <- c(MAF.range.high, MAF.high[af])
      impact.range <- c(impact.range,"All")
      gene_group.range <- c(gene_group.range,g)
      zygosity.range <- c(zygosity.range,z)
    }
    for(i in is){
      site_count <- rbind(site_count,
                      table(factor(unlist(df.all[df.all$IMPACT == i & 
                                                   variant_gene_group[,g] & 
                                                   min_criteria &
                                                   df.all$MAF <= MAF.high[af] & 
                                                   df.all$MAF > MAF.low[af]
                                                 , zs]), levels = samples)))
      MAF.range.low <- c(MAF.range.low, MAF.low[af])
      MAF.range.high <- c(MAF.range.high, MAF.high[af])
      impact.range <- c(impact.range,i)
      gene_group.range <- c(gene_group.range,g)
      zygosity.range <- c(zygosity.range,"Both")
    }
    site_count <- rbind(site_count,
                    table(factor(unlist(df.all[variant_gene_group[,g] & 
                                                 min_criteria &
                                                 df.all$MAF <= MAF.high[af] & 
                                                 df.all$MAF > MAF.low[af]
                                               , zs]), levels = samples)))
    MAF.range.low <- c(MAF.range.low, MAF.low[af])
    MAF.range.high <- c(MAF.range.high, MAF.high[af])
    impact.range <- c(impact.range,"All")
    gene_group.range <- c(gene_group.range,g)
    zygosity.range <- c(zygosity.range,"Both")
  }
}


df.counts <- data.frame(
  MAF.range.low,
  MAF.range.high,
  impact.range,
  gene_group.range,
  zygosity.range,
  mean_sites = rowMeans(site_count),
  sd_sites = apply(site_count,1, sd),
  site_count
)

write.table(df.counts,
            file = paste(results_dir,"/MAF_gene_groups_indv.tsv", sep = ""),
            sep = "\t", quote = FALSE, row.names = FALSE)







