rm(list = ls())

options(stringsAsFactors = FALSE)


ortho <- read.table("~/Documents/projects/project_files/find_orthologs/Felis_catus-Homo_sapian/Felis_catus-Homo_sapiens.orthologs.txt", header = TRUE)
ortho$Felis_catus <- sapply(strsplit(x = ortho$Felis_catus,"\\."),FUN = "[[", 1)
ortho$Homo_sapiens <- sapply(strsplit(x = ortho$Homo_sapiens,"\\."),FUN = "[[", 1)


humanFile = "/home/buckley/Documents/projects/project_files/vcf_impact_analysis_V9/gene_names/Homo_sapiens.GRCh38.98.ena.tsv"
humanNames <- read.table(humanFile, sep="\t", header = TRUE)
humanNames <- unique(humanNames[,c("protein_stable_id","transcript_stable_id","gene_stable_id")])
humanNames <- humanNames[humanNames$protein_stable_id != "",]
rownames(humanNames) <- humanNames$protein_stable_id

catFile = "/home/buckley/Documents/projects/project_files/vcf_impact_analysis_V9/gene_names/Felis_catus.Felis_catus_9.0.98.ena.tsv"
catNames <- read.table(catFile, sep="\t", header = TRUE)
catNames <- unique(catNames[,c("protein_stable_id","transcript_stable_id","gene_stable_id")])
catNames <- catNames[catNames$protein_stable_id != "",]
all_cat_genes <- catNames
rownames(catNames) <- catNames$protein_stable_id


ortho <- ortho[ortho$Homo_sapiens %in% humanNames$protein_stable_id,]
ortho <- ortho[ortho$Felis_catus %in% catNames$protein_stable_id,]

humanNames <- humanNames[humanNames$protein_stable_id %in% ortho$Homo_sapiens,]
catNames <- catNames[catNames$protein_stable_id %in% ortho$Felis_catus,]

colnames(humanNames) <- paste("human", colnames(humanNames), sep = "_")
colnames(catNames) <- paste("cat", colnames(catNames), sep = "_")


orthoName <- data.frame(ortho, catNames[ortho$Felis_catus,], humanNames[ortho$Homo_sapiens,])

pLI <- read.table("~/Documents/projects/project_files/vcf_impact_analysis_V9/GNOMAD/gnomad.v2.1.1.lof_metrics.by_transcript.txt", header = TRUE, fill = TRUE, sep = "\t")
rownames(pLI) <- pLI$transcript

orthoName$pLI <- pLI[orthoName$human_transcript_stable_id,"pLI"]


gtf <- read.table("/home/buckley/Documents/projects/project_files/vcf_impact_analysis_V9/gtf/Felis_catus.Felis_catus_9.0.98.gtf", sep = "\t")
gtf.cds <- gtf[gtf$V3 == "CDS",]
meta_data <- strsplit(gtf.cds$V9, "; ")
gene.id <- gsub("gene_id ", "",sapply(meta_data,"[[", 1))
tx.id <- gsub("transcript_id ", "",sapply(meta_data,"[[", 3))

df.cds <- data.frame(chr = paste("chr",gtf.cds[,1], sep = ""),
                     start = gtf.cds[,4],
                     end = gtf.cds[,5],
                     strand = gtf.cds[,7],
                     tx_id = tx.id,
                     gene_id = gene.id
                     )
                     
rownames(orthoName) <- orthoName$cat_transcript_stable_id

df.cds <- data.frame(df.cds, orthoName[df.cds$tx_id, c("cat_transcript_stable_id","cat_gene_stable_id", "human_transcript_stable_id", "human_gene_stable_id", "pLI")])

length(unique(df.cds$gene_id))
length(unique(df.cds$gene_id[!is.na(df.cds$human_transcript_stable_id)]))
length(unique(df.cds$gene_id[!is.na(df.cds$pLI)]))
length(unique(df.cds$gene_id[!is.na(df.cds$pLI) & df.cds$pLI > .9]))

all(df.cds$tx_id[!is.na(df.cds$human_transcript_stable_id)] == df.cds$cat_transcript_stable_id[!is.na(df.cds$human_transcript_stable_id)])
all(df.cds$gene_id[!is.na(df.cds$human_transcript_stable_id)] == df.cds$cat_gene_stable_id[!is.na(df.cds$human_transcript_stable_id)])

df.cds <- df.cds[,c("chr","start","end","strand","tx_id","gene_id","human_transcript_stable_id","human_gene_stable_id","pLI")]
rownames(df.cds) <- NULL

write.table(df.cds, file = "~/Documents/projects/project_files/vcf_impact_analysis_V9/human_tx_ortho_pLI.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

