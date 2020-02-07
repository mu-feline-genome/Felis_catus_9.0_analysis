#!/usr/bin/env Rscript

library(optparse)

options(stringsAsFactors = FALSE)

option_list <- list( 
  make_option(c("-s", "--species1"), type="character", default = "NA", action = "store",
              help="tab blast out file for species 1"),
  make_option(c("-S", "--species2"), type="character", default = "NA", action = "store",
              help="tab blast out file for species 2"),
  make_option(c("-n", "--species1-name"), type="character", default = "NA", action = "store",
              help="name for species 1"),
  make_option(c("-N", "--species2-name"), type="character", default = "NA", action = "store",
              help="name for species 2")
)

opt <- parse_args(OptionParser(option_list = option_list))



# get reciprocal best hits

# sort each and remove duplicates



s1File <- opt$species1
s1_s2 <- read.table(s1File, sep = "\t",
                    col.names = c("qaccver", "saccver",
                                  "pident","length", 
                                  "mismatch", "gapopen",
                                  "qstart", "qend", 
                                  "sstart", "send", 
                                  "evalue", "bitscore")
                    )

s1_s2 <- s1_s2[order(s1_s2$bitscore, decreasing = TRUE),]
s1_s2 <- s1_s2[!duplicated(s1_s2$qaccver),]

s2File <- opt$species2
s2_s1 <- read.table(s2File, sep = "\t",
                    col.names = c("qaccver", "saccver",
                                  "pident","length", 
                                  "mismatch", "gapopen",
                                  "qstart", "qend", 
                                  "sstart", "send", 
                                  "evalue", "bitscore")
)

s2_s1 <- s2_s1[order(s2_s1$bitscore, decreasing = TRUE),]
s2_s1 <- s2_s1[!duplicated(s2_s1$qaccver),]


s1_s2 <- s1_s2[paste(s1_s2$saccver,s1_s2$qaccver) %in% paste(s2_s1$qaccver,s2_s1$saccver), ]

s2_s1 <- s2_s1[paste(s2_s1$saccver,s2_s1$qaccver) %in% paste(s1_s2$qaccver,s1_s2$saccver), ]

s1_s2 <- s1_s2[order(s1_s2$qaccver),]

write(paste(c(opt$`species1-name`,s1_s2$qaccver), 
            c(opt$`species2-name`,s1_s2$saccver), 
            sep = "\t"), 
      stdout())
