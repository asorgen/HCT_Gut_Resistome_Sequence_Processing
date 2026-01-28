#Author: Alicia Sorgen
#Date: 10-12-21
#Description: Builds raw taxa count tables from MetaPhlAn2 output.

# Libraries -----
library(data.table); message("data.table:", packageVersion("data.table"))
library(stringr); message("stringr:", packageVersion("stringr"))

rm(list=ls())

# Prep -----
params <- commandArgs(trailingOnly = TRUE)
rel_abun_file <- params[1]
outputDir <- params[2]
funcScript <- params[3]
profileDir <- params[4]

# rel_abun_file <- "~/git/Heston_Reproduction/data/processed/UNC/UNC_short_tables/UNC_short_metaphlan4_counts.tsv"
# outputDir <- "~/git/Heston_Reproduction/data/processed/UNC/UNC_short_tables"
# funcScript <- "~/git/Heston_Reproduction/pipelineScripts/scripts/post_scripts/functions.R"
# profileDir <- "~/git/Heston_Reproduction/data/processed/UNC/UNC_short_tables"

source(funcScript)


# Relative abundance -----
table <- read.table(rel_abun_file, skip = 1, header = T, sep = "\t")

for (i in c("k", "p", "c", "o", "f", "g", "s", "t")) {
   if (i == 'k') {
      unclassified_delim <- paste0(i, '__UNCLASSIFIED')
   } else {
      unclassified_delim <- paste0(unclassified_delim, "|", i, '__UNCLASSIFIED')
   }
}
table$clade_name[which(table$clade_name == "UNCLASSIFIED")] <- unclassified_delim

phylum <- getphylumTableMetaPhlAn2(table)
class <- getclassTableMetaPhlAn2(table)
order <- getorderTableMetaPhlAn2(table)
family <- getfamilyTableMetaPhlAn2(table)
genus <- getgenusTableMetaPhlAn2(table)
species <- getspeciesTableMetaPhlAn2(table)
strain <- getstrainTableMetaPhlAn2(table)

taxonomy <- strain$strain
taxonomy <- gsub("[|]", "_", taxonomy)
taxonomy <- str_split(taxonomy, pattern = "_[pcofgst]__")
taxonomy=do.call(rbind,taxonomy)
taxonomy <- as.data.frame(taxonomy)

colnames(taxonomy)[colnames(taxonomy)=="V1"] <- "Kingdom"
taxonomy$Kingdom <- gsub("k__", "", taxonomy$Kingdom)
colnames(taxonomy)[colnames(taxonomy)=="V2"] <- "Phylum"
colnames(taxonomy)[colnames(taxonomy)=="V3"] <- "Class"
colnames(taxonomy)[colnames(taxonomy)=="V4"] <- "Order"
colnames(taxonomy)[colnames(taxonomy)=="V5"] <- "Family"
colnames(taxonomy)[colnames(taxonomy)=="V6"] <- "Genus"
colnames(taxonomy)[colnames(taxonomy)=="V7"] <- "Species"
colnames(taxonomy)[colnames(taxonomy)=="V8"] <- "Strain"
write.table(taxonomy, paste0(outputDir, "_MetaPhlAn4_Full_Taxonomy.tsv"),sep="\t", quote = FALSE, row.names = F)


rel_tableList <- list(phylum, class, order, family, genus, species)

for (i in rel_tableList) {
   taxaLevel <- names(i)[1]
   
   df <- parseMetaPhlAn2Taxa(i)
   df$SampleID <- sapply(strsplit(df$SampleID, "_rel_abun"), "[", 1)
   # df$SampleID <- gsub(pattern = "_rel_abun_w_read_stats", replacement = "", df$SampleID)
   df$SampleID <- gsub(pattern = "[.]", replacement = "-", df$SampleID)
   
   write.table(df, paste0(outputDir, "_", taxaLevel, "_MetaPhlAn4_rel_abun.tsv"),sep="\t", quote = FALSE, row.names = F)
}
message("Relative abundance tables COMPLETE!\n")






# Raw Counts -----
# files <- list.files(profileDir)
files <- list.files(path = profileDir, pattern = "_rel_abun_w_read_stats\\.txt$", full.names = FALSE)

index <- 1

for (file in files) {

   filePath <- file.path(profileDir, file); #message("filePath = ", filePath)
   sampleID <- sapply(strsplit(file, "_rel_abun"), `[`, 1); #message("sampleID = ",sampleID)
   table <- read.delim(filePath, skip = 5, sep = "\t", header = TRUE)
   colnames(table) <- gsub("^X\\.", "", colnames(table))
   names(table)[names(table) == "estimated_number_of_reads_from_the_clade"] <- sampleID
   table <- table[, c("clade_name", sampleID)]
   
   for (i in c("k", "p", "c", "o", "f", "g", "s", "t")) {
      if (i == 'k') {
         unclassified_delim <- paste0(i, '__UNCLASSIFIED')
      } else {
         unclassified_delim <- paste0(unclassified_delim, "|", i, '__UNCLASSIFIED')
      }
   }
   table$clade_name[which(table$clade_name == "UNCLASSIFIED")] <- unclassified_delim
   

   phylum.df <- getphylumTableMetaPhlAn2(table)
   class.df <- getclassTableMetaPhlAn2(table)
   order.df <- getorderTableMetaPhlAn2(table)
   family.df <- getfamilyTableMetaPhlAn2(table)
   genus.df <- getgenusTableMetaPhlAn2(table)
   species.df <- getspeciesTableMetaPhlAn2(table)
   # message("***** NO ERROR *****\n")

   if (index == 1) {
      phylum <- phylum.df
      class <- class.df
      order <- order.df
      family <- family.df
      genus <- genus.df
      species <- species.df
   } else {
      phylum <- merge(phylum, phylum.df, by = "phylum", all = TRUE)
      class <- merge(class, class.df, by = "class", all = TRUE)
      order <- merge(order, order.df, by = "order", all = TRUE)
      family <- merge(family, family.df, by = "family", all = TRUE)
      genus <- merge(genus, genus.df, by = "genus", all = TRUE)
      species <- merge(species, species.df, by = "species", all = TRUE)
   }

   index <- index + 1
} # for (file in files)


raw_tableList <- list(phylum, class, order, family, genus, species)

for (i in raw_tableList) {
   taxaLevel <- names(i)[1]
   
   df <- parseMetaPhlAn2Taxa(i)
   df$SampleID <- gsub(pattern = "[.]", replacement = "-", df$SampleID)
   
   write.table(df, paste0(outputDir, "_", taxaLevel, "_MetaPhlAn4_raw_abun.tsv"),sep="\t", quote = FALSE, row.names = F)
   
}
message("Raw count tables COMPLETE!\n")

