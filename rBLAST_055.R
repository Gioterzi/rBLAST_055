library(ggplot2)
library(rBLAST)
library(Biostrings)
library(patchwork)
library(dplyr)
library(msa)
library(ape)
library(DECIPHER)
library(stringr)

# replace_plus_with_minus <- function(seq) {
#   return(chartr("+", "-", seq))
# }

write_fasta <- function(sequences, names, filename) {
  con <- file(filename, "w")
  for (i in seq_along(sequences)) {
    writeLines(paste0(">", names[i]), con)
    if (!is.character(sequences[[i]])) {
      sequences[[i]] <- as.character(sequences[[i]])
    }
    writeLines(sequences[[i]], con)
  }
  close(con)
}

read_fasta_df <- function (file = "") {
  fasta <- readLines(file)
  ind <- grep(">", fasta)
  s <- data.frame(ind = ind, from = ind + 1, to = c((ind - 
                                                       1)[-1], length(fasta)))
  seqs <- rep(NA, length(ind))
  for (i in 1:length(ind)) {
    seqs[i] <- paste(fasta[s$from[i]:s$to[i]], collapse = "")
  }
  tib <- tibble(label = gsub(">", "", fasta[ind]), sequence = seqs)
  return(tib)
}

  
palette <- c("A" = "#46ff2d", 
             "G" = "#ffae01", 
             "C" = "#f24641", 
             "T" = "#4294fa", 
             "K" = "#8b4816",
             "M" = "#83831f",
             "R" = "#ffff81",
             "S" = "#ff9d80",
             "Y" = "#e381f2",
             "W" = "#80fff2",
             "V" = "#fde4b8",
             "B" = "#f9c1bf",
             "H" = "#c0d9f9",
             "D" = "#c7ffba",
             "U" = "#8989fb",
             "N" = "black", 
             "-" = "white",
             "+" = "deeppink")


pal_df <- data.frame(names = names(palette), col = palette)

setwd("/Users/bross/Desktop/AIMS/Analysis/rBLAST")

genome <- readDNAStringSet("SCF055_genomic.fna")
mitochondrial <- readDNAStringSet("SCF055_mitochondrial_fa.gz")
plastid <- readDNAStringSet("SCF055_plastid.gz")
writeXStringSet(mitochondrial, "mito_055.fasta")
writeXStringSet(plastid, "plastid_055.fasta")


LSU_pro <- readDNAStringSet("LSU_proliferum.fasta")
LSU_055 <- readDNAStringSet("LSU_055.fasta")
LSU_049 <- readDNAStringSet("LSU_049.fasta")
Cob_pro <- readDNAStringSet("cob_proliferum.fasta")
Cob_055 <- readDNAStringSet("cob_055.fasta")
Cox_pro <- readDNAStringSet("cox_proliferum.fasta")
Cox_055 <- readDNAStringSet("cox_055.fasta")
cp23s_pro <- readDNAStringSet("cp23s_proliferum.fasta")
cp23s_055 <- readDNAStringSet("cp23s_055.fasta")

makeblastdb("SCF055_genomic.fna", dbtype = "nucl")
makeblastdb("mito_055.fasta", dbtype = "nucl")
makeblastdb("plastid_055.fasta", dbtype = "nucl")



bl <- blast(db="SCF055_genomic.fna", type = "blastn")
bl2 <- blast(db="mito_055.fasta", type = "blastn")
bl3 <- blast(db="plastid_055.fasta", type = "blastn")


run_LSU_pro <- predict(bl, LSU_pro)
run_LSU_055 <- predict(bl, LSU_055)

run_cob_pro <- predict(bl2, Cob_pro)
run_cob_055 <- predict(bl2, Cob_055)

run_cox_pro <- predict(bl2, Cox_pro)
run_cox_055 <- predict(bl2, Cox_055)

run_cp23s_pro <- predict(bl3, cp23s_pro)
run_cp23s_055 <- predict(bl3, cp23s_055)


df <- rbind(run_LSU_pro, run_LSU_055, run_cob_pro, run_cob_055, run_cox_pro, run_cox_055) %>% 
  subset(mismatch <= 4) %>%  #exclude the cases in which the base missmatches are more than 4 
  filter(
    !(grepl("LSU", qseqid) & length < 540) &
      !(grepl("Cob", qseqid) & length < 850) &
      !(grepl("Cox1", qseqid) & length < 930)
  )    %>% #exclude partial sequence matches
  subset(send >= sstart) #exclude 3'-5' sequences

extract_sequences <- function(genome_type, df) {
  sequences <- DNAStringSet()  # Initialize an empty DNAStringSet
  
  for (i in seq_len(nrow(df))) {
    qseqid <- df$qseqid[i]
    sseqid <- df$sseqid[i]
    sstart <- df$sstart[i]
    send <- df$send[i]
    
    if (grepl("LSU", qseqid)) {
      genome <- genome_type$genomic
    } else if (grepl("Cob", qseqid) || grepl("Cox", qseqid)) {
      genome <- genome_type$mitochondrial
    } else if (grepl("cp23s", qseqid)) {
      genome <- genome_type$plastid
    } else {
      warning(paste("Unrecognized qseqid:", qseqid))
      next
    }
    
    # Find matching contig in the selected genome
    matching_contigs <- grep(sseqid, names(genome), value = TRUE, ignore.case = TRUE)
    
    if (length(matching_contigs) == 0) {
      warning(paste("No contig matching", sseqid, "found in genome. Skipping..."))
      next
    }
    
    # Take the first matching contig (assuming there's only one match)
    contig_name <- matching_contigs[1]
    
    # Get the contig sequence from the genome
    contig <- genome[[contig_name]]
    
    # Extract the sequence from the contig
    seq <- subseq(contig, start=sstart, end=send)
    
    # Append the sequence to the DNAStringSet
    sequences <- append(sequences, DNAStringSet(seq))
  }
  
  # Set names for the sequences
  names(sequences) <- paste(df$qseqid, df$sseqid, sep="_")
  
  return(sequences)
}

genome_type <- list(
  genomic = genome,
  mitochondrial = mitochondrial,
  plastid = plastid
)

sequences <- extract_sequences(genome_type, df)

subset_LSU <- sequences[grepl("LSU", names(sequences))]
subset_Cob <- sequences[grepl("Cob", names(sequences))]
subset_Cox1 <- sequences[grepl("Cox1", names(sequences))]

aligned_LSU <- AlignSeqs(subset_LSU)
aligned_Cob <- AlignSeqs(subset_Cob)
aligned_Cox1 <- AlignSeqs(subset_Cox1)

aligned_sequences <- c(aligned_LSU, aligned_Cob, aligned_Cox1)

# Write sequences to a FASTA file
writeXStringSet(sequences, "BLAST_sequences_PROvs055.fasta")
writeXStringSet(aligned_sequences, "BLAST_sequences_PROvs055_ALIGNED.fasta")
writeXStringSet(aligned_LSU, "BLAST_aligned_LSU.fasta")
writeXStringSet(aligned_Cob, "BLAST_aligned_Cob.fasta")
writeXStringSet(aligned_Cox1, "BLAST_aligned_Cox1.fasta")


seq_df <- read_fasta_df("BLAST_sequences_PROvs055_ALIGNED.fasta")

aligned_plotting <- seq_df %>%
  mutate(sample_id = case_when(
    str_detect(label, "cp23S") ~ "cp23S",
    str_detect(label, "Cob") ~ "Cob",
    str_detect(label, "Cox1") ~ "Cox1",
    str_detect(label, "LSU") ~ "LSU"
  ))

key <- aligned_plotting %>%
  tibble::rownames_to_column(var = "id")

# Create long dataframe for ggplot
long_sequences <- str_split(aligned_plotting$sequence, "") %>%
  reshape2::melt() %>%
  group_by(L1) %>%
  mutate(x = row_number(),
         L1 = as.character(L1)) %>%
  left_join(., key, by = c("L1" = "id")) %>%
  ungroup()

# Plot alignment
p1 <- ggplot(long_sequences, aes(y = label, x = x)) +
  geom_tile(aes(fill = value), size = 1, name = "base") +
  facet_wrap(~ sample_id, nrow = 3, scales = "free") +
  geom_vline(xintercept = seq(50, max(long_sequences$x), by = 50), color = "black", linetype = "dashed", size = 0.05) +
  scale_fill_manual(values = palette) +
  theme(aspect.ratio = 0.3,
        axis.title.y = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 10, 0)), # Add space below strip text
        strip.background = element_blank(), # Remove strip background
        strip.placement = "outside",
        axis.text.x = element_text(size = 0.5),
        axis.text.y = element_text(size = 1),
        axis.ticks.width = unit(0.1, "cm"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(min(long_sequences$x), max(long_sequences$x), by = 100)) +
  xlab("Position")
p1

ggsave("BLASTn_sequences_aligned.pdf", plot = p1, device="pdf", width = 10, height = 4)


LSU_dna <- read.dna("BLAST_aligned_LSU.fasta", format = "fasta")
distance_matrix_LSU <- dist.dna(LSU_dna, model = "raw", as.matrix = TRUE) 
Cob_dna <- read.dna("BLAST_aligned_Cob.fasta", format = "fasta")
distance_matrix_Cob <- dist.dna(Cob_dna, model = "raw", as.matrix = TRUE) 
Cox1_dna <- read.dna("BLAST_aligned_Cox1.fasta", format = "fasta")
distance_matrix_Cox1 <- dist.dna(Cox1_dna, model = "raw", as.matrix = TRUE) 



nj_tree_LSU <- nj(distance_matrix_LSU)
nj_tree_Cob <- nj(distance_matrix_Cob)
nj_tree_Cox1 <- nj(distance_matrix_Cox1)

pdf("Blastn_NJ_Trees.pdf", width = 20, height = 10)

par(mfrow = c(1, 3))  # Set up a 1x3 plotting layout

plot(nj_tree_LSU, main = "LSU Tree", cex = 0.8)
plot(nj_tree_Cob, main = "Cob Tree", cex = 0.8)
plot(nj_tree_Cox1, main = "Cox1 Tree", cex = 0.8)

dev.off()


pca_LSU <- prcomp(distance_matrix_LSU)
pca_scores_LSU <- pca_LSU$x
sample_names_LSU <- rownames(distance_matrix_LSU)
# Plot PCA
p4 <- ggplot(data = as.data.frame(pca_scores_LSU), aes(x = PC1, y = PC2, label=sample_names_LSU)) +
  geom_point(size = 3, shape = 1) +
  geom_text(size = 3, hjust = -0.2, vjust = -0.2) +
  labs(title = "PCA LSU",
       x = "PC1",
       y = "PC2") +
  theme_minimal()

pca_Cob <- prcomp(distance_matrix_Cob)
pca_scores_Cob <- pca_Cob$x
sample_names_Cob <- rownames(distance_matrix_Cob)
# Plot PCA
p5 <- ggplot(data = as.data.frame(pca_scores_Cob), aes(x = PC1, y = PC2, label=sample_names_Cob)) +
  geom_point(size = 3, shape = 1) +
  geom_text(size = 3, hjust = -0.2, vjust = -0.2) +
  labs(title = "PCA Cob",
       x = "PC1",
       y = "PC2") +
  theme_minimal()

pca_Cox1 <- prcomp(distance_matrix_Cox1)
pca_scores_Cox1 <- pca_Cox1$x
sample_names_Cox1 <- rownames(distance_matrix_Cox1)
# Plot PCA
p6 <- ggplot(data = as.data.frame(pca_scores_Cox1), aes(x = PC1, y = PC2, label=sample_names_Cox1)) +
  geom_point(size = 3, shape = 1) +
  geom_text(size = 3, hjust = -0.2, vjust = -0.2) +
  labs(title = "PCA Cox1",
       x = "PC1",
       y = "PC2") +
  theme_minimal()

p7 <- p4+p5+p6

ggsave("BLASTn_PCA.pdf", plot = p7, device="pdf", width = 40, height = 20)

  