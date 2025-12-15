# List of cancers to work with: OV, UCEC, UCS
rm(list=ls())
setwd("U:/Sparse_Covariance_with_BLOC/Real data analysis")

library(ggplot2)
library(reshape2)
library(magick)
################################################################################
### Reading dataset ############################################################
################################################################################

names_vector <- read.csv("NExUS data/cancer_serial_no.csv", header = FALSE)[[1]]
target_terms <- c("BRCA", "CESC", "OV", "UCEC", "UCS")
match_indices <- which(names_vector %in% target_terms)
matched_names <- names_vector[match_indices]
matrix_list <- lapply(match_indices, function(i) {
  file_path <- file.path("NExUS data", paste0(i, ".csv"))
  as.matrix(read.csv(file_path, header = FALSE))
})
names(matrix_list) <- matched_names


BRCA <- matrix_list$BRCA
CESC <- matrix_list$CESC
OV <- matrix_list$OV
UCEC <- matrix_list$UCEC
UCS <- matrix_list$UCS


################################################################################
### Selecting variables ########################################################
################################################################################

### Full pathway information ###################################################)
name_list <- vector("list", 12)
array_names_FULL <- c("APOPTOSIS", "CELL CYCLE", "DNA DMG RSPNS", "EMT", 
                      "HORMONE RECPTR", "HORMONE SIG BRST", "PI3K/AKT",
                      "RAS/MAPK", "RTK", "TSC/mTOR", "BREAST REACTIVE", 
                      "CORE REACTIVE") 
Apoptosis <- c("BAK", "BAX", "BID", "BIM", "CASPASE7CLEAVEDD198", "BAD_pS112", 
               "BCL2", "BCLXL", "CIAP")
Cell_cycle <- c("CDK1", "CYCLINB1", "CYCLINE2", "P27_pT157", "P27_pT198", "PCNA",
                "FOXM1")
DNA_damage_response <- c("X53BP1", "ATM", "CHK1_pS345", "CHK2_pT68", "KU80", 
                         "MRE11", "P53", "RAD50", "RAD51", "XRCC1")
EMT <- c("FIBRONECTIN", "NCADHERIN", "COLLAGENVI", "CLAUDIN7", "ECADHERIN", 
         "BETACATENIN", "PAI1")
Hormone_receptor <- c("ERALPHA", "ERALPHA_pS118", "PR", "AR")
Hormone_signaling_Breast <- c("BCL2", "INPP4B", "GATA3")
PI3K_AKT <- c("P27_pT157", "P27_pT198", "INPP4B", "AKT_pS473", "AKT_pT308", 
              "GSK3ALPHABETA_pS21S9",
              "GSK3_pS9", "PRAS40_pT246", "TUBERIN_pT1462", "PTEN")
RAS_MAPK <- c("ARAF_pS299", "CJUN_pS73", "CRAF_pS338", "JNK_pT183Y185", 
              "MAPK_pT202Y204", "MEK1_pS217S221", "P38_pT180Y182", 
              "P90RSK_pT359S363", "YB1_pS102")
RTK <- c("EGFR_pY1068", "EGFR_pY1173", "HER2_pY1248", "HER3_pY1289", "SHC_pY317",
         "SRC_pY416", "SRC_pY527")
TSC_mTOR <- c("X4EBP1_pS65", "X4EBP1_pT37T46", "X4EBP1_pT70", "P70S6K_pT389", 
              "MTOR_pS2448", "S6_pS235S236", "S6_pS240S244", "RB_pS807S811")
Breast_reactive <- c("BETACATENIN", "CAVEOLIN1", "MYH11", "RAB11", "GAPDH", "RBM15")
Core_reactive <- c("CLAUDIN7", "ECADHERIN", "BETACATENIN", "CAVEOLIN1", "RBM15")

PI3K_AKT_unique <- setdiff(
  PI3K_AKT,
  union(
    union(Breast_reactive, Cell_cycle),
    union(Hormone_receptor, Hormone_signaling_Breast)
  )
)
################################################################################

proteins_here <- unique(c(Breast_reactive, Cell_cycle, Hormone_receptor, Hormone_signaling_Breast, PI3K_AKT_unique))
selected_variables <- read.csv("NExUS data/selected_variables.csv", header = FALSE)[[1]]
match_indices <- match(proteins_here,selected_variables)


ALL_samples <- list()
ALL_samples[[1]] <- BRCA[, match_indices]
ALL_samples[[2]] <- CESC[, match_indices]
ALL_samples[[3]] <- OV[, match_indices]
ALL_samples[[4]] <- UCEC[, match_indices]
ALL_samples[[5]] <- UCS[, match_indices]

C <- 5
p <- length(match_indices)
sample_sizes <- c(dim(BRCA)[1], dim(CESC)[1], dim(OV)[1], dim(UCEC)[1], dim(UCS)[1])


write.table(as.matrix(ALL_samples[[1]]), "proteomics_BRCA.csv", row.names = FALSE,
          col.names = FALSE, sep = ",")
write.table(as.matrix(ALL_samples[[2]]), "proteomics_CESC.csv", row.names = FALSE,
            col.names = FALSE, sep = ",")
write.table(as.matrix(ALL_samples[[3]]), "proteomics_OV.csv", row.names = FALSE,
            col.names = FALSE, sep = ",")
write.table(as.matrix(ALL_samples[[4]]), "proteomics_UCEC.csv", row.names = FALSE,
            col.names = FALSE, sep = ",")
write.table(as.matrix(ALL_samples[[5]]), "proteomics_UCS.csv", row.names = FALSE,
            col.names = FALSE, sep = ",")

################################################################################
### Saving protein-pathway map #################################################
################################################################################

df <- rbind(
  data.frame(protein = Breast_reactive, pathway = "BREAST_REACTIVE"),
  data.frame(protein = Cell_cycle, pathway = "CELL_CYCLE"),
  data.frame(protein = Hormone_receptor, pathway = "HORMONE_RECPTR"),
  data.frame(protein = Hormone_signaling_Breast, pathway = "HORMONE_SIG_BRST"),
  data.frame(protein = PI3K_AKT_unique, pathway = "PI3K_AKT")
)

df <- df[df$protein %in% proteins_here, ]

write.csv(df, "protein_pathways_selected.csv", row.names = FALSE)

