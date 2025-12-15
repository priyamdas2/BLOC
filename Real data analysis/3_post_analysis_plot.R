rm(list=ls())
setwd("U:/Sparse_Covariance_with_BLOC/Real data analysis")
library(cowplot)
library(ggtext)
library(ggplot2)
library(reshape2)
library(magick)
library(gridExtra)
library(dplyr) 

################################################################################
## Pathways and proteins #######################################################
################################################################################
Cell_cycle <- c("CDK1", "CYCLINB1", "CYCLINE2", "P27_pT157", "P27_pT198", "PCNA",
                "FOXM1")
Hormone_receptor <- c("ERALPHA", "ERALPHA_pS118", "PR", "AR")
Hormone_signaling_Breast <- c("BCL2", "INPP4B", "GATA3")
PI3K_AKT <- c("P27_pT157", "P27_pT198", "INPP4B", "AKT_pS473", "AKT_pT308", 
              "GSK3ALPHABETA_pS21S9",
              "GSK3_pS9", "PRAS40_pT246", "TUBERIN_pT1462", "PTEN")
Breast_reactive <- c("BETACATENIN", "CAVEOLIN1", "MYH11", "RAB11", "GAPDH", "RBM15")

PI3K_AKT_unique <- setdiff(
  PI3K_AKT,
  union(
    union(Breast_reactive, Cell_cycle),
    union(Hormone_receptor, Hormone_signaling_Breast)
  )
)

proteins_here <- unique(c(Breast_reactive, Cell_cycle, Hormone_receptor, Hormone_signaling_Breast, PI3K_AKT_unique))


case_id <- 5                    # 1=BRCA, 2=CESC, 3=OV, 4=UCEC, 5=UCS
case_names <- c("BRCA","CESC","OV","UCEC","UCS")
tol <- 1e-6                        # threshold for "Zero" tiles

# --- read correlation matrix saved from MATLAB ---
filename <- sprintf("corr_%s.csv", case_names[case_id])
C_est <- as.matrix(read.csv(filename, header = FALSE))
stopifnot(nrow(C_est) == 27, ncol(C_est) == 27)
rownames(C_est) <- colnames(C_est) <- proteins_here
################################################################################

################################################################################

# Read matrix
cor_mat <- as.matrix(C_est)
cor_melt <- melt(cor_mat)

# Label colors
label_colors <- setNames(rep("black", length(proteins_here)), proteins_here)
label_colors[Breast_reactive] <- "forestgreen"
label_colors[Cell_cycle] <- "firebrick"
label_colors[Hormone_receptor] <- "royalblue"
label_colors[Hormone_signaling_Breast] <- "orange3"
label_colors[PI3K_AKT_unique] <- "magenta"

colored_labels <- sapply(names(label_colors), function(prot) {
  sprintf("<span style='color:%s'>%s</span>", label_colors[prot], prot)
})

# Heatmap with built-in correlation color scale
heatmap_plot <- ggplot(cor_melt, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1),
                       name = "Correlation") +
  coord_fixed() +
  theme_minimal() +
  labs(title = case_names[case_id], x = "", y = "") +
  # labs(title = "Correlation Heatmap", x = "", y = "") +
  scale_x_discrete(labels = colored_labels) +
  scale_y_discrete(labels = colored_labels) +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1),
    axis.text.y = element_markdown(),
    legend.position = "right",
    legend.justification = c(0, 1),
    legend.box.just = "left",
    legend.margin = margin(l = -0),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_colorbar(barheight = 8))


# Extract the legend from the heatmap
legend_only <- get_legend(heatmap_plot)

# Create custom pathway legend as vertical boxes + text
legend_df <- data.frame(
  label = c("Breast Reactive", "Cell Cycle", "Hormone Receptor", "Hormone Signaling", "PI3K/AKT"),
  color = c("forestgreen", "firebrick", "royalblue", "orange3", "magenta")
)

pathway_legend <- ggplot(legend_df, aes(x = 1, y = reorder(label, desc(label)), fill = color)) +
  geom_tile(width = 0.2, height = 0.5, show.legend = FALSE) +
  geom_text(aes(label = label), hjust = c(-0.25,-0.38,-0.2,-0.2, -0.4), size = 5) +
  xlim(0.8, 2.5) +
  scale_fill_identity() +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  coord_fixed(ratio = 0.2)

# Combine colorbar + pathway legend vertically
right_stack <- plot_grid(legend_only, pathway_legend +
                           theme(plot.margin = margin(0, 0, 0, -30)),  # <--- shifted closer
                         ncol = 1, rel_heights = c(5, 5))

# Final layout: heatmap left, legends right (tighter)
final_plot <- plot_grid(heatmap_plot + theme(legend.position = "none"),
                        right_stack,
                        ncol = 2,
                        rel_widths = c(1, 0.42),  # <--- narrower gap
                        align = "h")
# Save and show
print(final_plot)
ggsave(
  filename = paste0("plot_corr_", case_names[case_id], ".jpg"),
  plot = final_plot,
  width = 12, height = 8, dpi = 300, units = "in"
)