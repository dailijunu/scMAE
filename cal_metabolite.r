library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(DoubletFinder)

library(ggplot2)
library(ggsignif)

all_df <- readRDS("./rds/dm_obj_in_expT.rds")

all_cell_type <- table( all_df@meta.data$cell_types )
rownames(all_cell_type)
all_cell_type[1]

cell_type <- rownames(all_cell_type)[9]

DefaultAssay(all_df) <- "SCT"
df <- subset(all_df, cell_types == cell_type)

expr <- GetAssayData(df, assay = "SCT", layer = "data")[c(product_genes, substrate_genes), ]

fdps_val <- expr[product_genes, ]

substrate_mean <- colMeans(expr[substrate_genes, ])

fpp_balance <- adj_coef_p * log2(fdps_val) - adj_coef_s * log2(substrate_mean)

df$fpp_balance <- as.vector(fpp_balance)

library(ggplot2)
library(ggsignif)

meta <- df@meta.data[, c("group","fpp_balance")]

my_comparisons <- list(
  c("Control","Diabetic"),
  c("Control","DR"),
  c("Diabetic","DR")
)

options(repr.plot.width = 5, repr.plot.height = 5)
ggplot(meta, aes(x = group, y = fpp_balance, fill = group)) +
  geom_boxplot(outlier.size = 0.3) +
  geom_signif(
    comparisons = my_comparisons,
    annotate = "p-value",
    tip_length = 0.01,
    vjust = 0,
    y_position = y_pos 
  ) +
  ggtitle(cell_type) +
  theme_bw() +
  labs(x = "Group", y = "Prediction level of FPP") +
  theme(legend.position = "none")

expr_matrix <- AverageExpression(
  object = df,
  assays = "SCT",
  features = c(product_genes, substrate_genes),
  group.by = "group", 
  verbose = TRUE
)

expr_matrix
