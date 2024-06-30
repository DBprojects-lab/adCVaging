library(ggplot2)
library(tidyverse)
library(ggrepel)
setwd("")
load("Ctx_F_new.Rdata")
load("Ctx_M_new.Rdata")

Ctx_F[, 3] <- rownames(Ctx_F)
Ctx_M[, 3] <- rownames(Ctx_M)

colnames(Ctx_F)[c(3, 4)] <- c("geneID", "Cluster")
Ctx_F$Cluster <- "Endo_F"

colnames(Ctx_M)[c(3, 4)] <- c("geneID", "Cluster")
Ctx_M$Cluster <- "Endo_M"

df <- rbind(Ctx_F, Ctx_M)

df$label <-
  ifelse(
    df$p_val < 0.05 & df$avg_log2FC < (-0.1),
    "Down-Significant",
    ifelse(
      df$p_val < 0.05 & df$avg_log2FC > 0.1,
      "UP-Significant",
      "Unconspicuous"
    )
  )

df1 <- df[df$p_val < 0.05 & abs(df$avg_log2FC) > 0.1, ]

result <- subset(df1, label == "UP-Significant") %>%
  group_by(Cluster) %>%
  top_n(5, avg_log2FC * (-log10(p_val)))

result2 <- subset(df1, label == "Down-Significant") %>%
  group_by(Cluster) %>%
  top_n(-5, avg_log2FC * (-log10(p_val)))

top10sig <- df1[df1$geneID %in% unique(
  c(
    result$geneID,
    result2$geneID,
    "ADM",
    "JUNB",
    "SLCO4A1",
    "SLCO4A2",
    "RAPGEF5",
    "ANGPT2",
    "NR4A1",
    "PIK3C2A",
    "NFKBIZ",
    "IRAK3",
    "CCN2",
    "PIK3CA",
    "GNAI2",
    "LITAF",
    "MYH9",
    "COL4A2",
    "PLOD2",
    "LDHA",
    "MAP2K1",
    "HIF1A",
    "ZFPM2",
    "NR3C1",
    "COL4A2",
    "MBD5",
    "PLCB1",
    "PIK3CA",
    "ADM",
    "INPP5D",
    "SOCS3",
    "HIF1A",
    "LAMA5",
    "DDIT4",
    "CP",
    "SGPP2"
  )
), ]

gene <- unique(top10sig$geneID)

ggplot(df, aes(
  x = -log10(p_val),
  y = avg_log2FC,
  color = label
)) +
  geom_point(alpha = 2 , size = 2) +
  facet_wrap( ~ Cluster, ncol = 2, scale = 'free') +
  scale_color_manual(name = NULL, values = c("blue", "grey", "red")) +
  theme_minimal() +
  theme(
    axis.title = element_text(
      size = 5,
      color = "black",
      face = "bold"
    ),
    axis.line.y = element_line(color = "black",
                               size = 1),
    panel.grid = element_blank(),
    legend.position = "",
    legend.direction = "vertical",
    legend.justification = c(1, 0),
    legend.text = element_text(size = 15)
  ) +
  geom_label_repel(
    data = top10sig,
    aes(
      x = -log10(p_val),
      y = avg_log2FC,
      label = geneID
    ),
    fill = "#FFDDAA",
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    force = 5,
    size = 4.5,
    max.overlaps = 100,
    arrow = arrow(
      length = unit(0.01, "npc")
      ,
      type = "open",
      ends = "last"
    )
  )


#
load("Ctx_F_new.Rdata")
load("Ctx_M_new.Rdata")
Ctx_F <- Ctx_F[Ctx_F$p_val < 0.05, ]
Ctx_M <- Ctx_M[Ctx_M$p_val < 0.05, ]
Ctx_F <- Ctx_F[abs(Ctx_F$avg_log2FC) > 0.1, ]
Ctx_M <- Ctx_M[abs(Ctx_M$avg_log2FC) > 0.1, ]
Ctx_F[, 3] <- rownames(Ctx_F)
Ctx_M[, 3] <- rownames(Ctx_M)
colnames(Ctx_F)[c(3, 4)] <- c("geneID", "Cluster")
Ctx_F$Cluster <- "Female"
colnames(Ctx_M)[c(3, 4)] <- c("geneID", "Cluster")
Ctx_M$Cluster <- "Male"

F_UP <- Ctx_F[Ctx_F$avg_log2FC > 0.1, ]
F_DOWN <- Ctx_F[Ctx_F$avg_log2FC < (-0.1), ]

M_UP <- Ctx_M[Ctx_M$avg_log2FC > 0.1, ]
M_DOWN <- Ctx_M[Ctx_M$avg_log2FC < (-0.1), ]
library(UpSetR)
upset(
  fromList(
    list(
      F_UP = F_UP$geneID,
      F_DOWN = F_DOWN$geneID,
      M_UP = M_UP$geneID,
      M_DOWN = M_DOWN$geneID
      
    )
  ),
  order.by = "freq",
  keep.order = F,
  mb.ratio = c(0.6, 0.4),
  text.scale = 2,
  queries = list(
    list(
      query = intersects,
      params = list("F_UP"),
      color = "#60281E",
      active = T
    ),
    list(
      query = intersects,
      params = list("F_DOWN"),
      color = "#8A9FD1",
      active = T
    ),
    list(
      query = intersects,
      params = list("M_UP"),
      color = "#D51F26",
      active = T
    ),
    list(
      query = intersects,
      params = list("M_DOWN"),
      color = "#272E6A",
      active = T
    ),
    
    list(
      query = intersects,
      params = list("F_UP", "M_UP"),
      color = "#208A42",
      active = T
    ),
    list(
      query = intersects,
      params = list("F_DOWN", "M_DOWN"),
      color = "#F47D2B",
      active = T
    ),
    list(
      query = intersects,
      params = list("F_DOWN", "M_UP"),
      color = "#89288F",
      active = T
    ),
    list(
      query = intersects,
      params = list("F_UP", "M_DOWN"),
      color = "#89288F",
      active = T
    )
  )
)
