library(viridis)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggpubr)


meta <- readRDS("ROSMAP.VascularCells.meta_full.rds")

library(readxl)
cellfraction <-
  read_excel("cellfraction_cellsubtype_vsALL.xlsx", sheet = "cellfraction_cellsubtype_vsALL")
cellfraction$Allcellnum <- 0

#Fig 1E : Calculate the total number of cells provided by each sample in the prefrontal cortex
projidfraction <-
  as.data.frame(table(meta$projid, meta$cellsubtype))

for (i in 1:length(cellfraction$projid)) {
  cellfraction[i, ]$Allcellnum <-
    mean(as.numeric(projidfraction[projidfraction$Var1 == cellfraction$projid[i], ]$Freq /
                      cellfraction[cellfraction$projid ==
                                     cellfraction$projid[i],
                                   projidfraction[projidfraction$Var1 ==
                                                    cellfraction$projid[i], ]$Var2]),
         na.rm = TRUE)
}

meta <- meta[meta$brain_region == "Prefrontal_cortex", ]
table(meta$msex, meta$ADdiag2types)

data <- meta[, c("projid", "msex", "ADdiag2types", "celltype")]
data$class <- 0
data[data$msex == "0" & data$ADdiag2types == "AD", "class"] <-
  "F_AD"
data[data$msex == "0" &
       data$ADdiag2types == "nonAD", "class"] <- "F_nonAD"
data[data$msex == "1" & data$ADdiag2types == "AD", "class"] <-
  "M_AD"
data[data$msex == "1" &
       data$ADdiag2types == "nonAD", "class"] <- "M_nonAD"

projidcell <- as.data.frame(table(data$projid, data$celltype))

A <-
  data.frame(prop.table(table(data$projid, data$celltype), margin = 1))
colnames(A) <- c("projid", "celltype", "Freq")
A$msex <- 3
A$ADdiag2types <- 0
A$class <- 0
identical(A$projid, projidcell$Var1)
identical(A$celltype, projidcell$Var2)
A$Freq <- projidcell$Freq

identical(subset(A, celltype == "Endo")$projid,
          subset(A, celltype == "Fib")$projid)
identical(subset(A, celltype == "Endo")$projid,
          subset(A, celltype == "Ependymal")$projid)
subset(A, celltype == "Endo")$projid

data1 <- data[!duplicated(data$projid), ]

msex <- c()
ADdiag2types <- c()
class <- c()
for (i in 1:length(subset(A, celltype == "Endo")$projid)) {
  msex <-
    c(msex, data1[data1$projid == subset(A, celltype == "Endo")$projid[i], "msex"])
  ADdiag2types <-
    c(ADdiag2types, data1[data1$projid == subset(A, celltype == "Endo")$projid[i], "ADdiag2types"])
  class <-
    c(class, data1[data1$projid == subset(A, celltype == "Endo")$projid[i], "class"])
}

A$msex <- rep(msex, 5)
A$ADdiag2types <- rep(ADdiag2types, 5)
A$class <- rep(class, 5)
A$celltype <- as.factor(A$celltype)
A$class <- as.factor(A$class)
A$msex <- as.factor(A$msex)
A$class <- as.factor(A$class)

A$Allcellnum <- 0
rownames(cellfraction) <- cellfraction$projid
cellfraction1 <- cellfraction[unique(as.character(A$projid)), ]
identical(as.character(A$projid), as.character(rep(cellfraction1$projid, 5)))
A$Allcellnum <- rep(cellfraction1$Allcellnum, 5)


meta11 <- readRDS("ROSMAP.VascularCells.meta_full.rds")
A1 <-
  prop.table(table(meta11$projid, meta11$brain_region), margin = 1)
A1 <- as.data.frame(A1)
A1 <- subset(A1, subset = Var2 == "Prefrontal_cortex")
A1 <- A1[A1$Freq != 0, ]
identical(as.character(rep(A1$Var1, 5)), as.character(A$projid))

A$brain_region_num <- A$Allcellnum * c(rep(A1$Freq, 5))

A$cellfraction_celltype_vsALL <- A$Freq / A$brain_region_num

A$celltype <-
  factor(A$celltype, levels = c("Endo", "Per", "Fib", "SMC", "Ependymal"))
ggboxplot(
  A,
  x = "celltype",
  y = "cellfraction_celltype_vsALL",
  color = "class",
  outline = FALSE,
  add = "jitter",
  short.panel.labs = FALSE
) +
  stat_compare_means(aes(group = class), label = "p.format") +
  scale_y_continuous(limits = c(0, 0.02))

#Fig 1A
meta$msex
meta$ADdiag2types
data1 <- subset(meta, subset = msex == "1" & ADdiag2types == "AD")
data2 <-
  subset(meta, subset = msex == "1" & ADdiag2types == "nonAD")
data3 <- subset(meta, subset = msex == "0" & ADdiag2types == "AD")
data4 <-
  subset(meta, subset = msex == "0" & ADdiag2types == "nonAD")

#Calculate the proportion of brain area data provided by the sample
A <-
  data.frame(
    name = c(
      "Prefrontal_cortex",
      "Hippocampus",
      "Angular_gyrus",
      "Anterior_thalamus",
      "Midtemporal_cortex",
      "Entorhinal_cortex"
    ),
    percent = c(
      length(table(meta[meta$brain_region == "Prefrontal_cortex",]$projid)) /
        length(table(meta$projid))
      ,
      length(table(meta[meta$brain_region == "Hippocampus",]$projid)) / length(table(meta$projid))
      ,
      length(table(meta[meta$brain_region == "Angular_gyrus",]$projid)) / length(table(meta$projid))
      ,
      length(table(meta[meta$brain_region == "Anterior_thalamus",]$projid)) /
        length(table(meta$projid))
      ,
      length(table(meta[meta$brain_region == "Midtemporal_cortex",]$projid)) /
        length(table(meta$projid))
      ,
      length(table(meta[meta$brain_region == "Entorhinal_cortex",]$projid)) /
        length(table(meta$projid))
    )
  )
A$percent <- round(A$percent, 3)
A$name <- factor(
  A$name,
  levels = c(
    "Entorhinal_cortex",
    "Midtemporal_cortex",
    "Anterior_thalamus",
    "Angular_gyrus",
    "Hippocampus",
    "Prefrontal_cortex"
  )
)
dfbar <-
  data.frame(
    x = c(
      "Entorhinal_cortex",
      "Midtemporal_cortex",
      "Anterior_thalamus",
      "Angular_gyrus",
      "Hippocampus",
      "Prefrontal_cortex"
    ),
    y = c(1, 1, 1, 1, 1, 1)
  )
dfbar$x <- factor(
  dfbar$x,
  levels = c(
    "Entorhinal_cortex",
    "Midtemporal_cortex",
    "Anterior_thalamus",
    "Angular_gyrus",
    "Hippocampus",
    "Prefrontal_cortex"
  )
)

ggplot(A, aes(name, percent, fill = name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "#E6D933",
    "#B19B19",
    "#2CAD3F",
    "#1EB5B8",
    "#6792CD",
    "#F96D00"
  )) +
  geom_col(
    data = dfbar,
    mapping = aes(x = x, y = y),
    fill = "#999999",
    alpha = 0.2
  ) +
  theme_test(base_size = 20) +
  geom_text(aes(name, percent, label = percent)) +
  theme(axis.text.x = element_text(
    vjust = 0.5,
    hjust = 0.5,
    angle = 45
  ))




A <-
  data.frame(
    class = c("AD", "AD", "nonAD", "nonAD"),
    num = c(95, 113, 114, 106),
    sex = c("male", "female", "male", "female"),
    AAA = c("1", "1", "1", "1")
  )
#Sample number, male and female, and normal AD number drawing
ggplot(A, aes(class, num, fill = sex)) +
  theme_test(base_size = 30) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(class, num, label = num), position = "fill") +
  scale_fill_manual(
    name = "tuli",
    labels = c("female", "male"),
    values = c("#F672A2", "#0091FF")
  )

ggplot(A, aes(AAA, num, fill = class)) +
  theme_test(base_size = 30) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(AAA, num, label = num), position = "stack") +
  scale_fill_manual(
    name = "tuli",
    labels = c("AD", "nonAD"),
    values = c("#804D9B", "#00D1CD")
  )


#
mydata <- as.data.frame(table(meta$projid, meta$celltype))
mydata <- mydata[mydata$Freq != 0,]
colnames(mydata) <- c("projid", "celltype", "Freq")

cellnum <-
  aggregate(mydata$Freq, by = list(tpye = mydata$projid), sum)
cellnum <- cellnum[order(cellnum$x),]
cellnum$tpye <- as.character(cellnum$tpye)
mydata$projid <- factor(mydata$projid, levels = cellnum$tpye)

ggplot(mydata, aes(projid, Freq, fill = celltype)) +
  labs(x = 'high', cex = 4, y = 'stack') + theme_test(base_size = 10) +
  geom_bar(stat = "identity", position = "stack")

#Sample provides a bar graph of brain regions
mydata2 <- as.data.frame(table(meta$projid, meta$brain_region))
mydata2 <- mydata2[mydata2$Freq != 0,]
colnames(mydata2) <- c("projid", "brain_region", "Freq")
mydata2$projid <- factor(mydata2$projid, levels = cellnum$tpye)

AAA <-
  ggplot(mydata2, aes(projid, brain_region, fill = brain_region)) +
  theme_test(base_size = 10) +
  geom_bar(stat = "identity", position = "fill")
AAA
table(ggplot_build(AAA)$data[[1]]$fill)

#

need <-
  meta[, c('projid',
           'ADdiag2types',
           'ADdiag3types',
           'braaksc',
           'cogdx',
           'ceradsc',
           'msex')]
need <- need[!duplicated(need$projid),]
need$projid <- as.character(need$projid)
rownames(need) <- need$projid
need <- need[cellnum$tpye,]
identical(need$projid, cellnum$tpye)
need$cellnum <- cellnum$x

data <- as.matrix(need$cellnum)
rownames(data) <- rownames(need)
colnames(data) <- c("cellnum")

annotation_col <- need[, c(-1,-8)]
annotation_col$braaksc <- as.character(annotation_col$braaksc)
annotation_col$cogdx <- as.character(annotation_col$cogdx)
annotation_col$ceradsc <- as.character(annotation_col$ceradsc)
annotation_col$msex <- as.character(annotation_col$msex)

colnames(annotation_col)

annColors <- list(
  "ADdiag2types" = c("AD" = "#804D9B", "nonAD" = "#00D1CD"),
  "ADdiag3types" = c(
    "earlyAD" = "#A86A7D",
    "lateAD" = "#FFBED5",
    "nonAD" = "#00D1CD"
  ),
  
  "braaksc" = c(
    "0" = "#F6BCC6",
    "1" = "#F1891A",
    "2" = "#EE7B51",
    "3" = "#E47C7B",
    "4" = "#EB6248",
    "5" = "#E9491A",
    "6" = "#E7211A"
  ),
  "cogdx" = c(
    "1" = "#C48ABB",
    "2" = "#C486B8",
    "3" = "#AB549C",
    "4" = "#C12787",
    "5" = "#A8315C",
    "6" = "#BA1D7C"
  ),
  
  "ceradsc" = c(
    "1" = "#F9F8CC",
    "2" = "#CEE9EE",
    "3" = "#F8C9C9",
    "4" = "#F9C898"
  ),
  "msex" = c("0" = "#F672A2", "1" = "#0091FF")
)

library(viridis)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
annCol = annotation_col
bk <- c(seq(0, 500, by = 20))
pheatmap(
  mat = t(data),
  
  border_color = NA,
  
  cellwidth = 1,
  cellheight = 30,
  color = c(colorRampPalette(colors = c("white", "white"))(length(bk))),
  cluster_rows = F,
  
  cluster_cols = F,
  
  show_rownames = F,
  
  show_colnames = F,
  
  annotation_col = annCol,
  
  annotation_colors = annColors,
  
  fontsize = 6.3
)
