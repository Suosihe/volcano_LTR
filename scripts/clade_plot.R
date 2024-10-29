#!/bin/env Rscript
args = commandArgs(T)
treefile = args[1]
annofile = args[2]
outfig = args[3]
if (is.na(outfig)) {outfig = paste0(treefile, '.png')}


library(ggtree)
library(treeio)
library(ape)
library(dplyr)

# Read the tree file
tree <- read.tree(treefile)

# Read the annotation file
anno <- read.table(annofile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Create a data frame for annotations
anno_df <- data.frame(label = as.character(anno$ID), clade = anno$Clade)

# Create a function to clean the labels
clean_labels <- function(label) {
  label <- gsub("[':]", "", label)
  return(label)
}

# Clean the labels in the tree
tree$tip.label <- sapply(tree$tip.label, clean_labels)


p0 <- ggtree(tree, layout="equal_angle") %<+% anno_df + 
  #  geom_tiplab(aes(label = label)) +
  geom_tippoint(aes(color = clade), size = 1) +
  #  geom_tree(aes(color = clade)) +
  scale_color_manual(values = c(
    'Ivana' = 'firebrick', 
    'Tork' = 'lightblue', 
    'CRM' = 'green', 
    'SIRE' = 'pink', 
    'Reina' = 'orange', 
    'Ale' = 'navy',
    'Athila' = 'tomato',
    'Ikeros' = 'darkorange',
    'Tekay' = 'darkgreen',
    'TAR' = 'gold'
  ), na.translate = FALSE)  +
  theme(legend.position = "right")

ggsave(outfig, p0, width=13.5, height=8.4, dpi=300, units="in")
