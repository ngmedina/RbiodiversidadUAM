library(dplyr)
library(tidyr)
library(stringr)
library(vegan)
library(taxizedb)

setwd("path")

# Descargar base datos AFLIBER
aflclim <- read.csv("/afliber_clim.csv")

# Paper
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0006362#pone-0006362-t005

# Los nombres de las especies tienen subespecies
sps_mat <- str_split_fixed(aflclim$taxon, " ", 3)
sps_bin <- paste(sps_mat[,1], sps_mat[,2])
afliber2 <- cbind(aflclim, sps = sps_bin)
afliber2 <- cbind(afliber2, gen = sps_mat[,2])

# extraer familia
sps_unique <- unique(afliber2$sps)
ids <- name2taxid(sps_unique, out_type="summary")
taxonomy <- classification(ids$id)
# i <- 1
family <- lapply(1:length(ids$id), function(i){
  names <- taxonomy[[ids$id[i]]][,"name"]
  tf <- unlist(lapply(names, function(x) str_detect(x, "ceae")))
  taxonomy[[ids$id[i]]][tf,"name"]
})

family <- do.call(rbind, family)

bclass <- cbind(family = family[,1], ids )
bclass$gen <- str_split_fixed(bclass$name, " ", 2)[,1]
colnames(bclass)[2] <- "sps"

afliber3 <- left_join(afliber2, bclass[,-4], by = "sps")

afliber3 <- afliber3[complete.cases(afliber3),]
write.csv(afliber3, "alfclim_UAM.csv")

# generate the tree ####
# load the package

library("V.PhyloMaker2")

# input the sample species list

### make the example file
c1 <- bclass$Taxon
c2 <- bclass$gen
c3 <- bclass$family

example <- data.frame(species = c1, genus = c2, family = c3)

# generate a phylogeny for the sample species list

tree <- phylo.maker(sp.list = example, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL, scenarios = "S3")

write.tree(tree$scenario.3, "afliber.tre")

# more about trees
# https://pedrohbraga.github.io/PhyloCompMethods-in-R-workshop/PhyloCompMethodsMaterial.html
# https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.05904

# get trait data ####
# from here 
# https://doi.org/10.6084/m9.figshare.c.3843841.v1
# original paper here
# https://www.nature.com/articles/sdata2018135

BROT2_dat <- read.csv("C:/Users/NG.5027073/Dropbox (SCENIC MNCN CSIC)/práce/Nagore uam/Clases/Master Biodiversidad/MODIFICA/Analisis de la biodiversidad en R y SIG/brot/BROT2_dat.csv")

# select a set of traits
target_trait <- c("SLA", "LeafArea", "Height", "SeedMass", "LNCm", "GrowthForm")
brot_sel <- BROT2_dat %>% 
  filter(Trait %in% target_trait)

# Los nombres de las especies tienen subespecies
sps_trait <- str_split_fixed(brot_sel$Taxon, " ", 3)
sps_bin <- paste(sps_trait[,1], sps_trait[,2])
brot_sel<- cbind(brot_sel, sps = sps_bin)

# write.csv(brot_sel, "brot.csv")

