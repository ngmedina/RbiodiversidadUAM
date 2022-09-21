library(dplyr)
library(tidyr)
library(stringr)
library(picante)
library(vegan)

setwd("path")

# Datos de comunidad ####
com <- read.csv("alfclim_UAM.csv")

# long to wide
com$count <- 1
com_wide <- com %>%
  select(sps, utm.cell, count) %>%
  pivot_wider(names_from = sps,
              values_from = count,
              values_fn = mean)

com_wide[is.na(com_wide)] <- 0

# la base de datos es enorme. Para simplificar cogemos solo las cuadrículas en 
# las que hay especies de Cistus
cis_utm <- com %>%
  filter(str_detect(sps, "Cistus")) %>% # seleccionar las filas que tienen plots que tienen Quercus
  distinct(utm.cell) %>% # extraer las cuadriculas que tienen al menos una especie de Cistus
  pull(utm.cell) # seleccionar en la base de datos los plots que tienen Cistus

# generamos la lista de especies
cis_sps <- com %>%
  filter(str_detect(sps, "Cistus")) %>% # seleccionar las filas que tienen plots que tienen Quercus
  distinct(sps) %>% # extraer las cuadriculas que tienen al menos una especie de Cistus
  pull(sps) # seleccionar en la base de datos los plots que tienen Cistus

# todas las especies que están en cuadrículas en las que está cistus
cis_tbl <- com %>%
  filter(utm.cell %in% cis_utm)
# lista de especies acompañantes
acomp <- cis_tbl %>% 
  distinct(sps)
acomp <- acomp[!acomp %in% cis_sps]

sp <- cis_sps[[1]]
hi_cor <- lapply(cis_sps, function(sp){
  # seleccionamos las 10 especies más correlacionadas con cada especie
  # las utm en las que está la especie de cistus
  cis_utm_tmp <- com %>%
    filter(str_detect(sps,sp)) %>% 
    distinct(utm.cell) %>% # extraer las cuadriculas que tienen al menos una especie de Cistus
    pull(utm.cell)
  
  # tabla con la UTM que tienen la especie de cistus
  cis_tbl_tmp <- com %>%
    filter(utm.cell %in% cis_utm_tmp)
    
  acomp_tmp <- cis_tbl_tmp %>% 
    distinct(sps)
  acomp_tmp <-pull(acomp_tmp[!acomp_tmp %in% sp])
  
  # sustituir los NA y otros valores raros por 0
  
  # spcol <- acomp_tmp[[2]]
  pairs <- lapply(acomp_tmp, function(spcol){
    x <- cbind(com_wide[, sp], com_wide[,spcol])
    # tryCatch(1-bipartite::C.score(x, normalize = T), 
    #          error=function(e) NA)
    tryCatch(cor(x, method = "spearman")[1,2], 
             error=function(e) NA)
  })
  names(pairs) <- acomp_tmp   
  pairs <- unlist(pairs)
  pairs_sort <- sort(pairs, decreasing = T)
  names(pairs_sort)[1:11]
})

hi_cor <- unlist(hi_cor)

# creamos la base de datos recortada
com_wide_sel <- com_wide %>%
  select(utm.cell, all_of(hi_cor))

com_wide_sel <- com_wide_sel %>%
  filter(rowSums(com_wide_sel[,-1]) !=0)

# guardamos la base de datos
write.csv(com_wide_sel, "alfiber_com_sel.csv" )

# Datos ambientales ####
# creamos la tabla de variables ambientales

utm_sel <- com_wide_sel %>%
  distinct(utm.cell) %>%
  pull(utm.cell)

clim_sel <- com %>%
  select(utm.cell, x, y, bio19:bio01) %>%
  filter(utm.cell %in% utm_sel) %>%
  distinct()

# guardar la matriz de características ambientales
write.csv(clim_sel, "afiber_clim_sel.csv")

# Datos de rasgos ####

# cargar datos de rasgos 
trait <- read.csv("brot.csv")

traitGF <- trait %>% filter(Trait == "GrowthForm")
traitCuant <- trait%>% filter(Trait != "GrowthForm")
traitCuant$Data <- as.numeric(traitCuant$Data)

traitCuant <- traitCuant %>% select(sps, Trait, Data) %>%
  pivot_wider(names_from = Trait,
              values_from = Data,
              values_fn = mean)

traitGF <- traitGF %>% select(sps, Trait, Data) %>%
  pivot_wider(names_from = Trait,
             values_from = Data,
             values_fn = first)

trait_wide <- as.data.frame(left_join(traitGF, traitCuant, by = "sps"))
rownames(trait_wide) <- trait_wide$sps

# Datos filogenéticos ####

# prune and check integrity
tree <- read.tree("afliber.tre")
tree$tip.label <- str_replace(tree$tip.label, "_", " ")

# Match data in phylogeny, community and trait #####
matchPC <- match.phylo.comm(tree, com_wide_sel)
matchPT <- match.phylo.data(matchPC$phy, trait_wide)
match <- match.phylo.comm(matchPT$phy, com_wide_sel)

comm <- data.frame(utm.cell = com_wide_sel$utm.cell, match$comm)
phy <- match$phy
trait <- matchPT$data

# guardar 
write.csv(comm, "match_comm.csv", row.names = F)
write.csv(trait, "match_trait.csv", row.names = F)
write.tree(phy, "match_phylo.tre")
