setwd("~/repos/Test_mappoly/potato/")
require(mappoly)
library(ggplot2)
require(tidyverse)
require(reshape2)
source("utils.R")
ploidy <- 4

## Formatting pedigree file
pedigree.in <- read.csv("~/repos/Test_mappoly/potato/TableS3_ped.csv")
pedigree <- pedigree.in[,c(1,3,4)]
pedigree[1:3, 2:3][]<-NA
colnames(pedigree) <- c("Name", "Parent1", "Parent2")

#### Loop over all chromosomes
DAT <- MAPs <- vector("list", 12)
for(ch in 1:12){
  #### Reading multi-parental data ####
  input.data <- read.csv("~/repos/Test_mappoly/potato/TableS2_dose.csv")
  ## Selecting chromosome and converting data set to MAPpoly format
  input.data.temp <- input.data[input.data$Chrom == ch, ]
  colnames(input.data.temp) <- stringr::str_replace_all(string = colnames(input.data.temp), pattern = "\\.", replacement = "-")
  colnames(input.data.temp)[1:3] <- c("marker", "chromosome", "position")
  dat <- full_pop_to_full_sib(input.data.temp, pedigree, ploidy)
  ## Full-sib maps ##
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n~~~~~~~~~~~~~~~~~~ ", ch, " ~~~~~~~~~~~~~~~~~~~")
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  input.data = dat[[1]]
  map1 <- mapping_poptato(input.data)
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n~~~~~~~~~~~~~~~~~~ ", ch, " ~~~~~~~~~~~~~~~~~~~")
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  input.data = dat[[2]]
  map2 <- mapping_poptato(input.data)
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n~~~~~~~~~~~~~~~~~~ ", ch, " ~~~~~~~~~~~~~~~~~~~")
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  input.data = dat[[3]]
  map3 <- mapping_poptato(input.data)
  DAT[[ch]] <- dat
  MAPs[[ch]] <- list(Pop1 = map1, Pop2 = map2, Pop3 = map3)
}
##save.image(file = "~/repos/Test_mappoly/potato/poptato_maps.rda")

load("~/repos/Test_mappoly/potato/poptato_maps.rda")
#### Time in minutes (considering two-points) ####
all.times <- round(get_time(MAPs), 2)
formattable::formattable(as.data.frame(all.times))
#### ~7.4 hours ####
sum(all.times)/60

#### Get map lengths ####
all.map.lengths <- get_length(MAPs)
apply(all.map.lengths, 2, sum)

#### Map list for each population ####
l1<-lapply(MAPs, function(x) x$Pop1$final.map)
l2<-lapply(MAPs, function(x) x$Pop2$final.map)
l3<-lapply(MAPs, function(x) x$Pop3$final.map)

#### Plot map list ####
plot_map_list(l1, col = "ggstyle")
plot_map_list(l2, col = "ggstyle")
plot_map_list(l3, col = "ggstyle")

#### Genome vs. map ####
plot_genome_vs_map(l1, same.ch.lg = TRUE)
plot_genome_vs_map(l2, same.ch.lg = TRUE)
plot_genome_vs_map(l3, same.ch.lg = TRUE)

## Homolog probability for individual "W15268.27R", ch5
input.data <- DAT[[5]][[1]]
genoprob.ch5.pop1 <- calc_genoprob_error(MAPs[[5]][[1]]$final.map, error = 0.05)
h.porb <- calc_homoprob(genoprob.ch5.pop1)
H <- h.porb$homoprob %>% filter(individual == "W15268.27R")
I <- acast(H, map.position~homolog, value.var="probability")
image(I, col = blues9)

## Homolog probability for individual "W15268.27R", ch6
input.data <- DAT[[6]][[1]]
genoprob.ch6.pop1 <- calc_genoprob_error(MAPs[[6]][[1]]$final.map, error = 0.05)
h.porb <- calc_homoprob(genoprob.ch6.pop1)
H <- h.porb$homoprob %>% filter(individual == "W15268.27R")
I <- acast(H, map.position~homolog, value.var="probability")
image(I, col = blues9)


