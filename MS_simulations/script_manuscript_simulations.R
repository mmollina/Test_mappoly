install.packages("mappoly")
## Loading MAPpoly
require(mappoly)
setwd("~/repos/Test_mappoly/MS_simulations/")
## Declaring arguments
mrk.tail <- c(60, 60, 60, 100, 100, 150, 200, 150)
n.ph <- c(20, 20, 40, 60, 60, 100, 300, 100)
## List to store results
RES <- vector("list", 8)
## Naming according to number of individuals
names(RES) <- names(mrk.tail) <- names(n.ph) <- c('200', '100', '50', '30', '20', '15', '10', '15.DR')
for(i in names(RES)) ## Loop over all sample sizes
{
  reps <- vector("list", 3)
  for(j in 1:3){ ## Over 3 reps
    if(i != '15.DR'){ ## if different from DR0.5 data set
      ## Read data
      dat <- read_geno_csv(file.in = paste0("~/repos/Test_mappoly/MS_simulations/data/F1_noff",
                                            i,"_DR0_rep",
                                            j,"_polyorigin_geno_snparray_mappoly.csv"),
                           ploidy = 4)
    } else {  ## if DR0.5 data set
      ## Read data
      dat <- read_geno_csv(file.in = paste0("~/repos/Test_mappoly/MS_simulations/data/F1_noff15_DR0.5_rep",
                                            j,"_polyorigin_geno_snparray_mappoly.csv"),
                           ploidy = 4)
    }
    ## make sequence
    s <- make_seq_mappoly(dat, "all")
    ## two-point computation
    time.tpt <- system.time(tpt <- est_pairwise_rf(s))
    ## map/phasing estimation for a given order
    time.map <- system.time(map <- est_rf_hmm_sequential(input.seq = s,
                                                         start.set = 4,
                                                         thres.twopt = 5,
                                                         thres.hmm = 5,
                                                         extend.tail = mrk.tail[i],
                                                         twopt = tpt,
                                                         sub.map.size.diff.limit = 10,
                                                         phase.number.limit = n.ph[i]))
    ## Re-estimating recombination fractions with error
    time.map.err <- system.time(map.err <- est_full_hmm_with_global_error(map, error = 0.05))
    ## Gathering results
    reps[[j]] <- list(time.tpt = time.tpt[3],
                      tpt = tpt,
                      time.map = time.map[3],
                      map = map,
                      time.map.err = time.map.err[3],
                      map.err = map.err)
    ## Temporary save
    save(reps, RES, file = "~/repos/Test_mappoly/out.rda")
  }
  ## Saving in the final list
  RES[[i]] <- reps
}
## Saving results
save(reps, RES, file = "~/repos/Test_mappoly/MS_simulations/mappoly_manuscript_simulations.rda")

## Results
require(mappoly)
require(ggplot)
require(tidyverse)
load(file = "~/repos/Test_mappoly/MS_simulations/mappoly_manuscript_simulations.rda")
source("~/repos/Test_mappoly/utils.R")
## Re-arranging results
a3 <- a2 <- a1 <- vector("list", length(RES))
names(a3) <- names(a2) <- names(a1) <- names(RES)
for(i in names(RES)){
  a1[[i]] <- RES[[i]][[1]]
  a2[[i]] <- RES[[i]][[2]]
  a3[[i]] <- RES[[i]][[3]]
}
RES2 <- list(a1, a2, a3)
Z <- parse.results(RES2)
fact.ord <- c("10", "15.DR", "15",  "20",  "30", "50", "100", "200")
## Time (s) distribution for all sample sizes
a1<-Z %>%
  mutate(n.ind = fct_relevel(n.ind, fact.ord)) %>%
  group_by(n.ind) %>%
  summarise(mean = round(mean(time.sec),2), sd = round(sd(time.sec),2), n = n())
colnames(a1) <- c("N.ind.", "Mean (min)", "Std. Dev", "N.Rep")
formattable(a1)
Z %>%
  mutate(n.ind = fct_relevel(n.ind, fact.ord)) %>%
  ggplot( aes(x=n.ind, y=time.sec, fill=n.ind)) +
  geom_boxplot() + ylab("Time (min)") +
  xlab("Number of individuals")

## Distribution of phased markers (%)
a2<-Z %>%
  mutate(n.ind = fct_relevel(n.ind, fact.ord)) %>%
  group_by(n.ind) %>%
  summarise(mean = round(mean(perc.phased),2), sd = round(sd(perc.phased),2), n = n())
colnames(a2) <- c("N.ind.", "Mean (perc)", "Std. Dev", "N.Rep")
formattable(a2)
Z %>%
  mutate(n.ind = fct_relevel(n.ind, fact.ord)) %>%
  ggplot(aes(x=n.ind, y=perc.phased, fill=n.ind)) +
  geom_boxplot() + ylab("Percentage phased") +
  xlab("Number of individuals") + ylim(0, 100)

## Distribution of correct phased markers among the phased ones (% - P1)
p1 <- Z %>%
  mutate(n.ind = fct_relevel(n.ind, fact.ord)) %>%
  ggplot(aes(x=n.ind, y=perc.corect.phase.p1, fill=n.ind)) +
  geom_boxplot() + ylab("percentage phased") +
  xlab("number of individuals") + ylim(0, 100)

## Distribution of correct phased markers among the phased ones (% - P2)
p2 <- Z %>%
  mutate(n.ind = fct_relevel(n.ind, fact.ord)) %>%
  ggplot(aes(x=n.ind, y=perc.corect.phase.p2, fill=n.ind)) +
  geom_boxplot() + ylab("percentage phased") +
  xlab("number of individuals") + ylim(0, 100)

gridExtra::grid.arrange(p1,p2, ncol = 2)
