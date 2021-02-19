## Results
require(mappoly)
require(ggplot2)
require(tidyverse)
require(formattable)
source("~/repos/Test_mappoly/utils.R")
## listing all result files
fls <- list.files(path = "~/repos/Test_mappoly/1200_reps/out", pattern = "out_", full.names = T)
## reading/parsing files
Z <- NULL
for(i in fls){
  cat(i, "\n")
  load(i)
  X <- parse.results(RES)
  Z <- rbind(Z, X)
}
head(Z, 30)
Z<-Z[!is.na(Z$perc.phased),]
## Time (s) distribution for all sample sizes
fact.ord <- c("20",  "30", "50", "100", "200")
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
  xlab("Number of individuals") + ylim(0, 50)

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
  geom_boxplot() + ylab("percentage phased") +
  xlab("number of individuals") + ylim(0, 100)

## Distribution of correct phased markers among the phased ones (% - P1)
p1<-Z %>%
  mutate(n.ind = fct_relevel(n.ind, fact.ord)) %>%
  ggplot(aes(x=n.ind, y=perc.corect.phase.p1, fill=n.ind)) +
  geom_boxplot() + ylab("percentage phased") +
  xlab("number of individuals") + ylim(0, 100) + ggtitle("P1")

## Distribution of correct phased markers among the phased ones (% - P2)
p2<-Z %>%
  mutate(n.ind = fct_relevel(n.ind, fact.ord)) %>%
  ggplot(aes(x=n.ind, y=perc.corect.phase.p2, fill=n.ind)) +
  geom_boxplot() + ylab("percentage phased") +
  xlab("number of individuals") + ylim(0, 100) + ggtitle("P2")

gridExtra::grid.arrange(p1,p2, ncol = 2)
