## Run with
## R CMD BATCH --no-save --no-restore '--args seed=1' multiple_simulations.R  res_1.Rout &

## Reading arguments
options(echo=TRUE)
arg <- commandArgs(trailingOnly = TRUE)
print(arg)
for(j in 1:length(arg)){
  eval(parse(text=arg[[j]]))
}
rm(arg)
# seed = 1
require(mappoly)
setwd("~/repos/Test_mappoly/scripts/")
mrk.tail <- c(60, 60, 60, 100, 100)
    n.ph <- c(20, 20, 40, 60,  60)
temp.res <- vector("list", 5)
names(temp.res) <- names(mrk.tail) <- names(n.ph) <- c(200, 100, 50, 30, 20)
RES <-  vector("list", 10)
save(RES, file = paste0("~/repos/Test_mappoly/1200_reps/out/out_", seed, ".rda"))
set.seed(seed)
for(j in 1:10){
  for(i in c(200, 100, 50, 30, 20)){
    cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    cat("~~~~~~~~~~~~~~~~~ rep: ", j, " n.ind:", i, " ~~~~~~~~~~~~~~~~~~~\n")
    ## Reading first rep, 200 individuals
    dat <- read_geno_csv(file.in = paste0("~/repos/Test_mappoly/MS_simulations/data/F1_noff200_DR0_rep1_polyorigin_geno_snparray_mappoly.csv"),
                         ploidy = 4, elim.redundant = FALSE, verbose = FALSE)
    ## Sampling 100, 50,30, and 20 individuals (10 and 15 would take to much time, I do not recommend)
    if(i != 200)
      dat <- sample_data(dat, n = i, type = "individual")
    ## Making a sequence with all markers
    s <- make_seq_mappoly(dat, "all")
    ## Two-point analysis
    time.tpt <- system.time(tpt <- est_pairwise_rf(s, verbose = FALSE))
    ## Map and phase estimation
    time.map <- system.time(
      map <- tryCatch({
        est_rf_hmm_sequential(input.seq = s,
                              start.set = 4,
                              thres.twopt = 5,
                              thres.hmm = 5,
                              extend.tail = mrk.tail[as.character(i)],
                              twopt = tpt,
                              sub.map.size.diff.limit = 8,
                              phase.number.limit = n.ph[as.character(i)],
                              verbose = FALSE)
      }, error = function(e) {NA}))

    if(is.na(map))
    {
      temp.res[[as.character(i)]] <- list(time.tpt = time.tpt[3],
                            tpt = tpt,
                            time.map = time.map[3],
                            map = NA,
                            time.map.err = NA,
                            map.err = NA)
    } else {
      ## Modeling a global error of 5% in the final map
      time.map.err <- system.time(map.err <- est_full_hmm_with_global_error(map, error = 0.05))
      ## Gathering results
      temp.res[[as.character(i)]] <- list(time.tpt = time.tpt[3],
                            tpt = tpt,
                            time.map = time.map[3],
                            map = map,
                            time.map.err = time.map.err[3],
                            map.err = map.err)
      cat("~~~~~~~~~~~~~~~~~~~ n mrk map: ", map$info$n.mrk , " ~~~~~~~~~~~~~~~~~~~~~\n")
      cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    }
  }
  RES[[j]] <- temp.res
}
## save simulation for a given seed
save(RES, file = paste0("~/repos/Test_mappoly/1200_reps/out/out_", seed, ".rda"))
