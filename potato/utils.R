## Wraper for map construction 
mapping_poptato<-function(input.data){
  s <- make_seq_mappoly(input.data, "all")
  time.tpt<-system.time(tpt <- est_pairwise_rf(s))
  s <- rf_snp_filter(tpt)
  time.map<-system.time(map <- est_rf_hmm_sequential(input.seq = s,
                                                     start.set = 3,
                                                     thres.twopt = 10,
                                                     thres.hmm = 10,
                                                     extend.tail = 50,
                                                     twopt = tpt,
                                                     verbose = TRUE,
                                                     tol = 10e-2,
                                                     tol.final = 10e-3,
                                                     phase.number.limit = 40,
                                                     sub.map.size.diff.limit = 5,
                                                     info.tail = TRUE,
                                                     reestimate.single.ph.configuration = TRUE))
  map<-filter_map_at_hmm_thres(map, thres.hmm = 0.001)
  time.maperr<-system.time(final.map <- est_full_hmm_with_global_error(map, error = 0.05))
  return(list(time.tpt = time.tpt, 
              time.map = time.map, 
              time.maperr = time.maperr, 
              map = map, 
              final.map  = final.map))
}

## PolyOrigin to MAPpoly data set
full_pop_to_full_sib <- function(input.data, pedigree, ploidy){
  id <- which(apply(pedigree, 1, function(x) any(is.na(x))))
  w <- apply(pedigree[-id, 2:3], 1, paste, collapse = "x")
  names(w) <- pedigree[-id,1]
  output.data <- vector("list", length(unique(w)))
  names(output.data) <- unique(w)
  for(i1 in unique(w))
  {
    p <- unlist(strsplit(i1, "x"))
    dat <- data.frame(snp_name = input.data[,"marker"],
                      P1 = input.data[,p[1]],
                      P2 = input.data[,p[2]],
                      sequence = input.data[, "chromosome"],
                      sequence.position = input.data[, "position"],
                      input.data[,names(which(w==i1))])
    dat <- dat[!apply(dat[,-c(1:6)], 1, function(x) all(is.na(x))),]
    m <- ploidy
    ## get number of individuals -------------
    n.ind <- ncol(dat) - 5
    ## get number of markers -----------------
    n.mrk <- nrow(dat)
    ## get marker names ----------------------
    mrk.names <- dat[,1]
    ## get individual names ------------------
    ind.names <- colnames(dat)[-c(1:5)]
    ## get dosage in parent P ----------------
    dosage.p <- as.integer(dat[,2])
    ## get dosage in parent Q ----------------
    dosage.q <- as.integer(dat[,3])
    ## monomorphic markers
    dp<-abs(abs(dosage.p-(m/2))-(m/2))
    dq<-abs(abs(dosage.q-(m/2))-(m/2))
    id<-dp+dq!=0
    ## get sequence info ---------------------
    sequence <- as.character(dat[,4])
    ## get sequence position info ------------
    sequencepos <- as.numeric(dat[,5])
    names(sequencepos) <- names(sequence) <- names(dosage.q) <- names(dosage.p) <-  mrk.names
    nphen <- 0
    phen <- NULL
    cat("Reading the following data:")
    cat("\n    Ploidy level:", m)
    cat("\n    No. individuals: ", n.ind)
    cat("\n    No. markers: ", n.mrk)
    cat("\n    No. informative markers:  ", sum(id), " (", round(100*sum(id)/n.mrk,1), "%)", sep = "")
    if (all(unique(nphen) != 0))
      cat("\n    This dataset contains phenotypic information.")
    
    if (length(sequence) > 1)
      cat("\n    This dataset contains sequence information.")
    cat("\n    ...")
    ## get genotypic info --------------------
    geno.dose <- dat[,-c(1:5)]
    dimnames(geno.dose)<-list(mrk.names, ind.names)
    geno.dose[is.na(geno.dose)] <- m + 1
    ## returning the 'mappoly.data' object
    cat("\n    Done with reading.\n")
    geno.dose<-geno.dose[id,]
    names(dosage.p) <- names(dosage.q) <- names(sequence) <- names(sequencepos) <- mrk.names
    res <- structure(list(m = m,
                          n.ind = n.ind,
                          n.mrk = sum(id),
                          ind.names = ind.names,
                          mrk.names = mrk.names[id],
                          dosage.p = dosage.p[id],
                          dosage.q = dosage.q[id],
                          sequence = sequence[id],
                          sequence.pos = sequencepos[id],
                          seq.ref = NULL,
                          seq.alt = NULL,
                          all.mrk.depth = NULL,
                          prob.thres = NULL,
                          geno.dose = geno.dose,
                          nphen = nphen,
                          phen = phen,
                          kept = NULL,
                          elim.correspondence = NULL),
                     class = "mappoly.data")
    ## Filtering non-conforming markers.
    cat("    Filtering non-conforming markers.\n    ...")
    res<-filter_non_conforming_classes(res)
    ##Computing chi-square p.values
    Ds <- array(NA, dim = c(m+1, m+1, m+1))
    for(i in 0:m)
      for(j in 0:m)
        Ds[i+1,j+1,] <- segreg_poly(m = m, dP = i, dQ = j)
    Dpop<-cbind(res$dosage.p, res$dosage.q)
    M<-t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
    dimnames(M)<-list(res$mrk.names, c(0:m))
    M<-cbind(M, res$geno.dose)
    res$chisq.pval<-apply(M, 1, mappoly:::mrk_chisq_test, m = m)
    cat("\n    Done with filtering.\n")
    output.data[[i1]] <- res
  }
  return(output.data)
}

## Get time 
get_time<-function(MAPs){
  tm<-matrix(NA, 12, 3, dimnames = list(paste0("ch", 1:12), c("Pop1", "Pop2", "Pop3")))
  for(j in 1:12){
    for(i in 1:3){
      tm[j,i] <- sum(sapply(MAPs[[j]][[i]][1:3], function(x) x[3]))/60    
    }
  }
  tm
}

## Get map length
get_length<-function(MAPs){
  len<-matrix(NA, 12, 3, dimnames = list(paste0("ch", 1:12), c("Pop1", "Pop2", "Pop3")))
  for(j in 1:12){
    for(i in 1:3){
      len[j,i] <- sum(imf_h(MAPs[[j]][[i]]$final.map$maps[[1]]$seq.rf))
    }
  }
  len
}
  