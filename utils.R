## Parse results
parse.results <- function(RES, ploidy = 4){
  ## First, match the data:
  U <- NULL
  for(k in 1:length(RES)){
    X <- RES[[k]]
    ## Though I don't have the simulated phases, I am assuming 200 individual phase as gold standard
    true.phase <- cbind(ph_list_to_matrix(X[[1]]$map.err$maps[[1]]$seq.ph$P, ploidy),
                        ph_list_to_matrix(X[[1]]$map.err$maps[[1]]$seq.ph$Q, ploidy))
    dimnames(true.phase) <- list(X[[1]]$map.err$info$mrk.names, paste0("h", 1:(2*ploidy)))
    V<-NULL
    for(i in 1:length(X)){
      if(any(unlist(sapply(X[[i]], is.na)))){
        V <- rbind(V, rep(NA, 5))
      } else {
        est.phase <- cbind(ph_list_to_matrix(X[[i]]$map.err$maps[[1]]$seq.ph$P, ploidy),
                           ph_list_to_matrix(X[[i]]$map.err$maps[[1]]$seq.ph$Q, ploidy))
        dimnames(est.phase) <- list(X[[i]]$map.err$info$mrk.names,
                                    paste0("h", 1:(2*ploidy)))
        mrk.id <- intersect(rownames(true.phase), rownames(est.phase))
        true.mapped <- true.phase[mrk.id,]
        est.mapped <- as.matrix(est.phase[mrk.id,])

        ##Generate all permutations of the phasing order
        perms <- mappoly::perm_tot(1:ploidy)

        p1.cors <- apply(perms,1,function(x) sum(abs(true.mapped[,1:ploidy] - est.mapped[,1:ploidy][,x])))
        best.p1 <- which.min(p1.cors)

        p2.cors <- apply(perms,1,function(x) sum(abs(true.mapped[,(ploidy + 1):(2*ploidy)] - est.mapped[,(ploidy + 1):(2*ploidy)][,x])))
        best.p2 <- which.min(p2.cors)

        nm <- nrow(est.mapped)

        est.mapped2 <- cbind(est.mapped[,1:ploidy][,perms[best.p1,]],
                             est.mapped[,(ploidy + 1):(2*ploidy)][,perms[best.p2,]])

        errors <- abs(true.mapped - est.mapped2)

        ## Return the percentage and number of correctly phased, also per parent:
        np1 <- sum(rowSums(errors[,1:ploidy]) > 0)
        np2 <- sum(rowSums(errors[,(ploidy + 1):(2*ploidy)]) > 0)

        v <- c(100*nm/nrow(true.phase), nm - np1, nm - np2, 100*(nm - np1)/nm, 100*(nm - np2)/nm)
        V <- rbind(V, v)
      }
    }
    V <- data.frame(names(X),
                    V,
                    sapply(X, function(X) sum(unlist(X[c(1,3,5)]), na.rm = TRUE)/60))

    dimnames(V) <- list(names(X),
                        c("n.ind",
                          "perc.phased",
                          "n.corect.phase.p1",
                          "n.corect.phase.p2",
                          "perc.corect.phase.p1",
                          "perc.corect.phase.p2",
                          "time.sec"))
    U <- rbind(U, V)
  }
  return(U)
}
