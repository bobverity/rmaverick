
#------------------------------------------------
# simulate data

sim_data <- function(n = 10, ploidy = 2, loci = 10, allele_num = 5, lambda = 1.0, admix_on = FALSE, alpha = 1, K = 3, pop_col_on = TRUE, ploidy_col_on = TRUE) {
  
  # check inputs
  check_sim_data(n, ploidy, loci, allele_num, lambda, admix_on, alpha, K)
  
  # expand ploidy to vector
  if (length(ploidy)==1) {
    ploidy <- rep(ploidy, n)
  }
  
  # generate names
  ind_names <- paste0("ind", 1:n)
  pop_names <- paste0("pop", 1:K)
  locus_names <- paste0("locus", 1:loci)
  allele_names <- paste0("allele", 1:allele_num)
  
  # generate allele frequencies
  lambda <- rep(lambda, allele_num)
  allele_freqs <- list()
  for (k in 1:K) {
    m <- t(replicate(loci, rdirichlet(lambda)))
    row.names(m) <- locus_names
    colnames(m) <- allele_names
    allele_freqs[[k]] <- m
  }
  names(allele_freqs) <- pop_names
  
  # generate admixture matrix
  alpha <- rep(alpha, K)
  if (!admix_on) {  # no admixture
    admix_freqs <- matrix(0, n, K)
    group <- sort(sample(1:K, n, replace = TRUE))
    admix_freqs[cbind(1:n, group)] <- 1
    
  } else {  # independent admixture for each sample
    admix_freqs = t(replicate(n, rdirichlet(alpha)))
    group <- order(apply(admix_freqs, 1, which.max))
    admix_freqs <- admix_freqs[group,]
  }
  colnames(admix_freqs) <- pop_names
  row.names(admix_freqs) <- ind_names
  
  # generate data
  dat <- matrix(NA, sum(ploidy), loci)
  i2 <- 1
  for (i in 1:n) {
    for (c in 1:ploidy[i]) {
      for (j in 1:loci) {
        g <- sample(1:K, 1, prob=admix_freqs[i,])
        p <- allele_freqs[[g]][j,]
        dat[i2, j] <- sample(1:allele_num, 1, prob = p)
      }
      i2 <- i2 + 1
    }
  }
  colnames(dat) <- locus_names
  
  # create data frame
  df <- data.frame(ID = rep(ind_names, times = ploidy), stringsAsFactors = FALSE)
  if (pop_col_on) {
    df <- cbind(df, pop = rep(group, times = ploidy))
  }
  if (ploidy_col_on) {
    df <- cbind(df, ploidy = rep(ploidy, times = ploidy))
  }
  df <- cbind(df, dat)
  
  # return list
  ret <- list(data = df, allele_freqs = allele_freqs, admix_freqs = admix_freqs, group = group)
  return(ret)
}
