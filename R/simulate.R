
pick_snp = function(maf, snpdat, ldscore=NULL, sample_proportion=0.1){
#  if(ldscore!=NULL){ ; }
  snpdat$FREQ1[snpdat$FREQ1>0.5] = 1- snpdat$FREQ1[snpdat$FREQ1>0.5]
  snpdat$diffs = abs(maf - snpdat$FREQ1)
  snpdat_2choose = order(snpdat$diffs)[1:max(1,dim(snpdat)[1]*sample_proportion)]
#  hist(snpdat$FREQ1[snpdat_2choose], xlim=c(0,0.5))
  idx = sample(snpdat_2choose, size=1)
  return(idx)
}

simulate_unscaled_betas = function(effect, idx, ld_matrix){
  n_snps = dim(ld_matrix)[1]
  effects = rep(0, n_snps)
  effects[idx] = effect

  noises = rnorm(n_snps,0,1)

  print(system.time(chol_ld_matrix <- chol(ld_matrix)))
  print( dim(chol_ld_matrix) )
#  C = sqrt(1/N)*chol_ld_matrix
  C = chol_ld_matrix
  print(system.time( D_I <- chol2inv(chol_ld_matrix) ))
  print( dim(D_I) )
  print(system.time( betas_ld <- ld_matrix %*% effects ))
  print(system.time( noises_ld <- t(C) %*% noises ))

  betas = betas_ld + noises_ld
#  z = betas * sqrt(N)

  return(list(true_betas=betas_ld, betas=betas))
}


