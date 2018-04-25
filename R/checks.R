
#------------------------------------------------
# checks on inputs to sim_data function
# (not exported)

check_sim_data <- function(n, ploidy, loci, allele_num, lambda, admix_on, alpha, K) {
  
  assert_that(is.int(n))
  assert_that(length(n)==1)
  assert_that(n>0)
  
  assert_that(is.int(ploidy))
  assert_that(length(ploidy)==1 | length(ploidy)==n)
  assert_that(all(ploidy>0))
  
  assert_that(is.int(loci))
  assert_that(length(loci)==1)
  assert_that(loci>0)
  
  assert_that(is.int(allele_num))
  assert_that(length(allele_num)==1)
  assert_that(allele_num>0)
  
  assert_that(is.numeric(lambda))
  assert_that(length(lambda)==1)
  assert_that(lambda>0)
  
  assert_that(is.logical(admix_on))
  assert_that(length(admix_on)==1)
  
  assert_that(is.numeric(alpha))
  assert_that(length(alpha)==1)
  assert_that(alpha>0)
  
  assert_that(is.int(K))
  assert_that(length(K)==1)
  assert_that(K>0)
}

#------------------------------------------------
# check data

check_data <- function() {
  
  print("foo")
  
}
