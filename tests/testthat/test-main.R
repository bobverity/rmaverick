#------------------------------------------------
test_that("run_mcmc fails with incorrect input", {
  
  # simulate some data
  K_sim <- 3
  true_alpha <- 0.1
  n <- 100
  loci <- 10
  
  set.seed(1)
  d <- sim_data(n = n, loci = loci, ploidy = 2, ploidy_col_on = FALSE,
                K = K_sim, admix_on = FALSE)
  
  # test that project fails with no data or parameter sets
  mavproject() %>%
    run_mcmc(K = 2, burnin = 1e3, samples = 1e3, rungs = 5) %>%
    expect_error()
  
  # test that project fails with parameter set but no data
  mavproject() %>%
    new_set(admix_on = FALSE) %>%
    run_mcmc(K = 2, burnin = 1e3, samples = 1e3, rungs = 5) %>%
    expect_error()
  
  # test that project fails with data but no active parameter set
  mavproject() %>%
    bind_data(df = d$data, ID_col = 1, pop_col = 2, ploidy_col = NULL, ploidy = 2) %>%
    run_mcmc(K = 2, burnin = 1e3, samples = 1e3, rungs = 5) %>%
    expect_error()
  
  # define working project with bound data and defined parameters (no-admixture model)
  p <- mavproject() %>%
    bind_data(df = d$data, ID_col = 1, pop_col = 2, ploidy_col = NULL, ploidy = 2) %>%
    new_set(admix_on = FALSE, name = "no-admixture")
  
  # check that MCMC runs with single rung
  run_mcmc(p, K = 2, burnin = 1e3, samples = 1e3, rungs = 1) %>%
    expect_error(NA)
  
  # check that MCMC runs with multiple rungs
  run_mcmc(p, K = 2, burnin = 1e3, samples = 1e3, rungs = 5) %>%
    expect_error(NA)
  
  # equivalent checks for admixture model with alpha fixed
  p <- p %>%
    new_set(admix_on = TRUE, alpha = 1.0, estimate_alpha = FALSE, name = "admixture, alpha fixed")
  
  run_mcmc(p, K = 2, burnin = 1e3, samples = 1e3, rungs = 1) %>%
    expect_error(NA)
  run_mcmc(p, K = 2, burnin = 1e3, samples = 1e3, rungs = 5) %>%
    expect_error(NA)
  
  # equivalent checks for admixture model with alpha free
  p <- p %>%
    new_set(admix_on = TRUE, alpha = 1.0, estimate_alpha = TRUE, name = "admixture, alpha free")
  
  run_mcmc(p, K = 2, burnin = 1e3, samples = 1e3, rungs = 1) %>%
    expect_error(NA)
  run_mcmc(p, K = 2, burnin = 1e3, samples = 1e3, rungs = 5) %>%
    expect_error(NA)
  
})

#------------------------------------------------
test_that("plotting functions run without error", {
  
  # simulate some data
  K_sim <- 3
  true_alpha <- 0.1
  n <- 100
  loci <- 10
  
  set.seed(1)
  d <- sim_data(n = n, loci = loci, ploidy = 2, ploidy_col_on = FALSE,
                K = K_sim, admix_on = FALSE)
  
  # define project, bind data, define parameters (no-admixture model), run MCMC
  p <- mavproject() %>%
    bind_data(df = d$data, ID_col = 1, pop_col = 2, ploidy_col = NULL, ploidy = 2) %>%
    new_set(admix_on = FALSE) %>%
    run_mcmc(K = 2, burnin = 1e3, samples = 1e3, rungs = 5)
  
  # expect error plots as no-admixture model
  plot_trace(p) %>% expect_error()
  plot_acf(p) %>% expect_error()
  plot_alpha(p) %>% expect_error()
  
  # expect no error plots
  plot_GTI_path(p) %>% expect_error(NA)
  plot_trace(p, param = "loglike") %>% expect_error(NA)
  plot_acf(p, param = "loglike") %>% expect_error(NA)
  plot_density(p, param = "loglike") %>% expect_error(NA)
  plot_logevidence_K(p) %>% expect_error(NA)
  plot_logevidence_model(p) %>% expect_error(NA)
  plot_loglike(p) %>% expect_error(NA)
  plot_loglike_quantiles(p) %>% expect_error(NA)
  plot_posterior_K(p) %>% expect_error(NA)
  plot_posterior_model(p) %>% expect_error(NA)
  plot_qmatrix(p) %>% expect_error(NA)
  plot_mc_acceptance(p) %>% expect_error(NA)
  
  # same again for admixture model (i.e. alpha in output)
  p <- mavproject() %>%
    bind_data(df = d$data, ID_col = 1, pop_col = 2, ploidy_col = NULL, ploidy = 2) %>%
    new_set(admix_on = TRUE, estimate_alpha = TRUE) %>%
    run_mcmc(K = 2, burnin = 1e3, samples = 1e3, rungs = 5)
  
  # expect no error plots
  plot_trace(p) %>% expect_error(NA)
  plot_acf(p) %>% expect_error(NA)
  plot_alpha(p) %>% expect_error(NA)
  plot_GTI_path(p) %>% expect_error(NA)
  plot_trace(p, param = "loglike") %>% expect_error(NA)
  plot_acf(p, param = "loglike") %>% expect_error(NA)
  plot_density(p, param = "loglike") %>% expect_error(NA)
  plot_logevidence_K(p) %>% expect_error(NA)
  plot_logevidence_model(p) %>% expect_error(NA)
  plot_loglike(p) %>% expect_error(NA)
  plot_loglike_quantiles(p) %>% expect_error(NA)
  plot_posterior_K(p) %>% expect_error(NA)
  plot_posterior_model(p) %>% expect_error(NA)
  plot_qmatrix(p) %>% expect_error(NA)
  plot_mc_acceptance(p) %>% expect_error(NA)
  
})