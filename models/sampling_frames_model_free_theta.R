## This model is adapted from Danielle Navarro's sampling_frames_model found
## in https://github.com/djnavarro/sampling-differences . 99.9% of this model code
## is her work, and it has been adapted here to allow for relatively easy fitting
## procedures.

sampling_frames_model <- function(
  tau = 3,     # gaussian process prior: baseline smoothness constraint
  rho = 1,     # gaussian process prior: decay rate for smoothness constraint
  sigma = .5,  # gaussian process prior: inherent noise level in the data
  mu = .05,   # prior mean (on the probability scale, not logit scale)
  theta = 1.0,  # weight of property sampling, with 1-theta being the weight of category,
  x_grid=1:6, # location of test items--feed in outside function to improve speed.
  hypotheses, # hypothesis space
  conditions # sample size conditions [e.g., c(2,8,20)]
) {

  `%>%` <- magrittr::`%>%` # use the pipe in the local environment

  # convenience functions ---------------------------------------------------

  # logit function: phi is a logistic function of f, so f is a logit of phi
  logit <- function(phi) {
    log(phi/(1-phi))
  }

  # the psychological distance between values is just the difference between
  # coordinates in the stimulus space. this is a little unrealistic, but good
  # enough as an approximation in the current setting
  distance <- function(xi, xj) {
    abs(xi - xj)
  }

  # the kernel function
  kernel <- function(dist) {
    tau^2 * exp(-rho * dist)
  }

  # inherent noise
  noise <- function() {
    diag(sigma^2, length(x_grid), length(x_grid))
  }

  # normalise to 1
  normalise <- function(x) {x / sum(x)}

  # prior mean vector is constant
  prior_mean <- function() {
    logit(rep.int(x = mu, times = length(x_grid)))
  }

  # covariance matrix
  covariance_matrix <- function() {
    distance_matrix <- outer(x_grid, x_grid, distance)
    kernel(distance_matrix) + noise()
  }

  # transform the data frame of hypotheses on the probability
  # scale into a matrix of hypotheses on the logit scale
  logit_scale_hypotheses <- function() {
    logit(as.matrix(hypotheses))
  }

  logit_scale_hypotheses_data <- function(data) {
    logit(data)
  }


  # gaussian process prior --------------------------------------------------

  model_prior <- function() {
    mvtnorm::dmvnorm(
      x = logit_scale_hypotheses(),
      mean = prior_mean(),
      sigma = covariance_matrix()
    ) %>%
      normalise()
  }

  # the three possible sampling models --------------------------------------

  # sampling models when the items are selected *from* the extension
  # of the category, or the extension of the property... for property
  # sampling the item can be sampled at any location in the stimulus
  # space so long as it is property-positive...
  property_sampling <- function(stimulus) {
    numerator <- eval(
      expr = substitute(stimulus),
      envir = hypotheses
    )
    denominator <- eval(
      expr = quote(test1 + test2 + test3 + test4 + test5 + test6),
      envir = hypotheses
    )
    return(numerator / denominator)
  }


  # category sampling is also a strong sampling model, but where we
  # are restricted to sampling items that belong to some subset of the
  # stimulus space (regardless of whether they are property positive).
  # In this context, that subset is simply the first two test items.
  category_sampling <- function(stimulus) {
    numerator <- eval(
      expr = substitute(stimulus),
      envir = hypotheses
    )
    denominator <- 2
    return(numerator / denominator)
  }

  # weak sampling assumes the items are sampled uniformly at random from
  # all items, and the fact that the observed item turned out to be property
  # positive is coincidental. in this experimental design there is not much
  # of a difference between weak sampling and category sampling... NOTE:
  # the weak_sampling() function does not get called anymore, because I have
  # fixed theta = 1. See comments below.
  weak_sampling <- function(stimulus) {
    numerator <- eval(
      expr = substitute(stimulus),
      envir = hypotheses
    )
    denominator <- 6
    return(numerator / denominator)
  }


  # sampling frame imposes likelihood ---------------------------------------

  # the samples observed by the participant occur only in
  # the locations corresponding to test1 and test2, and are
  # always property positive. depending on condition the
  # participant may have seen 2, 8 or 20 items. for simplicity
  # we assume that they are evenly split between test1 and
  # test 2. conditional on knowing the probabilities of at
  # each test point (and assuming the base rates for entities
  # at each location are uniform) the likelihood is proportional
  # to the function value at the location, where the normalisation
  # depends on the admissability rule specified by the frame


  model_likelihood <- function(cond) {

    # read off condition parameters
    n <- cond$n

    # the learner is assumed to apply a mixture of strong and weak sampling,
    # per the model in Navarro et al (2012). in practice, because weak and
    # category sampling are indistinguishable in this design, the role of
    # the theta parameter is to shift the generalisation gradients in the
    # property sampling condition only (in the direction of the category/weak
    # sampling gradients)...
    #
    # ...but, one thing that emerged from the exploratory analysis is that
    # the Navarro et al (2012) formalism isn't buying us much here, so I have
    # fixed the value of theta to be 1 rather than treat it as a free parameter.
    # this choice (and its implications) are noted in the scratchpad log.
    # SXL: Attempting to rewrite theta to essentially "measure" the effect of
    # sampling frames. May be theoretically problematic, as indicated in DJN's
    # scratchpad... but may be useful for individual fitting?
    # Might just assume that weak sampling doesn't exist here--only property or category sampling.

    prob1 <- theta * property_sampling(test1) + (1 - theta) * category_sampling(test1)
    prob2 <- theta * property_sampling(test2) + (1 - theta) * category_sampling(test2)

    # observations are conditionally iid so the probability of the sample
    # is the product of the individual observations:
    sample_likelihood <- prob1^(cond$n/2) * prob2^(cond$n/2)

    return(sample_likelihood)
  }


  # compute generalisation gradient -----------------------------------------

  expected_value <- function(probability, value) {
    sum(probability * value)
  }

  # to compute the generalisation gradient we first compute the posterior
  # distribtution over hypotheses, given the observed data and the sampling
  # condition. then, for each test item, the posterior expected probability
  # that a new item possesses the property is just the expected value of
  # the unknown function at that location in the stimulus space:

  generalisation_gradient <- function(cond, prior) {

    likelihood <- model_likelihood(cond)
    posterior <- normalise(likelihood * prior)
    generalisation <- tibble::tibble(
      test1 = expected_value(posterior, hypotheses$test1),
      test2 = expected_value(posterior, hypotheses$test2),
      test3 = expected_value(posterior, hypotheses$test3),
      test4 = expected_value(posterior, hypotheses$test4),
      test5 = expected_value(posterior, hypotheses$test5),
      test6 = expected_value(posterior, hypotheses$test6)
    )

    return(generalisation)
  }

  # do the work and return to user ------------------------------------------

  # this is the important part of the computation
  result <- conditions %>%
    purrr::transpose() %>%
    purrr::map_dfr(generalisation_gradient, prior = model_prior()) %>%
    dplyr::bind_cols(conditions, .)

  # clean it up and make it nice and pretty
  output <- result %>%
    tidyr::pivot_longer(
      cols = tidyselect::starts_with("test"),
      names_to = "test_item",
      values_to = "response"
    ) %>%
    dplyr::rename(
      #sampling_frame = frame,
      sample_size = n
    ) %>%
    dplyr::mutate(
      test_item = as.numeric(stringr::str_remove_all(test_item, "test")),
      source = "model"
    )

  return(output)
}



#----------------------------------------
#### Same model but used for ML fitting
sampling_frames_model_likelihood <- function(data, n,
  tau = 3,     # gaussian process prior: baseline smoothness constraint
  rho = 1,     # gaussian process prior: decay rate for smoothness constraint
  sigma = .5,  # gaussian process prior: inherent noise level in the data
  mu = .05,   # prior mean (on the probability scale, not logit scale)
  theta = 1.0  # weight of property sampling, with 1-theta being the weight of category
) {

  `%>%` <- magrittr::`%>%` # use the pipe in the local environment

  # useful quantities -------------------------------------------------------

  x_grid = 1:6                # locations of the test items
  #y_grid <- seq(.1, .9, .2)   # possible values for the probability


  # convenience functions ---------------------------------------------------

  # logit function: phi is a logistic function of f, so f is a logit of phi
  logit <- function(phi) {
    log(phi/(1-phi))
  }

  # the psychological distance between values is just the difference between
  # coordinates in the stimulus space. this is a little unrealistic, but good
  # enough as an approximation in the current setting
  distance <- function(xi, xj) {
    abs(xi - xj)
  }

  # the kernel function
  kernel <- function(dist) {
    tau^2 * exp(-rho * dist)
  }

  # inherent noise
  noise <- function() {
    diag(sigma^2, length(x_grid), length(x_grid))
  }

  # normalise to 1
  normalise <- function(x) {x / sum(x)}

  # prior mean vector is constant
  prior_mean <- function() {
    logit(rep.int(x = mu, times = length(x_grid)))
  }

  # covariance matrix
  covariance_matrix <- function() {
    distance_matrix <- outer(x_grid, x_grid, distance)
    kernel(distance_matrix) + noise()
  }

  # transform the data on the probability
  # scale into a matrix of hypotheses on the logit scale
  logit_scale_hypotheses_data <- function(data) {
    logit(data)
  }


  # gaussian process prior --------------------------------------------------

  model_prior_data <- function(data) {
    mvtnorm::dmvnorm(
      x = logit_scale_hypotheses_data(data),
      mean = prior_mean(),
      sigma = covariance_matrix()
    )
  }

  # the three possible sampling models --------------------------------------

  # sampling models when the items are selected *from* the extension
  # of the category, or the extension of the property... for property
  # sampling the item can be sampled at any location in the stimulus
  # space so long as it is property-positive...

  property_sampling_data <- function(stimulus,data) {
    numerator <- stimulus
    denominator <- sum(data)
    return(numerator / denominator)
  }


  # category sampling is also a strong sampling model, but where we
  # are restricted to sampling items that belong to some subset of the
  # stimulus space (regardless of whether they are property positive).
  # In this context, that subset is simply the first two test items.

  category_sampling_data <- function(stimulus,data) {
    numerator <- stimulus
    denominator <- 2
    return(numerator / denominator)
  }


  # sampling frame imposes likelihood ---------------------------------------

  model_likelihood_data <- function(data,n) {
    #model likelihood given some specific data
    # read off condition parameters

    prob1 <- theta * property_sampling_data(data[1],data) + (1 - theta) * category_sampling_data(data[1],data)
    prob2 <- theta * property_sampling_data(data[2],data) + (1 - theta) * category_sampling_data(data[2],data)

    # observations are conditionally iid so the probability of the sample
    # is the product of the individual observations:
    sample_likelihood <- prob1^(n/2) * prob2^(n/2)

    return(sample_likelihood)
  }

  # compute generalisation gradient -----------------------------------------

  expected_value <- function(probability, value) {
    sum(probability * value)
  }


  #Xian: We can calculate (proportional) likelihoods using model_likelihood and prior, given the data
  generalisation_likelihood <- function(data, n){
    likelihood <- model_likelihood_data(data, n) * model_prior_data(data)
    return(likelihood)
  }

  output <- generalisation_likelihood(data, n)
  return(output)
}
