BrainFeatureAnalysis <- function(feature_data, including_equal_mixture_model = T)
{
  hypothesis_results <- applyHypothesisOverResiduals(feature_data,
                                                     including_equal_mixture_model)
  return(list(
    hypothesis_results = hypothesis_results
  ))
}

#' Apply a set of hypothesis test over each brain feature
#' @param feature_data a tibble of scaled residuals for brain feature
#' @return test results
applyHypothesisOverResiduals <- function(feature_data,
                                         including_equal_mixture_model = T) {

  t_test_for_difference_between_genders <-
    t.test(feature_data$value[feature_data$sex == 0],
           feature_data$value[feature_data$sex == 1])

  cohen_d_test <- cohen.d(feature_data$value[feature_data$sex == 0],
                          feature_data$value[feature_data$sex == 1])

  single_population_model <- singleDistributionHypothesisLoglikelihood(feature_data$value)
  pure_types_model <- pureTypeMixtureModelHypothesis(feature_data)
  mixture_model <- tryCatch(
    {
      doubleDoubleEM(feature_data,T)
    },
    error=function(cond)
    {
      return(data.frame(llk = -Inf))
    }
  )
  if(including_equal_mixture_model){
    equal_proportions_model <- computeEqualProprtionsEM(feature_data)
    equal_proportions_vs_mixture_model_llr <- mixture_model$llk - equal_proportions_model$loglik
  }
  else{
    equal_proportions_model <- NULL
    equal_proportions_vs_mixture_model_llr <- NULL
  }
  single_population_vs_mixture_model_llr <- mixture_model$llk - single_population_model
  pure_types_vs_mixture_model_llr <- mixture_model$llk - pure_types_model$llk

  return(list(
    t_test_for_difference_between_genders =
      t_test_for_difference_between_genders,
    cohen_d_test = cohen_d_test,
    simple_vs_compostie_llrt = single_population_vs_mixture_model_llr,
    pure_types_vs_mixture_model_llr = pure_types_vs_mixture_model_llr,
    equal_proportions_vs_mixture_model_llr = equal_proportions_vs_mixture_model_llr,
    pure_types_model = pure_types_model,
    equal_proportions_model = equal_proportions_model,
    mixture_model = mixture_model
  ))
}

# Pure type analysis ------------------------------------------------------

#' Computes the MLE for the single distribution model
#' @param feature_data a tibble of scaled residuals for brain feature
#' @return loglikelihood for the aromentioned model

singleDistributionHypothesisLoglikelihood <- function(feature_data)
{
  null_mean <- mean(feature_data)
  null_sd <- sd(feature_data)
  null_hypothesis_llk <- sum(dnorm(feature_data,null_mean,null_sd,T))
  return(null_hypothesis_llk)
}

#' Computes the MLE for the pure-types model
#' @param feature_data a tibble of scaled residuals for brain feature
#' @return loglikelihood for the aromentioned model

pureTypeMixtureModelHypothesis <- function(feature_data)
{
  pure_type_model <- logNormalNullLlkComputation(split_to_gender,T)
  return(pure_type_model)
}

# Computing the -llk under the null hypothesis
#' @param feature_data a tibble of scaled residuals for brain feature
#' @return loglikelihood for the aromentioned model

logNormalNullLlkComputation <- function(feature_data)
{
  men <- feature_data$value[feature_data$bio_sex == 1]
  women <- feature_data$value[feature_data$bio_sex == 0]
  theta_mas <- mean(men)
  theta_fem <- mean(women)
  p <- 1
  q <- 0
  sigma_2_men = sd(men)
  sigma_2_women = sd(women)
  parameter_row <- c(theta_mas,theta_fem,p,q,sigma_2_men,sigma_2_women)
  llk <- logLikelihoodComputationForFeatureSet(feature_data,
                                               parameter_row,
                                               is_bio_bank_data)
  return(list(llk = llk,
              men_mean = theta_mas,
              women_mean = theta_fem,
              sigma_2_men = sigma_2_men,
              sigma_2_women = sigma_2_women))
}

#' Computes the loglikelihood for each hypothesis
#' @param feature_data a tibble of scaled residuals for brain feature
#' @param parameter_row a vector with at most 6 different parameters
#' @return loglikelihood for the aromentioned model

logLikelihoodComputationForFeatureSet <- function(feature_data,
                                                  parameter_row,
                                                  ...)
{
  men_data <- feature_data$value[feature_data$bio_sex == 1]
  women_data <- feature_data$value[feature_data$bio_sex == 0]
  theta_mas <- parameter_row[1]
  theta_fem <- parameter_row[2]
  p <- parameter_row[3]
  q <- parameter_row[4]
  sigma_2_men <- parameter_row[5]
  sigma_2_women <- parameter_row[6]
  llk <- computeLogLikelihood(men_data,
                              women_data,
                              theta_mas,
                              theta_fem,
                              p,
                              q,
                              sigma_2_men,
                              sigma_2_women)
  return(llk)
}

computeLogLikelihood <- function(men,women,theta_mas, theta_fem, p, q, sigma_2_men, sigma_2_women)
{

  men_llk <- sum(log(p * dnorm(men, theta_mas, sqrt(sigma_2_men), log = F) +
                       (1 - p) * dnorm(men, theta_fem, sqrt(sigma_2_women), log = F)))
  women_llk <- sum(log(q * dnorm(women, theta_mas, sqrt(sigma_2_men), log = F) +
                         (1-q) * dnorm(women, theta_fem, sqrt(sigma_2_women), log = F)))

  return(men_llk + women_llk)
}

computeLogLikelihoodFull <- function(observations,mle_parameters)
{
  men <- observations$men
  women <- observations$women
  mu_1 <- mle_parameters$mu_1
  mu_2 <- mle_parameters$mu_2
  p <-  mle_parameters$p
  q <-  mle_parameters$q
  sigma_2_men <- mle_parameters$sigma_2_men
  sigma_2_women <- mle_parameters$sigma_2_women

  men_llk <- sum(log(p * dnorm(men, mu_1, sqrt(sigma_2_men), log = F) +
                       (1 - p) * dnorm(men, mu_2, sqrt(sigma_2_women), log = F)))
  women_llk <- sum(log(q * dnorm(women, mu_1, sqrt(sigma_2_men), log = F) +
                         (1-q) * dnorm(women, mu_2, sqrt(sigma_2_women), log = F)))

  return(men_llk + women_llk)
}

# Mixed-types analysis ------------------------------------------------------

#' Computes the MLE for the mixed-types hypothesis using EM algorithm
#' @param feature_data a tibble of scaled residuals for brain feature
#' @param max_iters maximal iterations
#' @param init_parameters initial values for the EM parameters
#' @return EM results a list of the E-step parameters, the M-step parameters, the loglikelihood and the responsebilities for each sex

doubleDoubleEM <- function(feature_data
                           ,max_iter = 3000
                           ,init_parameters = NULL, ...)
{
  full_data <- feature_data
  if(!setequal(names(feature_data),c("men","women")))
  {
    observations <- prepareDataBioBank(feature_data)
  }

  if(is.null(init_parameters) || !validateInitEMParameters(init_parameters))
  {
    p <- 0.51
    q <- 0.49
    men_bar <- mean(observations$men)
    women_bar <- mean(observations$women)
    distributions_parameters <- computeInitParameters(p,q,men_bar,women_bar)
    m_parameters <- list(p = p,
                         q = q,
                         mu_1 = distributions_parameters$mu_1,
                         mu_2 = distributions_parameters$mu_2,
                         sigma_2_men = distributions_parameters$sigma_2,
                         sigma_2_women = distributions_parameters$sigma_2)
  }
  sample_list <- c(1000 + 1:5,2000 + 1:5 ,3000 + 1:5)
  llk <- NULL
  stopping_condition <- F
  m_samples <- c()
  for(i in 1:max_iter)
  {
    # Estep
    e_parameters <- EStep(m_parameters,observations)
    # Mstep
    m_parameters <-
      tryCatch(
        {
          MStep(e_parameters,observations,do_not_ignore_NA = T)
        },
        error=function(cond)
        {
	  stderr('EM encountered a problem with posterior computation')
          MStep(e_parameters,observations,do_not_ignore_NA = F)
        }
        )
    if(i %in% sample_list)
    {
      m_samples <- rbind(m_samples,unlist(m_parameters))
    }
    # llk
    llk[i] <- computeLogLikelihoodFull(observations,m_parameters)
  }
  men_responsebilities <- tibble(eid = full_data[full_data$sex == 1,]$eid,
                                  responsebility = e_parameters$I)
  women_responsebilities <- tibble(eid = full_data[full_data$sex == 0,]$eid,
                                  responsebility = e_parameters$J)
  em_results <- list(e_parameters = e_parameters,
                     m_parameters = m_parameters,
                     m_samples = m_samples,
                     all.llk = llk,
                     llk = tail(llk,1),
                     men_responsebilities = men_responsebilities,
                     women_responsebilities = women_responsebilities)
  return(em_results)
}

prepareDataBioBank <- function(observations)
{
  men <- observations$value[observations$sex == 1]
  women <- observations$value[observations$sex == 0]
  obs <- list(men = men, women = women)
  return(obs)
}

validateInitEMParameters <- function(init_parameters)
{
  if(names(init_parameters) == c("p","q","mu_1","mu_2","sigma_2"))
  {
    if(init_parameters$q < 0 || init_parameters$q > 1)
    {
      return(F)
    }
    if(init_parameters$p < 0 || init_parameters$p > 1)
    {
      return(F)
    }
    if(init_parameters$sigma_2 < 0)
    {
      return(F)
    }
    return(T)
  }
  return(F)
}

computeInitParameters <- function(p,q,men_bar,women_bar)
{
  if((p == q) & (p == 0.5))
  {
    warning("Cannot start for 50-50 mixture model")
    p = 0.51
    q = 0.49
  }
  mu_1 = ((1-q) * men_bar - (1-p) * women_bar) / (p*(1-q) - q*(1-p))
  mu_2 = (q * men_bar - p * women_bar) / (q*(1-p) - p*(1-q))
  sigma_2 = max((mu_1 - mu_2)^2*(p+q)/2*(2-p-q)/2,2)
  init_paramters <- list(mu_1 = mu_1, mu_2 = mu_2, sigma_2 = sigma_2)
  return(init_paramters)
}

#### E-Step #####

EStep <- function(parameters, observations, warn = NULL, ...)
{
  # the parameters from iteration (t-1) are stored in a list
  stopifnot(class(parameters) == "list")
  # the parameters from iteration (t-1) are 5
  stopifnot(length(names(parameters)) == 6)
  # Validate that all the parameters are in the parameter list and are valid
  stopifnot(validateEStepParameters(parameters))

  p <- parameters$p; q <- parameters$q; mu_1 <- parameters$mu_1;
  mu_2 <- parameters$mu_2; sigma_men <- sqrt(parameters$sigma_2_men);
  sigma_women <- sqrt(parameters$sigma_2_women);

  men <- observations$men
  women <- observations$women

  I <- computeConditionalExpectation(men,p,mu_1,mu_2,sigma_men,sigma_women)
  J <- computeConditionalExpectation(women,q,mu_1,mu_2,sigma_men,sigma_women)

  e_step <- list(I = I, J = J)
  return(e_step)
}

validateEStepParameters <- function(parameters)
{
  if(parameters$p < 0 || parameters$p > 1)
  {
    return(F)
  }
  if(parameters$q < 0 || parameters$q > 1)
  {
    return(F)
  }
  if(parameters$sigma_2_men < 0)
  {
    return(F)
  }
  if(parameters$sigma_2_women < 0)
  {
    return(F)
  }
  return(T)
}

computeConditionalExpectation <- function(observations, prop, mean_1, mean_2, std_1, std_2)
{
  cond_exp <- prop * dnorm(observations,mean_1,std_1) /
    (prop * dnorm(observations,mean_1,std_1) + (1-prop)*dnorm(observations,mean_2,std_2))
  return(cond_exp)
}

#### M-Step #####
MStep <- function(parameters
                  ,observations
                  ,do_not_ignore_NA = T,
                  warn = NULL, ...)
{
  # the parameters from iteration (t-1) are stored in a list
  stopifnot(class(parameters) == "list")
  # the parameters from iteration (t-1) are 2
  stopifnot(length(names(parameters)) == 2)
  # Validate that all the parameters are in the parameter list and are valid
  if(do_not_ignore_NA)
  {
    stopifnot(validateMStepParameters(parameters))
    I <- parameters$I;
    J <- parameters$J
    men <- observations$men
    women <- observations$women
  }
  else
  {
    I <- parameters$I[!is.na(parameters$I)]
    J <- parameters$J[!is.na(parameters$J)]
    men <- observations$men[!is.na(parameters$I)]
    women <- observations$women[!is.na(parameters$J)]
  }

  m <- length(men); n <- length(women)
  # M-Step
  p <- mean(I)
  q <- mean(J)
  mu_1 <- as.numeric((I %*% men + J %*% women) / (sum(I) + sum(J)))
  mu_2 <- as.numeric((sum(men) + sum(women) - (I %*% men + J %*% women)) /
    (m + n - (sum(I) + sum(J))))
  sigma_2 <- as.numeric((I %*% (men - mu_1)^2 + (1-I) %*% (men - mu_2)^2 +
                J %*% (women - mu_1)^2 + (1-J) %*% (women - mu_2)^2) /
    (m + n))
  sigma_2_men <- as.numeric((I %*% (men - mu_1)^2  + J %*% (women - mu_1)^2) /
                          (sum(I) + sum(J)))
  sigma_2_women <- as.numeric(
    ((1-I) %*% (men - mu_2)^2 +  (1-J) %*% (women - mu_2)^2) /
                            (m + n - (sum(I) + sum(J)))
    )
  if(sigma_2_men < sigma_2_women / 3)
  {
    sigma_2_men <- sigma_2_women / 3
  }
  m_step <- list(p = p, q = q, mu_1 = mu_1, mu_2 = mu_2,
                 sigma_2_men = sigma_2_men, sigma_2_women = sigma_2_women)
  return(m_step)
}

validateMStepParameters <- function(parameters)
{
  if(any(parameters$I < 0) || any(parameters$I > 1))
  {
    return(F)
  }
  if(any(parameters$J < 0) || any(parameters$J > 1))
  {
    return(F)
  }
  return(T)
}

# Mixed-types equal proportions analysis ------------------------------------------------------

computeEqualProprtionsEM <- function(feature){
  centers_for_em <- feature %>%
    group_by(sex) %>%
    select(sex, value) %>%
    summarise(avg = mean(value), std = sd(value))
  equal_proportion_em <- normalmixEM(feature$value, lambda = 0.5, mu = centers_for_em$avg, sigma = centers_for_em$std, maxit = 10000)
  return(equal_proportion_em)
}

