#' Wrapper function for applyLinearModelOverBrainFeature
#' @param feature_name - brain region name
#' @param subject_data - subject's sex (data.frame with two columns - id and sex)
#' @param total_brain_volume_data - subjects total brain volume
#' @return tibble with three columns - subject's id, subject's sex, subject's residual volume

ComputeLogMeasureResiduals <- function(feature_name, subject_data, brain_volume,
                                              winsorized = F,normalized = F)
{
  feature_residuals <- applyLinearModelOverBrainFeature(feature_name = feature_name,
                                                        subject_data = subject_data,
                                                        total_brain_volume_data = brain_volume,
                                                        winsorized = winsorized,
                                                        normalized = normalized)  
  return(feature_residuals)
}

#' Compute the residuals of log-log regression of log feature volume vs log total brain volume
#' @param feature_name - brain region name
#' @param subject_data - subject's sex (data.frame with two columns - id and sex)
#' @param total_brain_volume_data - subjects total brain volume
#' @return tibble with three columns - subject's id, subject's sex, subject's residual volume

applyLinearModelOverBrainFeature <- function(feature_name,
                                             subject_data,
                                             total_brain_volume_data,
                                             winsorized = F,
                                             normalized = F)
{
  filtering_criteria <- !is.na(relevant_data[,feature_name]) & 
    relevant_data[,feature_name] > 0
  
  feature_data <- tibble(
    eid = subject_data$eid[filtering_criteria],
    sex = subject_data$Sex[filtering_criteria],
    value = relevant_data[,feature_name][filtering_criteria])
  
  data_for_lm <- tibble(y = feature_data$value,
                        x1 = total_brain_volume_data[filtering_criteria,1],
                        x2 = total_brain_volume_data[filtering_criteria,2])
  # Log
  data_for_lm_log_scale <- data_for_lm %>% 
    filter(y > 0) %>% 
    mutate(log_y = log(y), log_x1 = log(x1), log_x2 = log(x2))
  
  # Removing outliers - optional
  if(winsorized)
  {
    upper_and_lower <- quantile(data_for_lm_log_scale$log_y, c(0.01,0.99))
    data_for_lm_log_scale$trimmed_log_y <- Winsorize(data_for_lm_log_scale$log_y, 
                               probs = c(0.01,0.99))
  }
  else
  {
    data_for_lm_log_scale$trimmed_log_y <- data_for_lm_log_scale$log_y
  }

  residual_model <- lm(trimmed_log_y ~ log_x1 + log_x2, 
                             data = data_for_lm_log_scale)
  
  residuals <- tibble(
    eid = feature_data$eid,
    sex = feature_data$sex,
    value = residual_model$residuals)
  
  standardized_residuals <- 
    residuals %>% 
    normalizeResiduales(.,use_standard_data = normalized)
  
  return(standardized_residuals)
}

#' Forcing residuals to follow a noraml distirbution by replacing their values with the
#' proper standard normal quantile
#' @param residuals - linear model residuals
#' @return standard normally distributed residuals

normalizeResiduales <- function(residuals_tibble, use_standard_data = F)
{
  residuals <- residuals_tibble$value
  
  residuals_location <- order(residuals)
  
  n <- length(residuals)
  
  standard_normal <- qnorm(residuals_location / (n+1))
  
  if(use_standard_data)
  {
    residuals_tibble$value <- standard_normal
  }
  else
  {
    residuals_tibble$value <- (residuals_tibble$value - mean(residuals_tibble$value)) / 
      sd(residuals_tibble$value)
  }
  return(residuals_tibble)
}
