#' Dataset contains survival outcomes and quality of life for head and neck cancer patients
#'
#' Calculates the domain-wise relative hazard ratio (95% CI) using the arm-wise data from QLQ-HN35
#'
#' @description Creates a dataset containing the domain-based relative hazard ratio (95% CI) using
#' the arm-wise data from QLQ-HN35
#'
#' @details surv_hn35 function inputs either a dataset containing missing information, represented as,
#' 9 or 99 or NA or a data not containing any missing information. It passes the data to the hnc_qol()
#' function, which in turn gives the domain-wise scale scores. These domain-wise scale scores are used
#' for calculating the relative hazard ratio (95% CI) by first computing the hazard ratios by dividing
#' the data arm-wise.
#'
#' The surv_hn35 function includes the hnc_qol() function which will consider the arm-wise data and
#' calculate the domain-wise scale scores. Hence, two set of domain-wise scale scores will be obtained,
#' one for each arm.
#'
#' Each of the domain-wise scales are considered as the covariates. Using these
#' columns, Cox-Proportional model will be used for univariate analysis for each of the covariates. The
#' hazard ratio (95% CI) obtained for each arm is used to find out the relative hazard ratio (95% CI).
#'
#' Thus, the output will contain three columns, Hazard Ratio(HR), Lower 95% CI and Upper 95% CI, for each of the covariates.
#'
#' surv_hn35(x)
#'
#' 1) Subject ID column should be named as 'ID'.
#'
#' 2) Each question column should be named as 'HN_Q31' for data from question 31, HN_Q32' for data from question 32, and so on until 'HN_Q65' for data from question 65.
#'
#' 3) Data must contain columns for 'time', 'event' and 'arm'.
#'
#' 4) Data may contain more variables, such as, Age, Gender, etc.
#'
#' x - A data frame with ID, time, event, arm, HN_Q31,HN_Q32,...,HN_Q65 columns along with other columns if data is available.
#'
#' @param x A data frame with ID, time, event, arm, HN_Q31,HN_Q32,...,HN_Q65 columns along with other columns if data is available.
#'
#' @import utils
#' @import dplyr
#' @import survival
#' @return A data frame containing the Hazard Ratio(HR), Lower 95% CI and Upper 95% CI, for each of the covariates.
#'
#' @references QoLMiss: Package for Repeatedly measured Quality of Life of Cancer Patients Data
#'
#' @examples
#' ##
#' data(hnc_df)
#' surv_hn35(hnc_df)
#' ##
#'
#' @export
#' @author Atanu Bhattacharjee and Ankita Pal
#' @seealso https://github.com/apstat/QoLMiss-Package

surv_hn35 <- function(x){
  covariates <- c('HNPA','HNSW','HNSE','HNSP','HNSO','HNSC','HNSX','HNTE','HNOM','HNDR','HNSS','HNCO','HNFI','HNPK','HNNU','HNFE','HNWL','HNWG')

  arm1_qol <- hnc_qol(x[x$arm==1,])
  arm2_qol <- hnc_qol(x[x$arm==2,])

  univ_formulas <- sapply(covariates,function(y) as.formula(paste('Surv(time, event)~', y)))
  univ_arm1 <- lapply(univ_formulas,function(y){coxph(y, data = arm1_qol)})
  univ_arm2 <- lapply(univ_formulas,function(y){coxph(y, data = arm2_qol)})
  univ_results <- mapply(function(x1,x2){
    x1 <- summary(x1)
    x2 <- summary(x2)
    HR <- signif(x2$coef[2]/x1$coef[2], digits=3);
    HR.confint.lower <- signif(x2$conf.int[,"lower .95"]/x1$conf.int[,"lower .95"], 3)
    HR.confint.upper <- signif(x2$conf.int[,"upper .95"]/x1$conf.int[,"upper .95"],3)
    res <- c(HR,HR.confint.lower,HR.confint.upper)
    names(res) <- c('HR','Lower 95% CI','Upper 95% CI')
    return(res)},
    univ_arm1,univ_arm2)
  relative.HR <- t(as.data.frame(univ_results, check.names = FALSE))
  return(relative.HR)
}
utils::globalVariables(c("as.data.frame","as.formula"))
