#' Dataset contains survival outcomes and quality of life for cancer patients with missing observation
#'
#' Calculates the domain-wise relative hazard ratio (95% CI) using the arm-wise data from QLQ-C30
#'
#' @description Creates a dataset containing the domain-based relative hazard ratio (95% CI) using
#' the arm-wise data from QLQ-C30
#'
#' @details surv_c30_miss function inputs a dataset where information of some patients
#' are completely missing, that is, some rows contain only NA. It passes the data to the qol_miss()
#' function, which in turn gives the domain-wise scale scores. These domain-wise scale scores are used
#' for calculating the relative hazard ratio (95% CI) by first computing the hazard ratios by dividing
#' the data arm-wise.
#'
#' The surv_c30_miss function includes the qol_miss() function which will consider the arm-wise data and
#' calculate the domain-wise scale scores. Hence, two set of domain-wise scale scores will be obtained,
#' one for each arm.
#'
#' Each of the domain-wise scales, 'QL','PF','RF','EF','CF','SF','FA','NV','PA','DY','SL','AP','CO','DI','FI', are considered
#' as the covariates. Using these columns, Cox-Proportional model will be used for univariate analysis for
#' each of the covariates. The hazard ratio (95% CI) obtained for each arm is used to find out the relative hazard ratio (95% CI).
#'
#' Thus, the output will contain three columns, Hazard Ratio(HR), Lower 95% CI and Upper 95% CI, for each of the covariates.
#'
#' surv_c30_miss(x)
#'
#' 1) Subject ID column should be named as 'ID'.
#'
#' 2) Each question column should be named as 'Q1' for data from question 1,'Q2' for data from question 2, and so on until 'Q30' for data from question 30.
#'
#' 3) Only those data can be used which contains no information for some patients, that is, some rows contain only NA.
#'
#' 4) Data must contain columns for 'time', 'event' and 'arm'.
#'
#' 5) Data may contain more variables, such as, Age, Gender, etc.
#'
#' x - A data frame with ID, time, event, arm, Q1,Q2,...,Q30 columns along with other columns if data is available.
#'
#' @param x A data frame with ID, time, event, arm, Q1,Q2,...,Q30 columns along with other columns if data is available.
#'
#' @import dplyr
#' @import survival
#' @return A data frame containing the Hazard Ratio(HR), Lower 95% CI and Upper 95% CI, for each of the covariates.
#'
#' @references QoLMiss: Package for Repeatedly measured Quality of Life of Cancer Patients Data
#'
#' @examples
#' ##
#' data(patient_miss)
#' surv_c30_miss(patient_miss)
#' ##
#'
#' @export
#' @author Atanu Bhattacharjee and Ankita Pal
#' @seealso https://github.com/apstat/QoLMiss-Package

surv_c30_miss <- function(x){
  covariates <- c('QL','PF','RF','EF','CF','SF','FA','NV','PA','DY','SL','AP','CO','DI','FI')

  arm1_qol <- qol_miss(x[x$arm==1,])
  arm2_qol <- qol_miss(x[x$arm==2,])

  univ_formulas <- sapply(covariates,function(y) as.formula(paste('Surv(time, event)~', y)))
  univ_arm1 <- lapply(univ_formulas,function(y){coxph(y, data = arm1_qol)})
  univ_arm2 <- lapply(univ_formulas,function(y){coxph(y, data = arm2_qol)})
  univ_results <- mapply(function(x1,x2){
    x1 <- summary(x1)
    x2 <- summary(x2)
    HR <- signif(x2$coef[2]/x1$coef[2], digits=3);
    HR.confint.lower <- signif(x2$conf.int[,"lower .95"]/x1$conf.int[,"lower .95"], 3)
    HR.confint.upper <- signif(x2$conf.int[,"upper .95"]/x1$conf.int[,"upper .95"],3)
    res<-c(HR,HR.confint.lower,HR.confint.upper)
    names(res) <- c('HR','Lower 95% CI','Upper 95% CI')
    return(res)},
    univ_arm1,univ_arm2)
  relative.HR <- t(as.data.frame(univ_results, check.names = FALSE))
  return(relative.HR)
}
