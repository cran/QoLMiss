#' Calculates the domain-based scale scores using the data of QLQ-HN35
#'
#' Calculates the domain-based scale scores using the data of QLQ-HN35
#'
#' @description Creates a dataset containing the domain-based scale scores using
#' the data from QLQ-HN35
#'
#' @details hn_miss function inputs either a dataset containing missing information, represented as,
#' 9 or 99 or NA or a data not containing any missing information. It extracts only the columns
#' named 'HN_Q31','HN_Q32',...,'HN_Q65' and replaces the missing data with the minimum value of the particular question.
#'
#' Using each of the 30 columns, the Raw Score is computed, and one column is obtained containing
#' the Raw Score for each patient.
#'
#' Further, using each of the Raw Scores, three domain-based Scale Scores are computed,
#' they are, Global Scales Score, Functional Scales Score and Symptoms Scales Score.
#'
#' Thus, the columns 'HN_Q31','HN_Q32',...,'HN_Q65' are replaced by the domain-based scale scores,
#' which is obtained as the output.
#'
#' hnc_qol(x)
#'
#' 1) Subject ID column should be named as 'ID'.
#'
#' 2) Each question column should be named as 'HN_Q31' for data from question 31,
#' 'HN_Q32' for data from question 32, and so on until 'HN_Q65' for data from question 65.
#'
#' 3) Data may contain more variables, such as, Age, Gender, etc.
#'
#' x - A data frame with ID, HN_Q31,HN_Q32,...,HN_Q65 columns along with other columns if data
#' is available.
#'
#' rs - A matrix containing the Raw Score computed using all HN_Q31 to HN_Q65 data for each
#' patient. The RS(a) function is used in this case.
#'
#' ss - A matrix containing the Global Scale Scores computed using all HN_Q31 to HN_Q65
#' data for each patient. The SS(a,b) function is used in this case.
#'
#' final_data - A data frame formed by replacing the columns 'HN_Q31','HN_Q32',...,'HN_Q65' by
#' the domain-based scale scores.
#'
#' @param x A data frame with ID, HN_Q31,HN_Q32,...,HN_Q65 columns along with other columns if data is available.
#'
#' @import dplyr
#'
#' @return A data frame by replacing the columns 'HN_Q31','HN_Q32',...,'HN_Q65' by the domain-based scale scores.
#'
#' @references QoLMiss: Package for Repeatedly measured Quality of Life of Cancer Patients Data
#'
#' @examples
#' ##
#' data(hnc_df)
#' hnc_qol(hnc_df)
#' data(hnc_df_miss)
#' hnc_qol(hnc_df_miss)
#' ##
#'
#' @export
#' @author Atanu Bhattacharjee and Ankita Pal
#' @seealso https://github.com/apstat/QoLMiss-Package

hnc_qol <- function(x){
  d <- as.matrix(select(x,'HN_Q31':'HN_Q65'))

  # Imputing missing values with minimum value of respective question
  for(j in 1:ncol(d)){
    for(i in 1:nrow(d)){
      if(is.na(d[i,j])==TRUE || d[i,j]==9 || d[i,j]==99){
        d[i,j] <- min(d[,j],na.rm = TRUE)
      }
    }
  }

  # Raw Score
  RS <- function(a){
    nr <- nrow(a)
    rs <- rep(0, nr)
    for(i in 1:nr){
      rs[i] <- mean(a[i,])
    }
    return(rs)
  }

  # Symptoms Scales Score
  SS <- function(a,b){
    nr <- length(a)
    ss <- rep(0, nr)
    for(i in 1:nr){
      ss[i] <- ((a[i]-1)/diff(range(b)))*100
    }
    return(ss)
  }

  # Dataset with Raw Scores
  RS_data <- data.frame(RS_HNPA = RS(d[,1:4]),
                        RS_HNSW = RS(d[,5:8]),
                        RS_HNSE = RS(d[,c(13,14)]),
                        RS_HNSP = RS(d[,c(16,23,24)]),
                        RS_HNSO = RS(d[,19:22]),
                        RS_HNSC = RS(d[,c(18,25,26,27,28)]),
                        RS_HNSX = RS(d[,c(29,30)]),
                        RS_HNTE = d[,9],
                        RS_HNOM = d[,10],
                        RS_HNDR = d[,11],
                        RS_HNSS = d[,12],
                        RS_HNCO = d[,15],
                        RS_HNFI = d[,17],
                        RS_HNPK = d[,31],
                        RS_HNNU = d[,32],
                        RS_HNFE = d[,33],
                        RS_HNWL = d[,34],
                        RS_HNWG = d[,35])

  # Dataset with Score Values
  score_data <- data.frame(HNPA = SS(RS_data$RS_HNPA,d[,1:4]),
                           HNSW = SS(RS_data$RS_HNSW,d[,5:8]),
                           HNSE = SS(RS_data$RS_HNSE,d[,c(13,14)]),
                           HNSP = SS(RS_data$RS_HNSP,d[,c(16,23,24)]),
                           HNSO = SS(RS_data$RS_HNSO,d[,19:22]),
                           HNSC = SS(RS_data$RS_HNSC,d[,c(18,25,26,27,28)]),
                           HNSX = SS(RS_data$RS_HNSX,d[,c(29,30)]),
                           HNTE = SS(RS_data$RS_HNTE,d[,9]),
                           HNOM = SS(RS_data$RS_HNOM,d[,10]),
                           HNDR = SS(RS_data$RS_HNDR,d[,11]),
                           HNSS = SS(RS_data$RS_HNSS,d[,12]),
                           HNCO = SS(RS_data$RS_HNCO,d[,15]),
                           HNFI = SS(RS_data$RS_HNFI,d[,17]),
                           HNPK = SS(RS_data$RS_HNPK,d[,31]),
                           HNNU = SS(RS_data$RS_HNNU,d[,32]),
                           HNFE = SS(RS_data$RS_HNFE,d[,33]),
                           HNWL = SS(RS_data$RS_HNWL,d[,34]),
                           HNWG = SS(RS_data$RS_HNWG,d[,35]))
  new_data <- select(x,-('HN_Q31':'HN_Q65'))
  final_data <- data.frame(new_data,score_data)
  return(final_data)
}
