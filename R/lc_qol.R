#' Calculates the domain-based scale scores using the data of QLQ-LC13.
#'
#' Calculates the domain-based scale scores using the data of QLQ-LC13
#'
#' @description Creates a dataset containing the domain-based scale scores using
#' the data from QLQ-LC13
#'
#' @details lc_miss function inputs either a dataset containing missing information, represented as,
#' 9 or 99 or NA or a data not containing any missing information. It extracts only the columns
#' named 'LC_Q31','LC_Q32',...,'LC_Q42' and replaces the missing data with the minimum value of the particular question.
#'
#' Using each of the 30 columns, the Raw Score is computed, and one column is obtained containing
#' the Raw Score for each patient.
#'
#' Further, using each of the Raw Scores, three domain-based Scale Scores are computed,
#' they are, Global Scales Score, Functional Scales Score and Symptoms Scales Score.
#'
#' Thus, the columns 'LC_Q31','LC_Q32',...,'LC_Q42' are replaced by the domain-based scale scores,
#' which is obtained as the output.
#'
#' lc_qol(x)
#'
#' 1) Subject ID column should be named as 'ID'.
#'
#' 2) Each question column should be named as 'LC_Q31' for data from question 31,
#' 'LC_Q32' for data from question 32, and so on until 'LC_Q42' for data from question 42.
#'
#' 3) Data may contain more variables, such as, Age, Gender, etc.
#'
#' x - A data frame with ID, LC_Q31,LC_Q32,...,LC_Q42 columns along with other columns if data
#' is available.
#'
#' rs - A matrix containing the Raw Score computed using all LC_Q31 to LC_Q42 data for each
#' patient. The RS(a) function is used in this case.
#'
#' ss - A matrix containing the Global Scale Scores computed using all LC_Q31 to LC_Q42
#' data for each patient. The SS(a,b) function is used in this case.
#'
#' final_data - A data frame formed by replacing the columns 'LC_Q31','LC_Q32',...,'LC_Q42' by
#' the domain-based scale scores.
#'
#' @param x A data frame with ID, LC_Q31,LC_Q32,...,LC_Q42 columns along with other columns if data is available.
#'
#' @import dplyr
#'
#' @return A data frame by replacing the columns 'LC_Q31','LC_Q32',...,'LC_Q42' by the domain-based scale scores.
#'
#' @references QoLMiss: Package for Repeatedly measured Quality of Life of Cancer Patients Data
#'
#' @examples
#' ##
#' data(lc_df)
#' lc_qol(lc_df)
#' data(lc_df_miss)
#' lc_qol(lc_df_miss)
#' ##
#'
#' @export
#' @author Atanu Bhattacharjee and Ankita Pal
#' @seealso https://github.com/apstat/QoLMiss-Package

lc_qol <- function(x){
  d <- as.matrix(select(x,'LC_Q31':'LC_Q42'))

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
  RS_data <- data.frame(RS_LCDY = RS(d[,3:5]),
                        RS_LCCO = d[,1],
                        RS_LCHA = d[,2],
                        RS_LCSM = d[,6],
                        RS_LCDS = d[,7],
                        RS_LCPN = d[,8],
                        RS_LCHR = d[,9],
                        RS_LCPC = d[,10],
                        RS_LCPA = d[,11],
                        RS_LCPO = d[,12])

  # Dataset with Score Values
  score_data <- data.frame(LCDY = SS(RS_data$RS_LCDY,d[,3:5]),
                           LCCO = SS(RS_data$RS_LCCO,d[,1]),
                           LCHA = SS(RS_data$RS_LCHA,d[,2]),
                           LCSM = SS(RS_data$RS_LCSM,d[,6]),
                           LCDS = SS(RS_data$RS_LCDS,d[,7]),
                           LCPN = SS(RS_data$RS_LCPN,d[,8]),
                           LCHR = SS(RS_data$RS_LCHR,d[,9]),
                           LCPC = SS(RS_data$RS_LCPC,d[,10]),
                           LCPA = SS(RS_data$RS_LCPA,d[,11]),
                           LCPO = SS(RS_data$RS_LCPO,d[,12]))
  new_data <- select(x,-('LC_Q31':'LC_Q42'))
  final_data <- data.frame(new_data,score_data)
  return(final_data)
}
