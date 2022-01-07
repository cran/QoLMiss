#' Calculates the domain-based scale scores using the data of QLQ-OV28.
#'
#' Calculates the domain-based scale scores using the data of QLQ-OV28
#'
#' @description Creates a dataset containing the domain-based scale scores using
#' the data from QLQ-OV28
#'
#' @details brc_miss function inputs either a dataset containing missing information, represented as,
#' 9 or 99 or NA or a data not containing any missing information. It extracts only the columns
#' named 'OV_Q31','OV_Q32',...,'OV_Q58' and replaces the missing data with the minimum value of the particular question.
#'
#' Using each of the 30 columns, the Raw Score is computed, and one column is obtained containing
#' the Raw Score for each patient.
#'
#' Further, using each of the Raw Scores, three domain-based Scale Scores are computed,
#' they are, Global Scales Score, Functional Scales Score and Symptoms Scales Score.
#'
#' Thus, the columns 'OV_Q31','OV_Q32',...,'OV_Q58' are replaced by the domain-based scale scores,
#' which is obtained as the output.
#'
#' ovc_qol(x)
#'
#' 1) Subject ID column should be named as 'ID'.
#'
#' 2) Each question column should be named as 'OV_Q31' for data from question 31,
#' 'OV_Q32' for data from question 32, and so on until 'OV_Q58' for data from question 58
#'
#' 3) Data may contain more variables, such as, Age, Gender, etc.
#'
#' x - A data frame with ID, OV_Q31,OV_Q32,...,OV_Q58 columns along with other columns if data
#' is available.
#'
#' rs - A matrix containing the Raw Score computed using all OV_Q31 to OV_Q58 data for each
#' patient. The RS(a) function is used in this case.
#'
#' ss - A matrix containing the Global Scale Scores computed using all OV_Q31 to OV_Q58
#' data for each patient. The SS(a,b) function is used in this case.
#'
#' final_data - A data frame formed by replacing the columns 'OV_Q31','OV_Q32',...,'OV_Q58' by
#' the domain-based scale scores.
#'
#' @param x A data frame with ID, OV_Q31,OV_Q32,...,OV_Q58 columns along with other columns if data is available.
#'
#' @import dplyr
#'
#' @return A data frame by replacing the columns 'OV_Q31','OV_Q32',...,'OV_Q58' by the domain-based scale scores.
#'
#' @references QoLMiss: Package for Repeatedly measured Quality of Life of Cancer Patients Data
#'
#' @examples
#' ##
#' data(ovc_df)
#' ovc_qol(ovc_df)
#' data(ovc_df_miss)
#' ovc_qol(ovc_df_miss)
#' ##
#'
#' @export
#' @author Atanu Bhattacharjee and Ankita Pal
#' @seealso https://github.com/apstat/QoLMiss-Package

ovc_qol <- function(x){
  d <- as.matrix(select(x,'OV_Q31':'OV_Q58'))

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
  RS_data <- data.frame(RS_GI = RS(d[,1:6]),
                        RS_PN = RS(d[,c(11,12)]),
                        RS_HOR = RS(d[,c(18,19)]),
                        RS_BI = RS(d[,c(20,21)]),
                        RS_AD = RS(d[,c(22,23,24)]),
                        RS_CSE = RS(d[,13:17]),
                        RS_SI = RS(d[,7:10]),
                        RS_SX = RS(d[,25:28]))

  # Dataset with Score Values
  score_data <- data.frame(Abdominal_GI = SS(RS_data$RS_GI,d[,1:6]),
                           Peripheral_Neuropathy = SS(RS_data$RS_PN,d[,c(11,12)]),
                           Hormonal = SS(RS_data$RS_HOR,d[,c(18,19)]),
                           Body_Image = SS(RS_data$RS_BI,d[,c(20,21)]),
                           Attitude_to_Disease = SS(RS_data$RS_AD,d[,c(22,23,24)]),
                           Chemotherapy_side_effects = SS(RS_data$RS_CSE,d[,13:17]),
                           Other_Single_Items = SS(RS_data$RS_SI,d[,7:10]),
                           Sexuality = SS(RS_data$RS_SX,d[,25:28]))
  new_data <- select(x,-('OV_Q31':'OV_Q58'))
  final_data <- data.frame(new_data,score_data)
  return(final_data)
}
