#' Calculates the domain-based scale scores using the data of QLQ-BR23
#' @description Creates a dataset containing the domain-based scale scores using
#' the data from QLQ-BR23
#'
#' @details brc_miss function inputs either a dataset containing missing information, represented as,
#' 9 or 99 or NA or a data not containing any missing information. It extracts only the columns
#' named 'BR_Q31','BR_Q32',...,'BR_Q53' and replaces the missing data with the minimum value of the particular question.
#'
#' Using each of the 30 columns, the Raw Score is computed, and one column is obtained containing
#' the Raw Score for each patient.
#'
#' Further, using each of the Raw Scores, three domain-based Scale Scores are computed,
#' they are, Global Scales Score, Functional Scales Score and Symptoms Scales Score.
#'
#' Thus, the columns 'BR_Q31','BR_Q32',...,'BR_Q53' are replaced by the domain-based scale scores,
#' which is obtained as the output.
#'
#' brc_qol(x)
#'
#' 1) Subject ID column should be named as 'ID'.
#'
#' 2) Each question column should be named as 'BR_Q31' for data from question 31,
#' 'BR_Q32' for data from question 32, and so on until 'BR_Q53' for data from question 53
#'
#' 3) Data may contain more variables, such as, Age, Gender, etc.
#'
#' x - A data frame with ID, BR_Q31,BR_Q32,...,BR_Q53 columns along with other columns if data
#' is available.
#'
#' rs - A matrix containing the Raw Score computed using all BR_Q31 to BR_Q53 data for each
#' patient. The RS(a) function is used in this case.
#'
#' fs - A matrix containing the Functional Scale Scores computed using all BR_Q31 to BR_Q53
#' data for each patient. The FS(a,b) function is used in this case.
#'
#' ss - A matrix containing the Global Scale Scores computed using all BR_Q31 to BR_Q53
#' data for each patient. The SS(a,b) function is used in this case.
#'
#' final_data - A data frame formed by replacing the columns 'BR_Q31','BR_Q32',...,'BR_Q53' by
#' the domain-based scale scores.
#'
#' @param x A data frame with ID, BR_Q31,BR_Q32,...,BR_Q53 columns along with other columns if data is available.
#'
#' @import dplyr
#' @import utils
#' @return A data frame by replacing the columns 'BR_Q31','BR_Q32',...,'BR_Q53' by the domain-based scale scores.
#'
#' @references QoLMiss: Package for Repeatedly measured Quality of Life of Cancer Patients Data
#'
#' @examples
#' ##
#' data(brc_df)
#' brc_qol(brc_df)
#' data(brc_df_miss)
#' brc_qol(brc_df_miss)
#' ##
#'
#' @export
#' @author Atanu Bhattacharjee and Ankita Pal
#' @seealso https://github.com/apstat/QoLMiss-Package

brc_qol <- function(x){
  d <- as.matrix(select(x,'BR_Q31':'BR_Q53'))

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

  # Functional Scales Score
  FS <- function(a,b){
    nr <- length(a)
    fs <- rep(0, nr)
    for(i in 1:nr){
      z <- (a[i]-1)/diff(range(b))
      fs[i] <- (1-z)*100
    }
    return(fs)
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
  RS_data <- data.frame(RS_BRBI = RS(d[,9:12]),
                        RS_BRSEF = RS(d[,c(14,15)]),
                        RS_BRSEE = d[,16],
                        RS_BRFU = d[,13],
                        RS_BRST = RS(d[,c(1,2,3,4,6,7,8)]),
                        RS_BRBS = RS(d[,20:23]),
                        RS_BRAS = RS(d[,17:19]),
                        RS_BRHL = d[,5])

  # Dataset with Score Values
  score_data <- data.frame(BRBI = FS(RS_data$RS_BRBI,d[,9:12]),
                           BRSEF = SS(RS_data$RS_BRSEF,d[,c(14,15)]),
                           BRSEE = SS(RS_data$RS_BRSEE,d[,16]),
                           BRFU = FS(RS_data$RS_BRFU,d[,13]),
                           BRST = SS(RS_data$RS_BRST,d[,c(1,2,3,4,6,7,8)]),
                           BRBS = SS(RS_data$RS_BRBS,d[,20:23]),
                           BRAS = SS(RS_data$RS_BRAS,d[,17:19]),
                           BRHL = SS(RS_data$RS_BRHL,d[,5]))
  new_data <- select(x,-('BR_Q31':'BR_Q53'))
  final_data <- data.frame(new_data,score_data)
  return(final_data)
}
utils::globalVariables(c("select"))
