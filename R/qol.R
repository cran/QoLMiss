#' Calculates the domain-based scale scores using the data from
#' Quality of Life questionnaire
#'
#' Calculates the domain-based scale scores using the data from
#' Quality of Life questionnaire
#'
#' @description Creates a dataset containing the domain-based scale scores using
#' the data from Quality of Life questionnaire
#'
#' @details qol function inputs either a dataset containing missing information, represented as,
#' 9 or 99 or NA or a data not containing any missing information. It extracts only the columns
#' named 'Q1','Q2',...,'Q30' and replaces the missing data with the minimum value of the particular question.
#'
#' Using each of the 30 columns, the Raw Score is computed, and one column is obtained containing
#' the Raw Score for each patient.
#'
#' Further, using each of the Raw Scores, three domain-based Scale Scores are computed,
#' they are, Global Scales Score, Functional Scales Score and Symptoms Scales Score.
#'
#' Thus, the columns 'Q1','Q2',...,'Q30' are replaced by the domain-based scale scores,
#' which is obtained as the output.
#'
#' qol(x)
#'
#' 1) Subject ID column should be named as 'ID'.
#'
#' 2) Each question column should be named as 'Q1' for data from question 1,
#' 'Q2' for data from question 2, and so on until 'Q30' for data from question 30.
#'
#' 3) Data may contain more variables, such as, Age, Gender, etc.
#'
#' x - A data frame with ID, Q1, Q2,..., Q30 columns along with other columns if data
#' is available.
#'
#' rs - A matrix containing the Raw Score computed using all Q1 to Q30 data for each
#' patient. The RS(a) function is used in this case.
#'
#' fs - A matrix containing the Functional Scale Scores computed using all Q1 to Q30
#' data for each patient. The FS(a,b) function is used in this case.
#'
#' ss_gs - A matrix containing the Global Scale Scores computed using all Q1 to Q30
#' data for each patient. The SS_GS(a,b) function is used in this case.
#'
#' final_data - A data frame formed by replacing the columns 'Q1','Q2',...,'Q30' by
#' the domain-based scale scores.
#'
#' @param x A data frame with ID, Q1, Q2,..., Q30 columns along with other columns if data is available.
#'
#' @import dplyr
#'
#' @return A data frame by replacing the columns 'Q1','Q2',...,'Q30' by the domain-based scale scores.
#'
#' @references QoLMiss: Package for Repeatedly measured Quality of Life of Cancer Patients Data
#'
#' @examples
#' ##
#' data(c30_df)
#' qol(c30_df)
#' data(c30_df_miss)
#' qol(c30_df_miss)
#' ##
#'
#' @export
#' @author Atanu Bhattacharjee and Ankita Pal
#' @seealso https://github.com/apstat/QoLMiss-Package

qol <- function(x){
  d <- as.matrix(select(x,'Q1':'Q30'))

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
  SS_GS <- function(a,b){
    nr <- length(a)
    ss_gs <- rep(0, nr)
    for(i in 1:nr){
      ss_gs[i] <- ((a[i]-1)/diff(range(b)))*100
    }
    return(ss_gs)
  }

  # Dataset with Raw Scores
  RS_data <- data.frame(RS_QL = RS(d[,c(29,30)]),
                        RS_PF = RS(d[,1:5]),
                        RS_RF = RS(d[,c(6,7)]),
                        RS_EF = RS(d[,21:24]),
                        RS_CF = RS(d[,c(20,25)]),
                        RS_SF = RS(d[,c(26,27)]),
                        RS_FA = RS(d[,c(10,12,18)]),
                        RS_NV = RS(d[,c(14,15)]),
                        RS_PA = RS(d[,c(9,19)]),
                        RS_DY = d[,8],
                        RS_SL = d[,11],
                        RS_AP = d[,13],
                        RS_CO = d[,16],
                        RS_DI = d[,17],
                        RS_FI = d[,28])

  # Dataset with Score Values
  score_data <- data.frame(QL = SS_GS(RS_data$RS_QL,d[,c(29,30)]),
                           PF = FS(RS_data$RS_PF,d[,1:5]),
                           RF = FS(RS_data$RS_RF,d[,c(6,7)]),
                           EF = FS(RS_data$RS_EF,d[,21:24]),
                           CF = FS(RS_data$RS_CF,d[,c(20,25)]),
                           SF = FS(RS_data$RS_SF,d[,c(26,27)]),
                           FA = SS_GS(RS_data$RS_FA,d[,c(10,12,18)]),
                           NV = SS_GS(RS_data$RS_NV,d[,c(14,15)]),
                           PA = SS_GS(RS_data$RS_PA,d[,c(9,19)]),
                           DY = SS_GS(RS_data$RS_DY,d[,8]),
                           SL = SS_GS(RS_data$RS_SL,d[,11]),
                           AP = SS_GS(RS_data$RS_AP,d[,13]),
                           CO = SS_GS(RS_data$RS_CO,d[,16]),
                           DI = SS_GS(RS_data$RS_DI,d[,17]),
                           FI = SS_GS(RS_data$RS_FI,d[,28]))
  new_data <- select(x,-('Q1':'Q30'))
  final_data <- data.frame(new_data,score_data)
  return(final_data)
}
