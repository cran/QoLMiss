#' Calculates the domain-based scale scores of Thyroid cancer using the data of QLQ-THY34
#'
#' @description Creates a dataset containing the domain-based scale scores using
#' the data from QLQ-THY34
#'
#' @details brc_miss function inputs either a dataset containing missing information, represented as,
#' 9 or 99 or NA or a data not containing any missing information. It extracts only the columns
#' named 'THY_Q31','THY_Q32',...,'THY_Q64' and replaces the missing data with the minimum value of the particular question.
#'
#' Using each of the 30 columns, the Raw Score is computed, and one column is obtained containing
#' the Raw Score for each patient.
#'
#' Further, using each of the Raw Scores, three domain-based Scale Scores are computed,
#' they are, Functional Scales Score and Symptoms Scales Score.
#'
#' Thus, the columns 'THY_Q31','THY_Q32',...,'THY_Q64' are replaced by the domain-based scale scores,
#' which is obtained as the output.
#'
#' thyc_qol(x)
#'
#' 1) Subject ID column should be named as 'ID'.
#'
#' 2) Each question column should be named as 'THY_Q31' for data from question 31,
#' 'THY_Q32' for data from question 32, and so on until 'THY_Q64' for data from question 64
#'
#' 3) Data may contain more variables, such as, Age, Gender, etc.
#'
#' x - A data frame with ID, THY_Q31,THY_Q32,...,THY_Q64 columns along with other columns if data
#' is available.
#'
#' rs - A matrix containing the Raw Score computed using all THY_Q31 to THY_Q64 data for each
#' patient. The RS(a) function is used in this case.
#'
#' ss - A matrix containing the Global Scale Scores computed using all THY_Q31 to THY_Q64
#' data for each patient. The SS(a,b) function is used in this case.
#'
#' final_data - A data frame formed by replacing the columns 'THY_Q31','THY_Q32',...,'THY_Q64' by
#' the domain-based scale scores.
#'
#' @param x A data frame with ID, THY_Q31,THY_Q32,...,THY_Q64 columns along with other columns if data is available.
#'
#' @import dplyr
#'
#' @return A data frame by replacing the columns 'THY_Q31','THY_Q32',...,'THY_Q64' by the domain-based scale scores.
#'
#' @references QoLMiss: Package for Repeatedly measured Quality of Life of Cancer Patients Data
#'
#' @examples
#' ##
#' data(thyc_df)
#' thyc_qol(thyc_df)
#' data(thyc_df_miss)
#' thyc_qol(thyc_df_miss)
#' ##
#'
#' @export
#' @author Atanu Bhattacharjee and Ankita Pal
#' @seealso https://github.com/apstat/QoLMiss-Package

thyc_qol <- function(x){
  d <- as.matrix(select(x,'THY_Q31':'THY_Q64'))

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
  RS_data <- data.frame(RS_THYFA = RS(d[,1:3]),
                        RS_THYDI = RS(d[,c(4,5,19)]),
                        RS_THYVO = RS(d[,6:8]),
                        RS_THYHA = RS(d[,9:10]),
                        RS_THYSW = RS(d[,11:12]),
                        RS_THYDM = d[,13],
                        RS_THYTO = d[,14],
                        RS_THYBI = d[,15],
                        RS_THYRE = RS(d[,16:17]),
                        RS_THYSH = d[,18],
                        RS_THYFE = RS(d[,20:22]),
                        RS_THYJP = d[,23],
                        RS_THYTI = RS(d[,24:25]),
                        RS_THYCR = d[,26],
                        RS_THYWO = RS(d[,27:30]),
                        RS_THYJE = d[,31],
                        RS_THYSO = RS(d[,32:34]))

  # Dataset with Score Values
  score_data <- data.frame(THYFA = SS(RS_data$RS_THYFA,d[,1:3]),
                           THYDI = SS(RS_data$RS_THYDI,d[,c(4,5,19)]),
                           THYVO = SS(RS_data$RS_THYVO,d[,6:8]),
                           THYHA = SS(RS_data$RS_THYHA,d[,9:10]),
                           THYSW = SS(RS_data$RS_THYSW,d[,11:12]),
                           THYDM = SS(RS_data$RS_THYDM,d[,13]),
                           THYTO = SS(RS_data$RS_THYTO,d[,14]),
                           THYBI = SS(RS_data$RS_THYBI,d[,15]),
                           THYRE = SS(RS_data$RS_THYRE,d[,16:17]),
                           THYSH = SS(RS_data$RS_THYSH,d[,18]),
                           THYFE = SS(RS_data$RS_THYFE,d[,20:22]),
                           THYJP = SS(RS_data$RS_THYJP,d[,23]),
                           THYTI = SS(RS_data$RS_THYTI,d[,24:25]),
                           THYCR = SS(RS_data$RS_THYCR,d[,26]),
                           THYWO = SS(RS_data$RS_THYWO,d[,27:30]),
                           THYJE = SS(RS_data$RS_THYJE,d[,31]),
                           THYSO = FS(RS_data$RS_THYSO,d[,32:34]))
  new_data <- select(x,-('THY_Q31':'THY_Q64'))
  final_data <- data.frame(new_data,score_data)
  return(final_data)
}
