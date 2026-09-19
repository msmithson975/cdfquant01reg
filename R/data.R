#' Yoon Data Set
#'
#' The data are records of 894 patients whose times were recorded for five stages of ED assessment and treatment from a study consisting of recorded information on patients admitted to the emergency department of the University of Alberta Hospital between midnight January 23 and midnight January 29, 1999 (Yoon, Steiner, and Reinhardt 2003).
#'
#' @format A data frame with these columns:
#' \describe{
#'   \item{pphysdecis}{The proportion of total time spent by each patient in the ED under assessment and treatment by a physician}
#'   \item{Ambulance}{A binary variable indicating how each patient arrived at the ED, 0 = walk-in, 1 = ambulance arrival}
#' }
#' @source Yoon, P., I. Steiner, and G. Reinhardt. 2003. Analysis of factors influencing length of stay in the emergency department. CJEM 5 (3):155–61.
"yoon"

#' Blame Data Set
#'
#' The data for this example are from a replication and extension of Halevy et al. (2022), the replication conducted by Smithson.  An Israeli adult sample was recruited using the Qualtrics panel platform during April 2021 (N = 462; mean age = 39.9, age SD = 15.5; 245 females, 214 males, 3 other).  The participants were asked to assign percentages of blame to Israelis and Palestinians for the unresolved Israeli-Palestinian conflict.  Two experimental conditions varied whether the Palestinians were listed as one population or as three separate movements, i.e., listed as Palestinians vs listed as The Palestinian Hamas movement, The Palestinian Fatah movement, The Palestinian Islamic Jihad movement.  The hypothesis to be tested here was that unpacking the Palestinians into three groups would increase the percentage of blame that the Israeli participants assigned to Palestinians.
#'
#' @format A data frame with these columns:
#' \describe{
#'   \item{ptotpalest1}{The proportion of blame assigned to Palestinians}
#'   \item{palestpk}{A binary variable indicating whether the Palestinians were listed as three groups or one, 0 = unpacked, 1 = packed}
#' }
#' @source Halevy, N., Maoz, I., Vani, P. and Reit, E.S., 2022. Where the blame lies: Unpacking groups into their constituent subgroups shifts judgments of blame in intergroup conflict. Psychological Science, 33(1), pp.76-89.
"blame"

#' Gun Ownership Data Set
#'
#' The data for this example are from a 2015 survey of 321 American voters conducted by Smithson.  Participants were asked to assign a number representing their feeling about gun ownership from 0 (very negative) to 100 (very positive).  They were also asked to choose their political affiliation from four categories: Democrat, Republican, independent, and no preference.
#'
#' @format A data frame with these columns:
#' \describe{
#'   \item{poliorient}{Political affiliation by assignment to one of four categories}
#'   \item{demo}{A binary variable, coded 1 = Democrat, 0 = any other category}
#'   \item{indep}{A binary variable, coded 1 = independent, 0 = any other category}
#'   \item{nopref}{A binary variable, coded 1 = no preference, 0 = any other category}
#'   \item{repub}{A binary variable, coded 1 = Republican, 0 = any other category}
#'   \item{gun01}{A numeric variable, linearly transformed from the 0-100 range to 0-1}
#' }
#' @source Smithson, M. 2015. Unpublished dataset.
"gunown"
