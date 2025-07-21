#' Subject pronouns 
#'
#' @description
#' Usage of subject pronouns and its predictors.
#' 
#' @format
#' A data frame with 1024 observations and 7 columns
#' \describe{
#'   \item{real}{subject pronoun realized? Factor (2 levels "zero","realized")}
#'   \item{AGE}{age: Numerical (in months)}
#'   \item{LiBa}{linguistic background: Factor (3 levels "mono","multi", NA)}
#'   \item{ETH}{ethnic group: Factor (6 levels "C1a","C1b","C1c","C2a","C2b","C2c") (anonymized)}
#'   \item{SEX}{gender: Factor (2 levels "female","male")}
#'   \item{MLU}{mean length of utterance: Factor (4 levels "1","2","3","OL")}
#'   \item{PRN}{pronoun: Factor (7 levels "I","you_s","he","she","it","we","they")}
#' }
#' @source
#' Sarah Buschfeld, TU Dortmund
"data_zero"
#' Participants of subject pronoun study
#'
#' @description
#' Participants (random excerpt according to data_zero)
#'
#' @format
#' A dataset with 1024 observations and 1 column
#' \describe{
#'   \item{participant}{participant of study (P1 - P51)}
#' }
#' @source
#' Sarah Buschfeld, TU Dortmund
"participant_zero"
#' Subject pronouns and a predictor with one very frequent level
#'
#' @description
#' Usage of subject pronouns and its predictors; speaker level "adult" very frequent.
#' 
#' @format
#' A data frame with 3370 observations and 6 columns
#' \describe{
#'   \item{class}{subject pronoun realized? Factor (2 levels "zero","realized")}
#'   \item{AGE}{age: Numerical (in months)}
#'   \item{ETHN_GROUP}{ethnic group: Factor (3 levels "C","I","n_a") (anonymized)}
#'   \item{MLU}{mean length of utterance: Factor (5 levels "1","2","3","adult","OL")}
#'   \item{PRN_TYPE}{pronoun type: Factor (5 levels "dem","it_con","it_ex","it_ref","refer")}
#'   \item{SPEAKER}{speaker: Factor (2 levels "adult","child")}
#' }
#' @source
#' Sarah Buschfeld, TU Dortmund
"data_speaker"
#' Vowel length
#'
#' @description
#' Vowel length and categorization criteria.
#' 
#' @format
#' A data frame with 82 observations and 22 columns
#' \describe{
#'   \item{Nickname}{nickname: Factor (43 levels "Nick1","Nick2",...,"Nick43") (anonymized)}
#'   \item{LiBa}{linguistic background: Factor (2 levels "mono","multi")}
#'   \item{MLU}{mean length of utterance: Factor (3 levels "1","2","3")}
#'   \item{phone_label}{phone label: Factor (2 levels "fleece","kit")}
#'   \item{lexeme}{lexeme: Factor (14 levels "bee","cheek","cheese","chicken", ...)}
#'   \item{phone_left_1_duration}{duration of phone to the left of the vowel: Numerical (in msec)}
#'   \item{phone_right_1_duration}{duration of phone to the right of the vowel: Numerical (in msec)}
#'   \item{word_duration}{duration of word: Numerical (in msec)}
#'   \item{vowel_minimum_pitch}{minimum pitch of vowel: Numerical (in Hertz)}
#'   \item{vowel_maximum_pitch}{maximum pitch of vowel: Numerical (in Hertz)}
#'   \item{vowel_intensity_mean}{mean intensity of vowel: Numerical (in decibel)}
#'   \item{f1_fifty}{first formant F1 at midpoint of vowel (50\%): Numerical (in Hertz)}
#'   \item{f2_fifty}{second formant F2 at midpoint of vowel (50\%): numercial (in Hertz)}
#'   \item{target}{vowel length: Numerical (in msec)}
#'   \item{cons_class_l}{class of consonant to the left of the vowel: Factor (6 levels "l","r","tsh",...)}
#'   \item{cons_class_r}{class of consonant to the right of the vowel: Factor (7 levels "?"(glottal stop),"empty","nas",...)}
#'   \item{ETH}{ethnic group: Factor (6 levels "C1a","C1b","C1c","C2a","C2b","C2c") (anonymized)}
#'   \item{SEX}{gender: Factor (2 levels "female","male)}
#'   \item{AGE}{age: Numerical (in months)}
#'   \item{syllables}{number of syllables in lexeme: integer (1,2)}
#'   \item{speed}{speed of speech: Numerical (word duration / syllables; in msec)}
#'   \item{country}{country: Factor (2 levels "E","S")  (anonymized)}
#' } 
#' @source
#' Sarah Buschfeld, TU Dortmund
"data_vowel"
#' Landscape analysis
#'
#' @description
#' The use of language(s) on public signs and categorization criteria.
#' 
#' @format
#' A data frame with 149 observations and 28 columns
#' \describe{
#'   \item{coder}{who coded the sign: Factor (3 levels "C","E","S")  (anonymized)}
#'   \item{researcher}{who took a photograph of the sign: Factor (2 levels "L","S")  (anonymized)}
#'   \item{sign}{where the sign was found: Factor (11 levels "digi","door","graf",...)}
#'   \item{type.of.sign}{kind of sign: Factor (5 levels "com","commem","infra","reg","trans")}
#'   \item{permanent}{was the sign permanent? Factor (2 levels "no","yes")}
#'   \item{proper.noun}{kind of proper noun on the sign: Factor (9 levels "bn","bn+","cn",...)}
#'   \item{no.languages}{number of languages on sign: Factor (4 levels "1","2","3","4+")}
#'   \item{French}{French on sign? Factor (2 levels "0","1")}
#'   \item{Dutch}{Dutch on sign? Factor (2 levels "0","1")}
#'   \item{English}{English on sign? Factor (2 levels "0","1")}
#'   \item{Italian}{Italian on sign? Factor (2 levels "0","1")}
#'   \item{Spanish}{Spanish on sign? Factor (2 levels "0","1")}
#'   \item{German}{German on sign? Factor (2 levels "0","1")}
#'   \item{Indian,Mandarin,Greek,Indonesian,Portuguese,Libanese,Japanese,Russian,Danish,
#'             Norwegian,Hebrew,Catalan}{12 infrequent languages: Factor (2 levels "0","1")}
#'   \item{bn.unclear}{brand name unclear: Factor (2 levels "0","1")}
#'   \item{multilingual.type}{type of multilingualism on sign: Factor (6 levels "0","1","2","3","4","uc")}
#'   \item{location}{location of sign: Factor (2 levels "M","P")  (anonymized)}
#' }
#' @source
#' Sarah Buschfeld, TU Dortmund
"data_land"