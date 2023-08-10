
#' Polarity index calculation
#'
#' This function calculates the polarity index based on organic modifier and organic modifier percentage
#' @param organic_modifier The organic modifier for which polarity index is calculated (either "MeOH" or "MeCN")
#' @param organic_percent The organic modifier percentage at a time of elution or the isocratic percentage value
#' @return Polarity index (units?)
#' @examples 
#' polarity_index_suspect <- polarity_index(organic_percent = 80, 
#'                                          organic_modifier = "MeOH");
#' @export
polarity_index = function(organic_percent = numeric(), 
                          organic_modifier = character()){
  polarity_index = case_when(
    organic_modifier == "MeCN" ~ (organic_percent/100)*5.1+((100-organic_percent)/100)*10.2,
    organic_modifier == "MeOH" ~ (organic_percent/100)*5.1+((100-organic_percent)/100)*10.2)
  return(polarity_index)
}