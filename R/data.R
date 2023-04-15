.toWide <- function(data, RAM){

  manifestNames <- rownames(RAM@F)
  manifestNamesReduced <- unique(stringr::str_remove(string = manifestNames,
                                                     pattern = "_u[0-9]*"))

  dataWide <- data |>
    tidyr::pivot_wider(
      names_from = .data$occasion,
      values_from = tidyr::all_of(manifestNamesReduced),
      values_fill = NA,
      names_sep = "_u",
      names_prefix = ifelse(length(manifestNamesReduced) == 1,
                            paste0(manifestNamesReduced,"_u"),
                            "")
    )

  return(dataWide)
}
