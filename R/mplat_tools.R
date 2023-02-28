# mplat tools

#' Change mplat to a long format dataframe
#'
#' @param mplat A matrix
#' @return A dataframe with columns lexeme, morphosite, allele
pivot_mplat_longer = function(mplat) {
  
  colnames(mplat) <- str_c("c", 1:ncol(mplat))
  as.data.frame(mplat) %>%
    mutate("lexeme" = 1:nrow(mplat)) %>%
    pivot_longer(
      cols = starts_with("c"),
      names_to = "morphosite",
      names_prefix = "c",
      values_to = "allele"
    ) %>%
    mutate(
      morphosite = as.integer(morphosite),
      allele = as.integer(allele)
      ) %>%
    arrange(morphosite, lexeme)
}


#' Unpack evolution
#'
#' @param evolution A list, the output from evolve_mplat.
#' @param repetition An integer. Which repetition to unpack.
#' @param step_method A string, specifying the kind of steps to take through the
#'   evolution: by generation or by change.
#' @param skip_steps An integer, specifying how many steps to skip ahead each
#'   time.
#' @return A dataframe of the evolution in long form.
unpack_evolution = function(
  evolution,
  repetition = 1, 
  step_method = c("generation", "change"),
  skip_steps = 10
) {
  
  mplat <- evolution$mplat_0
  changes <- evolution$changes
  steps <- get_steps(evolution, repetition, step_method, skip_steps)
  
  # Initialise evo_df
  evo_df <- long_mplat <-
    pivot_mplat_longer(mplat) %>%
    mutate(generation = 0, step = 1)
  
  ## Cycle through all steps
  
  for(i in 2:length(steps$rows)) {
    
    this_row <- steps$rows[i]
    prev_row <- steps$rows[i - 1]
    
    if (this_row > prev_row) {
      # mplat needs updating via all intervening changes:
      chg <- changes[(prev_row + 1):this_row, ]
      long_mplat <- update_long_mplat(long_mplat, chg)
    }
    
    # Add the current mplat to evo_df
    new_rows <- 
      long_mplat %>%
      mutate(generation = steps$gens[i], step = i)
    evo_df <- bind_rows(evo_df, new_rows)
  }
  
  evo_df %>% select(step, generation, lexeme, morphosite, allele)
}


#' Get step rows and generations
#' 
#' @param evolution A list, the output from evolve_mplat.
#' @param repetition An integer. Which repetition to unpack.
#' @param step_method A string, specifying the kind of steps to take through the
#'   evolution: by generation or by change.
#' @param skip_steps An integer, specifying how many steps to skip ahead each
#'   time.
#' @return A list.
get_steps = function(
  evolution,
  repetition = 1,
  step_method = c("generation", "change"),
  skip_steps
) {
  
  s_method <- substr(step_method[1], 1, 1)
  changes <- evolution$changes %>% filter(repetition == !!repetition)
  change_gens <- changes$generation
  final_chg <- nrow(changes)
  final_gen <- change_gens[final_chg]
  
  if (s_method == "g") {
    # For plotting evenly spaced generations
    step_gens <- seq(1, final_gen, by = skip_steps)
    step_rows <- 
      step_gens %>% 
      sapply(function(g) { c(rev(which(change_gens <= g)), 0)[1] })
    
  } else if (s_method == "c") {
    # For plotting evenly spaced changes
    step_rows <- seq(1, final_chg, by = skip_steps)
    step_gens <- changes$generation[step_rows]
  }
  
  list(
    rows = c(0, step_rows, final_chg),
    gens = c(0, step_gens, final_gen)
  )
}


#' Update mplat by one or more changes
#' @param mplat A matrix
#' @param changes A dataframe
#' @return A list, the mplat updated by the changes; a vector of changed rows
#'   and a vector of changed columns
update_mplat = function(
  mplat,
  changes
) {
  
  chg_rows <- NULL
  chg_cols <- NULL
  
  for(i in 1:nrow(changes)) {
    change <- changes[i, ]
    lex <- str2ints(change$foci_lex)
    morphosite <- str2ints(change$foci_morphosite)
    mplat[lex, morphosite] <- str2ints(change$a_new)
    chg_rows <- c(chg_rows, lex)
    chg_cols <- c(chg_cols, morphosite)
  }
  
  list(
    mplat = mplat,
    rows = unique(chg_rows),
    cols = unique(chg_cols)
  )
}


#' Update longform mplat by one or more changes
#' @param mplat A dataframe
#' @param changes A dataframe
#' @return A dataframe, the mplat updated by the changes
update_long_mplat = function(
  mplat,
  changes
) {
  nr <- nrow(mplat)
  nlex <- mplat$lexeme[nr]
  for(i in 1:nrow(changes)) {
    change <- changes[i, ]
    lex <- str2ints(change$foci_lex)
    morphosite <- str2ints(change$foci_morphosite)
    chg_row <- (morphosite - 1) * nlex + lex
    mplat$allele[chg_row] <- str2ints(change$a_new)
  }
  mplat
}