# Evolution


#' Simulate mplat evolution
#'
#' @param x Either a list or a matrix. If a matrix, x is interpreted as an mplat
#'   representing the initial state of the inflectional system; if a list, then
#'   a set of parameters to generate that state using initialise_mplat().
#' @param pivot_type A string, giving the model type: morphosite pivot (as in
#'   A&M) or lexeme pivot (as in Esher 2015)
#' @param foci_lex_selector A string, giving the method for selecting the foci
#'   lexemes, which change.
#' @param foci_morphosite_selector A string, giving the method for selecting the
#'   foci morphosites, which change.
#' @param n_foci An integer, giving the number of foci to change in a single
#'   generation.
#' @param pivot_weighting A string, giving the method for weighting pivots; used
#'   either for sampling them, or weighting them when all are taken into
#'   account.
#' @param pivot_use_proportion An numeric, giving the proportion of the
#'   available pivots to use. If set to 0, then just one pivot is used.
#' @param evidence_weighting A string, giving the method for weighting evidence
#'   morphosites; used either for sampling them, or weighting them when all are
#'   taken into account.
#' @param evidence_use_proportion An numeric, giving the proportion of the
#'   available evidence morphosites to use. If set to 0, just one is used.
#' @param evidence_balance A numeric, between 0 and 1, giving the balance of
#'   negative evidence (max when 0) and positive evidence (max when 1).
#' @param neg_evidence_coeff A non-negative numeric, an alternative
#'   parameterisation of evidence balance, giving negative evidence strength as
#'   a multiple of positive evidence.
#' @param change_classwise A logical, whether to change whole classes at once.
#' @param halt_method A string, giving the method for choosing when to halt the
#'   simulation.
#' @param halt_n An integer, indicating either the total number of generations
#'   to evolve through or the number of consecutive generations of
#'   quasi-stability before halting.
#' @param n_rep An integer, the number of repetitions of the simulation to
#'   perform
#' @return A list, containing: containing: mplat_0; total_generations; a list of
#'   the model parameters; a list of final mplats; a dataframe of changes (with
#'   columns generation, foci_lex, foci_morphosite, a_new); a dataframe of
#'   statistics.

evolve_mplat = function(
  x = NULL,
  param_dict = mplat_params(),
  pivot_type = c("morphosite", "lexeme"),
  foci_lex_selector = c("uniform", "zipf"),
  foci_morphosite_selector = c("uniform", "zipf"),
  n_foci = 1,
  pivot_weighting = c("uniform", "zipf", "classwise"),
  pivot_use_proportion = 0,
  evidence_weighting = c("uniform", "zipf", "classwise"),
  evidence_use_proportion = 1,
  evidence_balance = 1 / (neg_evidence_coeff + 1),
  neg_evidence_coeff = 0,
  change_classwise = FALSE,
  halt_method = c("n_generations", "quasi_stability"),
  halt_n = 5000,
  n_rep = 1
) {

    
  ## Initialise some values
  
  # Get mplat_0
  if (is.list(x)) {
    param_dict <- x
    mplat_0 <- initialise_mplat(x)
  } else {
    param_dict <- NULL
    mplat_0 <- x
  }
  
  # Parameters
  halt_method <- substr(halt_method[1], 1, 1)
  pvt_type <- substr(pivot_type[1], 1, 1)
  is_rotated <- pvt_type == "l"
  foc_lex_sel  <- substr(foci_lex_selector[1], 1, 1)
  foc_morphosite_sel  <- substr(foci_morphosite_selector[1], 1, 1)
  
  # Label the 'evidence axis' e, and the 'base axis' b:
  nr <- nrow(mplat_0)
  nc <- ncol(mplat_0)
  e_length  <- if (is_rotated) nc else nr # Size on the evidence axis.
  b_length  <- if (is_rotated) nr else nc # Size on the base axis.
  foc_e_sel <- if (is_rotated) foc_morphosite_sel else foc_lex_sel
  foc_b_sel <- if (is_rotated) foc_lex_sel else foc_morphosite_sel

  # Evolution data records
  mplat_finals <- vector("list", n_rep)
  changes_df <- data.frame(
    stringsAsFactors = FALSE,
    repetition = integer(0), generation = integer(0), 
    foci_lex = character(0), foci_morphosite = character(0), 
    a_new = character(0)
    )

  ## Start simulations
  
  for (rep_i in 1:n_rep) {

    ## Initialise
    
    mplat <- if (is_rotated) t(mplat_0) else mplat_0 # mplat, rotated if required
    generation <- 0
    n_gen_stability <- 0
    is_halt_generation <- FALSE
  
    ## Evolve
    
    while (!is_halt_generation) {
      
      ## Focus location(s)
      
      # For weighting, use inverse of zipfian distribution, so less common
      # items are more likely to change
      wt_b <- 1 / get_weights(foc_b_sel, b_length)
      wt_e <- 1 / get_weights(foc_e_sel, e_length)
      foci_loc_b <- sample(b_length, size = n_foci, prob = wt_b)
      focus_loc_e  <- sample(e_length, size = 1, prob = wt_e)
  
      # If changes are classwise, selection of the focus should also be
      # classwise on the e axis:
      if (change_classwise) {
        class_rows <- which(!duplicated(mplat)) # rows the represent a class
        focus_loc_e <- class_rows[sample(length(class_rows), size = 1)]
      }
      
      
      ## Select the new alleles
  
      a_new <- 
        select_replacement_alleles(
          mplat, 
          foci_loc_b, 
          focus_loc_e, 
          pivot_weighting, 
          pivot_use_proportion, 
          evidence_weighting, 
          evidence_use_proportion, 
          evidence_balance
        )
      
      
      ## Update
      
      generation <- generation + 1
      a_old <- mplat[focus_loc_e, foci_loc_b]
      is_changed <- any(a_old != a_new)
      
      if (!is_changed) {
        
        # No change in this step
        
        n_gen_stability <- n_gen_stability + 1
        
      } else {
        
        # Change in this step
        
        n_gen_stability <- 0
        
        if (change_classwise) {
          # Change all rows that are identical to the focus row
          chg_loc_e <- match_rows(mplat, mplat[focus_loc_e, ])
          n_chg_loc <- length(chg_loc_e)
          a_new <- rep(a_new, each = n_chg_loc)
        } else {
          chg_loc_e <- focus_loc_e
        }
        
        # Update mplat
        mplat[chg_loc_e, foci_loc_b] <- a_new
        
        # Update changes dataframe
        foci_str <- 
          c(paste(chg_loc_e, collapse = ","),
            paste(foci_loc_b, collapse = ","))
        changes_df <-
          bind_rows(
            changes_df,
            data.frame(
              stringsAsFactors = FALSE,
              repetition = rep_i,
              generation, 
              foci_lex = foci_str[1 + is_rotated], 
              foci_morphosite = foci_str[1 + !is_rotated],
              a_new = paste(a_new, collapse = ","))
          )
      }
      
      ## Check whether to stop
      
      is_halt_generation <-
        (halt_method == "n" & generation == halt_n) | 
        (halt_method == "q" & n_gen_stability == halt_n) 
      
    }
    
    # Add a row for the final generation
    changes_df[nrow(changes_df) + 1, ] <-
        changes_df[nrow(changes_df), ] %>% 
          mutate(generation = !!generation)
    
    # Add final mplat to list
    mplat_finals[[rep_i]] <- if (is_rotated) t(mplat) else mplat
  }
  
  ## Return results
  list(
    mplat_0 = mplat_0,
    total_generations = generation,
    param_dict = param_dict,
    n_rep = n_rep,
    model = list(
      pivot_type = pivot_type,
      foci_lex_selector = foci_lex_selector,
      foci_morphosite_selector = foci_morphosite_selector,
      n_foci = n_foci,
      pivot_weighting = pivot_weighting,
      pivot_use_proportion = pivot_use_proportion,
      evidence_weighting = evidence_weighting,
      evidence_use_proportion = evidence_use_proportion,
      evidence_balance = evidence_balance,
      change_classwise = change_classwise,
      halt_method = halt_method,
      halt_n = halt_n
    ),
    mplat_finals = mplat_finals,
    changes = changes_df,
    stats = NULL
  )
}


#' Select replacement alleles according to the generalised substance or identity
#' model
#'
#' @param mplat A matrix of integers representing the current state of the
#'   inflectional system, rotated so that pivots are in columns.
#' @param foci_loc_b A vector of integers, column indices for mplat giving the
#'   position of the foci.
#' @param focus_loc_e An integers, a row index for mplat giving the positions of
#'   the foci.
#' @param pivot_weighting A string, giving the method for weighting pivots; used
#'   either for sampling them, or weighting them when all are taken into
#'   account.
#' @param pivot_use_proportion A numeric, giving the proportion of the available
#'   pivots to use. If set to 0, then just one is used.
#' @param evidence_weighting A string, giving the method for weighting evidence
#'   morposites; used either for sampling them, or weighting them when all are taken
#'   into account.
#' @param evidence_use_proportion A numeric, giving the proportion of the
#'   available evidence morposites to use. If set to 0, then just one is used.
#' @param evidence_balance A numeric, between 0 and 1, giving the balance of
#'   negative evidence (max when 0) and positive evidence (max when 1).
#' @return A vector of integers, the replacements for foci
select_replacement_alleles = function(
  mplat,
  foci_loc_b,
  focus_loc_e,
  pivot_weighting,
  pivot_use_proportion,
  evidence_weighting,
  evidence_use_proportion,
  evidence_balance
) {
  
  ## If there is only one possible result, return it now
  
  n_potential_results <- 
    nrow(unique(mplat[, foci_loc_b, drop = FALSE]))
  if (n_potential_results == 1) { return(mplat[1, foci_loc_b]) }
  
  
  ## Initialise some values
  
  # Parameters
  pvt_wtng <- substr(pivot_weighting[1], 1, 1)
  ev_wtng  <- substr(evidence_weighting[1], 1, 1)
  
  # mplat
  b_length <- ncol(mplat)
  e_length <- nrow(mplat)
  
  # Potential pivot and evidence locations (= all but the focus)
  pvt_pot_loc_b <- (1:b_length)[-foci_loc_b]
  ev_pot_loc_e  <- (1:e_length)[-focus_loc_e]
  n_pvt_pot_locs <- length(pvt_pot_loc_b)
  n_ev_pot_locs  <- length(ev_pot_loc_e)
  
  # Calculate number of locations from proportions
  n_pivot <- 
    round(n_pvt_pot_locs * pivot_use_proportion) %>%
    max(1) %>% min(n_pvt_pot_locs)
  n_evidence <- 
    round(n_ev_pot_locs * evidence_use_proportion) %>%
    max(1) %>% min(n_ev_pot_locs)
 
  
  ## Sample and weight the pivot & evidence locations

  # Pivot locations
  if (pvt_wtng == "c") {
    # Sample classwise: one per class (ignoring the foci cols)
    first_in_class <- !duplicated(mplat[, pvt_pot_loc_b, drop = FALSE])
    pvt_loc_b <- pvt_pot_loc_b[first_in_class]
    pvt_wt <- 1
  } else if (pvt_wtng %in% c("u", "z")) {
    # Get uniform or zipfian weights
    pvt_wt <- get_weights(pvt_wtng, b_length)[pvt_pot_loc_b]
    if (pivot_use_proportion < 1) {
      # Sample n_pivot locations using the weights
      pvt_loc_b <- sample(pvt_pot_loc_b, size = n_pivot, prob = pvt_wt)
      # Having sampled, now weight evenly
      pvt_wt <- rep(1, n_pivot) 
    } else {
      # Use all locations, and below, use pvt_wt to weight them
      pvt_loc_b <- pvt_pot_loc_b
    } 
  }
  
  # Evidence locations
  if (ev_wtng == "c") {
    # Sample classwise: one per class. This needs to include the focus row as a
    # possible selection, otherwise, otherwise the evolutionary process will
    # reduce down to just two remaining classes, whose member lexemes always
    # change to the other class.
    first_in_class <- !duplicated(mplat)
    ev_loc_e <- which(first_in_class)
    # first_in_class <- !duplicated(mplat[ev_pot_loc_e, , drop = FALSE])  # original code here
    # ev_loc_e <- ev_pot_loc_e[first_in_class]
    ev_wt <- 1
  } else if (ev_wtng %in% c("u", "z")) {
    # Get uniform or zipfian weights
    ev_wt  <- get_weights(ev_wtng, e_length)[ev_pot_loc_e]
    if (evidence_use_proportion < 1) {
      # Sample n_evidence locations using the weights
      ev_loc_e <- sample(ev_pot_loc_e, size = n_evidence, prob = ev_wt)
      # Having sampled, now weight evenly
      ev_wt <- rep(1, n_evidence)
    } else {
      # Use all locations, and below, use ev_wt to weight them
      ev_loc_e <- ev_pot_loc_e
    } 
  }
  n_ev_loc  <- length(ev_loc_e)
  
  
  ## Weight the location contents to select a winner
  
  # Attend only to the evidence rows. Also define the subset with just pivot or
  # just foci columns
  ev_mplat <- mplat[ev_loc_e, , drop = FALSE]
  ev_mplat_pvt  <- ev_mplat[, pvt_loc_b, drop = FALSE]
  ev_mplat_foci <- ev_mplat[, foci_loc_b, drop = FALSE]
  
  # Get a matrix of (non)identity observations: Rows r are rows in ev_mplat,
  # columns are for pivots. For each pivot, morphosites are 0/1 for
  # (non)identity of the alleles in row r vs those in focal row.
  compare_alleles <- mplat[focus_loc_e, pvt_loc_b, drop = FALSE]
  compare_mat  <- repeat_row(compare_alleles, n_ev_loc)
  identity_mat <- ev_mplat_pvt == compare_mat
  
  # Positive and negative weights = weighting due to: overall evidence balance;
  # individual evidence row weights; and individual pivot column weights:
  pos_mat <- 
    (identity_mat * evidence_balance) %>%
    scalar_mult_rows(ev_wt) %>%
    scalar_mult_cols(pvt_wt)
  neg_mat <- 
    ((!identity_mat) * (1 - evidence_balance)) %>%
    scalar_mult_rows(ev_wt) %>%
    scalar_mult_cols(pvt_wt)
  total_wt <- (pos_mat - neg_mat) %>% rowSums()
  
  # Take the foci evidence, then conflate repeated allele sets and sum their
  # weights:
  foci_wt <- conflate_and_sum_wt(ev_mplat_foci, total_wt)
  
  # Choose the replacement alleles by sampling from the equally highest weighted
  # options
  the_winner <- which_is_sampled_max(foci_wt$weight)
  foci_wt$matrix[the_winner, ]
  
}


#' Evaluate an expression n times, returning results in a list.
#' @param n An integer
#' @param expr An expression
#' @return A list
replicate_evo = function(n, expr) { 
  lapply(1:n, eval.parent(substitute(function(...) expr))) 
}