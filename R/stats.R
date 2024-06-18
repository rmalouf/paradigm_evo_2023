# Analysis of evolution runs


#' Powerlaw distribution
#'
#' @param n An integer, the number of items
#' @param exponent A numeric
power_dist = function(n, exponent) {
  exponent <- -abs(exponent)
  x <- (1:n)^exponent
  x / sum(x)
}


#' Powerlaw sample
#'
#' @param n An interger, the number of levels to sample from, with replacement
#' @param size The sample size
#' @param exponent A numeric
sample_power_dist = function(n, size, exponent) {
  sample.int(n, size, replace = TRUE,
             prob = power_dist(n, exponent))
}


#' Uniform and zipf distribution
#' @param dist_type A string, giving the distribution type
#' @param n An integer, the number of data points
#' @return A vector of numerics of length n
get_weights = function(dist_type, n) {
  d <- substr(dist_type[1], 1, 1)
  if (d == "u") { rep(1, n) }
  else if (d == "z") { power_dist(n, 1) }
  else { stop("Unrecognised dist_type. Use 'u' or 'z'.") }
}


#' Which is the maximum value, sampling among equal top
#' @param x A vector of numerics
#' @return A numeric
which_is_sampled_max = function(x) {
  if (length(x) == 1) { return(1) }
  tops <- which(x == max(x))
  if (length(tops) == 1) { return(tops) }
  sample(tops, 1)
}


#' Add statistics to an evolution
#'
#' @param evolution A list, the output from evolve_mplat.
#' @param step_method A string, specifying the kind of steps to take through the
#'   evolution: by generation or by change.
#' @param skip_steps An integer, specifying how many steps to skip ahead each
#'   time.
#' @return A list, the evolution including the stats
add_stats = function(
    evolution,
    step_method = c("generation", "change"),
    skip_steps = 20,
    entropy_classwise = FALSE
) {
  
  ## Approach to stats
  
  # Our stats involve a lot of comparison of pairs of rows
  # and columns. Since individual generations alter only
  # certain rows and columns, it would be wasteful to
  # redo all pairwise computation. Instead, we update only
  # those that will have changed.
  
  # We examine:
  #   1. Mean conditional entropy across columns (for "m" pivots) or 
  #      rows (for "l" pivots).
  #   2. Mean uncertainty coefficient across columns (for "m" pivots) or 
  #      rows (for "l" pivots).
  #   3. The number of rowwise and columnwise classes
  #   4. The number of distinct exponents in the system
  #   5. The turnover of classes since the last step
  
  ## Initialise some values
  p_type  <- substr(evolution$model$pivot_type[1],1,1)
  n_rep   <- evolution$n_rep
  mplat_0 <- evolution$mplat_0
  stats_df <- NULL
  
  ## stats_df assembly function
  get_stats_df <- function(gen) {
    # Create dataframe containing one row of stats
    data.frame(
      rep = rep_i,
      generation = gen,
      n_classes = n_classes,
      turnover = turnover,
      size_class1 = size_class1,
      size_class2 = size_class2,
      mean_exponents = mean_exponents,
      mean_cond_H = mean_ignoring_diag(cond_H),
      mean_U = mean_ignoring_diag(U)
    )
  }
  
  ## Cycle through simulation repetitions
  
  for (rep_i in 1:n_rep){
    
    if (p_type == "l") {
      mplat      <- t(mplat_0)
    } else if (p_type == "m") {
      mplat      <- mplat_0
    }
    
    ## Initialise simulation stats
    turnover <- 0
    cplat <- unique(mplat) # plat of classes
    if (rep_i == 1) { cplat_prev <- cplat }
    n_classes <- nrow(cplat)
    class_sizes <- sort(table(do.call(paste, data.frame(mplat))), decreasing = T)
    size_class1 <- class_sizes[1]
    size_class2 <- if (n_classes == 1) 0 else class_sizes[2]
    mean_exponents <- mean(count_exponents(cplat))
    
    # Entropy
    joint_H <- # matrix of joint H of pairs of rows
      if (entropy_classwise) get_col_joint_H(cplat) else get_col_joint_H(mplat) 
    H_X         <- square_mat(cols = diag(joint_H)) # matrix, every col = H_x
    H_Y         <- square_mat(rows = diag(joint_H)) # matrix, every row = H_y
    cond_H      <- joint_H - H_Y # H(X|Y)
    I           <- H_X - cond_H # mutual information = H(X) - H(X|Y)
    U           <- I / H_X  # uncertainty coefficient U(X|Y)
    U[is.na(U)] <- 0 # Define U as 0 when H_X == 0
    
    # Bookkeeping 
    plat_stats_df <- get_stats_df(gen = NA)
    changes <- evolution$changes %>% filter(repetition == rep_i)
    steps   <- get_steps(evolution, rep_i, step_method, skip_steps)
    
    ## Cycle through simulation steps
    
    for(i in 2:length(steps$rows)) {
      
      this_row <- steps$rows[i]
      prev_row <- steps$rows[i - 1]
      is_change <- this_row > prev_row  # did a change happen between steps?
      
      if (is_change) {
        # A change happened. mplat needs updating. 
        # Get an update of it via all intervening changes:
        chg <- changes[(prev_row + 1):this_row, ]
        if (p_type == "m") {
          # Operating with an unrotated plat
          update <- update_mplat(mplat, chg)
          mplat <- update$mplat # updated mplat
          ro <- update$rows # Which rows and cols changed. Updates below target them.
          co <- update$cols
        } else if (p_type == "l") {
          # Operating with a rotated plat
          update <- update_mplat(t(mplat), chg)
          mplat <- t(update$mplat)
          ro <- update$cols 
          co <- update$rows
        }
        
        # Update n_classes, n_exponents and turnover
        cplat <- unique(mplat) # plat of classes
        n_classes <- nrow(cplat)
        class_sizes <- sort(table(do.call(paste, data.frame(mplat))), decreasing = T)
        size_class1 <- class_sizes[1]
        size_class2 <- if (n_classes == 1) 0 else class_sizes[2]
        mean_exponents <- mean(count_exponents(cplat))
        
        # Update entropy
        joint_H <- # matrix of joint H of pairs of rows
          if (entropy_classwise) get_col_joint_H(cplat) else get_col_joint_H(mplat) 
        H_X         <- square_mat(cols = diag(joint_H)) # matrix, every col = H_x
        H_Y         <- square_mat(rows = diag(joint_H)) # matrix, every row = H_y
        cond_H      <- joint_H - H_Y # H(X|Y)
        I           <- H_X - cond_H # mutual information = H(X) - H(X|Y)
        U           <- I / H_X  # uncertainty coefficient U(X|Y)
        U[is.na(U)] <- 0 # Define U as 0 when H_X == 0
        
        # Update stats df: will contain one row
         plat_stats_df <- get_stats_df(gen = steps$gens[i])
      }
      
      # Calculate class turnover
      cplat_and_prev <- rbind(cplat_prev, cplat)
      union_n <- sum(!duplicated(cplat_and_prev))
      intersection_n <- sum(duplicated(cplat_and_prev))
      plat_stats_df$turnover <- union_n - intersection_n 
      cplat_prev <- cplat
      
      # Append the current stats to stats_df
      plat_stats_df$generation <- steps$gens[i]
      stats_df <- bind_rows(stats_df, plat_stats_df)
    }
  }
  
  # Shift the turnover data back one time step:
  stats_df <-
    stats_df %>%
    group_by(rep) %>%
    arrange(generation) %>%
    mutate(turnover = lead(turnover)) %>%
    ungroup()
  
  evolution$stats <- stats_df
  evolution
}


#' Get pairwise joint entropies of rows
#' @param m A matrix
#' @return A matrix of numerics,
#'   the joint entropies.
get_row_joint_H = function(m) {
  nr = nrow(m)
  h_mat <- matrix(0, nrow = nr, ncol = nr)
  for (i in 1:nr) {
    for (j in i:nr) {
      h_mat[i,j] <- h_mat[j,i] <-
        # Entropy in bits
        entropy(data.frame(r1 = m[i, ], r2 = m[j, ])) / log(2)
    }
  }
  h_mat
}


get_col_joint_H = function(m) { get_row_joint_H(t(m)) }


#' Update pairwise joint entropies of rows
#' @param h_mat A matrix of H
#' @param m A matrix, whose contents to get the entropy of.
#' @param rows A vector of integers, the rows to update.
#' @return A matrix of numerics, the joint entropies.
update_row_joint_H = function(h_mat, m, rows) {
  nr = nrow(m)
  for (i in rows) {
    for (j in 1:nr) {
      h_mat[i,j] <- h_mat[j,i] <-
        # Entropy in bits
        entropy(data.frame(r1 = m[i, ], r2 = m[j, ])) / log(2)
    }
  }
  h_mat
}


update_col_joint_H = function(h_mat, m, cols) {
  update_row_joint_H(h_mat, t(m), cols)
}


#' Comma delimited string to vector of numerics
#' @param s A string
#' @return A vector of numerics
str2ints = function(s) {
  s %>% str_split(",") %>% unlist() %>% as.integer()
}


#' Count the number of unique values in rows or columns
#' @param m A matrix
#' @param p_type A string
#' @return A vector of counts
count_exponents = function(m) {
  apply(m, 2, function(v) { length(unique(v)) })
}


#' Count the number of ways to place p stones in an m x n grid, under the
#' constraints that no row is empty; no column is empty; only one stone per
#' cell; and cell (1,1) is filled.
#' @param m An integer. The number of grid rows.
#' @param n An integer. The number of grid columns.
#' @param p An integer. The number of stones
#' @return An integer. The number of solutions.
count_arrangements = function(m,n,p) {
  # Calculate an inclusion-exclusion sum
  summands <- matrix(0, m, n)
  for (n_0cols in 0:(n-1)) { # number of empty columns
    for (n_0rows in 0:(m-1)) { # number of empty rows
      choose_0cols <- choose(n - 1, n_0cols) # choices of the empty columns
      choose_0rows <- choose(m - 1, n_0rows) # choices of the empty rows
      n_available_cols <- n - n_0cols
      n_available_rows <- m - n_0rows
      n_available_cells <- n_available_cols * n_available_rows - 1
      choose_cells = choose(n_available_cells, p - 1)
      incl_excl_sign <- (-1)^(n_0cols + n_0rows)
      summands[n_0rows+1, n_0cols+1] <- incl_excl_sign * choose_0rows * choose_0cols * choose_cells
    }
  }
  print(summands)
  sum(summands)
}

#' Count the number of ways to place p stones in an m x n grid, under the
#' constraints that no row is empty; no column is empty; only one stone per
#' cell; cell (1,1) is filled; and cell (1,j) must also be filled, for some j >
#' 1. This indicates the how many times a cell (1,j) is filled, among all the
#' solutions to the grid-filling constraints described for
#' `count_arrangements()`.
#' @param m An integer. The number of grid rows.
#' @param n An integer. The number of grid columns.
#' @param p An integer. The number of stones
#' @return An integer. The number of solutions.
count_arrangements_1j = function(m, n, p) {
  # Calculate an inclusion-exclusion sum
  summands <- matrix(0, m, n-1)
  for (n_0cols in 0:(n-2)) { # number of empty columns
    for (n_0rows in 0:(m-1)) { # number of empty rows
      choose_0cols <- choose(n - 2, n_0cols) # choices of the empty columns
      choose_0rows <- choose(m - 1, n_0rows) # choices of the empty rows
      n_available_cols <- n - n_0cols
      n_available_rows <- m - n_0rows
      n_available_cells <- n_available_cols * n_available_rows - 1
      choose_cells = choose(n_available_cells, p - 2)
      incl_excl_sign <- (-1)^(n_0cols + n_0rows)
      summands[n_0rows+1, n_0cols+1] <- incl_excl_sign * choose_0rows * choose_0cols * choose_cells
    }
  }
  print(summands)
  sum(summands)
}

#' Count the number of ways to place p stones in an m x n grid, under the
#' constraints that no row is empty; no column is empty; only one stone per
#' cell; cell (1,1) is filled; and cell (i,1) must also be filled, for some i >
#' 1. This indicates the how many times a cell (i,1) is filled, among all the
#' solutions to the grid-filling constraints described for
#' `count_arrangements()`.
#' @param m An integer. The number of grid rows.
#' @param n An integer. The number of grid columns.
#' @param p An integer. The number of stones
#' @return An integer. The number of solutions.
count_arrangements_i1 = function(m, n, p) {
  # Calculate an inclusion-exclusion sum
  summands <- matrix(0, m-1, n)
  for (n_0cols in 0:(n-1)) { # number of empty columns
    for (n_0rows in 0:(m-2)) { # number of empty rows
      choose_0cols <- choose(n - 1, n_0cols) # choices of the empty columns
      choose_0rows <- choose(m - 2, n_0rows) # choices of the empty rows
      n_available_cols <- n - n_0cols
      n_available_rows <- m - n_0rows
      n_available_cells <- n_available_cols * n_available_rows - 1
      choose_cells = choose(n_available_cells, p - 2)
      incl_excl_sign <- (-1)^(n_0cols + n_0rows)
      summands[n_0rows+1, n_0cols+1] <- incl_excl_sign * choose_0rows * choose_0cols * choose_cells
    }
  }
  print(summands)
  sum(summands)
}

#' Count the number of ways to place p stones in an m x n grid, under the
#' constraints that no row is empty; no column is empty; only one stone per
#' cell; cell (1,1) is filled; and cell (i,j) must also be filled, for some i >
#' 1, j > 1. This indicates the how many times a cell (i,j) is filled, among all
#' the solutions to the grid-filling constraints described for
#' `count_arrangements()`.
#' @param m An integer. The number of grid rows.
#' @param n An integer. The number of grid columns.
#' @param p An integer. The number of stones
#' @return An integer. The number of solutions.
count_arrangements_ij = function(m, n, p) {
  # Calculate an inclusion-exclusion sum
  summands <- matrix(0, m-1, n-1)
  for (n_0cols in 0:(n-2)) { # number of empty columns
    for (n_0rows in 0:(m-2)) { # number of empty rows
      choose_0cols <- choose(n - 2, n_0cols) # choices of the empty columns
      choose_0rows <- choose(m - 2, n_0rows) # choices of the empty rows
      n_available_cols <- n - n_0cols
      n_available_rows <- m - n_0rows
      n_available_cells <- n_available_cols * n_available_rows - 1
      choose_cells = choose(n_available_cells, p - 2)
      incl_excl_sign <- (-1)^(n_0cols + n_0rows)
      summands[n_0rows+1, n_0cols+1] <- incl_excl_sign * choose_0rows * choose_0cols * choose_cells
    }
  }
  print(summands)
  sum(summands)
}

#' Generate the solutions to placing p stones in an m x n grid, under the
#' constraints that no row is empty; no column is empty; only one stone per
#' cell; and cell (1,1) is filled.
#' @param m An integer. The number of grid rows.
#' @param n An integer. The number of grid columns.
#' @param p An integer. The number of stones
#' @return A list of matrices. The solutions, where 1 = filled cell.
generate_arrangements = function(m, n, p){
  cell_indices <- 2:(m * n)
  combos <- combn(cell_indices, p - 1)
  n_combos <- ncol(combos)
  matdat <- c(1, rep(0, m * n - 1))
  solution_list <- list()
  blanks <- matrix(0, m, n)
  for (i in 1:n_combos) {
    matdati <- matdat
    matdati[combos[,i]] <- 1
    mati <- matrix(matdati, nrow = m)
    rblank <- sum(rowSums(mati) == 0)
    cblank <- sum(colSums(mati) == 0)
    blanks[rblank + 1, cblank + 1] <- blanks[rblank + 1, cblank + 1] + 1
    if (rblank == 0 | cblank == 0) { next }
    solution_list[[length(solution_list) +1]] <- mati
  }
  print(blanks)
  solution_list
}
