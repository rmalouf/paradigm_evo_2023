# Matrix tools

#' Sort a matrix by all columns, left to right
#' @param m A matrix
#' @return A matrix, same as m but with
#'   rows sorted.
sort_mat = function(m) {
  m[do.call(order, lapply(1:ncol(m), function(i) m[, i])), , drop = FALSE]
}


#' Which rows in a matrix are identical
#' to a given vector
#' 
#' @param m A matrix
#' @param v A vector, of length ncol(m)
#' @return A vector of row indices
match_rows = function(m, v) {
  which(apply(m, 1, function(r) all(r == v)))
}


#' Which rows in matrix 1 are in matrix 2
#' Only works for sorted matrices!
#' 
#' @param m1 A matrix
#' @param m2 A martix
#' @return A vector of row indices
match_sorted_rows = function(m1, m2) {
  key1 = setkey(data.table(m1))
  key2 = setkey(data.table(m2))
  r <- key1[key2, which=TRUE]
  r[!is.na(r)]
}


#' Replace m1 cells with non-na cells
#' of m2
#' @param m1 A matrix
#' @param m2 A matrix
#' @return A matrix
update_mat = function(m1, m2) {
  if (any(dim(m1) != dim(m2))) { 
    stop("Matrices must be same dimensions") 
    }
  is_update <- !is.na(m2)
  m1[is_update] <- m2[is_update]
  m1
}


#' Shuffle within rows
#' @param m A matrix.
#' @param rows A vector of integers,
#'   the rows to shuffle within.
#' @return A matrix, with each
#'   row's contents shuffled.
shuffle_within_rows = function(
  m, 
  rows = 1:nrow(m)
) {
  m[rows, ] <- apply(m[rows, , drop = FALSE], 1, sample)
  m
}


#' Shuffle within columns
#' @param m A matrix.
#' @param cols A vector of integers,
#'   the columns to shuffle within.
#' @return A matrix, with each.
#'   column's contents shuffled
shuffle_within_cols = function(
  m,
  cols = 1:ncol(m)
) {
  m[, cols] <- apply(m[, cols, drop = FALSE], 2, sample)
  m
}


#' Merge matrix rows that are identical outside
#' of the first column; sum the first column when 
#' rows are merged
#' @param m A matrix
#' @param w A vector of weights
#' @return A list, containing 
#'   $matrix and $weight
conflate_and_sum_wt = function(m, w) {
  
  if (!is.matrix(m)) { stop("m must be a matrix.") }
  
  nr <- nrow(m)
  if (nr == 1) { return(list(matrix = m, weight = 1)) }
  
  m <- sort_mat(cbind(m, w))
  nc <- ncol(m)
  r <- 1
  while (r < nr) {
    
    if (all(m[r, -nc] == m[r + 1, -nc])) {
      
      # Rows r and r+1 are identical: compress them 
      
      # Sum the first column of rows r and r+1, put
      # the total in row r
      m[r, nc] <- m[r, nc] + m[r + 1, nc]
      
      # Delete row r+1, update nr
      m <- m[-(r+1), , drop = FALSE]
      nr <- nr -1
      
    } else {
      
      # Rows r and r+1 are non-identical: progress to next
      r <- r + 1
      
    }
  }
  list(matrix = m[, -nc, drop = FALSE], weight = m[, nc])
}


#' Cluster-sort the rows of a matrix
#' @param m A matrix
#' @return A matrix
cluster_sort_mat = function(m) {
  m[cluster_sort_order(m), ]
}


#' Cluster order
#' @param m A matrix
#' @return A vector of integers,
#'   the sort order of rows.
cluster_sort_order = function(m) {
  m %>%
    dist() %>%
    hclust() %>%
    as.dendrogram() %>%
    order.dendrogram()
}


#' Multiply each column by an amount
#' @param m A matrix
#' @param x A vector length ncol(m)
#' @return A matrix
scalar_mult_cols = function(m, x) {
  t(t(m) * x)
}


#' Multiply each row by an amount
#' @param m A matrix
#' @param x A vector length nrow(m)
#' @return A matrix
scalar_mult_rows = function(m, x) {
  m * x
}


#' Make a matrix of a repeated row
#' @param row A vector
#' @param n An integer
#' @return A matrix
repeat_row = function(row, n) {
  do.call(rbind, rep(list(row), n))
}


#' Sum rows' minimums
#' @param m A matrix
sum_row_mins = function(m) {
  sum(apply(m, 1, min, na.rm = T))
}


#' Sum rows' minimums in lower triangle
#' @param m A matrix
sum_row_mins_in_lower_triangle = function(m) {
  m[!lower.tri(m)] <- NA
  sum(apply(m[-1, ], 1, min, na.rm = T))
}


#' Mean without diagonal
#' @param m A matrix
mean_ignoring_diag = function(m) {
  if (length(m) == 1) { return(m[1,1]) }
  diag(m) <- NA
  mean(m, na.rm = T)
}


#' Turn a vector into square matrix by repeating it
#' either in every row or every column
#' @param rows A vector
#' @param cols A vector
square_mat = function(
  rows, cols
) {
  if (missing(rows)) {
    n <- length(cols)
    dat <- rep(cols, n)
  } else {
    n <- length(rows)
    dat <- rep(rows, each = n)
  }
  matrix(dat, nrow = n)
}