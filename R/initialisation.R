# Generate initial states

#' Initialise a parametrically random mplat
#'
#' @param x Either an integer or a list. If a list, x is treated a list of
#'   parameters, else as the number of lexemes.
#' @param polymorphies A vector of integers: the number of alleles for each
#'   morphosite.
#' @param allele_distribution A numeric, giving the exponent of a powerlaw
#'   distribution for skewing the sampling from a morphosite's alleles. A value
#'   of 0 is uniform; a value of 1 is Zipfian.
#' @param class_distribution A numeric, giving the exponent of a powerlaw
#'   distribution for skewing the repetitive sampling of lexical classes.
#' @param col_class_distribution A numeric, giving the exponent of a powerlaw
#'   distribution for skewing the repetitive sampling of morphosites, i.e.,
#'   having multiple morphosites repeating the same pattern.
#' @param max_classes An integer, placing an upper bound on the number of
#'   distinct classes sampled when class_distribution is not zero.
#' @param max_col_classes An integer, placing an upper bound on the number of
#'   distinct morphosites sampled when col_class_distribution is not zero.
#' @return A matrix, the mplat.
initialise_mplat = function(
  x = 50,
  polymorphies = rep(6, times = 8),
  allele_distribution = c("uniform", "zipfian"),
  class_distribution = c("uniform", "zipfian"),
  col_class_distribution = c("uniform", "zipfian"),
  max_classes = NULL,
  max_col_classes = NULL
) {

  if (is.list(x)) {
    # Assign variable values from param_dict
    for (i in 1:length(x)) {
      assign(names(x)[i], x[[i]])
    }
  } else {
    n_lexemes <- x
  }
    
  al_dist = substr(allele_distribution[1],1,1)
  cl_dist = substr(class_distribution[1],1,1)
  co_dist = substr(col_class_distribution[1],1,1)
  
  nr <- n_lexemes
  nc <- length(polymorphies)
  mplat <- matrix(0, nr, nc)
    
  for(i in 1:nc) {
    weights <- get_weights(al_dist, n = polymorphies[i])
    alleles <- sample(polymorphies[i], size = nr, 
                      replace = TRUE, prob = weights)
    mplat[, i] <- alleles
  }
  
  if (!is.null(max_classes)) {
    max_classes <- min(max_classes, nr)
    weights <- get_weights(cl_dist, n = max_classes)
    rows <- sample(max_classes, size = nr, 
                   replace = TRUE, prob = weights)
    mplat <- mplat[rows, ]
  }
  
  if (!is.null(max_col_classes)) {
    max_col_classes <- min(max_col_classes, nc)
    weights <- get_weights(co_dist, n = max_col_classes)
    cols <- sample(max_col_classes, size = nc, 
                   replace = TRUE, prob = weights)
    mplat <- mplat[, cols]
  }
  
  mplat
}


#' Make a list of mplat parameters
#' @return A list, of parameters for `initialise_mplat()`
mplat_params = function(
  n_lexemes = 50,
  polymorphies = rep(6, times = 8),
  allele_distribution = "uniform",
  class_distribution = "uniform",
  col_class_distribution = "uniform",
  max_classes = NULL,
  max_col_classes = NULL
){
  list(
    n_lexemes = n_lexemes,
    polymorphies = polymorphies,
    allele_distribution = allele_distribution,
    class_distribution = class_distribution,
    col_class_distribution = col_class_distribution,
    max_classes = max_classes,
    max_col_classes = max_col_classes
  )
}

