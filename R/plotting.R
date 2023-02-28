# Plotting

#' Plot an mplat
#'
#' @param x A matrix or a list. If a list, it is used as parameters for
#'   initialise_mplat(); if a matrix, then it is an mplot.
#' @param cluster_lexemes A logical. Whether to sort the lexemes by similarity.
#' @param cluster_morphosites A logical. Whether to sort the morphosites by
#'   similarity.
#' @param multicolour A logical. Whether to put column/rows in different
#'   colours.
#' @param multicolour_axis A string, which axis to apply multicolour to.
#' @return A ggplot raster plot.
plot_mplat = function(
  x, 
  cluster_lexemes = FALSE,
  cluster_morphosites = FALSE,
  multicolour = FALSE,
  multicolour_axis = c("columns", "rows")
)  {

  if (is.list(x)) {
    mplat <- initialise_mplat(x)
  } else {
    mplat <- x
  }
  
  if (cluster_lexemes) {
    mplat <- cluster_sort_mat(mplat)
  }
  
  if (cluster_morphosites) {
    mplat <- t(cluster_sort_mat(t(mplat)))
  }
  
  if (multicolour) {
    mc <- multicolourise(mplat, multicolour_axis)
    mplat <- mc$data
    colours <- mc$colours
    mplat %>%
      pivot_mplat_longer() %>%
      ggplot(aes(morphosite, -lexeme)) +
      geom_raster(aes(fill = as.factor(allele))) +
      theme_void() +
      theme(legend.position="none") +
      scale_fill_manual(values = colours)

  } else {
    
    mplat %>%
      pivot_mplat_longer() %>%
      ggplot(aes(morphosite, -lexeme)) +
      geom_raster(aes(fill = allele)) +
      theme_void() +
      theme(legend.position="none")
  }
}



#' Animate the evolution of a mplat
#'
#' @param evolution A list, the output from evolve_mplat.
#' @param which_rep An integer, specifying which repetition to plot.
#' @param step_method A string, specifying the kind of steps to take through the
#'   evolution: by generation or by change.
#' @param skip_steps An integer, specifying how many steps to skip ahead each
#'   time.
#' @param telic_order A logical, Whether to order columns and rows so they
#'   cluster at the final generation.
#' @param multicolour A logical. Whether to put column/rows in different
#'   colours.
#' @param multicolour_axis A string, which axis to apply multicolour to.
#' @param fps An integer
#' @return A dataframe of the evolution in long form.
plot_evolution = function(
  evolution,
  which_rep = 1,
  step_method = c("generation", "change"),
  skip_steps = 
    if (substr(step_method[1],1,1) == "g") { 
      evolution$total_generations / 100 
    } else {
      nrow(evolution$change) / 100 
    },
  telic_order = FALSE,
  multicolour = FALSE,
  multicolour_axis = c("columns", "rows"),
  fps = 10
) {
  
  m <- evolution$mplat_final[[which_rep]]
  if (telic_order) {
    row_order <- cluster_sort_order(m)
    # col_order <- cluster_sort_order(t(m))
    col_order <- 1:ncol(m)
  } else {
    row_order <- nrow(m):1
    col_order <- 1:ncol(m)
  }
  
  evo_df <- 
    unpack_evolution(
      evolution, repetition = which_rep, step_method, skip_steps
    )
  
  if (multicolour) {
    
    mc <- multicolourise(evo_df, multicolour_axis)
    evo_df <- mc$data
    colours <- mc$colours
    p <-
      evo_df %>%
      ggplot(aes(match(morphosite, col_order), 
                 match(lexeme, row_order))) +
      geom_tile(aes(fill = as.factor(allele))) +
      theme_void() +
      theme(legend.position="none") +
      scale_fill_manual(values = colours) +
      # animation specific part
      transition_manual(step)
    
  } else {
    
    p <-
      evo_df %>%
      ggplot(aes(match(morphosite, col_order), 
                 match(lexeme, row_order))) +
      geom_tile(aes(fill = allele)) +
      theme_void() +
      theme(legend.position="none") +
      # animation specific part
      transition_manual(step)
  }
  
  animate(p, 
          nframes = max(evo_df$step),
          # fps = fps, 
          renderer = magick_renderer(), 
          end_pause = 2, detail = 0)
}


#' Add multiple colours
#' @param data A matrix, giving an mplat or a dataframe giving a long form
#'   evo_df.
#' @param axis A string, indicating which axis to apply multicolour to.
#' @return A list, containing the modified mplat and a vector of colours.
multicolourise = function(
    data, 
    axis = c("columns", "rows")
) {
  
  data_type = class(data)
  
  ## Examine data contents relevant to
  #  colour choices
  
  if ("matrix" %in% data_type) {
    
    # Data is mplat
    m <- data
    if (substr(axis[1],1,1) == "c") {
      # temporarily rotate, and apply colouring
      # to rows
      m <- t(m)
    }
    m_size <- nrow(m)
    n_alleles <- apply(m, 1, function(v) length(unique(v)))
    
  } else {
    
    # Data is evo_df
    evo_df <- data
    if (substr(axis[1],1,1) == "c") {
      evo_df <- evo_df %>% rename(row = morphosite)
    } else {
      evo_df <- evo_df %>% rename(row = lexeme)
    }
    m_size <- max(evo_df$row)
    unique_alleles_df <- 
      evo_df %>% 
      select(row, allele) %>%
      distinct() %>%
      arrange(row, allele)
    n_alleles <- 
      unique_alleles_df %>% 
      group_by(row) %>%
      summarise(n = n()) %>%
      .$n
  }
  
  ## Choose colours
  
  # We use four colour types, unless there are
  # fewer than four rows/cols.
  n_types <- min(4, m_size)
  n_padding <- max(0, 4 - n_types)
  
  # Each type is used every four columns/rows, so 
  # needs enough variants to handle the most-diverse one
  n_colour_variants <-
    c(
      sapply(1:n_types, function(i) {
        rows_with_this_colour <- seq(i, m_size, by = 4)
        max(n_alleles[rows_with_this_colour])
      }),
      rep(0, n_padding)
    )
  
  colours <- c(
    get_colours(n_colour_variants[1], "Blues"), 
    get_colours(n_colour_variants[2], "Greys"), 
    get_colours(n_colour_variants[3], "Greens"), 
    get_colours(n_colour_variants[4], "Reds"))
  
  n_padding <- max(0, 5 - length(colours))
  colours <- c(colours, rep(0, n_padding))
  
  
  ## Modify data to work with colours
  
  # Now change the values in the cells, so
  # colours will match. This will mean, e.g.
  # increasing values in the green cells by
  # n, where n is the number of blue colour
  # variants
  value_increase <- 
    c(0, cumsum(n_colour_variants[1:3])) %>%
    rep_len(length.out = m_size)
  
  if ("matrix" %in% data_type) {
    
    # Data is mplat
    
    # Renumber alleles in each set to be 1,2,3...
    m <- t(apply(m, 1, function(v) { as.numeric(as.factor(v))}))
    # Apply increase
    m <- m + value_increase
    if (substr(axis[1],1,1) == "c") {
      m <- t(m) # rotate back again
    }
    data <- m
    
  } else {
    
    # Data is evo_df
    
    # Renumber alleles in each set to be 1,2,3...
    renum <-
      unique_alleles_df %>%
      group_by(row) %>%
      mutate(new_val = row_number())
    evo_df <-
      evo_df %>% 
      left_join(renum, by = c("row", "allele")) %>%
      select(-allele) %>%
      rename(allele = new_val)
    # Apply increase
    evo_df <-
      evo_df %>%
      mutate(allele = allele + value_increase[row])
    
    if (substr(axis[1],1,1) == "c") {
      evo_df <- evo_df %>% rename(morphosite = row)
    } else {
      evo_df <- evo_df %>% rename(lexeme = row)
    }
    data <- evo_df
    
  }
  
  ## Return modified data and colours
  
  list(
    data = data,
    colours = colours
  )
}


#' Choose n brewer.pal colours, but allow n < 4, in which case choose more
#' saturated end of the n=4 set
#' @param n An integer, how many colours
#' @param col A string, the palette to choose from
#' @return A vector of hex colour specifications
get_colours <- function(n, col) {
  if (n == 0) { character(0) }
  else if (n == 1 | n == 2) { rev(brewer.pal(4, col))[1:n] }
  else { rev(brewer.pal(n + 1, col)[-1]) }
}


#' Plot stats
#' @param evolution A list, the output from evolve_mplat.
#' @param suppress_top2 A logical, whether to suppress the plotting of the top 2
#'   classes
#' @param turnover_limits A vector of two numerics. The y axis limits for the
#'   plot of turnover
#' @param turnover_breaks A vector of numerics. The y axis tick mark locations
#'   for the plot of turnover
#' @param classes_limits A vector of two numerics. The y axis limits for the
#'   plot of classes
#' @param classes_breaks A vector of numerics. The y axis tick mark locations
#'   for the plot of classes
plot_stats = function(
    evolution,
    suppress_top2 = FALSE,
    turnover_limits = c(0, 200), 
    turnover_breaks = c(0, 2, 8, 20, 50, 100, 200),
    classes_limits = c(1, 100),
    classes_breaks = c(1, 2, 4, 8, 20, 50, 100)
) {
  
  p_H <- plot_H(evolution)
  p_U <- plot_U(evolution)
  p_c <- plot_class_contrasts(evolution, classes_limits, classes_breaks)
  p_t <- plot_turnover(evolution, turnover_limits, turnover_breaks)
  p_e <- plot_exp_contrasts(evolution)
  p_t2 <- plot_top2classes(evolution)

  plotlist <- list(
    p_H + theme_minimal() + rremove("xlab"),
    p_c + theme_minimal() + rremove("xlab"),
    p_U + theme_minimal() + rremove("xlab"),
    p_t + theme_minimal() + rremove("xlab"),
    p_e + theme_minimal(),
    p_t2 + theme_minimal()
    )
  
  if (suppress_top2) { plotlist <- plotlist[-6] }
  
  ggarrange(plotlist = plotlist, nrow = 3, ncol = 2)
}


#' Plot entropy stats
#' @param evolution A list, the output from evolve_mplat.
#' @param y_limits A vector of two numerics. The y axis limits. 
#' @param y_breaks A vector of numerics. The y axis tick mark locations.
#' @param step A string, specifying the kind of steps to take through the
#'   evolution: by generation or by change.
plot_H = function(
    evolution,
    y_limits = c(0, 3),
    y_breaks = c(0, 0.05, 0.2, 0.5, 1, 3),
    step = c("generation", "change")
) {
  
  dat <- evolution$stats
  
  if (substr(step[1],1,1) == "c") {
    # Use changes rather than generations
    dat <- dat %>% select(-generation) %>% distinct()
    dat$generation <- 1:nrow(dat)
    x_label <- "Change"
  } else {
    x_label <- "Cycle"
  }
  
  plot_dat <-
    dat %>% 
    rename(value = mean_cond_H) %>% 
    compile_mean_hi_lo() %>%
    rename(entropy = mean_value) %>%
    mutate(lower = lower + 0.0001) %>% # Prevent bug where 0 doesn't plot
    filter(entropy != -1)
  
  p <- ggplot(plot_dat)
    
  
  if (evolution$n_rep > 1) {
    p <- p +
      geom_ribbon(
        aes(x = generation, ymin = lower, ymax = upper),
        alpha = 0.3, colour = NA) 
  } 
  
  p + 
    geom_line(aes(x = generation, y = entropy), size = 1) +
    labs(x = x_label, y = "Mean H(X|Y)") +
    scale_y_continuous(
      trans = loglike_trans(),
      limits = y_limits,
      breaks = y_breaks,
      labels = scales::label_number(accuracy = 0.01)
    )
}


#' Plot entropy stats
#' @param evolution A list, the output from evolve_mplat.
#' @param y_limits A vector of two numerics. The y axis limits. 
#' @param y_breaks A vector of numerics. The y axis tick mark locations.
#' @param step A string, specifying the kind of steps to take through the
#'   evolution: by generation or by change.
plot_U = function(
  evolution,
  y_limits = c(0, 1),
  y_breaks = c(0, 0.05, 0.2, 0.5, 1),
  step = c("generation", "change")
) {
  
  dat <- evolution$stats
  
  if (substr(step[1],1,1) == "c") {
    # Use changes rather than generations
    dat <- dat %>% select(-generation) %>% distinct()
    dat$generation <- 1:nrow(dat)
    x_label <- "Change"
  } else {
    x_label <- "Cycle"
  }
  
  plot_dat <-
    dat %>% 
    rename(value = mean_U) %>% 
    compile_mean_hi_lo() %>%
    rename(ThielsU = mean_value) %>%
    filter(ThielsU != -1)
  
  p <- ggplot(plot_dat)
  
  if (evolution$n_rep > 1) {
    p <- p +
      geom_ribbon(
        aes(x = generation, ymin = lower, ymax = upper),
        alpha = 0.3, colour = NA) 
  } 
  
  p + 
    geom_line(aes(x = generation, y = ThielsU), size = 1) +
    scale_y_continuous(
      trans = loglike_trans(),
      limits = y_limits,
      breaks = y_breaks,
      labels = scales::label_number(accuracy = 0.01) 
    ) +
    scale_x_continuous(limits = c(1, max(dat$generation))) +
    labs(x = x_label, y = "Mean U(X|Y)")
}



#' Plot class contrast stats
#' @param evolution A list, the output from evolve_mplat.
#' @param y_limits A vector of two numerics. The y axis limits. 
#' @param y_breaks A vector of numerics. The y axis tick mark locations.
#' @param step A string, specifying the kind of steps to take through the
#'   evolution: by generation or by change.
plot_class_contrasts = function(
    evolution,
    y_limits = c(1, max(evolution$stats$n_classes)),
    y_breaks = c(1,5,20,50,100),
    step = c("generation", "change")
) {
  
  dat <- evolution$stats
  
  if (substr(step[1],1,1) == "c") {
    # Use changes rather than generations
    dat <- dat %>% select(-generation) %>% distinct()
    dat$generation <- 1:nrow(dat)
    x_label <- "Change"
  } else {
    x_label <- "Cycle"
  }
  
  pvt_type <- substr(evolution$model$pivot_type[1], 1, 1)
  if (pvt_type == "l") {
    y_label <- "Morphomic zones" 
  } else {
    y_label <- "Classes" 
  }

  plot_dat <-
    dat %>% 
    rename(value = n_classes) %>% 
    compile_mean_hi_lo() %>%
    rename(n_classes = mean_value)
  
  p <- ggplot(plot_dat)
  
  if (evolution$n_rep > 1) {
    p <- 
      p +
      geom_ribbon(
        aes(x = generation, ymin = lower, ymax = upper),
        alpha = 0.3, colour = NA) 
  } 
  
  p <- 
    p + 
    geom_line(aes(x = generation, y = n_classes), size = 1) +
    labs(linetype = "", x = x_label, y = y_label) +
    scale_y_continuous(
      trans = loglike_trans(),
      breaks = y_breaks,
      limits = y_limits,
      labels = scales::label_number(accuracy = 1)
    )

  p
  
}


#' Plot exponent contrast stats
#' @param evolution A list, the output from evolve_mplat.
#' @param y_limits A vector of two numerics. The y axis limits. 
#' @param y_breaks A vector of numerics. The y axis tick mark locations.
#' @param step A string, specifying the kind of steps to take through the
#'   evolution: by generation or by change.
plot_exp_contrasts = function(
    evolution,
    y_limits = c(1, 5),
    y_breaks = c(1, 1.5, 2, 3, 4, 5),
    step = c("generation", "change")
) {
  
  dat <- evolution$stats
  
  if (substr(step[1],1,1) == "c") {
    # Use changes rather than generations
    dat <- dat %>% select(-generation) %>% distinct()
    dat$generation <- 1:nrow(dat)
    x_label <- "Change"
  } else {
    x_label <- "Cycle"
  }
  
  pvt_type <- substr(evolution$model$pivot_type[1], 1, 1)
  if (pvt_type == "l") {
    y_label <- "Indices per lexeme" 
  } else {
    y_label <- "Exponents per cell" 
  }
  
  plot_dat <-
    dat %>% 
    rename(value = mean_exponents) %>% 
    compile_mean_hi_lo() %>%
    rename(mean_exponents = mean_value)
  
  p <- ggplot(plot_dat)
  
  if (evolution$n_rep > 1) {
    p <- 
      p +
      geom_ribbon(
        aes(x = generation, ymin = lower, ymax = upper),
        alpha = 0.3, colour = NA) 
  } 
  
  p <- 
    p + 
    geom_line(aes(x = generation, y = mean_exponents), size = 1) +
    labs(linetype = "", x = x_label, y = y_label) +
    scale_y_continuous(
      trans = loglike_trans(),
      breaks = y_breaks,
      limits = y_limits,
      labels = scales::label_number(accuracy = 0.01)
    )
  
  p
  
}


#' Plot class turnover stats
#' @param evolution A list, the output from evolve_mplat.
#' @param y_limits A vector of two numerics. The y axis limits. 
#' @param y_breaks A vector of numerics. The y axis tick mark locations.
#' @param step A string, specifying the kind of steps to take through the
#'   evolution: by generation or by change.
plot_turnover = function(
    evolution,
    y_limits = c(0, 200),
    y_breaks = c(0,2,5,20,50,100,200),
    step = c("generation", "change")
) {
  
  dat <- evolution$stats
  
  if (substr(step[1],1,1) == "c") {
    # Use changes rather than generations
    dat <- dat %>% select(-generation) %>% distinct()
    dat$generation <- 1:nrow(dat)
    x_label <- "Change"
  } else {
    x_label <- "Cycle"
  }
  
  pvt_type <- substr(evolution$model$pivot_type[1], 1, 1)
  if (pvt_type == "l") {
    y_label <- "Zone turnover" 
  } else {
    y_label <- "Class turnover" 
  }

  plot_dat <-
    dat %>% 
    rename(value = turnover) %>% 
    mutate(value = ifelse(is.na(value), lag(value), value)) %>%
    compile_mean_hi_lo() %>%
    rename(turnover = mean_value)

  p <- ggplot(plot_dat)
  
  if (evolution$n_rep > 1) {
    p <- 
      p +
      geom_ribbon(
        aes(x = generation, ymin = lower, ymax = upper),
        alpha = 0.3, colour = NA) 
  } 
  
  p <- 
    p + 
    geom_line(aes(x = generation, y = turnover), size = 1) +
    labs(linetype = "", x = x_label, y = y_label) +
    scale_y_continuous(
      trans = loglike2_trans(),
      limits = y_limits,
      breaks = y_breaks,
      labels = scales::label_number(accuracy = 1)
    )
  
  p
}


#' Plot top two classes stats
#' @param evolution A list, the output from evolve_mplat.
#' @param step A string, specifying the kind of steps to take through the
#'   evolution: by generation or by change.
plot_top2classes = function(
    evolution,
    step = c("generation", "change")
) {
  
  dat <- evolution$stats
  
  if (substr(step[1],1,1) == "c") {
    # Use changes rather than generations
    dat <- dat %>% select(-generation) %>% distinct()
    dat$generation <- 1:nrow(dat)
    x_label <- "Change"
  } else {
    x_label <- "Cycle"
  }
  
  pvt_type <- substr(evolution$model$pivot_type[1], 1, 1)
  if (pvt_type == "l") {
    y_label <- "Largest two zones" 
  } else {
    y_label <- "Largest two classes" 
  }
  
  # Handle different widths of the columns
  gens <- sort(unique(dat$generation))
  n_steps <- length(gens)
  gen_width <- gens[c(2,2:n_steps)] - gens[c(1,1:(n_steps-1))]
  final_width_diff <- gen_width[n_steps - 1] - gen_width[n_steps]
  gen_width[n_steps] <- gen_width[n_steps] - final_width_diff
  
  dat <-
    dat %>% 
    group_by(generation) %>%
    summarise(
      largest = mean(size_class1, na.rm = TRUE),
      second = mean(size_class2, na.rm = TRUE),
    ) %>%
    pivot_longer(cols = c(largest, second), 
                 names_to = "class", 
                 values_to = "size")
  
  plot_dat <- dat %>% select(generation, class, size)
  
  p <- 
    plot_dat %>%
    ggplot() +
    geom_col(aes(x = generation, y = size, fill = class), 
             width = rep(gen_width, each = 2), alpha = 0.6) +
    labs(linetype = "", x = x_label, y = y_label) +
    scale_fill_grey() + guides(fill="none")
  
  p
}


#' A log-like transformation for data that goes to 0
loglike_trans = function() {
  scales::trans_new(
    name = "loglike",
    transform = function(x) log(x + 0.1),
    inverse = function(x) exp(x) - 0.1)
}


#' A log-like transformation for data that goes to 0
loglike2_trans = function() {
  scales::trans_new(
    name = "loglike2",
    transform = function(x) log(x + 1),
    inverse = function(x) exp(x) - 1)
}


#' Get mean and ribbon edges for plotting, for the column "value"
#' @param df A dataframe with columns generation and value
#' @return A dataframe with columns generation, upper, lower & mean_value
compile_mean_hi_lo = function(df) {
  df %>% 
  group_by(generation) %>%
  # Get the overall mean, and for the ribbon, the means of the top and 
  # bottom 10%
  summarise(
    n_10pc = ceiling(n()/10),
    upper = mean(tail(sort(value), n_10pc)),
    lower = mean(head(sort(value), n_10pc)),
    top = max(value, na.rm = TRUE),
    bottom = min(value, na.rm = TRUE),
    mean_value = mean(value)
  )  %>%
  # Replace NAs
  mutate(
    upper = ifelse(is.na(upper), lag(upper), upper),
    lower = ifelse(is.na(lower), lag(lower), lower),
    mean_value = ifelse(is.na(mean_value), lag(mean_value), mean_value)
  ) %>%
  # Smooth
  mutate(
    upper =
      predict(loess(upper ~ generation, span=0.10)) %>%
      pmax(bottom) %>%
      pmin(top),
    lower =
      predict(loess(lower ~ generation, span=0.10)) %>%
      pmax(bottom) %>%
      pmin(top)
  ) %>%
  select(-top, -bottom, -n_10pc)
}
