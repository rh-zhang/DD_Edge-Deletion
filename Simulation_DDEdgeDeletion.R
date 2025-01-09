library(igraph)
library(parallel)
library(doParallel)
library(scales)
library(ggplot2)
library(patchwork)

# Model A
Model_consgrowth_g <- function(size, g0, p, q, r, directed = FALSE) {
  g1 <- graph_from_adjacency_matrix(g0, mode = "undirected")
  Nm <- gorder(g1)
  Expected <- vector()
  Expected[1] <- Nm
  Nm0 <- vector()
  Nm0[1] <- sum(degree(g1) == 0)
  Nm1 <- vector()
  Nm1[1] <- sum(degree(g1) == 1)
  Proportion <- vector()
  Proportion[1] <- Nm0[1] / Expected[1]
  ExpectedProportion <- vector()
  ExpectedProportion[1] <- Proportion[1]
  Proportionk <- vector("list", size)
  Nmk <- vector("list", size - 1)
  i <- 1
  vertex <- Nm + 1
  matrix <- cbind(g0, 0)
  matrix <- rbind(matrix, 0)
  
  while (i <= size) {
    prob <- runif(1, 0, 1)
    if (prob <= q) {
      duplication <- sample(1:(vertex - 1), 1)
      matrix[vertex, ] <- matrix[duplication, ]
      matrix[, vertex] <- matrix[vertex, ]
      edges <- which(matrix[vertex, ] != 0)
      for (e in edges) {
        if (runif(1) <= p) {
          matrix[vertex, e] <- 0
          if (!directed) {
            matrix[e, vertex] <- 0
          }
        }
      }
      if (runif(1) <= r) {
        matrix[vertex, duplication] <- 1
        matrix[duplication, vertex] <- 1
      }
    } else {
      deletion <- sample(1:(vertex - 1), 1)
      matrix[deletion, ] <- 0
      matrix[, deletion] <- 0
      matrix[vertex, ] <- 0
      matrix[, vertex] <- 0
    }
    
    g1 <- graph_from_adjacency_matrix(matrix, mode = "undirected")
    Nm0[i + 1] <- sum(degree(g1) == 0)
    Nm1[i + 1] <- sum(degree(g1) == 1)
    Nm <- gorder(g1)
    Proportion[i + 1] <- Nm0[i + 1] / Nm
    
    dupvalsum <- 0
    for (k in 2:size) {
      Nmk[[k - 1]][i] <- sum(degree(g1) == k)
      Proportionk[[k]][i] <- Nmk[[k - 1]][i] / Nm
      dupvalsum <- dupvalsum + Nmk[[k - 1]][i] * r * p^k + (1 - r) * p^(k - 1)
    }
    
    matrix <- cbind(matrix, 0)
    matrix <- rbind(matrix, 0)
    vertex <- nrow(matrix)
    i <- i + 1
  }
  
  return(list(Nm0, Nm1, Nmk, Proportion, Proportionk))
}

# Define parameters
q_values <- c(0.9, 0.7, 0.5, 0.2)
p_values_map <- list(
  `0.9` = c(0, 0.2, 2 - sqrt(22) / 3, 0.6, 0.8),
  `0.7` = c(0, 2 - sqrt(182) / 7, 0.4, 0.6, 0.8),
  `0.5` = c(0, 0.2, 0.4, 0.6, 0.8),
  `0.2` = c(0, 0.2, 0.4, 0.6, 0.8)
)
r_values <- c(0, 0.5, 1)
n_steps <- c(1500, 1000)

# Define matrices
matrix_edge <- matrix(0, nrow = 2, ncol = 2)
matrix_edge[1, 2] <- 1
matrix_edge[2, 1] <- 1

matrix_v1e1 <- matrix(0, nrow = 3, ncol = 3)
matrix_v1e1[1, 2] <- 1
matrix_v1e1[2, 1] <- 1

matrix_t1 <- matrix(1, nrow = 3, ncol = 3)
diag(matrix_t1) <- 0

# Function to run simulations
run_simulation <- function(n_steps, matrix_type, p, q, r) {
  Model_consgrowth_g(n_steps, matrix_type, p, q, r)
}

# Parallel setup
cl <- makeCluster(6)
registerDoParallel(cl)
clusterExport(cl, list("Model_consgrowth_g", "run_simulation", "matrix_edge", "matrix_v1e1", "matrix_t1"))

# Simulations
for (rep in 1:30) {
  for (q in q_values) {
    p_values <- p_values_map[[as.character(q)]]
    for (r in r_values) {
      for (p in p_values) {
        steps <- ifelse(r == 0 && q==0.9, n_steps[1], n_steps[2])
        # starting with an edge
        variable_name <- paste0("cg_e1_r", r, "q", gsub("\\.", "", as.character(q * 10)), "p", gsub("\\.", "", as.character(p)), "_", rep)
        result <- run_simulation(steps, matrix_edge, p, q, r)
        assign(variable_name, result, envir = .GlobalEnv)
        # starting with a vertex and edge
        variable_name <- paste0("cg_v1e1_r", r, "q", gsub("\\.", "", as.character(q * 10)), "p", gsub("\\.", "", as.character(p)), "_", rep)
        result <- run_simulation(steps, matrix_v1e1, p, q, r)
        assign(variable_name, result, envir = .GlobalEnv)
      }
    }
  }
}

stopCluster(cl)


# Calculate average and 2 standard errors across repetitions
calculate_summary <- function(..., level = 0.95) {
  data <- data.frame(...)
  avg <- rowMeans(data)
  n <- ncol(data)
  stderr <- apply(data, 1, sd) / sqrt(n)
  error <- 2 * stderr
  lower <- avg - error
  upper <- avg + error
  data.frame(avg, lower, upper)
}

# Function to convert group data into a data frame for ggplot
prepare_plot_data <- function(group_data) {
  do.call(rbind, lapply(names(group_data), function(group) {
    data.frame(
      x = 1:nrow(group_data[[group]]),
      avg = group_data[[group]]$avg,
      lower = group_data[[group]]$lower,
      upper = group_data[[group]]$upper,
      group = factor(group, levels = c("Red", "Orange", "Green", "Blue", "Purple"))
    )
  }))
}

# Create sub-plot
create_plot <- function(
    data,
    title,
    q_value,
    r_value,
    p_values_map,
    tol = 1e-7
) {
  ########################################
  # 1) LINES + RIBBONS (No Sorting)
  ########################################
  # data has columns: x, avg, lower, upper, group ∈ {Red,Orange,Green,Blue,Purple}
  # We'll color them in the same order (no final-avg sorting).
  
  # (A) Color + linetype logic for lines
  line_colors <- c(
    "Red"    = "red",
    "Orange" = "orange",
    "Green"  = "green",
    "Blue"   = "blue",
    "Purple" = "purple"
  )
  ########################################
  if (q_value == 0.9) {
    if (0.9 > 1/(2*(1 - r_value))){
      line_types <- c(
        "Purple"="dotted",
        "Blue"  ="dotted",
        "Green" ="dotted",
        "Orange"="dotted",
        "Red"   ="dotted"
      )}else {
        line_types <- c(
          "Purple"="solid",
          "Blue"  ="solid",
          "Green" ="dotted",
          "Orange"="dotted",
          "Red"   ="dotted"
        )}
  } else if (q_value == 0.7) {
    if (0.7 > 1/(2*(1 - r_value))){
      line_types <- c(
        "Purple"="dotted",
        "Blue"  ="dotted",
        "Green" ="dotted",
        "Orange"="dotted",
        "Red"   ="dotted"
      )
    } else {
      line_types <- c(
        "Purple"="solid",
        "Blue"  ="solid",
        "Green" ="solid",
        "Orange"="dotted",
        "Red"   ="dotted"
      )
    }} else {
    # q=0.5 or 0.2 => all solid
    line_types <- c(
      "Red"="solid", "Orange"="solid",
      "Green"="solid", "Blue"="solid", "Purple"="solid"
    )
  }
  
  # Build the basic plot
  plt <- ggplot(data, aes(x = x, y = avg, color = group)) +
    geom_line(aes(linetype = group), size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = group),
                alpha = 0.2, color = NA)
  
  ########################################
  # 2) BARS on the right side
  ########################################
  # Descending p => Purple, Blue, Green, Orange, Red
  p_vals <- p_values_map[[as.character(q_value)]]
  p_vals_sorted <- sort(p_vals, decreasing = TRUE)
  color_order_bars <- c("Purple","Blue","Green","Orange","Red")
  
  if (length(p_vals_sorted) != 5) {
    stop("Expected exactly 5 p-values for bars. Found: ", length(p_vals_sorted))
  }
  
  # We'll place bars beyond max(data$x)
  x_max_data <- max(data$x, na.rm=TRUE)
  XBAR_BASE  <- x_max_data + 15
  bar_width  <- 30
  gap        <- 10
  
  # Compute each bar's uncapped vs final (capped at 1) range
  bar_data <- do.call(rbind, lapply(seq_along(p_vals_sorted), function(i) {
    p_i   <- p_vals_sorted[i]
    col_b <- color_order_bars[i]
    
    rho_0 <- (1 - q_value) / (1 - q_value*(1 - r_value))
    
    uncapped_rho_1 <- (3*(1 - q_value) + p_i*q_value*(1 - r_value)) /
      (2*(1 - q_value*(1 - r_value)))
    
    # Capped at 1
    final_rho_1 <- min(uncapped_rho_1, 1)
    
    uncapped_ymin <- min(rho_0, uncapped_rho_1)
    uncapped_ymax <- max(rho_0, uncapped_rho_1)
    final_ymin    <- min(rho_0, final_rho_1)
    final_ymax    <- max(rho_0, final_rho_1)
    
    xmin_i <- XBAR_BASE + (i - 1)*(bar_width + gap)
    xmax_i <- xmin_i + bar_width
    
    data.frame(
      bar_color     = col_b,
      uncapped_ymin = uncapped_ymin,
      uncapped_ymax = uncapped_ymax,
      final_ymin    = final_ymin,
      final_ymax    = final_ymax,
      xmin          = xmin_i,
      xmax          = xmax_i
    )
  }))
  
  # Helper for tolerance-based equality
  all_equal_with_tol <- function(x, tol) all(abs(x - x[1]) < tol)
  
  # Check if all bars are equal BEFORE capping
  uncapped_lens     <- bar_data$uncapped_ymax - bar_data$uncapped_ymin
  same_len_uncapped <- all_equal_with_tol(uncapped_lens,  tol)
  same_ymin_uncapped<- all_equal_with_tol(bar_data$uncapped_ymin, tol)
  same_ymax_uncapped<- all_equal_with_tol(bar_data$uncapped_ymax, tol)
  all_equal_uncapped<- (same_len_uncapped && same_ymin_uncapped && same_ymax_uncapped)
  
  # Check if all bars are equal AFTER capping
  final_lens   <- bar_data$final_ymax - bar_data$final_ymin
  same_len_final<- all_equal_with_tol(final_lens, tol)
  same_ymin_final<- all_equal_with_tol(bar_data$final_ymin, tol)
  same_ymax_final<- all_equal_with_tol(bar_data$final_ymax, tol)
  all_equal_final<- (same_len_final && same_ymin_final && same_ymax_final)
  
  # Decide how to plot
  if (all_equal_uncapped) {
    # Single magenta bar
    single_bar <- data.frame(
      bar_color = "MagentaSingle",
      xmin      = XBAR_BASE,
      xmax      = XBAR_BASE + bar_width,
      ymin      = bar_data$final_ymin[1],
      ymax      = bar_data$final_ymax[1]
    )
    plt <- plt + geom_rect(
      data = single_bar,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill  = "magenta",
      alpha = 0.6,
      inherit.aes = FALSE
    )
    
  } else if (all_equal_final) {
    # Single cyan bar
    single_bar <- data.frame(
      bar_color = "CyanSingle",
      xmin      = XBAR_BASE,
      xmax      = XBAR_BASE + bar_width,
      ymin      = bar_data$final_ymin[1],
      ymax      = bar_data$final_ymax[1]
    )
    plt <- plt + geom_rect(
      data = single_bar,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill  = "cyan",
      alpha = 0.6,
      inherit.aes = FALSE
    )
    
  } else {
    # Side-by-side bars in their final (capped) range
    bar_data$ymin <- bar_data$final_ymin
    bar_data$ymax <- bar_data$final_ymax
    
    plt <- plt + geom_rect(
      data = bar_data,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = bar_color),
      alpha = 0.6,
      inherit.aes = FALSE
    )
  }
  
  ########################################
  # 3) Manual color/fill/linetype scales
  ########################################
  # lines: group ∈ {Red,Orange,Green,Blue,Purple}
  # bars:  bar_color ∈ {Purple,Blue,Green,Orange,Red} or Single(Magenta/Cyan)
  plt <- plt +
    scale_color_manual(values = c(
      "Red"="red",    "Orange"="orange", "Green"="green",
      "Blue"="blue",  "Purple"="purple",
      "MagentaSingle"="magenta",
      "CyanSingle"   ="cyan"
    )) +
    scale_fill_manual(values = c(
      "Red"="red",    "Orange"="orange", "Green"="green",
      "Blue"="blue",  "Purple"="purple",
      "MagentaSingle"="magenta",
      "CyanSingle"   ="cyan"
    )) +
    scale_linetype_manual(values = line_types)
  
  ########################################
  # 4) Hide x-axis ticks, show title,
  #    "cap" y at [min(data), 1] with coord_cartesian(clip="off")
  ########################################
  x_min_val <- min(data$x, na.rm = TRUE)
  y_min_val <- min(c(data$avg, data$lower, data$upper), na.rm = TRUE)
  
  plt <- plt +
    coord_cartesian(
      xlim = c(x_min_val, x_max_data),
      ylim = c(y_min_val, 1),
      clip = "off"  # so if a line goes above 1, it remains visible
    ) +
    labs(
      title = title,
      x = NULL,  # or "Number of Steps" if you want a label
      y = "Proportion of Isolated Vertices"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      # Hide x-axis ticks and labels
      axis.ticks.x    = element_blank(),
      axis.text.x     = element_blank()
    )
  
  return(plt)
}


      
# Function to generate combined plot dynamically
generate_combined_plot_dynamic <- function(
    prefix,
    q_value,
    r_values = c(0, 0.5, 1),
    p_values_map,
    num_reps = 30
) {
  # Same as before
  generate_group_data_inline <- function(prefix, q_value, r_value, p_values_map, num_reps) {
    q_str <- gsub("\\.", "", as.character(q_value * 10))
    colors <- setNames(
      p_values_map[[as.character(q_value)]],
      c("Red","Orange","Green","Blue","Purple")
    )
    
    group_data <- list()
    for (color in names(colors)) {
      p <- colors[[color]]
      base_name <- paste0(prefix, "_r", r_value,
                          "q", q_str,
                          "p", gsub("\\.", "", as.character(p)))
      suffixes  <- paste0("_", seq_len(num_reps))
      variable_names <- paste0(base_name, suffixes)
      
      vectors <- mget(variable_names, envir = .GlobalEnv, ifnotfound = NA)
      if (any(sapply(vectors, is.na))) {
        stop("Missing variables for: ", base_name)
      }
      vectors <- lapply(vectors, function(x) {
        if (!is.null(x) && all(!is.na(x))) x[[4]] else NULL
      })
      vectors <- Filter(Negate(is.null), vectors)
      
      if (length(vectors) >= 2) {
        group_data[[color]] <- do.call(calculate_summary, vectors)
      } else {
        message("Insufficient data for ", base_name)
      }
    }
    return(group_data)
  }
  
  # gather data for each r
  group_data_list <- lapply(r_values, function(r_val) {
    generate_group_data_inline(prefix, q_value, r_val, p_values_map, num_reps)
  })
  
  # Title
  title_prefix <- ifelse(
    prefix == "cg_v1e1", "Model A starting with 1 edge and 1 vertex",
    ifelse(prefix == "cg_t1", "Model A starting with 1 triangle",
           "Model A starting with 1 edge")
  )
  title_prefix <- paste0(title_prefix, ", q=", q_value)
  
  # Create subplots
  plots <- lapply(seq_along(group_data_list), function(i) {
    create_plot(
      data         = prepare_plot_data(group_data_list[[i]]),
      title        = paste0("r=", r_values[i]),
      q_value      = q_value,
      r_value      = r_values[i],
      p_values_map = p_values_map,
      tol          = 1e-7
    )
  })
  
  # Combine via patchwork
  plots <- lapply(plots, function(p) {
    p + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"))
  })
  
  combined_plot <- Reduce(`|`, plots) + 
    plot_annotation(
      title = title_prefix,
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
  
  print(combined_plot)
  return(combined_plot)
}



# Output plot
prefixes = c("cg_e1","cg_v1e1")
for (prefix in prefixes) {
  for (q in q_values) {
    generate_combined_plot_dynamic(
      prefix       = prefix,
      q_value      = q,
      r_values     = r_values,
      p_values_map = p_values_map,
      num_reps     = 30
    )
  }
}
  
