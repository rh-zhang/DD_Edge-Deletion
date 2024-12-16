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
      if (runif(1) <= r) {
        matrix[vertex, duplication] <- 1
        matrix[duplication, vertex] <- 1
      }
      edges <- which(matrix[vertex, ] != 0)
      for (e in edges) {
        if (runif(1) <= p) {
          matrix[vertex, e] <- 0
          if (!directed) {
            matrix[e, vertex] <- 0
          }
        }
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
  `0.2` = c(0, 0.2, 0.6, 0.6, 0.8)
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
run_simulation <- function(n_steps, matrix_type, q, p, r) {
  Model_consgrowth_g(n_steps, matrix_type, p, q, r)
}

# Parallel setup
cl <- makeCluster(6)
registerDoParallel(cl)
clusterExport(cl, list("Model_consgrowth_g", "run_simulation", "matrix_edge", "matrix_v1e1", "matrix_t1"))

# Simulations
for (rep in 1:10) {
  for (q in q_values) {
    p_values <- p_values_map[[as.character(q)]]
    for (r in r_values) {
      for (p in p_values) {
        steps <- ifelse(r == 0, n_steps[1], n_steps[2])
        variable_name <- paste0("cg_e1_r", r, "q", gsub("\\.", "", as.character(q * 10)), "p", gsub("\\.", "", as.character(p)), "_", rep)
        result <- run_simulation(steps, matrix_edge, q, p, r)
        assign(variable_name, result, envir = .GlobalEnv)
      }
    }
  }
}

stopCluster(cl)

# Group data by q and r values
generate_group_data <- function(prefix, q_values, r_values, p_values_map, num_reps = 10) {
  all_group_data <- list()
  
  for (q in q_values) {
    q_str <- gsub("\\.", "", as.character(q * 10))
    colors <- setNames(p_values_map[[as.character(q)]], c("Red", "Orange", "Green", "Blue", "Purple"))
    
    for (r in r_values) {
      r_group_data <- list()
      for (color in names(colors)) {
        p <- colors[[color]]
        base_name <- paste0(prefix, "_r", r, "q", q_str, "p", gsub("\\.", "", as.character(p)))
        suffixes <- paste0("_", seq_len(num_reps))
        variable_names <- paste0(base_name, suffixes)
        vectors <- mget(variable_names, envir = .GlobalEnv, ifnotfound = NA)
        vectors <- Filter(Negate(is.null), lapply(vectors, function(x) if (!is.null(x) && !is.na(x)) x[[3]] else NULL))
        if (length(vectors) >= 3) {
          r_group_data[[color]] <- do.call(calculate_summary, vectors)
        }
      }
      group_name <- paste0("group_data_", prefix, "_r", r, "q", q_str)
      all_group_data[[group_name]] <- r_group_data
    }
  }
  
  return(all_group_data)
}

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

# Prepare plot data
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
create_plot <- function(data, title, pho_0, q_value) {
  plot_config <- list(
    `0.9` = list(
      p_labels = c("p=0", "p=0.2", "p=2 - sqrt(22)/3", "p=0.6", "p=0.8"),
      linetypes = c("Red" = "dotted", "Orange" = "dotted", "Green" = "dotted", "Blue" = "solid", "Purple" = "solid")
    ),
    `0.7` = list(
      p_labels = c("p=0", "p=2-sqrt(182)/7", "p=0.4", "p=0.6", "p=0.8"),
      linetypes = c("Red" = "dotted", "Orange" = "dotted", "Green" = "solid", "Blue" = "solid", "Purple" = "solid")
    ),
    `0.5` = list(
      p_labels = c("p=0", "p=0.2", "p=0.4", "p=0.6", "p=0.8"),
      linetypes = c("Red" = "solid", "Orange" = "solid", "Green" = "solid", "Blue" = "solid", "Purple" = "solid")
    )
  )
  
  config <- plot_config[[as.character(q_value)]]
  
  ggplot(data, aes(x = x, y = avg, color = group, fill = group)) +
    geom_line(aes(linetype = group), size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    scale_color_manual(values = c("Red" = "red", "Orange" = "orange", "Green" = "green", "Blue" = "blue", "Purple" = "purple"), labels = config$p_labels) +
    scale_fill_manual(values = c("Red" = "red", "Orange" = "orange", "Green" = "green", "Blue" = "blue", "Purple" = "purple"), labels = config$p_labels) +
    scale_linetype_manual(values = config$linetypes, labels = config$p_labels) +
    geom_hline(yintercept = pho_0, color = "black", size = 0.8) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = title, x = "Number of Steps", y = "Proportion of Isolated Vertices") +
    theme_minimal()
}

# Generated combined plot
generate_combined_plot <- function(data_names, pho_values, q_value) {
  title_prefix <- ifelse(grepl("v1e1", data_names[[1]]), "Model A starting with 1 edge and 1 vertex", "Model A starting with 1 edge")
  title_prefix <- paste0(title_prefix, ", q=", q_value)
  
  data_groups <- lapply(data_names, function(name) mget(name, envir = .GlobalEnv))
  
  plots <- lapply(seq_along(data_groups), function(i) {
    create_plot(data = prepare_plot_data(data_groups[[i]]), title = paste0("r=", c(0, 0.5, 1)[i]), pho_0 = pho_values[i], q_value = q_value)
  })
  
  combined_plot <- Reduce(`|`, plots) + plot_annotation(title = title_prefix)
  print(combined_plot)
  return(combined_plot)
}

# Output when q=0.9
data_names_q9 <- c("group_data_e1_r0q9", "group_data_e1_r05q9", "group_data_e1_r1q9")
pho_values_q9 <- c(1, 2/11, 0.1)
generate_combined_plot(data_names_q9, pho_values_q9, q_value = 0.9)

#Output when q=0.7
data_names_q7 <- c("group_data_e1_r0q9", "group_data_e1_r05q9", "group_data_e1_r1q9")
pho_values_q9 <- c(1, 6/13, 0.3)
generate_combined_plot(data_names_q9, pho_values_q9, q_value = 0.9)

# Output when q=0.5
data_names_q9 <- c("group_data_e1_r0q9", "group_data_e1_r05q9", "group_data_e1_r1q9")
pho_values_q9 <- c(1, 2/3, 0.5)
generate_combined_plot(data_names_q9, pho_values_q9, q_value = 0.9)

# Output when q=0.2
data_names_q9 <- c("group_data_e1_r0q9", "group_data_e1_r05q9", "group_data_e1_r1q9")
pho_values_q9 <- c(1, 8/9, 0.8)
generate_combined_plot(data_names_q9, pho_values_q9, q_value = 0.9)
