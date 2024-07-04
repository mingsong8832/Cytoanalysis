#' Helper Function
#'
#' This is a helper function used by the main function.
#' @export
downsample <- function(data, n = 1e6, samples = names(data)) {
  tmp_data <- data[sapply(names(data), function(x) x %in% samples)]
  n <- n / length(tmp_data) # draw n random cells per sample
  tmp_sample <- lapply(tmp_data, function(x) {
    tmp_mat <- x[["data"]]
    if (nrow(tmp_mat) < n) {
      n <- nrow(tmp_mat)
    }
    rows <- sample(seq_len(nrow(tmp_mat)), size = n)
    tmp_mat <- tmp_mat[rows, ]
  })
  tmp_sample <- do.call(rbind, tmp_sample)
  return(tmp_sample)
}

#' @export
compute_som <- function(data, n_cells = 100) {
  require(kohonen)
  tmp_data <- data
  tmp_mat <- tmp_data$data
  # to define grid:each node to contain c. n_cells
  dimno <- round(sqrt(nrow(tmp_mat)) / sqrt(n_cells), 0)
  grd <- somgrid(xdim = dimno, ydim = dimno, topo = "hexagonal")
  tmp_data[['som']] <- som(tmp_data$data, grid = grd, keep.data = TRUE)
  tmp_data$data <- NULL
  return(tmp_data)
}

#' @export
assign_clusters <- function(data, clusters) {
  tmp_data <- data
  if  (exists("som", where = tmp_data)) {
    tmp_df <- tmp_data$som$unit.classif
  } else if (exists("classes", where = tmp_data)) {
    tmp_df <- tmp_data$classes
  } else {
    stop("data do not contain classification information\n")
  }
  tmp_name <- nrow(clusters) # get no of clusters
  tmp_name
  tmp_df <- pbapply::pblapply(tmp_df, function(x) {
    tmp <- clusters[x, , drop = FALSE]
    tmp <- cbind(x, tmp)
  }
  )
  tmp_df <- do.call(rbind, tmp_df)
  colnames(tmp_df) <- c(tmp_name, colnames(clusters))
  tmp_data$classes <- tmp_df
  return(tmp_data)
}
map_som <- function(data, trained, n_subset) {
  tmp_data <- data
  if (n_subset < nrow(tmp_data$data)) {
    tmp_mat <- tmp_data$data[sample.int(nrow(tmp_data$data), size = n_subset), ]
  } else {
    tmp_mat <- tmp_data$data
  }
  tmp_som <- kohonen::map(trained, tmp_mat) # map data to trained SOM
  tmp_data[["data"]] <- tmp_mat
  tmp_data[["classes"]] <- tmp_som$unit.classif
  return(tmp_data)
}

#' @export
assign_clusters <- function(data, clusters) {
  tmp_data <- data
  if  (exists("som", where = tmp_data)) {
    tmp_df <- tmp_data$som$unit.classif
  } else if (exists("classes", where = tmp_data)) {
    tmp_df <- tmp_data$classes
  } else {
    stop("data do not contain classification information\n")
  }
  tmp_name <- nrow(clusters) # get no of clusters
  tmp_name
  tmp_df <- pbapply::pblapply(tmp_df, function(x) {
    tmp <- clusters[x, , drop = FALSE]
    tmp <- cbind(x, tmp)
  }
  )
  tmp_df <- do.call(rbind, tmp_df)
  colnames(tmp_df) <- c(tmp_name, colnames(clusters))
  tmp_data$classes <- tmp_df
  return(tmp_data)
}

#' @export
count_observations <- function(data, clusters) {
  tmp_data <- data
  tmp_classes <- data$classes
  tmp_classes <- tmp_classes[, c(colnames(tmp_classes)[1], clusters)]
  tmp_counts <- lapply(colnames(tmp_classes), function(n) {
    tmp_cl <- tmp_data$classes[, n]
    tmp_mat <- sapply(seq_len(as.numeric(n)), function(x) {
      ct <- length(tmp_cl[tmp_cl == x])
    }
    )
    tmp_mat <- as(tmp_mat, "matrix")
    colnames(tmp_mat) <- data$name
    return(tmp_mat)
  }
  )
  names(tmp_counts) <- colnames(tmp_classes)
  tmp_data$counts <- tmp_counts
  return(tmp_data)
}

# function: combine count tables for all samples
#' @export
get_counts <- function(data) {
  tmp_cl_counts <- lapply(data, "[[", "counts")
  tmp_cl <- unique(unlist(lapply(tmp_cl_counts, names)))
  tmp_counts <- lapply(tmp_cl, function(n) {
    tmp_mat <- lapply(tmp_cl_counts, "[[", n)
    tmp_mat <- do.call(cbind, tmp_mat)
  }
  )
  names(tmp_counts) <- tmp_cl
  return(tmp_counts)
}
