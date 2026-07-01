# Generalized independence test
#
# Example:
# set.seed(1)
# X <- matrix(rnorm(100), ncol = 2)
# Y <- matrix(rnorm(100), ncol = 2)
# GIT_test(X, Y)
#
# Parameters:
# X, Y: numeric vectors, matrices, or data frames with the same number of rows.
#       Each row is one observation.
# k: number of nearest/farthest neighbors used to build the graph.
#    Default is floor(sqrt(nrow(X))).
# robustgraph: if TRUE, use the penalized robust graph construction.
#              If FALSE, use the unpenalized graph-rank construction.
# lambda: degree-regularization strength for the robust graph.
#         Larger values penalize graphs with uneven node degrees more strongly.
# symmetric: if TRUE, symmetrize the directed graph ranks by adding their
#            transpose. This matches the default method used in the paper code.
#
# Returns:
# A named numeric vector with the overall generalized p-value and four component
# p-values: overall, G1, G2, G3, and G4.

GIT_test <- function(X,
                     Y,
                     k = floor(sqrt(nrow(X))),
                     robustgraph = TRUE,
                     lambda = 0.3,
                     symmetric = TRUE) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  # Basic input checks. The test compares paired observations, so X and Y must
  # have the same sample size.
  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows.", call. = FALSE)
  }
  if (nrow(X) < 4) {
    stop("X and Y must contain at least 4 observations.", call. = FALSE)
  }
  if (!all(is.finite(X)) || !all(is.finite(Y))) {
    stop("X and Y must contain only finite values.", call. = FALSE)
  }

  n <- nrow(X)
  k <- as.integer(k)
  if (length(k) != 1 || is.na(k) || k < 1 || k >= n) {
    stop("k must be a single integer between 1 and nrow(X) - 1.", call. = FALSE)
  }
  if (length(lambda) != 1 || is.na(lambda) || lambda < 0) {
    stop("lambda must be a single non-negative number.", call. = FALSE)
  }

  # Build four graph-rank matrices:
  # RdX/RdY use farthest-neighbor ranks, and RsX/RsY use nearest-neighbor ranks.
  graphrank <- if (isTRUE(robustgraph)) {
    .obtain_robust_graph_rank(X, Y, k, lambda, symmetric = symmetric)
  } else {
    .obtain_graph_rank(X, Y, k, symmetric = symmetric)
  }

  Rxd0 <- graphrank$RdX
  Ryd0 <- graphrank$RdY
  Rxs0 <- graphrank$RsX
  Rys0 <- graphrank$RsY

  # Four component statistics combining nearest/farthest graph ranks from X
  # and Y. These are the quantities whose marginal p-values are reported.
  G1 <- sum(Rxd0 * Ryd0)
  G2 <- sum(Rxd0 * Rys0)
  G3 <- sum(Rxs0 * Ryd0)
  G4 <- sum(Rxs0 * Rys0)

  # Null means under random permutation of Y labels.
  G1_mean <- sum(Rxd0) * sum(Ryd0) / (n * (n - 1))
  G2_mean <- sum(Rxd0) * sum(Rys0) / (n * (n - 1))
  G3_mean <- sum(Rxs0) * sum(Ryd0) / (n * (n - 1))
  G4_mean <- sum(Rxs0) * sum(Rys0) / (n * (n - 1))

  # Reusable sums for closed-form variances and covariances.
  Rxd1 <- sum(Rxd0)
  Ryd1 <- sum(Ryd0)
  Rxd2 <- sum(Rxd0^2)
  Ryd2 <- sum(Ryd0^2)
  Rxd3 <- sum(colSums(Rxd0)^2)
  Ryd3 <- sum(colSums(Ryd0)^2)

  Rxs1 <- sum(Rxs0)
  Rys1 <- sum(Rys0)
  Rxs2 <- sum(Rxs0^2)
  Rys2 <- sum(Rys0^2)
  Rxs3 <- sum(colSums(Rxs0)^2)
  Rys3 <- sum(colSums(Rys0)^2)

  Rxds2 <- sum(Rxd0 * Rxs0)
  Ryds2 <- sum(Ryd0 * Rys0)
  Rxds3 <- sum(colSums(Rxd0) * colSums(Rxs0))
  Ryds3 <- sum(colSums(Ryd0) * colSums(Rys0))

  # Approximate variances of the four component statistics.
  G1_var <- .graph_stat_variance(Rxd1, Rxd2, Rxd3, Ryd1, Ryd2, Ryd3, G1_mean, n)
  G2_var <- .graph_stat_variance(Rxd1, Rxd2, Rxd3, Rys1, Rys2, Rys3, G2_mean, n)
  G3_var <- .graph_stat_variance(Rxs1, Rxs2, Rxs3, Ryd1, Ryd2, Ryd3, G3_mean, n)
  G4_var <- .graph_stat_variance(Rxs1, Rxs2, Rxs3, Rys1, Rys2, Rys3, G4_mean, n)

  # Approximate covariance matrix used by the generalized chi-square statistic.
  G12_cov <- .graph_stat_covariance(Rxd1, Rxd2, Rxd3, Ryd1, Rys1, Ryds2, Ryds3, n)
  G13_cov <- .graph_stat_covariance(Ryd1, Ryd2, Ryd3, Rxd1, Rxs1, Rxds2, Rxds3, n)
  G14_cov <- .graph_stat_covariance(Rxd1, Rxd2, Rxd3, Ryd1, Rys1, Ryds2, Ryds3, n,
                                    cross_x1 = Rxs1, cross_x2 = Rxds2, cross_x3 = Rxds3)
  G23_cov <- .graph_stat_covariance(Rxd1, Rxd2, Rxd3, Ryd1, Rys1, Ryds2, Ryds3, n,
                                    cross_x1 = Rxs1, cross_x2 = Rxds2, cross_x3 = Rxds3)
  G24_cov <- .graph_stat_covariance(Rys1, Rys2, Rys3, Rxd1, Rxs1, Rxds2, Rxds3, n)
  G34_cov <- .graph_stat_covariance(Rxs1, Rxs2, Rxs3, Ryd1, Rys1, Ryds2, Ryds3, n)

  covariance <- matrix(
    c(
      G1_var, G12_cov, G13_cov, G14_cov,
      G12_cov, G2_var, G23_cov, G24_cov,
      G13_cov, G23_cov, G3_var, G34_cov,
      G14_cov, G24_cov, G34_cov, G4_var
    ),
    byrow = TRUE,
    nrow = 4
  )

  centered_stats <- c(G1 - G1_mean, G2 - G2_mean, G3 - G3_mean, G4 - G4_mean)
  overall_stat <- t(centered_stats) %*% solve(covariance) %*% centered_stats

  # Report the overall generalized p-value first, followed by the four
  # component p-values.
  c(
    overall = as.numeric(1 - pchisq(overall_stat, df = 4)),
    G1 = .two_sided_normal_pvalue(G1, G1_mean, G1_var),
    G2 = .two_sided_normal_pvalue(G2, G2_mean, G2_var),
    G3 = .two_sided_normal_pvalue(G3, G3_mean, G3_var),
    G4 = .two_sided_normal_pvalue(G4, G4_mean, G4_var)
  )
}

.two_sided_normal_pvalue <- function(statistic, mean, variance) {
  as.numeric(2 * (1 - pnorm(abs((statistic - mean) / sqrt(variance)))))
}

.graph_stat_variance <- function(A1, A2, A3, B1, B2, B3, mean, n) {
  part1 <- 2 * A2 * B2 / (n * (n - 1))
  part2 <- 4 * (A3 - A2) * (B3 - B2) / (n * (n - 1) * (n - 2))
  part3 <- (A1^2 - 4 * A3 + 2 * A2) *
    (B1^2 - 4 * B3 + 2 * B2) /
    (n * (n - 1) * (n - 2) * (n - 3))
  part1 + part2 + part3 - mean^2
}

.graph_stat_covariance <- function(A1,
                                   A2,
                                   A3,
                                   B1,
                                   C1,
                                   BC2,
                                   BC3,
                                   n,
                                   cross_x1 = A1,
                                   cross_x2 = A2,
                                   cross_x3 = A3) {
  part1 <- 2 * cross_x2 * BC2 / (n * (n - 1))
  part2 <- 4 * (cross_x3 - cross_x2) * (BC3 - BC2) / (n * (n - 1) * (n - 2))
  part3 <- (A1 * cross_x1 - 4 * cross_x3 + 2 * cross_x2) *
    (B1 * C1 - 4 * BC3 + 2 * BC2) /
    (n * (n - 1) * (n - 2) * (n - 3))
  part1 + part2 + part3 - A1 * cross_x1 * B1 * C1 / (n^2 * (n - 1)^2)
}

.obtain_row_rank <- function(X, Y) {
  n <- nrow(X)

  # Pairwise Euclidean distances. Negative distances are used so that ordinary
  # rank() calls can represent nearest-neighbor orderings.
  DX <- as.matrix(dist(X))
  DY <- as.matrix(dist(Y))
  SX <- -DX
  SY <- -DY

  # Shift diagonals so self-distances do not enter the nearest/farthest
  # neighbor ranks.
  DX <- DX - 10 * diag(1, n, n)
  DY <- DY - 10 * diag(1, n, n)
  SX <- SX + 10 * diag(1, n, n)
  SY <- SY + 10 * diag(1, n, n)

  list(
    RdX = t(apply(DX, 1, rank)) - 1,
    RdY = t(apply(DY, 1, rank)) - 1,
    RsX = t(apply(SX, 1, rank)) - n * diag(1, n, n),
    RsY = t(apply(SY, 1, rank)) - n * diag(1, n, n)
  )
}

.obtain_graph_rank <- function(X, Y, k, symmetric = TRUE) {
  n <- nrow(X)
  ranks <- .obtain_row_rank(X, Y)

  # Keep only the top k nearest/farthest ranks; all other entries are set to 0.
  RdX <- ranks$RdX - n + 1 + k
  RdY <- ranks$RdY - n + 1 + k
  RsX <- ranks$RsX - n + 1 + k
  RsY <- ranks$RsY - n + 1 + k

  RdX[RdX < 0] <- 0
  RdY[RdY < 0] <- 0
  RsX[RsX < 0] <- 0
  RsY[RsY < 0] <- 0

  if (isTRUE(symmetric)) {
    RdX <- RdX + t(RdX)
    RdY <- RdY + t(RdY)
    RsX <- RsX + t(RsX)
    RsY <- RsY + t(RsY)
  }

  list(RdX = RdX, RdY = RdY, RsX = RsX, RsY = RsY)
}

.obtain_robust_graph_rank <- function(X, Y, k, lambda, symmetric = TRUE) {
  DX <- as.matrix(dist(X))
  DY <- as.matrix(dist(Y))

  # Convert distances so the same robust graph builder can be used for
  # farthest-neighbor and nearest-neighbor graphs.
  FX <- -DX + 2 * max(DX)
  FY <- -DY + 2 * max(DY)

  # Robust directed k-neighbor graphs for farthest and nearest relationships.
  XFNN <- .penalized_k_near_rank(FX, K = k, lambda = lambda)$trun_KNN
  YFNN <- .penalized_k_near_rank(FY, K = k, lambda = lambda)$trun_KNN
  XKNN <- .penalized_k_near_rank(DX, K = k, lambda = lambda)$trun_KNN
  YKNN <- .penalized_k_near_rank(DY, K = k, lambda = lambda)$trun_KNN

  # Mask distances by robust graph edges before converting them to graph ranks.
  XFNN_dist <- .rank_from_robust(XFNN) * FX
  YFNN_dist <- .rank_from_robust(YFNN) * FY
  XKNN_dist <- .rank_from_robust(XKNN) * DX
  YKNN_dist <- .rank_from_robust(YKNN) * DY

  XFNN_rank <- t(apply(-XFNN_dist, 1, rank))
  YFNN_rank <- t(apply(-YFNN_dist, 1, rank))
  XKNN_rank <- t(apply(-XKNN_dist, 1, rank))
  YKNN_rank <- t(apply(-YKNN_dist, 1, rank))

  XFNN_rank[XFNN_rank > k] <- 0
  YFNN_rank[YFNN_rank > k] <- 0
  XKNN_rank[XKNN_rank > k] <- 0
  YKNN_rank[YKNN_rank > k] <- 0

  if (isTRUE(symmetric)) {
    list(
      RdX = XFNN_rank + t(XFNN_rank),
      RdY = YFNN_rank + t(YFNN_rank),
      RsX = XKNN_rank + t(XKNN_rank),
      RsY = YKNN_rank + t(YKNN_rank)
    )
  } else {
    list(RdX = XFNN_rank, RdY = YFNN_rank, RsX = XKNN_rank, RsY = YKNN_rank)
  }
}

.rank_from_robust <- function(edge_connect) {
  # Convert an edge list into a 0/1 adjacency matrix.
  n <- max(edge_connect[, 1])
  robust_rank_matrix <- matrix(0, n, n)
  for (i in seq_len(nrow(edge_connect))) {
    robust_rank_matrix[edge_connect[i, 1], edge_connect[i, 2]] <- 1
  }
  robust_rank_matrix
}

.optimal_with_rank_current_node <- function(k, current_neighbors, neighbor, degree, lambda, row_rank) {
  # For one node, choose k outgoing neighbors that balance rank closeness and
  # degree regularity across the graph.
  degree[current_neighbors] <- degree[current_neighbors] - 1
  loss <- numeric(length(neighbor))

  for (i in seq_along(neighbor)) {
    new_connect <- rep(0, length(degree))
    new_connect[neighbor[i]] <- 1
    loss[i] <- row_rank[neighbor[i]] + lambda * sum((degree + new_connect - 2 * k)^2)
  }

  new_neighbors <- neighbor[order(loss)[seq_len(k)]]
  degree[new_neighbors] <- degree[new_neighbors] + 1

  list(new_neighbors = new_neighbors, degree = degree)
}

.degree_distribution <- function(G, sample_ids) {
  # Total degree counts both incoming and outgoing directed edges.
  degrees <- numeric(length(sample_ids))
  for (i in seq_along(sample_ids)) {
    degrees[i] <- sum(G[, 1] == sample_ids[i]) + sum(G[, 2] == sample_ids[i])
  }
  degrees
}

.out_direct <- function(K, nodes) {
  # Outgoing neighbors for each node in an edge list.
  out <- vector("list", length(nodes))
  for (i in nodes) {
    out[[i]] <- K[K[, 1] == i, 2]
  }
  out
}

.penalized_k_near_rank <- function(M, K = 5, lambda = 0.3) {
  # Start from the ordinary directed K-nearest-neighbor graph, then repeatedly
  # update each node's outgoing neighbors using the penalized rank objective.
  diag(M) <- NA
  Morder <- apply(M, 1, order)
  n <- nrow(M)
  Morder <- Morder[-n, ]
  KNN <- cbind(rep(seq_len(n), each = K), c(Morder[seq_len(K), ]))

  rank_matrix <- apply(M, 1, rank)
  out_nodes <- .out_direct(KNN, seq_len(n))
  degree <- .degree_distribution(KNN, seq_len(n))

  is_loop <- 1
  nodes <- sample(seq_len(n), n)
  id <- 1

  # Stop after a full pass through all nodes without changing any outgoing
  # neighbor set.
  while (is_loop <= n) {
    if (id > n) {
      id <- 1
      nodes <- sample(seq_len(n), n)
    }

    neighbor <- Morder[, nodes[id]]
    current_neighbors <- out_nodes[[nodes[id]]]
    op <- .optimal_with_rank_current_node(
      K,
      current_neighbors,
      neighbor,
      degree,
      lambda,
      rank_matrix[, nodes[id]]
    )

    if (setequal(current_neighbors, op$new_neighbors)) {
      is_loop <- is_loop + 1
      id <- id + 1
    } else {
      degree <- op$degree
      out_nodes[[nodes[id]]] <- op$new_neighbors
      is_loop <- 1
      id <- id + 1
    }
  }

  trun_KNN <- NULL
  for (i in seq_len(n)) {
    trun_KNN <- c(trun_KNN, rep(i, length(out_nodes[[i]])))
  }
  trun_KNN <- cbind(trun_KNN, unlist(out_nodes))
  colnames(trun_KNN) <- c("V1", "V2")

  list(trun_KNN = trun_KNN, degree = degree)
}
