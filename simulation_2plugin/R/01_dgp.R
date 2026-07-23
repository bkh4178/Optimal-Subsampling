# 01_dgp.R
# DGP 생성: p=9 혼합 공변량, 등/이분산 지원
#
# X = cbind(1, Z1,...,Z9),  col 인덱스:
#   1=intercept, 2=Z1, 3=Z2, 4=Z3, 5=Z4, 6=Z5, 7=Z6, 8=Z7, 9=Z8, 10=Z9
#
# 이분산: sigma(Z) = 0.5 + |0.5*Z1 + 0.5*Z2 + Z6 + Z8|
# 등분산: sigma(Z) = 1

  generate_data_final <- function(N, beta,
                                  heteroscedastic = TRUE, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    if (length(beta) != 10L) stop("beta must have length 10 (intercept + 9).")

    Z1 <- rnorm(N)
    Z2 <- rnorm(N)
    Z3 <- rnorm(N)
    Z4 <- rbinom(N, 1L, 0.5) - 0.5
    Z5 <- rbinom(N, 1L, 0.5) - 0.5
    Z6 <- rbinom(N, 1L, 0.5) - 0.5
    Z7 <- runif(N, -1, 1)
    Z8 <- runif(N, -1, 1)
    Z9 <- runif(N, -1, 1)

    X <- cbind(1, Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9)
    colnames(X) <- c("intercept", paste0("Z", 1:9))

    if (isTRUE(heteroscedastic)) {
      s_var     <- 0.5 * Z1 + 0.5 * Z2 + Z6 + Z8
      sigma_vec <- 0.5 + abs(s_var)
    } else {
      s_var     <- rep(0, N)
      sigma_vec <- rep(1, N)
    }
    sigma2_vec <- sigma_vec^2

    mu <- as.numeric(X %*% beta)
    y  <- mu + rnorm(N, 0, sigma_vec)

    A_hat     <- crossprod(X) / N
    A_hat_inv <- solve(A_hat)
    ell       <- rowSums((X %*% A_hat_inv) * X)

    list(
      X = X, y = y, mu = mu,
      s_var      = as.numeric(s_var),
      sigma_vec  = as.numeric(sigma_vec),
      sigma2_vec = as.numeric(sigma2_vec),
      ell        = as.numeric(ell),
      A_hat = A_hat, A_hat_inv = A_hat_inv,
      beta = beta, heteroscedastic = isTRUE(heteroscedastic)
    )
  }
