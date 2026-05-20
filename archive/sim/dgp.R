# dgp.R
# Data Generating Process — Phase 1 (저차원 검증)
#
# 지원 차원: p = 2, 3
#   p=2: X1 ~ N(0,1), X2 ~ Bernoulli(0.5) centered
#   p=3: X1 ~ N(0,1), X2 ~ Bernoulli(0.5) centered, X3 ~ N(0,1)
#
# Cases:
#   Case 1 — 등분산         : sigma2 = 1
#   Case 2 — 연속형 이분산  : sigma2_raw = exp(gamma * X1)
#   Case 3 — 이진형 이분산  : sigma2_raw = 1 + 2*G  (G = X2 uncentered binary)
#
# 공통: sigma2 = sigma2_raw / mean(sigma2_raw)  → 평균 분산 = 1로 정규화
# 분산 구조는 leverage와 무관 (X 자체 값에만 의존)

generate_data <- function(N = 10000, p = 2, case = 1, gamma = 1.0, beta0 = NULL) {

  if (!p %in% c(2, 3)) stop("p must be 2 or 3.")

  case <- as.character(case)
  if (!case %in% c("1", "2", "3")) stop("case must be 1, 2, or 3.")

  # beta0 기본값: p=2이면 (1.0, 2.0), p=3이면 (0.5, 1.0, 2.0)
  if (is.null(beta0)) {
    beta0 <- if (p == 2) c(1.0, 2.0) else c(0.5, 1.0, 2.0)
  }
  if (length(beta0) != p) stop("length(beta0) must equal p.")

  # ── Step 1. 공변량 생성 ────────────────────────────────────────────────────

  X1 <- rnorm(N)

  # Bernoulli: raw(uncentered)는 G로 저장, centered는 회귀에 사용
  G       <- rbinom(N, 1, 0.5)        # uncentered: Case 3 분산 구조에 사용
  X2      <- G - mean(G)              # centered: 회귀 공변량

  if (p == 2) {
    X <- cbind(X1, X2)
  } else {
    X3 <- rnorm(N)
    X  <- cbind(X1, X2, X3)
  }

  # ── Step 2. A_hat 및 leverage score ───────────────────────────────────────

  A_hat     <- crossprod(X) / N       # p x p
  A_hat_inv <- solve(A_hat)

  # ell_i = x_i' A_hat^{-1} x_i  (methods.R에서 샘플링 확률 계산에 사용)
  ell <- rowSums((X %*% A_hat_inv) * X)

  # ── Step 3. 오차 분산 설정 ────────────────────────────────────────────────

  sigma2_raw <- switch(case,
    "1" = rep(1, N),
    "2" = 1 + gamma * abs(X1),            # X1(Normal) 기준 연속형 이분산
    "3" = 1 + 2 * G                   # G(binary) 기준 group-wise 이분산
  )

  # 케이스 간 평균 분산을 1로 통일 → 비교 가능
  sigma2 <- sigma2_raw / mean(sigma2_raw)

  # ── Step 4. 반응변수 생성 ─────────────────────────────────────────────────

  mu      <- as.vector(X %*% beta0)
  epsilon <- rnorm(N, mean = 0, sd = sqrt(sigma2))
  y       <- mu + epsilon

  # ── 반환 ──────────────────────────────────────────────────────────────────

  list(
    X         = X,
    y         = y,
    mu        = mu,           # true mean (excess risk 계산용)
    sigma2    = sigma2,       # true sigma2 (oracle method, metrics용)
    ell       = ell,          # leverage score (methods.R 샘플링 확률용)
    A_hat     = A_hat,        # 기록용
    beta0     = beta0,        # true coefficient (metrics용)
    G         = G,            # X2 uncentered binary (diagnostic용)
    p         = p,
    case      = case,
    gamma     = gamma
  )
}