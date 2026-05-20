res <- readRDS("sim/results_p2.rds")
subset(res, case == 2 & k == 500 & method == "OPT-hetero-oracle")[1:5, ]

source("sim/dgp.R")
source("sim/methods.R")

set.seed(1)
dat <- generate_data(N = 10000, p = 2, case = 2, gamma = 1.0, beta0 = c(1.0, 2.0))
res_oracle <- run_opt_hetero_oracle(dat, k = 500)

summary(res_oracle$pi)
quantile(res_oracle$pi, c(0.9, 0.95, 0.99, 1.0))

res_homo <- run_opt_homo(dat, k = 500)
summary(res_homo$pi)
quantile(res_homo$pi, c(0.9, 0.95, 0.99, 1.0))




res <- readRDS("sim/results_p2.rds")
sub <- res[res$case == 2 & res$k == 500, ]
agg_sd <- aggregate(excess_risk ~ method, data = sub, FUN = sd)
agg_n  <- aggregate(excess_risk ~ method, data = sub, FUN = length)
agg_sd$se <- agg_sd$excess_risk / sqrt(agg_n$excess_risk)
agg_sd



source("sim/dgp4.R")
source("sim/methods.R")

set.seed(1)
dat <- generate_data(N = 10000, p = 4, case = 1, gamma = 1.0, beta0 = c(0.5, 0.5, 1.0, 1.0))

# sigma2 확인
summary(dat$sigma2)

# score 비교
score_homo   <- sqrt(dat$ell)
score_oracle <- sqrt(dat$sigma2 * dat$ell)

all.equal(score_homo, score_oracle)




set.seed(1)
dat <- generate_data(N=10000, p=4, case=2, gamma=1.0, beta0=c(0.5,0.5,1.0,1.0))
cor(dat$sigma2, dat$ell)
var(dat$sigma2)



res <- readRDS("sim/results_p4.rds")
sub <- res[res$gamma == 1 & res$case == 2 & res$k == 500, ]
aggregate(c1_pi ~ method, data = sub, FUN = mean)
cor(dat$sigma2, dat$ell)
cor(sqrt(dat$sigma2 * dat$ell), sqrt(dat$ell))


sub <- res[res$gamma == 1 & res$case == 2 & res$k == 500, ]

aggregate(n_sub ~ method, data = sub, FUN = function(x) c(mean=mean(x), sd=sd(x)))
aggregate(clip_rate ~ method, data = sub, FUN = mean)





res <- readRDS("sim/results_lowdim_oracle.rds")
sub <- res[res$p==1 & res$gamma==1 & res$case==1 & res$k==500, ]
homo   <- sub[sub$method=="OPT-homo", "excess_risk"]
oracle <- sub[sub$method=="OPT-hetero-oracle", "excess_risk"]
mean(oracle) / mean(homo)   # 이렇게 해야 맞아
mean(oracle / homo)         # 현재 코드는 이렇게 하고 있어




res <- readRDS("sim/results_lowdim.rds")
sub <- res[res$k == 1000 & !res$fail, ]
homo   <- sub$excess_risk[sub$method == "OPT-homo"]
oracle <- sub$excess_risk[sub$method == "OPT-hetero-oracle"]
cat(sd(homo)/sqrt(length(homo)), "\n")
cat(sd(oracle)/sqrt(length(oracle)), "\n")
cat(mean(oracle) - mean(homo), "\n")
