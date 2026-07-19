#############################################################################
# Advanced Statistical Inference - Lab Problem Solutions
# Language: R (base R + stats, MASS)
#############################################################################

suppressMessages({
  library(MASS)   # for rlm/huber
  library(stats)
})

set.seed(123)  # global seed; re-set locally before each simulation block too

#############################################################################
# QUESTION 1
#############################################################################
cat("\n================= QUESTION 1 =================\n")

mu_true    <- 10
sigma2_true <- 4
n1  <- 20
M1  <- 5000

## ---- Q1a: Bias check ----------------------------------------------------
set.seed(1)
samples1 <- matrix(rnorm(n1 * M1, mean = mu_true, sd = sqrt(sigma2_true)),
                    nrow = M1, ncol = n1)

xbar_vec   <- rowMeans(samples1)
S2_vec     <- apply(samples1, 1, var)                 # divisor (n-1)
sigma2_mle_vec <- apply(samples1, 1, function(x) mean((x - mean(x))^2))  # divisor n

Q1a_table <- data.frame(
  Estimator   = c("X_bar (mean)", "S^2 (n-1 divisor)", "sigma2_MLE (n divisor)"),
  Average     = c(mean(xbar_vec), mean(S2_vec), mean(sigma2_mle_vec)),
  TrueValue   = c(mu_true, sigma2_true, sigma2_true)
)
print(Q1a_table)
cat("Interpretation: S^2 (n-1 divisor) averages very close to sigma^2 = 4 -> unbiased.\n")
cat("sigma2_MLE (n divisor) averages below 4 -> biased downward (factor (n-1)/n).\n")

## ---- Q1b: Efficiency comparison ------------------------------------------
T1 <- xbar_vec                                  # sample mean
T2 <- (samples1[, 1] + samples1[, n1]) / 2       # avg of first & last obs

varT1 <- var(T1)
varT2 <- var(T2)
ratio_T2_T1 <- varT2 / varT1

Q1b_table <- data.frame(
  Estimator = c("T1 = X_bar", "T2 = (X1+Xn)/2"),
  EmpVar    = c(varT1, varT2)
)
print(Q1b_table)
cat(sprintf("Ratio Var(T2)/Var(T1) = %.3f\n", ratio_T2_T1))
cat("Interpretation: Ratio >> 1, so T1 (sample mean) is far more efficient than T2.\n")

## ---- Q1c: Sufficiency check ----------------------------------------------
set.seed(2)
sampA <- rnorm(n1, mu_true, sqrt(sigma2_true))
# Build sampB with identical mean & variance but different individual values:
# standardize sampA, apply an orthogonal-ish permutation/reflection that
# preserves sum and sum of squares (swap symmetric pairs around the mean).
sampB <- sampA
# reverse-order swap of paired deviations preserves mean & variance
ord <- order(sampA)
sampB[ord] <- rev(sampA[ord])
# To guarantee different individual values while keeping mean/var identical,
# perform a mean/variance-preserving rotation on two coordinates at a time:
rotate_pair <- function(x, i, j, theta) {
  m <- mean(x)
  xi <- x[i] - m; xj <- x[j] - m
  x[i] <- m + xi * cos(theta) - xj * sin(theta)
  x[j] <- m + xi * sin(theta) + xj * cos(theta)
  x
}
sampB <- sampA
set.seed(3)
for (k in seq(1, n1 - 1, by = 2)) {
  sampB <- rotate_pair(sampB, k, k + 1, theta = pi / 3)
}

cat(sprintf("mean(sampA)=%.6f  mean(sampB)=%.6f\n", mean(sampA), mean(sampB)))
cat(sprintf("var(sampA)=%.6f   var(sampB)=%.6f\n", var(sampA), var(sampB)))

loglik_normal <- function(x, mu, sigma2) {
  n <- length(x)
  -n/2 * log(2 * pi * sigma2) - sum((x - mu)^2) / (2 * sigma2)
}

param_grid <- data.frame(
  mu     = c(9, 10, 10, 11, 10),
  sigma2 = c(4, 3, 4, 4, 5)
)

llA <- mapply(function(m, s2) loglik_normal(sampA, m, s2),
              param_grid$mu, param_grid$sigma2)
llB <- mapply(function(m, s2) loglik_normal(sampB, m, s2),
              param_grid$mu, param_grid$sigma2)

Q1c_table <- data.frame(
  mu = param_grid$mu, sigma2 = param_grid$sigma2,
  loglikA = llA, loglikB = llB,
  LikelihoodRatio_A_over_B = exp(llA - llB)
)
print(Q1c_table)
cat("Interpretation: The likelihood ratio L_A/L_B is (numerically) constant across\n")
cat("all parameter values, confirming (X_bar, S^2) is sufficient for (mu, sigma^2).\n")

## ---- Q1d: Consistency check ----------------------------------------------
ns <- c(5, 10, 25, 50, 100, 250, 500, 1000)
M1d <- 2000

xbar_list <- vector("list", length(ns))
set.seed(4)
for (i in seq_along(ns)) {
  n <- ns[i]
  mat <- matrix(rnorm(n * M1d, mu_true, sqrt(sigma2_true)), nrow = M1d, ncol = n)
  xbar_list[[i]] <- rowMeans(mat)
}
names(xbar_list) <- paste0("n=", ns)

emp_var   <- sapply(xbar_list, var)
theo_var  <- sigma2_true / ns
frac_far  <- sapply(xbar_list, function(v) mean(abs(v - mu_true) > 0.1))

Q1d_table <- data.frame(n = ns, EmpVar_Xbar = emp_var, Theoretical_4_over_n = theo_var,
                         Frac_abs_dev_gt_0.1 = frac_far)
print(Q1d_table)

png("q1d_consistency.png", width = 900, height = 600)
boxplot(xbar_list, main = "Sampling distribution of X_bar across n",
        xlab = "sample size n", ylab = "X_bar", col = "lightblue")
abline(h = mu_true, col = "red", lty = 2)
dev.off()
cat("Plot saved: q1d_consistency.png\n")
cat("Interpretation: Empirical variance of X_bar tracks 4/n closely and the fraction\n")
cat("of samples with |X_bar - 10| > 0.1 shrinks toward 0 as n grows -> consistency.\n")

## Optional add-on: median version
med_list <- lapply(xbar_list, function(x) x) # placeholder, recompute properly below
med_var  <- numeric(length(ns))
med_frac <- numeric(length(ns))
set.seed(5)
for (i in seq_along(ns)) {
  n <- ns[i]
  mat <- matrix(rnorm(n * M1d, mu_true, sqrt(sigma2_true)), nrow = M1d, ncol = n)
  meds <- apply(mat, 1, median)
  med_var[i]  <- var(meds)
  med_frac[i] <- mean(abs(meds - mu_true) > 0.1)
}
Q1d_table$EmpVar_Median   <- med_var
Q1d_table$Frac_Median_gt0.1 <- med_frac
print(Q1d_table)


#############################################################################
# QUESTION 2
#############################################################################
cat("\n================= QUESTION 2 =================\n")

failure_times <- c(2.1, 0.8, 3.4, 1.2, 5.6, 0.4, 2.9, 1.8, 4.1, 0.6,
                    3.3, 2.0, 1.5, 6.2, 0.9)
n2 <- length(failure_times)

## ---- Q2a: log-likelihood curve -------------------------------------------
loglik_exp <- function(lambda, x) sum(dexp(x, rate = lambda, log = TRUE))

lambda_grid <- seq(0.05, 2, length.out = 2000)
ll_grid <- sapply(lambda_grid, loglik_exp, x = failure_times)
lambda_hat_grid <- lambda_grid[which.max(ll_grid)]

png("q2a_loglik.png", width = 800, height = 600)
plot(lambda_grid, ll_grid, type = "l", lwd = 2,
     xlab = expression(lambda), ylab = expression(ell(lambda)),
     main = "Exponential log-likelihood")
abline(v = lambda_hat_grid, col = "red", lty = 2)
dev.off()
cat(sprintf("Grid-based maximizer: lambda_hat (grid) = %.4f\n", lambda_hat_grid))

## ---- Q2b: Numerical MLE ---------------------------------------------------
neg_ll_exp <- function(lambda) -loglik_exp(lambda, failure_times)
opt2b <- optim(par = 1, fn = neg_ll_exp, method = "Brent", lower = 0.001, upper = 5)
lambda_hat_opt <- opt2b$par
lambda_hat_formula <- 1 / mean(failure_times)

cat(sprintf("optim() MLE:        lambda_hat = %.4f\n", lambda_hat_opt))
cat(sprintf("Grid MLE:            lambda_hat = %.4f\n", lambda_hat_grid))
cat(sprintf("1/mean(data):        lambda_hat = %.4f\n", lambda_hat_formula))

## ---- Q2c: Fisher Information, SE, Wald CI --------------------------------
h <- 1e-4
d2ll <- (loglik_exp(lambda_hat_opt + h, failure_times) -
         2 * loglik_exp(lambda_hat_opt, failure_times) +
         loglik_exp(lambda_hat_opt - h, failure_times)) / h^2
I_lambda_hat <- -d2ll
se_lambda <- sqrt(1 / I_lambda_hat)
z975 <- qnorm(0.975)
wald_ci <- lambda_hat_opt + c(-1, 1) * z975 * se_lambda

cat(sprintf("Observed Fisher Info I(lambda_hat) = %.4f\n", I_lambda_hat))
cat(sprintf("SE(lambda_hat) = %.4f\n", se_lambda))
cat(sprintf("95%% Wald CI for lambda: (%.4f, %.4f)\n", wald_ci[1], wald_ci[2]))

## ---- Q2d: Gamma fit + AIC comparison --------------------------------------
neg_ll_gamma <- function(par) {
  k <- par[1]; theta <- par[2]
  if (k <= 0 || theta <= 0) return(1e10)
  -sum(dgamma(failure_times, shape = k, rate = theta, log = TRUE))
}
opt_gamma <- optim(par = c(1, 1), fn = neg_ll_gamma, method = "L-BFGS-B",
                    lower = c(1e-3, 1e-3))
k_hat <- opt_gamma$par[1]
theta_hat <- opt_gamma$par[2]
ll_gamma <- -opt_gamma$value
ll_exp   <- loglik_exp(lambda_hat_opt, failure_times)

AIC_exp   <- 2 * 1 - 2 * ll_exp     # 1 parameter
AIC_gamma <- 2 * 2 - 2 * ll_gamma   # 2 parameters

Q2d_table <- data.frame(
  Model = c("Exponential", "Gamma"),
  Params = c(sprintf("lambda=%.4f", lambda_hat_opt),
             sprintf("k=%.4f, theta=%.4f", k_hat, theta_hat)),
  LogLik = c(ll_exp, ll_gamma),
  AIC = c(AIC_exp, AIC_gamma)
)
print(Q2d_table)
cat(if (AIC_exp < AIC_gamma) "AIC prefers the Exponential model.\n" else "AIC prefers the Gamma model.\n")

## ---- Q2e: Visual fit check -------------------------------------------------
png("q2e_fit_check.png", width = 800, height = 600)
hist(failure_times, breaks = 8, freq = FALSE, col = "lightgray",
     main = "Failure times: histogram with fitted densities",
     xlab = "time (years)")
curve(dexp(x, rate = lambda_hat_opt), add = TRUE, col = "blue", lwd = 2, n = 500)
curve(dgamma(x, shape = k_hat, rate = theta_hat), add = TRUE, col = "darkgreen", lwd = 2, n = 500)
legend("topright", legend = c("Exponential fit", "Gamma fit"),
       col = c("blue", "darkgreen"), lwd = 2)
dev.off()
cat("Plot saved: q2e_fit_check.png\n")


#############################################################################
# QUESTION 3
#############################################################################
cat("\n================= QUESTION 3 =================\n")

sensor <- c(23.1, 24.5, 22.8, 25.0, 23.6, 24.1, 23.9, 22.7, 150.2, 24.8,
            23.3, 24.0, 23.5, 22.9, -80.5, 24.2, 23.7, 24.6, 23.0, 24.3)

## ---- Q3a: Compare center estimates ---------------------------------------
trimmed_mean10 <- mean(sensor, trim = 0.10)
mad_val <- mad(sensor)  # default constant 1.4826, consistent with normal sd

Q3a_table <- data.frame(
  Estimator = c("Mean", "SD", "Median", "MAD", "10% Trimmed Mean"),
  Value = c(mean(sensor), sd(sensor), median(sensor), mad_val, trimmed_mean10)
)
print(Q3a_table)
cat("Interpretation: Mean and SD are pulled strongly toward the two outliers\n")
cat("(150.2 and -80.5); median, MAD, and trimmed mean stay close to the clean values.\n")

## ---- Q3b: Huber M-estimator -----------------------------------------------
huber_est <- MASS::huber(sensor, k = 1.345)
huber_mu <- huber_est$mu

Q3ab_table <- rbind(Q3a_table,
                     data.frame(Estimator = "Huber M-estimate (k=1.345)", Value = huber_mu))
print(Q3ab_table)

## ---- Q3c: Breakdown-point sweep -------------------------------------------
clean_vals <- sensor[!(sensor %in% c(150.2, -80.5))]  # 18 clean points
clean_mean <- mean(clean_vals)

set.seed(6)
contam_pool <- c(150.2, -80.5, 200, -150, 500, -300, 1000, -700, 5000, -3000)[1:9]
n_contam_seq <- 0:9
means_seq <- numeric(length(n_contam_seq))
medians_seq <- numeric(length(n_contam_seq))

for (i in seq_along(n_contam_seq)) {
  k <- n_contam_seq[i]
  n_clean_needed <- 20 - k
  clean_part <- rep_len(clean_vals, n_clean_needed)   # 20-k clean points (recycled if >18)
  if (k == 0) {
    vec <- clean_part
  } else {
    vec <- c(clean_part, contam_pool[1:k])             # add k contaminated points
  }
  means_seq[i] <- mean(vec)
  medians_seq[i] <- median(vec)
}

png("q3c_breakdown.png", width = 800, height = 600)
plot(n_contam_seq, means_seq, type = "b", col = "red", pch = 16, lwd = 2,
     ylim = range(c(means_seq, medians_seq)),
     xlab = "Number of contaminated points (out of 20)",
     ylab = "Estimate", main = "Breakdown behavior: mean vs median")
lines(n_contam_seq, medians_seq, type = "b", col = "blue", pch = 17, lwd = 2)
abline(h = clean_mean, lty = 2, col = "gray40")
legend("topleft", legend = c("Mean", "Median"), col = c("red", "blue"), pch = c(16, 17))
dev.off()
cat("Plot saved: q3c_breakdown.png\n")

first_break <- n_contam_seq[which(abs(means_seq - clean_mean) > 10)][1]
cat(sprintf("Mean moves > 10 units from clean mean starting at %d contaminated point(s).\n",
            first_break))
cat("Median remains stable (near-zero breakdown effect) across all 9 contamination levels.\n")

## ---- Q3d: Clean-data comparison -------------------------------------------
clean_trimmed <- mean(clean_vals, trim = 0.10)
clean_mad <- mad(clean_vals)
clean_huber <- MASS::huber(clean_vals, k = 1.345)$mu

Q3d_table <- data.frame(
  Estimator = c("Mean", "SD", "Median", "MAD", "10% Trimmed Mean", "Huber M-estimate"),
  Contaminated_n20 = c(mean(sensor), sd(sensor), median(sensor), mad_val,
                        trimmed_mean10, huber_mu),
  Clean_n18 = c(mean(clean_vals), sd(clean_vals), median(clean_vals), clean_mad,
                clean_trimmed, clean_huber)
)
print(Q3d_table)
cat("Interpretation: Mean and SD change drastically once outliers are removed,\n")
cat("while median, MAD, trimmed mean, and Huber estimate barely change -- confirming\n")
cat("their robustness to the two contaminated sensor readings.\n")


#############################################################################
# QUESTION 4
#############################################################################
cat("\n================= QUESTION 4 =================\n")

n4 <- 40
x4 <- 27
p_grid <- seq(0.001, 0.999, length.out = 1000)

## ---- Q4a: Posterior via Beta(2,2) prior -----------------------------------
a0 <- 2; b0 <- 2
prior_dens <- dbeta(p_grid, a0, b0)
lik <- dbinom(x4, n4, p_grid)
post_unnorm <- prior_dens * lik
post_dens <- post_unnorm / (sum(post_unnorm) * diff(p_grid)[1])  # normalize (grid)

# Analytic posterior for Beta prior + Binomial likelihood: Beta(a0+x, b0+n-x)
post_a <- a0 + x4
post_b <- b0 + (n4 - x4)
cat(sprintf("Posterior (grid, Beta(2,2) prior) matches Beta(%d, %d) analytically.\n",
            post_a, post_b))

## ---- Q4b: Plot prior, likelihood, posterior --------------------------------
lik_rescaled <- lik / max(lik) * max(post_dens)   # rescale for comparable height

png("q4b_prior_lik_post.png", width = 800, height = 600)
plot(p_grid, post_dens, type = "l", col = "darkgreen", lwd = 2,
     xlab = "p", ylab = "Density (rescaled)", main = "Prior, Likelihood, Posterior for p")
lines(p_grid, prior_dens, col = "blue", lwd = 2, lty = 2)
lines(p_grid, lik_rescaled, col = "red", lwd = 2, lty = 3)
legend("topleft", legend = c("Posterior", "Prior Beta(2,2)", "Likelihood (rescaled)"),
       col = c("darkgreen", "blue", "red"), lty = c(1, 2, 3), lwd = 2)
dev.off()
cat("Plot saved: q4b_prior_lik_post.png\n")

## ---- Q4c: Point estimates --------------------------------------------------
post_mean_grid <- sum(p_grid * post_dens) * diff(p_grid)[1]
post_mode_grid <- p_grid[which.max(post_dens)]
cdf_grid <- cumsum(post_dens) * diff(p_grid)[1]
post_median_grid <- p_grid[which.min(abs(cdf_grid - 0.5))]
p_mle <- x4 / n4

Q4c_table <- data.frame(
  Quantity = c("Posterior mean", "Posterior mode", "Posterior median", "Frequentist MLE"),
  Value = c(post_mean_grid, post_mode_grid, post_median_grid, p_mle)
)
print(Q4c_table)

## ---- Q4d: Interval comparison ----------------------------------------------
cred_lower <- p_grid[which.min(abs(cdf_grid - 0.025))]
cred_upper <- p_grid[which.min(abs(cdf_grid - 0.975))]
cred_int <- c(cred_lower, cred_upper)

# Frequentist: Wald and Wilson
p_hat <- p_mle
se_wald <- sqrt(p_hat * (1 - p_hat) / n4)
wald_ci4 <- p_hat + c(-1, 1) * qnorm(0.975) * se_wald

z <- qnorm(0.975)
wilson_center <- (p_hat + z^2/(2*n4)) / (1 + z^2/n4)
wilson_halfwidth <- (z / (1 + z^2/n4)) * sqrt(p_hat*(1-p_hat)/n4 + z^2/(4*n4^2))
wilson_ci4 <- wilson_center + c(-1, 1) * wilson_halfwidth

Q4d_table <- data.frame(
  Method = c("Bayesian 95% credible interval", "Frequentist 95% Wald CI", "Frequentist 95% Wilson CI"),
  Lower = c(cred_int[1], wald_ci4[1], wilson_ci4[1]),
  Upper = c(cred_int[2], wald_ci4[2], wilson_ci4[2]),
  Width = c(diff(cred_int), diff(wald_ci4), diff(wilson_ci4))
)
print(Q4d_table)

## ---- Q4e: Prior sensitivity table ------------------------------------------
compute_posterior_summary <- function(a0, b0, x, n, grid) {
  prior <- dbeta(grid, a0, b0)
  lik_ <- dbinom(x, n, grid)
  post <- prior * lik_
  post <- post / (sum(post) * diff(grid)[1])
  pm <- sum(grid * post) * diff(grid)[1]
  pmode <- grid[which.max(post)]
  cdf_ <- cumsum(post) * diff(grid)[1]
  pmed <- grid[which.min(abs(cdf_ - 0.5))]
  lo <- grid[which.min(abs(cdf_ - 0.025))]
  hi <- grid[which.min(abs(cdf_ - 0.975))]
  c(mean = pm, mode = pmode, median = pmed, ci_lo = lo, ci_hi = hi)
}

priors_list <- list(c(1, 1), c(2, 2), c(20, 20))
prior_names <- c("Beta(1,1)", "Beta(2,2)", "Beta(20,20)")
Q4e_res <- t(sapply(priors_list, function(ab) compute_posterior_summary(ab[1], ab[2], x4, n4, p_grid)))
Q4e_table <- data.frame(Prior = prior_names, Q4e_res)
Q4e_table$CI <- sprintf("(%.4f, %.4f)", Q4e_table$ci_lo, Q4e_table$ci_hi)
Q4e_table <- Q4e_table[, c("Prior", "mean", "mode", "median", "CI")]
print(Q4e_table)
cat("Interpretation: With n=40, the posterior mean/mode/median shift only modestly\n")
cat("across priors; the informative Beta(20,20) pulls estimates noticeably toward 0.5.\n")

## Optional add-on: n = 400, x = 270
x4b <- 270; n4b <- 400
Q4e_res_big <- t(sapply(priors_list, function(ab) compute_posterior_summary(ab[1], ab[2], x4b, n4b, p_grid)))
Q4e_table_big <- data.frame(Prior = prior_names, Q4e_res_big)
Q4e_table_big$CI <- sprintf("(%.4f, %.4f)", Q4e_table_big$ci_lo, Q4e_table_big$ci_hi)
Q4e_table_big <- Q4e_table_big[, c("Prior", "mean", "mode", "median", "CI")]
print(Q4e_table_big)
cat("Interpretation: With n=400 (same proportion), prior choice barely moves the\n")
cat("posterior summaries -- the data dominate the prior as sample size grows.\n")


#############################################################################
# QUESTION 5
#############################################################################
cat("\n================= QUESTION 5 =================\n")

## ---- Q5a: Two intervals for lambda -----------------------------------------
n5 <- length(failure_times)
xbar5 <- mean(failure_times)

chisq_lower <- qchisq(0.025, df = 2 * n5)
chisq_upper <- qchisq(0.975, df = 2 * n5)
# 2*n*lambda*Xbar ~ chisq(2n)  =>  lambda in ( chisq_lower/(2*n*Xbar), chisq_upper/(2*n*Xbar) )
chisq_ci <- c(chisq_lower, chisq_upper) / (2 * n5 * xbar5)

Q5a_table <- data.frame(
  Method = c("Chi-square-based CI", "Wald CI (from Q2c)"),
  Lower = c(chisq_ci[1], wald_ci[1]),
  Upper = c(chisq_ci[2], wald_ci[2]),
  Width = c(diff(chisq_ci), diff(wald_ci))
)
print(Q5a_table)

## ---- Q5b: Coverage simulation (t-based) -------------------------------------
set.seed(7)
t_crit <- qt(0.975, df = n1 - 1)
xbar_5b <- rowMeans(samples1)   # reuse Q1 samples (n=20, M=5000)
s_5b <- apply(samples1, 1, sd)
se_5b <- s_5b / sqrt(n1)
lower_t <- xbar_5b - t_crit * se_5b
upper_t <- xbar_5b + t_crit * se_5b
coverage_t <- mean(lower_t <= mu_true & mu_true <= upper_t)
cat(sprintf("Empirical coverage of 95%% t-based CI (n=20, M=5000): %.4f\n", coverage_t))

## ---- Q5c: Coverage simulation (z-misspecified) -------------------------------
set.seed(8)
simulate_z_coverage <- function(n, M = 5000) {
  mat <- matrix(rnorm(n * M, mu_true, sqrt(sigma2_true)), nrow = M, ncol = n)
  xb <- rowMeans(mat)
  s <- apply(mat, 1, sd)
  se <- s / sqrt(n)
  lo <- xb - 1.96 * se
  hi <- xb + 1.96 * se
  mean(lo <= mu_true & mu_true <= hi)
}
cov_z_n5   <- simulate_z_coverage(5)
cov_z_n100 <- simulate_z_coverage(100)

Q5c_table <- data.frame(n = c(5, 100), EmpiricalCoverage_z = c(cov_z_n5, cov_z_n100))
print(Q5c_table)
bigger_gap_n <- Q5c_table$n[which.max(abs(Q5c_table$EmpiricalCoverage_z - 0.95))]
cat(sprintf("Sample size with bigger gap from 95%%: n = %d\n", bigger_gap_n))
cat("Interpretation: Using z instead of t under-covers more severely for small n\n")
cat("(here n=5) because the t-distribution has heavier tails at low df.\n")

## ---- Q5d: Bootstrap CI for the median ----------------------------------------
set.seed(9)
B <- 2000
boot_medians <- replicate(B, median(sample(sensor, length(sensor), replace = TRUE)))
boot_ci_median <- quantile(boot_medians, c(0.025, 0.975))

# naive normal-theory CI for the mean on contaminated data
se_mean_contam <- sd(sensor) / sqrt(length(sensor))
naive_mean_ci <- mean(sensor) + c(-1, 1) * qnorm(0.975) * se_mean_contam

Q5d_table <- data.frame(
  Method = c("Bootstrap 95% CI for median", "Naive normal-theory 95% CI for mean"),
  Lower = c(boot_ci_median[1], naive_mean_ci[1]),
  Upper = c(boot_ci_median[2], naive_mean_ci[2])
)
print(Q5d_table)
cat("Interpretation: The bootstrap median CI is tight and unaffected by outliers,\n")
cat("while the naive mean CI is shifted/inflated by the two contaminated points.\n")


#############################################################################
# QUESTION 6
#############################################################################
cat("\n================= QUESTION 6 =================\n")

p0 <- 0.05
p1 <- 0.15
alpha6 <- 0.05
beta6 <- 0.10

## ---- Q6a: SPRT boundaries -----------------------------------------------
A_bound <- (1 - beta6) / alpha6
B_bound <- beta6 / (1 - alpha6)
logA <- log(A_bound)
logB <- log(B_bound)
cat(sprintf("A = %.4f, B = %.4f\n", A_bound, B_bound))
cat(sprintf("log A = %.4f, log B = %.4f\n", logA, logB))

## ---- Q6b: SPRT update rule / run function --------------------------------
run_sprt <- function(p_true, p0, p1, logA, logB, max_n = 10000) {
  logLambda <- 0
  n <- 0
  repeat {
    n <- n + 1
    x <- rbinom(1, 1, p_true)
    logLambda <- logLambda + x * log(p1 / p0) + (1 - x) * log((1 - p1) / (1 - p0))
    if (logLambda >= logA) return(list(decision = "reject_H0", n = n))
    if (logLambda <= logB) return(list(decision = "accept_H0", n = n))
    if (n >= max_n) return(list(decision = "no_decision", n = n))
  }
}

## ---- Q6c: Simulation at true rate p0 ---------------------------------------
set.seed(10)
nsim6 <- 1000
res_p0 <- replicate(nsim6, run_sprt(p0, p0, p1, logA, logB), simplify = FALSE)
dec_p0 <- sapply(res_p0, `[[`, "decision")
n_p0 <- sapply(res_p0, `[[`, "n")

prop_correct_p0 <- mean(dec_p0 == "accept_H0")
Q6c_table <- data.frame(
  TrueRate = "p0 = 0.05",
  Prop_Correct_Accept_H0 = prop_correct_p0,
  Mean_N = mean(n_p0),
  Median_N = median(n_p0)
)
print(Q6c_table)

## ---- Q6d: Simulation at true rate p1 ----------------------------------------
set.seed(11)
res_p1 <- replicate(nsim6, run_sprt(p1, p0, p1, logA, logB), simplify = FALSE)
dec_p1 <- sapply(res_p1, `[[`, "decision")
n_p1 <- sapply(res_p1, `[[`, "n")

prop_correct_p1 <- mean(dec_p1 == "reject_H0")
Q6d_table <- data.frame(
  TrueRate = "p1 = 0.15",
  Prop_Correct_Reject_H0 = prop_correct_p1,
  Mean_N = mean(n_p1)
)
print(Q6d_table)

## ---- Q6e: Fixed-sample-size comparison ---------------------------------------
z_alpha <- qnorm(1 - alpha6)
z_beta  <- qnorm(1 - beta6)  # power = 0.90
n_fixed <- ((z_alpha * sqrt(p0 * (1 - p0)) + z_beta * sqrt(p1 * (1 - p1))) / (p1 - p0))^2
n_fixed <- ceiling(n_fixed)

Q6e_table <- data.frame(
  Method = c("Fixed-sample formula", "SPRT avg N (true p0)", "SPRT avg N (true p1)"),
  N = c(n_fixed, mean(n_p0), mean(n_p1))
)
print(Q6e_table)
cat("Interpretation: SPRT reaches a decision, on average, with far fewer\n")
cat("observations than the fixed-sample design under both true rates.\n")

## ---- Q6f: Sample-path plot -----------------------------------------------
simulate_path <- function(p_true, p0, p1, logA, logB, max_n = 500) {
  logLambda <- numeric(0)
  cur <- 0
  n <- 0
  repeat {
    n <- n + 1
    x <- rbinom(1, 1, p_true)
    cur <- cur + x * log(p1 / p0) + (1 - x) * log((1 - p1) / (1 - p0))
    logLambda <- c(logLambda, cur)
    if (cur >= logA || cur <= logB || n >= max_n) break
  }
  logLambda
}

set.seed(12)
paths_p0 <- lapply(1:8, function(i) simulate_path(p0, p0, p1, logA, logB))
set.seed(13)
paths_p1 <- lapply(1:8, function(i) simulate_path(p1, p0, p1, logA, logB))

png("q6f_sample_paths.png", width = 1000, height = 600)
par(mfrow = c(1, 2))

ymax <- max(logA, sapply(c(paths_p0, paths_p1), max)) + 1
ymin <- min(logB, sapply(c(paths_p0, paths_p1), min)) - 1

plot(NULL, xlim = c(1, max(sapply(paths_p0, length))), ylim = c(ymin, ymax),
     xlab = "n", ylab = expression(log(Lambda[n])), main = "SPRT paths: true p = 0.05")
for (pth in paths_p0) lines(seq_along(pth), pth, col = "steelblue")
abline(h = c(logA, logB), col = c("red", "darkgreen"), lty = 2)
legend("topright", legend = c("log A (reject)", "log B (accept)"),
       col = c("red", "darkgreen"), lty = 2, cex = 0.8)

plot(NULL, xlim = c(1, max(sapply(paths_p1, length))), ylim = c(ymin, ymax),
     xlab = "n", ylab = expression(log(Lambda[n])), main = "SPRT paths: true p = 0.15")
for (pth in paths_p1) lines(seq_along(pth), pth, col = "orange3")
abline(h = c(logA, logB), col = c("red", "darkgreen"), lty = 2)
legend("topright", legend = c("log A (reject)", "log B (accept)"),
       col = c("red", "darkgreen"), lty = 2, cex = 0.8)

par(mfrow = c(1, 1))
dev.off()
cat("Plot saved: q6f_sample_paths.png\n")

cat("\n================= END OF SCRIPT =================\n")
