#---------------------------------------------
# 1. Load libraries

library(randomizr)
library(estimatr)
library(ggplot2)
library(readxl)
library(dplyr)
library(ri2)

# remove all previous objects from the environment
rm(list = ls())

# set random seed
set.seed(080403)

# design constants from the PAP
n_accounts <- 30    # 30 TikTok accounts
n_treat    <- 15    # 15 treatment / 15 control in each randomization
n_fyp      <- 50    # outcome = political videos among first 50 "For You" videos
n_videos   <- 150   # 25 videos/hashtag x 3 hashtags x 2 randomizations
n_sims     <- 1000  # randomization inference simulations

control_tags <- "#cockapoo, #cockapoosoftiktok, #cockapoolife"
treat_tags   <- "#cleanbeauty, #tradwife, #softliving"

# helper: F-statistic from an OLS regression
get_f <- function(formula, data) {
  unname(summary(lm(formula, data = data))$fstatistic["value"])
}

#---------------------------------------------
# 2. Conduct first randomization (complete RA, 15/15)
Z1 <- complete_ra(N = n_accounts, m = n_treat)

#---------------------------------------------
# 3. Conduct second randomization (complete RA, 15/15, independent of the first)
Z2 <- complete_ra(N = n_accounts, m = n_treat)

# Merge both randomizations with the account IDs / usernames from Table 1
assignments <- data.frame(
  account  = 1:n_accounts,
  username = paste0("notoriousdpg", 1:n_accounts),
  Z1       = Z1,
  Z2       = Z2
) |>
  mutate(
    # number of rounds of conservative-coded content (0, 1, or 2) -- for H2
    randomize      = Z1 + Z2,
    hashtags1 = if_else(Z1 == 1, treat_tags, control_tags),
    hashtags2 = if_else(Z2 == 1, treat_tags, control_tags)
  )

# 15/15 in each round, plus the 2x2 cross-tab of the two assignments
# print(table(Z1 = assignments$Z1, Z2 = assignments$Z2))

# save the assignment sheet to follow while running the experiment
write.csv(assignments, "assignments.csv", row.names = FALSE)

#---------------------------------------------
# 4. Randomization check 1: are the two randomizations independent?

declaration_Z2 <- declare_ra(N = n_accounts, m = n_treat)

balance_fun_1 <- function(data) get_f(Z2 ~ Z1, data)

ri_balance_1 <- conduct_ri(
  test_function = balance_fun_1,
  declaration   = declaration_Z2,
  assignment    = "Z2",
  sims          = n_sims,
  p             = "upper",
  data          = assignments
)
# print(summary(ri_balance_1))

#---------------------------------------------
# 5. Read the data files and do the data transformations according to the PAP
#
# Expected structure of the real data:
#   outcomes_data (one row per account):
#     account   -- 1 to 30
#     pol_wave0 -- # explicitly political videos in first 50 FYP videos,
#                  measured right after account creation, BEFORE
#                  randomization 1 (baseline covariate)
#     pol_wave1 -- # explicitly political videos in first 50 FYP videos,
#                  measured after randomization 1 (covariate)
#     pol_wave2 -- # explicitly political videos in first 50 FYP videos,
#                  measured after randomization 2 (outcome)
#   video_data (one row per liked video, 150 per account):
#     account, round (1/2), hashtag, likes, comments, followers

# ----- REAL DATA (uncomment once collected; delete the sample data below) -----
# outcomes_data <- read_excel("practicum_data.xlsx", sheet = "outcomes")
# video_data    <- read_excel("practicum_data.xlsx", sheet = "videos")

# ----- SAMPLE DATA (simulated so the code runs before data collection) --------
outcomes_data <- assignments |>
  mutate(
    # account-specific baseline propensity to be shown political content
    base_p    = runif(n_accounts, min = 0.02, max = 0.08),
    pol_wave0 = rbinom(n_accounts, size = n_fyp, prob = base_p),
    pol_wave1 = rbinom(n_accounts, size = n_fyp, prob = base_p + 0.08 * Z1),
    pol_wave2 = rbinom(n_accounts, size = n_fyp,
                       prob = base_p + 0.05 * Z1 + 0.08 * Z2)
  ) |>
  select(account, pol_wave0, pol_wave1, pol_wave2)

video_data <- expand.grid(
  video   = 1:25,
  hashtag = 1:3,
  round   = 1:2,
  account = 1:n_accounts
) |>
  left_join(select(assignments, account, Z1, Z2), by = "account") |>
  mutate(
    treated   = if_else(round == 1, Z1, Z2),
    # log-normal prevents negative numbers and creates "influencer" spikes
    followers = round(rlnorm(n(), meanlog = 8 + 0.5 * treated, sdlog = 1.5)),
    # Poisson distribution for typical event counts
    likes     = rpois(n(), lambda = 250 + 100 * treated),
    comments  = rpois(n(), lambda = 15 + 10 * treated)
  ) |>
  select(account, round, hashtag, video, likes, comments, followers)
# ------------------------------------------------------------------------------

# merge assignments and outcomes into one analysis dataset
practicum_data <- assignments |>
  left_join(outcomes_data, by = "account")

#---------------------------------------------
# 6. Covariate balance checks (F-statistic + RI p-value)

# Check 2: is Z1 balanced on the baseline covariate pol_wave0?
declaration_Z1 <- declare_ra(N = n_accounts, m = n_treat)

balance_fun_z1 <- function(data) get_f(Z1 ~ pol_wave0, data)

ri_balance_z1 <- conduct_ri(
  test_function = balance_fun_z1,
  declaration   = declaration_Z1,
  assignment    = "Z1",
  sims          = n_sims,
  p             = "upper",
  data          = practicum_data
)
# print(summary(ri_balance_z1))

# Check 3: is Z2 balanced on pol_wave0 and pol_wave1?
balance_fun_2 <- function(data) get_f(Z2 ~ pol_wave0 + pol_wave1, data)

ri_balance_2 <- conduct_ri(
  test_function = balance_fun_2,
  declaration   = declaration_Z2,
  assignment    = "Z2",
  sims          = n_sims,
  p             = "upper",
  data          = practicum_data
)
# print(summary(ri_balance_2))

#---------------------------------------------
# 7. Descriptive results: mean (and SE) of the outcome by condition

tab_1 <- practicum_data |>
  group_by(Z2) |>
  summarise(
    n           = n(),
    outcome_mu  = mean(pol_wave2),
    outcome_se  = sd(pol_wave2) / sqrt(n()),
    wave0_mu    = mean(pol_wave0),
    wave1_mu    = mean(pol_wave1),
    .groups = "drop"
  )
# print(tab_1)

tab_2 <- practicum_data |>
  group_by(randomize) |>
  summarise(
    n          = n(),
    outcome_mu = mean(pol_wave2),
    outcome_se = sd(pol_wave2) / sqrt(n()),
    .groups = "drop"
  )
# print(tab_2)

# Descriptive statistics for videos: likes, comments, followers
tab_3 <- video_data |>
  left_join(select(assignments, account, Z1, Z2), by = "account") |>
  mutate(condition = if_else(if_else(round == 1, Z1, Z2) == 1,
                             "Conservative-coded", "Apolitical (cockapoo)")) |>
  group_by(round, condition) |>
  summarise(
    n_videos       = n(),
    mean_likes     = mean(likes),
    mean_comments  = mean(comments),
    mean_followers = mean(followers),
    # composite on the log scale per the PAP: log(x + 1) avoids -Inf at 0
    log_mean_likes     = log(mean_likes + 1),
    log_mean_comments  = log(mean_comments + 1),
    log_mean_followers = log(mean_followers + 1),
    .groups = "drop"
  )
# print(tab_3)

#---------------------------------------------
# 8. Regression analysis
# H1: difference in means by assignment in randomization 2
# H2: outcome on number of rounds of conservative-coded content (0/1/2)
# Covariate adjustment uses only pre-treatment measures:
#   - Z2 models: pol_wave1 (PAP), then pol_wave0 + pol_wave1
#   - randomize model: pol_wave0 only (pol_wave1 is affected by Z1, part of randomize)

f_h1      <- pol_wave2 ~ Z2
f_h2      <- pol_wave2 ~ randomize
f_cov     <- pol_wave2 ~ Z2 + pol_wave1
f_full    <- pol_wave2 ~ Z2 + pol_wave0 + pol_wave1
f_h2_cov  <- pol_wave2 ~ randomize + pol_wave0

m1 <- lm(f_h1,     data = practicum_data)
m2 <- lm(f_h2,     data = practicum_data)
m3 <- lm(f_cov,    data = practicum_data)
m4 <- lm(f_full,   data = practicum_data)
m5 <- lm(f_h2_cov, data = practicum_data)

# Generate table output (regular OLS standard errors; in the .Rmd, put this in a chunk with results='asis')
# stargazer::stargazer(m1, m2, m3, m4, m5, type = "latex", header = FALSE,
#                      no.space = TRUE, title = "Treatment Effect Estimates")

#---------------------------------------------
# 9. Randomization inference

# --- H1: one-tailed (upper) test of the sharp null of no effect
ri_h1 <- conduct_ri(
  f_h1,
  declaration      = declaration_Z2,
  assignment       = "Z2",
  sharp_hypothesis = 0,
  sims             = n_sims,
  p                = "upper",
  data             = practicum_data
)
# summary(ri_h1)
# plot(ri_h1)

# 95% interval for the H1 estimate, assuming a constant effect equal to the ATE
ate_h1 <- ri_h1$sims_df$est_obs[1]
ri_h1_ci <- conduct_ri(
  f_h1,
  declaration      = declaration_Z2,
  assignment       = "Z2",
  sharp_hypothesis = ate_h1,
  sims             = n_sims,
  data             = practicum_data
)
ci_h1 <- quantile(ri_h1_ci$sims_df$est_sim, c(0.025, 0.975))

# --- H2: one-tailed (upper) test for the randomize variable.
randomize_est <- function(formula, data) coef(lm(formula, data = data))[["randomize"]]

ri_randomize <- function(formula, data, p = c("upper", "two-tailed")) {
  p   <- match.arg(p)
  obs <- randomize_est(formula, data)
  sims <- replicate(n_sims, {
    data$randomize <- complete_ra(N = n_accounts, m = n_treat) +
      complete_ra(N = n_accounts, m = n_treat)
    randomize_est(formula, data)  # outcome fixed under the sharp null
  })
  pval <- if (p == "upper") mean(sims >= obs) else mean(abs(sims) >= abs(obs))
  list(estimate = obs, p_value = pval, sims = sims)
}

ri_h2  <- ri_randomize(f_h2, practicum_data, p = "upper")
est_h2 <- ri_h2$estimate
p_h2   <- ri_h2$p_value

# --- H2 adjusted for the baseline covariate pol_wave0: two-tailed
ri_h2_cov <- ri_randomize(f_h2_cov, practicum_data, p = "two-tailed")

# --- Covariate-adjusted (pol_wave1): two-tailed
ri_cov <- conduct_ri(
  f_cov,
  declaration      = declaration_Z2,
  assignment       = "Z2",
  sharp_hypothesis = 0,
  sims             = n_sims,
  p                = "two-tailed",
  data             = practicum_data
)

# --- Adjusted for pol_wave0 + pol_wave1: two-tailed
ri_full <- conduct_ri(
  f_full,
  declaration      = declaration_Z2,
  assignment       = "Z2",
  sharp_hypothesis = 0,
  sims             = n_sims,
  p                = "two-tailed",
  data             = practicum_data
)

# --- Results table: estimate, RI p-value, and the model F-statistic
results <- data.frame(
  model     = c("H1: difference in means (Z2)",
                "H2: randomize (Z1 + Z2)",
                "H1: Z2 + pol_wave1",
                "H1: Z2 + pol_wave0 + pol_wave1",
                "H2: randomize + pol_wave0"),
  test      = c("one-tailed", "one-tailed", "two-tailed", "two-tailed",
                "two-tailed"),
  estimate  = c(ate_h1, est_h2,
                ri_cov$sims_df$est_obs[1], ri_full$sims_df$est_obs[1],
                ri_h2_cov$estimate),
  ri_p      = c(summary(ri_h1)$upper_p_value,
                p_h2,
                summary(ri_cov)$two_tailed_p_value,
                summary(ri_full)$two_tailed_p_value,
                ri_h2_cov$p_value),
  f_stat    = c(get_f(f_h1,     practicum_data),
                get_f(f_h2,     practicum_data),
                get_f(f_cov,    practicum_data),
                get_f(f_full,   practicum_data),
                get_f(f_h2_cov, practicum_data))
)
# print(results)
# print(ci_h1)

#---------------------------------------------
# 10. A graphical depiction of experimental findings in style of Coppock (2021)

coppock_plot <- function(data, group, labels) {
  data <- data |> mutate(group = factor(.data[[group]], labels = labels))
  # group means with regular OLS 95% confidence intervals
  fit <- lm(pol_wave2 ~ 0 + group, data = data)
  ci  <- confint(fit)
  summary_df <- data.frame(
    group     = factor(labels, levels = labels),
    pol_wave2 = unname(coef(fit)),
    conf.low  = unname(ci[, 1]),
    conf.high = unname(ci[, 2])
  )
  
  ggplot(summary_df, aes(group, pol_wave2)) +
    geom_point(
      data = data,
      position = position_jitter(width = 0.15, height = 0),
      alpha = 0.3
    ) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
    theme_bw() +
    scale_y_continuous(limits = c(0, n_fyp), breaks = seq(0, n_fyp, by = 10)) +
    theme(axis.title.x = element_blank()) +
    ylab("Explicitly political videos\n(of first 50 'For You' videos)")
}

plot_h1 <- coppock_plot(practicum_data, "Z2",
                        c("Apolitical (cockapoo)", "Conservative-coded"))
plot_h2 <- coppock_plot(practicum_data, "randomize",
                        c("0 rounds", "1 round", "2 rounds"))
# print(plot_h1)
# print(plot_h2)
