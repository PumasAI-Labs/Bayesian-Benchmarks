# This is the future LOO implementation described in
# https://mc-stan.org/loo/articles/loo2-lfo.html
# But accounting by index subjects and observations
# using the NONMEM's stuff from Stan's data

library(loo)
library(cmdstanr)
library(dplyr)
library(tibble)

set_cmdstan_path("Torsten/cmdstan")

data <- rstan::read_rdump("Torsten/example-models/poppk2cpt/poppk2cpt.data.R")

nonmem_data <- tibble(
    c = NA_real_,
    time = data$time,
    amt = data$amt,
    evid = data$evid,
    rate = data$rate,
    addl = data$addl,
    ii = data$ii,
    cmt = data$cmt,
    bloq = 0,
    lloq = 0,
    ss = data$ss
) %>%
    mutate(
        i = 1:n(),
        ID = cut(1:n(),
            breaks = c(data$start, last(data$end)),
            labels = FALSE, right = FALSE, include.lowest = TRUE
        )
    ) %>%
    left_join(tibble(
        DV = data$cObs,
        i = data$iObs
    ),
    by = "i"
    ) %>%
    mutate(DV = if_else(is.na(DV),
        5555555,
        DV
    )) %>% # This value can be anything > 0. It'll be indexed away
    select(c, ID, everything(), -i)

n_subjects <- nonmem_data %>% # number of individuals
    distinct(ID) %>%
    count() %>%
    deframe()

n_total <- nrow(nonmem_data) # total number of records

i_obs <- nonmem_data %>%
    mutate(row_num = 1:n()) %>%
    filter(evid == 0) %>%
    select(row_num) %>%
    deframe()

n_obs <- length(i_obs)

subj_start <- nonmem_data %>%
    mutate(row_num = 1:n()) %>%
    group_by(ID) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(row_num) %>%
    deframe()

subj_end <- c(subj_start[-1] - 1, n_total)

## LOO - helper functions
# more stable than log(sum(exp(x)))
log_sum_exp <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
    log_sum_exp(x) - log(length(x))
}

# compute log of raw importance ratios
# sums over observations *not* over posterior samples
sum_log_ratios <- function(loglik, ids = NULL) {
    if (!is.null(ids)) loglik <- loglik[, ids, drop = FALSE]
    rowSums(loglik)
}

# for printing comparisons later
rbind_print <- function(...) {
    round(rbind(...), digits = 2)
}

# Get back fitted run with `.rds` file
fit_normal <- readRDS(
    file.path(
        "Torsten-Runs",
        "poppk2cpt_fit_normal.rds"
    )
)

# ELPD (which is the mean of the log_lik)
fit_normal_loo <- fit_normal$loo()
loo_df <- fit_normal_loo$pointwise %>%
    as_data_frame() %>%
    mutate(
        row_number = row_number()
    ) %>%
    left_join(
        nonmem_data %>%
            filter(evid == 0) %>%
            transmute(
                row_number = row_number(),
                ID, time, DV, bloq
            ),
        by = "row_number"
    ) %>%
    relocate(
        row_number,
        .before = elpd_loo
    )


# log_lik matrix
log_lik_matrix <- fit_normal$draws("log_lik", format = "draws_matrix")
loo(log_lik_matrix) # "Vanilla" LOO
dim(log_lik_matrix) # 4000 by 530
nrow(log_lik_matrix)

# Approximate M-step-ahead (SAP) predictions
# set `M` to have 1-SAP
M <- 1
N <- length(i_obs) # the same numbers of cols as `log_lik_matrix`
n_time <- unique(loo_df$time) %>% length()
loglikm <- matrix(nrow = n_time, ncol = N)

# L is the number of
future_loo <- function(log_lik_matrix,
                       n_time, subj_start, subj_end,
                       N = ncol(log_lik_matrix),
                       nsamples = nrow(log_lik_matrix),
                       M = 1, L = 20, n_subj = 10) {
    loglikm <- matrix(nrow = nsamples, ncol = N)
    for (s in 1:n_subj) { # loop over subjects
        for (i in L:n_time) { # loop over n_time
            past <- (subj_start[s]):(subj_start[s] + i)
            oos <- (subj_start[s] + i):(subj_start[s] + i + M)
        }
    }
}
# TODO: nested for loop for subjects
for (i in L:(N - M)) {
    past <- 1:i
    oos <- (i + 1):(i + M)
    df_past <- df[past, , drop = FALSE]
    df_oos <- df[c(past, oos), , drop = FALSE]
    fit_past <- update(fit, newdata = df_past, recompile = FALSE)
    loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
    loglikm[, i + 1] <- rowSums(loglik[, oos])
}

exact_elpds_sap <- apply(loglikm, 2, log_mean_exp)
exact_elpd_sap <- c(ELPD = sum(exact_elpds_4sap, na.rm = TRUE))
