cat("sample size calculations")

library(tidyverse)
library(furrr)

plan(multisession)

sample_size_dir <- "sample-size"

# SECTION Functions

#' Returns OG scale titres
censor_logtitres <- function(logtitres) {
  cuts <- cut(exp(logtitres), c(-Inf, 5 * 2^(1:10), Inf)) %>% as.integer()
  5 * 2^(cuts - 1)
}

#' Returns log titres
midpoint_titres <- function(titres) {
  if_else(titres == 5L, log(titres), log(titres) + log(2) / 2)
}

#' Outcome - postvax titre
one_study <- function(n_per_group = 50,
                      logbaseline_mean = 3.7,
                      logtitre_sd = 1,
                      b0 = 1.1, # Postvax titre in infreq + qiv for base 1
                      blogbaseline = 0.8, # Add to postvax for 1 logbase
                      bfreq = -0.2, # Add to postvax titre for freq
                      bflucellvax = 0.25, # Add to postvax titre for flucell
                      bflublok = 0.5, # Add to postvax titre for flublock
                      bflucellvax_freq_add = 0, # Add to add for flucell
                      bflublok_freq_add = 0.2) { # Add to add for flublock
  groups <- map_dfr(
    0:1,
    function(freq) {
      map_dfr(
        c("qiv", "flucellvax", "flublok"),
        ~ tibble(.rows = n_per_group, freq = freq, vac = .x)
      )
    }
  )
  study <- groups %>%
    mutate(
      id = row_number(),
      freq_lab = recode(freq, "0" = "Infrequent", "1" = "Frequent"),
      logbaseline = rnorm(n(), logbaseline_mean, logtitre_sd),
      logpostvax = rnorm(
        n(),
        mean = b0 + blogbaseline * logbaseline +
          bfreq * freq +
          bflucellvax * (vac == "flucellvax") +
          bflublok * (vac == "flublok") +
          bflucellvax_freq_add * (vac == "flucellvax" & freq == 1) +
          bflublok_freq_add * (vac == "flublok" & freq == 1),
        sd = logtitre_sd
      ),
      baseline = censor_logtitres(logbaseline),
      postvax = censor_logtitres(logpostvax),
      logbaseline_mid = midpoint_titres(baseline),
      logpostvax_mid = midpoint_titres(postvax),
      vac = factor(vac, levels = c("qiv", "flucellvax", "flublok"))
    )
  attr(study, "params") <- tibble(
    n_per_group,
    logbaseline_mean,
    logtitre_sd,
    b0,
    blogbaseline,
    bfreq,
    bflucellvax,
    bflublok,
    bflucellvax_freq_add,
    bflublok_freq_add
  )
  study
}

plot_one_study <- function(...) {
  one_study(...) %>%
    select(id, freq_lab, vac, baseline, postvax) %>%
    pivot_longer(
      c(baseline, postvax),
      names_to = "timepoint", values_to = "titre"
    ) %>%
    ggplot(aes(timepoint, titre)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    scale_y_log10("Titre", breaks = 5 * 2^(0:10)) +
    scale_x_discrete(expand = c(0.1, 0)) +
    facet_grid(freq_lab ~ vac) +
    geom_line(aes(group = id), alpha = 0.3) +
    geom_point()
}

plot_diffs_one_study <- function(...) {
  one_study(...) %>%
    mutate(diff = exp(logpostvax_mid - logbaseline_mid)) %>%
    ggplot(aes(vac, diff)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank()
    ) +
    scale_y_log10(
      "Fold-change postvax vs baseline",
      breaks = c(1, 2, 3, 5, 10, 20, 40)
    ) +
    facet_wrap(~freq_lab) +
    geom_boxplot()
}

save_plot <- function(plot, name, ...) {
  ggdark::ggsave_dark(
    file.path(sample_size_dir, paste0(name, ".pdf")), plot,
    units = "cm",
    ...
  )
}

fit_one_study <- function(...) {
  study <- one_study(...)
  fit <- lm(logpostvax_mid ~ logbaseline_mid + freq * vac, study)

  # NOTE(sen) Work out the difference between recombinant and cell since it's
  # not an explicit parameter
  diff <- coef(fit)[["vacflublok"]] - coef(fit)[["vacflucellvax"]]
  var_names_of_interest <- c("vacflublok", "vacflucellvax")
  diff_var <- sum(vcov(fit)[var_names_of_interest, var_names_of_interest])
  diff_sd <- sqrt(diff_var)
  diff_statistic <- diff / diff_sd
  df_residual <- df.residual(fit)
  diff_p <- pt(abs(diff_statistic), df_residual, lower.tail = FALSE) * 2

  fit %>%
    broom::tidy() %>%
    bind_rows(tibble(
      term = "diffblokcell",
      estimate = diff,
      std.error = diff_sd,
      statistic = diff_statistic,
      p.value = diff_p
    )) %>%
    bind_cols(attr(study, "params")) %>%
    mutate(df_residual = df_residual)
}

fit_many_studies <- function(nsim = 20, ...) {
  map_dfr(1:nsim, function(i) fit_one_study(...) %>% mutate(i = i))
}

summ_many_studies <- function(nsim = 20,
                              n_per_group = 50,
                              bflublok = 0.5,
                              ...) {
  detections <- fit_many_studies(nsim, n_per_group = n_per_group, bflublok = bflublok, ...) %>%
    filter(term %in% c("vacflublok", "vacflucellvax", "diffblokcell")) %>%
    # NOTE(sen) One-sided test can probably be justified as there is evidence
    # that recombinant is better
    mutate(
      detect_one = (estimate > 0) & ((p.value / 2) < 0.05),
      detect_two = (estimate > 0) & (p.value < 0.05)
    )
  detections_add <- detections %>%
    # NOTE(sen) What's after `i` has to be the same as what's after `term` in
    # `summarise` below
    group_by(i, n_per_group, bflublok, bflucellvax) %>%
    summarise(
      .groups = "drop",
      detect_one = detect_one[term == "vacflublok"] & detect_one[term == "diffblokcell"],
      detect_two = detect_two[term == "vacflublok"] & detect_two[term == "diffblokcell"],
      term = "vacflublok+diffblokcell"
    )
  detections %>%
    bind_rows(detections_add) %>%
    group_by(term, n_per_group, bflublok, bflucellvax) %>%
    summarise(across(contains("detect"), ~ sum(.x) / n()), .groups = "drop")
}

save_data <- function(data, name) {
  write_csv(data, file.path(sample_size_dir, paste0(name, ".csv")))
}

# SECTION Script

plot_one_study() %>%
  save_plot("one", width = 15, height = 10)

plot_diffs_one_study(n_per_group = 5000) %>%
  save_plot("one-diff", width = 15, height = 10)

params <- tribble(
  ~n_per_group, ~bflublok,
  # NOTE(sen) bflublok one-sided
  50, 0.5,
  55, 0.5,
  60, 0.5,
  # NOTE(sen) bflublok two-sided
  60, 0.5,
  70, 0.5,
  80, 0.5,
  # NOTE(sen) bflublok one-sided lower
  200, 0.25,
  220, 0.25,
  240, 0.25,
  # NOTE(sen) Diff for one-sided
  440, 0.5,
  460, 0.5,
  480, 0.5,
  # NOTE(sen) Diff for two-sided
  600, 0.5,
  620, 0.5,
  640, 0.5,
)

sims <- future_pmap_dfr(
  params,
  function(n_per_group, bflublok) {
    summ_many_studies(nsim = 5000, n_per_group, bflublok)
  },
  .options = furrr_options(seed = 1)
)

save_data(sims, "sims")
