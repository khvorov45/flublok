cat("sample size calculations")

library(tidyverse)
library(furrr)

plan(multiprocess)

sample_size_dir <- "sample-size"

# Functions ===================================================================

# Returns OG scale titres
censor_logtitres <- function(logtitres) {
  cuts <- cut(exp(logtitres), c(-Inf, 5 * 2^(1:10), Inf)) %>% as.integer()
  5 * 2^(cuts - 1)
}

# Returns log titres
midpoint_titres <- function(titres) {
  if_else(titres == 5L, log(titres), log(titres) + log(2) / 2)
}

# Outcome - postvax titre
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
  groups %>%
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
  lm(logpostvax_mid ~ logbaseline_mid + freq * vac, one_study(...)) %>%
    broom::tidy()
}

fit_many_studies <- function(nsim = 20, ...) {
  future_map_dfr(1:nsim, function(i) fit_one_study(...) %>% mutate(i = i))
}

summ_many_studies <- function(n_per_group = 50, nsim = 20, ...) {
  fit_many_studies(nsim, n_per_group = n_per_group, ...) %>%
    filter(term == "vacflublok") %>%
    # One-sided test can probably be justified as there is evidence that
    # recombinant is better
    mutate(detect_flublok = (estimate > 0) & ((p.value / 2) < 0.05)) %>%
    summarise(
      detect_flublok = sum(detect_flublok) / n(),
      n_per_group = n_per_group,
      .groups = "drop"
    )
}

save_data <- function(data, name) {
  write_csv(data, file.path(sample_size_dir, paste0(name, ".csv")))
}

# Script ======================================================================

plot_one_study() %>%
  save_plot("one", width = 15, height = 10)

plot_diffs_one_study(n_per_group = 5000) %>%
  save_plot("one-diff", width = 15, height = 10)

sims <- map_dfr(c(50, 60, 70), ~ summ_many_studies(.x, nsim = 2000))

save_data(sims, "sims")
