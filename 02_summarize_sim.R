## These functions provide ways to summarize simulation output in order to
## contextualize results in terms of e.g. number of fish handled or number of
## recaptures.

plot_obs <- function(res) {
  par(mfrow = c(1, 2))
  ## Collected in the river each observation event
  barplot(colSums(res$y == 2 | res$y == 3, na.rm = TRUE), main = "In river")
  ## Collected in the trap
  barplot(colSums(res$y == 4, na.rm = TRUE), main = "Trap")
}

number_handled <- function(res) {
  y <- res$y
  ## Total number of fish handled in river
  mark <- apply(y, 1, \(r) 2 %in% r | 3 %in% r) |> sum()
  ## Total number of fish trapped
  trap <- apply(y, 1, \(r) 4 %in% r) |> sum()
  ## Total number of *recaptures*
  recap <- apply(y, 1, \(r) sum(r > 1) > 1) |> sum()

  list(marked = mark, trapped = trap, recaptured = recap)
}

sim_summarize <- function(res) {
  tibble::tibble(
            origin = res$origin,
            sex = res$sex,
            resid = res$resid,
            final_state = res$zfull[, ncol(res$zfull)],
            n_capture = apply(res$y, 1, \(r) sum(r %in% c(2, 3))),
            captured = n_capture > 0,
            recaptured = n_capture > 1,
            trapped = final_state == 4
          )
  }

sim_handled <- function(res) {
  sim_summarize(res) |>
    dplyr::summarize(
      n = n(),
      n_cap = sum(captured),
      n_recap = sum(recaptured),
      n_trapped = sum(trapped),
      .by = c(final_state)
    )
}
