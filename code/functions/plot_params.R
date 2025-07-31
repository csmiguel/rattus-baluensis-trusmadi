plot_params <- function(params, rej, param_range, xlab, prob = 0.5, lw = 1) {
  cols <- RColorBrewer::brewer.pal(6, "Set1")
  op <- par(lwd = lw)  # set line width globally for the plot
  z <-
    lapply(seq_along(params), function(x) {
    param_i <- params[x]
    # regression
    abc_reg <-
      makepd4(target = target,
              x = rej[[param_i]],
              sumstat = as.matrix(rej_stats),
              tol = 1, 
              rejmethod = F,
              transf = "log",
              bb = param_range)
    # plot posterior
    if (x == 1) {
      plot(
        locfit(formula = ~abc_reg$x,
               weight = abc_reg$wt,
               xlim = param_range),
        xlab = xlab,
        col = cols[x],
        xlim = param_range)
    } else if (x > 1) {
      plot(
        locfit(formula = ~abc_reg$x,
               weight = abc_reg$wt,
               xlim = param_range),
        xlab = xlab,
        col = cols[x],
        xlim = param_range,
        add = T)
    }
   
    # plot prior
    plot(locfit(~priors[, which(dimnames(priors)[[2]] == param_i)],
                xlim = param_range),
         col = cols[x],
         lty = 3,
         add = T)
    param_stats <-
      loc1statsx(x = abc_reg$x,
                 wt = abc_reg$wt,
                 prob = prob)
    param_stats
    }) |>
    setNames(params)
  legend("topright",
         legend = params,
         col = cols[seq_along(params)],
         lty = 1, bty = "n")
  message("Mode and HPD for estimates:\n")
  on.exit(par(op))     # reset pars
  return(z)
}

