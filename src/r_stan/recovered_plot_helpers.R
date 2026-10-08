# Recovered with dput from stan_full_workspace0915.RData.
# DEFAULT_FIG_PATH is supplied by the portable caller; definitions are unchanged.

save_plot <- 
function (plot, filename, width = 8.5, height = 7, units = "cm", 
    dpi = 600, path = DEFAULT_FIG_PATH, bg = "white", save_png = TRUE, 
    save_pdf = FALSE) 
{
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }
    if (save_png) {
        ggplot2::ggsave(filename = paste0(filename, ".png"), 
            plot = plot, path = path, width = width, height = height, 
            units = units, dpi = dpi, bg = bg)
    }
    if (save_pdf) {
        ggplot2::ggsave(filename = paste0(filename, ".pdf"), 
            plot = plot, path = path, width = width, height = height, 
            units = units, device = grDevices::pdf, bg = bg)
    }
}

save_journal_plot <- 
function (plot, filename, width = 8.5, height = 7, units = "cm", 
    dpi = 600, path = DEFAULT_FIG_PATH, save_png = TRUE, save_pdf = FALSE) 
{
    save_plot(plot = plot, filename = filename, width = width, 
        height = height, units = units, dpi = dpi, path = path, 
        bg = "white", save_png = save_png, save_pdf = save_pdf)
}

