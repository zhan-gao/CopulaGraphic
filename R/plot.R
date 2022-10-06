#' Plot the estimated survival function
#'
#' @import reshape2
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import ggpmisc
#' @import tibble
#'
#' @param time
#' @param surv
#' @param x_interval Interval length on x-axis
#'
#'
#' @export
#'
surv_plot <- function(time, surv, x_interval = 100) {
    if(length(time) != length(surv)) {
        stop("survival time and survival probability should of the same length.")
    }
    df <- data.frame(
        time = sort(time),
        prob = surv
    )

    # prepare the risk table
    br_pt <- seq(0, max(time) + x_interval, x_interval)
    m <- length(br_pt)
    risk_num <- rep(0, m)
    for (j in 1:m){
        br <- br_pt[j]
        risk_num[j] <- sum(df[, "time"] > br)
    }
    df_tab <- tibble(time = br_pt,
                     `Number at risk` = risk_num)
    df_tab <- melt(df_tab, id = "time")

    p <- ggplot(data = df) +
        geom_line(mapping = aes(
            x = time, y = prob
        )) +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        scale_x_continuous(breaks = br_pt, limits = c(0, max(time)), expand = c(0, 30)) +
        labs(x = "Time", y = "Survival Probability") +
        theme(
            panel.background =  element_blank(),
            panel.grid.major = element_line(linetype = 2, color = "grey90"),
            legend.title = element_blank(),
            legend.key = element_rect(colour = NA, fill = NA),
            axis.line = element_line(colour = "grey50"),
            axis.text = element_text(size = 11)
        )

    tab <- ggplot(df_tab) +
        geom_text(mapping = aes(x = time,
                                y = factor(variable),
                                label = value), size = 3.5) +
        theme_bw() +
        scale_y_discrete(limits = c("Number at risk")) +
        scale_x_continuous(breaks = br_pt, limits = c(0, max(time)), expand = c(0, 30)) +
        ggtitle("Number at risk") +
        theme(
            panel.background =  element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 10.5),
            axis.title = element_blank(),
            plot.title = element_text(size = 11, hjust = 0)
        )

    Layout <- grid.layout(nrow = 2, ncol = 1,
                          heights = unit(c(2, 0.25), c("null", "null")))
    grid.show.layout(Layout)
    vplayout <- function() {
        grid.newpage()
        pushViewport(viewport(layout = Layout))
    }
    subplot <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    mmplot <- function(a, b) {
        vplayout()
        print(a, vp = subplot(1, 1))
        print(b, vp = subplot(2, 1))
    }
    mmplot(p, tab)
}
