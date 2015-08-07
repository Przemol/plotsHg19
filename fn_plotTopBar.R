plotTopBar <- function(M, titles=LETTERS, limits=NULL) {
    #pdf(out, 6 * 4, 2 * 4)
    f = 37
    layout(
        matrix(c(1:100), 2, 4, byrow = FALSE), respect = TRUE,
        heights = c(.3, rep(1,5)), widths = c(1)
    )
    pp <- function() {
        par(mgp = c(0, 2.5, 0))
        axis(
            1, at = c(min(M[[1]]$all_ind), 0, max(M[[1]]$all_ind)),
            labels = c('-1kb', 'TSS', '+1.5kb'), cex.axis = f / 12
        )
        par(mgp = c(0, 1, 0))
    }
    plt <- function(p, ylim = NULL, ...) {
        for (i in 1:p$npaires()) {
            plot(
                p[i], ylim = ylim, keepratio = TRUE, legend = FALSE, ln.h = 0,
                cex.axis = f, cex.lab = f, ln.v = TRUE, xlab = '',
                cex.main = .01, panel.first = pp(), xaxt = "n", col = 'black', ...
            )
        }
    }
    
    par(mgp = c(3, 0, 0))
    
    plot.new(); title(xlab = titles[[1]], col.main = 'red', cex.lab = 6)
    plt(unlist(M[1,1]), ylim = limits[[1]], main = titles[[1]])
    
    par(mgp = c(3, 0, 0))
    
    plot.new(); title(xlab = titles[[2]], col.main = 'red', cex.lab = 6)
    plt(unlist(M[1,2]), ylim = limits[[2]], main = titles[[2]])
    
    par(mgp = c(3, 0, 0))
    
    plot.new(); title(xlab = titles[[3]], col.main = 'red', cex.lab = 6)
    plt(unlist(M[1,3]), ylim = limits[[3]], main = titles[[3]])
    
    par(mgp = c(3, 0, 0))
    
    plot.new(); title(xlab = titles[[4]], col.main = 'red', cex.lab = 6)
    plt(unlist(M[1,4]), ylim = limits[[4]], main = titles[[4]])
    
    #dev.off()
    
}
