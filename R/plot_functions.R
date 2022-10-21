#' Barplot of Shapley values
#'
#' @param x A list of class \code{shapley}.
#' @param subset Either an integer, \code{"chi2"}, or \code{NULL} (default) to select which rows of \code{phi} should be displayed.
#' If \code{NULL}, all \eqn{n} rows are displayed,
#' for a single integer the \code{subset} rows with the highest Mahalanobis distance are displayed,
#' for an integer vector the \code{subset} selected rows are displayed,
#' and for \code{"chi2"} all outlying rows are displayed (Mahalanobis distance greater than \eqn{\sqrt{}}\code{qchisq(chi2.q,p)}).
#' @param chi2.q Quantile, only used if \code{subset == "chi2"}.
#' @param abbrev.var Integer. If \code{abbrev.var} \eqn{> 0}, column names are abbreviated using abbreviate with \code{minlenght = abrev.var}.
#' @param abbrev.obs Integer. If \code{abbrev.obs} \eqn{> 0}, row names are abbreviated using abbreviate with \code{minlenght = abrev.obs}.
#' @param sort.var Logical. If \code{TRUE} (default), variables are sorted according to the  distance
#' @param sort.obs Logical. If \code{TRUE} (default), observations are sorted according to their Mahalanobis distance.
#' @param plot_md Logical. If \code{TRUE} (default), the Mahalanobis distance will be included in the plot.
#' @param md_squared Logical. If \code{TRUE} (default), the squared Mahalanobis distance is plotted otherwise the (not-squared) Mahalanobis distance.
#' @param rotate_x Logical. If \code{TRUE} (default), the x-axis labels are rotated.
#' @param ... Optional arguments passed to methods.
#'
#' @return Returns a barplot that displays the Shapley values (\code{\link{shapley}})for each observation and optionally (\code{plot_md = TRUE})
#' includes the squared Mahalanobis distance (black bar) and the corresponding (non-)central chi-square quantile (dotted line).
#' @export
#'
#' @examples
#' library(MASS)
#' set.seed(1)
#' n <- 100; p <- 10
#' mu <- rep(0,p)
#' Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
#' X <- mvrnorm(n, mu, Sigma)
#' X_clean <- X
#' X[sample(1:(n*p), 100, FALSE)] <- rep(c(-5,5),50)
#' call_shapley <- shapley(X, mu, Sigma)
#' plot(call_shapley, subset = 1:20)
#' plot(call_shapley, subset = 5, rotate_x = FALSE)
#' plot(call_shapley, subset = 5, md_squared = FALSE, rotate_x = FALSE)
plot.shapley <- function(x, subset = NULL, chi2.q = 0.99,
                         abbrev.var = 3, abbrev.obs = 10,
                         sort.var = FALSE, sort.obs = FALSE,
                         plot_md = TRUE, md_squared = TRUE, rotate_x = TRUE,...){
  phi <- x$phi
  non_centrality <- x$non_centrality
  bar_order <- fill <- name <- rowname <- upper <- value <- x <- y <- crit <- NULL #avoid warnings

  if(abbrev.obs > 0 & !is.null(rownames(phi))){
    rownames(phi) <- abbreviate(rownames(phi), minlength = abbrev.obs)
  }
  if(abbrev.var > 0 & !is.null(colnames(phi))){
    colnames(phi) <- abbreviate(colnames(phi), minlength = abbrev.var)
  }

  if(is.null(rownames(phi))){
    rownames(phi) <- paste("Obs. ", formatC(1:nrow(phi), width=nchar(nrow(phi)), flag="0"), sep="")
  }
  if(is.null(colnames(phi))){
    colnames(phi) <- paste("X", formatC(1:ncol(phi), width=nchar(ncol(phi)), flag="0"), sep="")
  }

  n = nrow(phi); p = ncol(phi)
  if(md_squared){
    a <- data.frame(phi)
    crit <- qchisq(chi2.q,p)
  } else{
    phi_perc <- t(apply(phi,1,function(x) x/sum(x)))
    a <- data.frame(sqrt(apply(phi,1,sum))*phi_perc)
    crit <- sqrt(qchisq(chi2.q,p))
  }
  MD <- apply(a,1,sum)
  #choose outlyingness criteria
  if(is.null(subset)){
    isout <- 1:n
  } else if(is.numeric(subset)){
    if(length(subset) == 1){
      isout <- order(MD,decreasing = TRUE)[1:subset]
    } else{
      isout <- subset
    }
  } else if(subset == "chi2"){
    isout <- which(MD>crit)
  } else{
    isout <- 1:n
  }

  #sort features for plot
  if(sort.var){
    features_sort <- names(sort(colSums(a[isout,]),decreasing = TRUE))
    features_sort_inc <- names(sort(colSums(a[isout,])))
  } else{
    features_sort <- colnames(a[isout,])
    features_sort_inc <- colnames(a[isout,])
  }
  if(sort.obs){
    observation_sort <- names(sort(rowSums(a[isout,]), decreasing = TRUE))
  } else{
    observation_sort <- rownames(a[isout,])
  }

  colorCount = ncol(phi)
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  if(colorCount>12){
    colorPalette <- getPalette(colorCount)
  } else{
    colorPalette <- getPalette(12)
  }

  a_barplot <- a[isout,] %>%
    arrange(desc(across(matches(features_sort)))) %>%
    rownames_to_column() %>%
    pivot_longer(cols = !rowname) %>%
    transmute(x = factor(rowname, levels = observation_sort),
              fill = factor(name, levels = features_sort_inc),
              y = value,
              sign = sign(y)) %>%
    arrange(x, desc(sign), abs(y)) %>%
    mutate(bar_order = seq_along(y))

  if(any(levels(a_barplot$fill) != colnames(phi))){
    levels(a_barplot$fill) <- colnames(phi)
  }

  if(plot_md){
    a_arrows <- data.frame(t(apply(a[isout,],1,function(x) {c("lower" = sum(x[x<0]), "upper" = sum(x))}))) %>%
      rownames_to_column() %>%
      mutate(x = factor(rowname, levels = observation_sort))
    if(is.null(non_centrality)){
      md_color = "chi-square\nquantile"
      if(md_squared){
        quantile_chisq_df <- data.frame(x = rownames(a[isout,]), y = crit)
      } else{
        quantile_chisq_df <- data.frame(x = rownames(a[isout,]), y = crit)
      }
    } else{
      md_color = "non-central\nchi-square\nquantile"
      if(md_squared){
        quantile_chisq <- sapply(non_centrality[isout], function(x)qchisq(p = chi2.q, df = p, ncp = x))
      } else{
        quantile_chisq <- sapply(non_centrality[isout], function(x){sqrt(qchisq(p = chi2.q, df = p, ncp = x))})
      }
      quantile_chisq_df <- data.frame(x = rownames(a[isout,]), y = quantile_chisq)
    }
  }

  y_axis <- "Decomposition of MD"
  if(md_squared){
    y_axis <- expression("Decomposition of MD "^2)
  }

  plt <- ggplot(data = NULL) +
    geom_bar(data = a_barplot, mapping = aes(fill = fill, y=y, x=x, group = bar_order),
             position="stack", stat="identity", width = 0.9) +
    geom_hline(yintercept = 0, linetype = 1) +
    theme_light() +
    theme(text = element_text(size=22), panel.grid.major.x = element_blank()) +
    labs(x = "Observation", y = y_axis, fill = "Feature") +
    scale_fill_manual(values = colorPalette) +
    labs(linetype = NULL) +
    scale_color_manual(name = NULL, values = c("black"))

  if(rotate_x){
    plt <- plt +
      theme(text = element_text(size=22),axis.text.x = element_text(angle = 38, vjust = 1, hjust=1))
  }

  if(plot_md){
    plt <- plt +
      geom_errorbar(data = quantile_chisq_df, aes(x = x, ymin = y, ymax = y, color = md_color), linetype = 2) +
      geom_segment(data = a_arrows, aes(x = x, xend = x, y = 0, yend = upper, linetype = "MD"),
                   size = 1, color = "black",
                   arrow = arrow(angle = 90, length = unit(0.1, "cm"), ends = "last", type = "open"))
    if(md_squared){
      plt <- plt + scale_linetype_manual(values = 1, name = "", labels = expression("MD"^2))
    }
  }
  return(plt)
}


#' Barplot and tileplot of Shapley values.
#'
#'
#' @param x A list of class \code{shapley_algorithm}.
#' @param type Either \code{"both"} (default), \code{"bar"}, or \code{"cell"}. If \code{"both"} (default) a barplot and a tileplot are created, otherwise only the selected plot is created.
#' @param n_digits Integer. If \code{n_digits}\eqn{ > 0}, the original values of the variables are given in each cell with \code{n_digits} decimals places.
#' @param continuous_rowname Logical. If \code{TRUE}, the rownames are converted to a numeric vector.
#' @param ... Arguments passed on to \code{\link{plot.shapley}}.
#' @inheritParams plot.shapley
#'
#' @return Returns plots for a list of class \code{shapley_algorithm}.
#' If \code{type} is \code{"bar"}, a barplot is generated. It displays the Shapley values (\code{\link{shapley}})
#' for each observation and optionally (\code{plot_md = TRUE}) includes the squared Mahalanobis distance (black bar)
#' and the corresponding (non-)central chi-square quantile (dotted line).
#' If \code{type} is \code{"cell"} a tileplot is generated. It displays each cells of the dataset and shows the original value from the observations,
#' color coding indicates whether those values were higher (red) or lower (blue) than the imputed values,
#' and the color intensity is based on the magnitude of the Shapley value.
#' If \code{type} is \code{"both"}, the barplot and the tileplot are generated.
#'
#' @export
#'
#' @examples
#' library(MASS)
#' set.seed(1)
#' n <- 100; p <- 10
#' mu <- rep(0,p)
#' Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
#' X <- mvrnorm(n, mu, Sigma)
#' X[sample(1:(n*p), 100, FALSE)] <- rep(c(-5,5),50)
#' MOE_X <- MOE(X, mu, Sigma)
#' plot(MOE_X, subset = 20, n_digits = 0)
plot.shapley_algorithm <- function(x, type = "both", subset = NULL,
                                   abbrev.var = FALSE, abbrev.obs = FALSE,
                                   sort.var = FALSE, sort.obs = FALSE,
                                   n_digits = 2, rotate_x = TRUE, continuous_rowname = FALSE, ...){
  if(type %in% c("both", "bar")){
    tmp <- new_shapley(phi = x$phi, mu_tilde = x$mu_tilde, non_centrality = x$non_centrality)
    plot_bar <- plot.shapley(tmp, subset = subset, abbrev.var = abbrev.var, abbrev.obs = abbrev.obs,
                             sort.var = sort.var, sort.obs = sort.obs, rotate_x = rotate_x, ...)
  }

  if(type %in% c("both", "cell")){
    rowname <- cell_outlier <- variable <- phi <- actual <- crit <- NULL #avoid warnings

    phi <- x$phi
    n = nrow(phi); p = ncol(phi)
    x_new <- matrix(x$x, nrow = n)
    x_original <- matrix(x$x_original, nrow = n)

    if(abbrev.obs > 0 & !is.null(rownames(phi))){
      rownames(phi) <- abbreviate(rownames(phi), minlength = abbrev.obs)
    }
    if(abbrev.var > 0 & !is.null(colnames(phi))){
      colnames(phi) <- abbreviate(colnames(phi), minlength = abbrev.var)
    }

    if(is.null(rownames(phi))){
      rownames(phi) <- paste("Obs. ", formatC(1:n, width=nchar(n), flag="0"), sep="")
    }
    if(is.null(colnames(phi))){
      colnames(phi) <- paste("X", formatC(1:p, width=nchar(p), flag="0"), sep="")
    }

    dimnames(x_new) <- dimnames(x_original) <- dimnames(phi)

    if(n > 1){
      MD <- apply(phi,1,sum)
      if(is.null(subset)){
        isout <- 1:n
      } else if(is.numeric(subset)){
        if(length(subset) == 1){
          isout <- order(MD,decreasing = TRUE)[1:subset]
        } else{
          isout <- subset
        }
      } else if(subset == "chi2"){
        isout <- which(MD>crit)
      } else{
        isout <- 1:n
      }

      x_new <- x_new[isout,]
      x_original <- x_original[isout,]
      phi <- phi[isout,]
    }

    if(sort.var){
      features_sort <- names(sort(colSums(phi),decreasing = TRUE))
      features_sort_inc <- names(sort(colSums(phi)))
    } else{
      features_sort <- colnames(phi)
      features_sort_inc <- colnames(phi)
    }

    #sort features for plot
    if(sort.obs){
      observation_sort <- names(sort(rowSums(phi), decreasing = TRUE))
    } else{
      observation_sort <- rownames(phi)
    }

    A1 <- as_tibble(round(x = x_original, digits = n_digits), rownames = NA) %>%
      rownames_to_column(var = "rowname") %>%
      pivot_longer(cols = -rowname, names_to = "variable", values_to = "actual")

    A2 <- as_tibble(sign(x_original - x_new), rownames = NA) %>%
      rownames_to_column(var = "rowname") %>%
      pivot_longer(cols = -rowname, names_to = "variable", values_to = "cell_outlier") %>%
      mutate(cell_outlier = fct_recode(factor(cell_outlier, levels = c(-1,0,1)), `low` = "-1", regular = "0", `high` = "1"))

    phi_rescaled <- t(apply(phi,1,function(x) x/(sum(x))))
    dimnames(phi_rescaled) <- dimnames(phi)
    A3 <- as_tibble(phi_rescaled, rownames = NA) %>%
      rownames_to_column(var = "rowname") %>%
      pivot_longer(cols = -rowname, names_to = "variable", values_to = "phi")

    phi_tiles <- inner_join(inner_join(A1,A2, by = c("rowname","variable")),A3, by = c("rowname","variable")) %>%
      mutate(rowname = factor(rowname, levels = observation_sort),
             variable = factor(variable, levels = features_sort_inc))

    if(continuous_rowname){
      phi_tiles <- phi_tiles %>% mutate(rowname = as.numeric(as.character(rowname)))
    }

    plt <- ggplot(phi_tiles, aes(x = rowname, y = variable, alpha = phi, fill = cell_outlier, label = actual)) +
      geom_tile(color = "lightgray") +
      coord_equal() +
      theme_minimal() +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        text = element_text(size=22),
        axis.title = element_blank()
      ) +
      guides(fill = guide_legend(title = "Cellwise\noutlier"),alpha = guide_legend(title = "Shapley\nvalue")) +
      scale_alpha(breaks = NULL)

    if(rotate_x){
      plt <- plt +
        theme(text = element_text(size=22),axis.text.x = element_text(angle = 38, vjust = 1, hjust=1))
    }

    if(all(c("low", "regular", "high") %in% phi_tiles$cell_outlier)){
      plt <- plt + scale_fill_manual(values = c("#1F78B4", "white", "#E31A1C"), labels = c("lower than\nexpected", "regular\ncell", "higher than\nexpected"))
    } else if(all(c("low", "regular") %in% phi_tiles$cell_outlier)){
      plt <- plt + scale_fill_manual(values = c("#1F78B4", "white"), labels = c("lower than\nexpected", "regular\ncell"))
    } else if(all(c("regular", "high") %in% phi_tiles$cell_outlier)){
      plt <- plt + scale_fill_manual(values = c("white", "#E31A1C"), labels = c("regular\ncell", "higher than\nexpected"))
    } else if(all(c("low", "high") %in% phi_tiles$cell_outlier)){
      plt <- plt + scale_fill_manual(values = c("#1F78B4", "#E31A1C"), labels = c("lower than\nexpected", "higher than\nexpected"))
    } else if(c("low") %in% phi_tiles$cell_outlier){
      plt <- plt + scale_fill_manual(values = c("#1F78B4"), labels = c("lower than\nexpected"))
    } else if(c("high") %in% phi_tiles$cell_outlier){
      plt <- plt + scale_fill_manual(values = c("#E31A1C"), labels = c("higher than\nexpected"))
    } else{
      plt <- plt + scale_fill_manual(values = c("white"), labels = c("regular\ncell"))
    }

    if(n_digits > 0){
      plt <- plt + geom_text(alpha = 1)
    }

    plot_cell <- plt
  }
  if(type == "both"){
    gridExtra::grid.arrange(plot_bar, plot_cell, ncol = 1)
  } else if(type == "cell"){
    plot_cell
  } else{
    plot_bar
  }
}


#' Plot of Shapley interaction indices
#'
#' @param x A \eqn{p \times p} matrix containing the Shapley interaction indices (\code{\link{shapley_interaction}}) of a single observation.
#' @param abbrev Integer. If \code{abbrev.var} \eqn{> 0}, variable names are abbreviated using abbreviate with \code{minlenght = abrev}.
#' @param title Character. Title of the plot.
#' @param legend Logical. If TRUE (default), a legend is plotted.
#' @param text_size Integer. Size of the text in the plot
#' @param ... Optional arguments passed to methods.
#'
#' @return Returns a figure consisting of two panels. The upper panel shows the Shapley values, and the lower panel the Shapley interaction indices.
#' @export
#'
#' @examples
#' p <- 5
#' mu <- rep(0,p)
#' Sigma <- matrix(0.9, p, p); diag(Sigma) = 1
#' Sigma_inv <- solve(Sigma)
#' x <- c(0,1,2,2.3,2.5)
#' PHI <- shapley_interaction(x, mu, Sigma)
#' plot(PHI)
plot.shapley_interaction <- function(x, abbrev = 4, title = "Shapley Interaction", legend = TRUE, text_size = 22,...){
  PHI <- x
  p <- ncol(PHI)
  check_matrix(PHI, p)

  rowname <- name <- value <- NULL #avoid warnings
  if(abbrev > 0 & !is.null(rownames(PHI)) & !is.null(colnames(PHI))){
    rownames(PHI) <- abbreviate(rownames(PHI), minlength = abbrev)
    colnames(PHI) <- abbreviate(colnames(PHI), minlength = abbrev)
  }

  if(is.null(rownames(PHI))){
    rownames(PHI) <- paste("X", formatC(1:nrow(PHI), width=nchar(nrow(PHI)), flag="0"), sep="")
  }
  if(is.null(colnames(PHI))){
    colnames(PHI) <- paste("X", formatC(1:ncol(PHI), width=nchar(ncol(PHI)), flag="0"), sep="")
  }



  plot_PHI <-  as_tibble(PHI, rownames = NA) %>%
    rownames_to_column() %>%
    pivot_longer(cols = -c(rowname)) %>%
    mutate(value = as.numeric(value))

  min_value <- min(PHI)
  max_value <- max(PHI)
  pp <- ggplot(plot_PHI, aes(x = name, y = rowname, fill = value)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue",mid = "white", high = "red", breaks=c(min_value,max_value),labels=c("low","high")) +
    theme_light() +
    theme(text = element_text(size=text_size),
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.x = element_text(angle = 45, vjust = 1.05, hjust=1),
          legend.position = "right") +
    labs(x = element_blank(), y = element_blank(), fill = "Shapley\ninteraction")

  if(!legend){
    pp <- pp + theme(legend.position = "none")
  }

  shval <- rowSums(PHI) #shapley values
  if(is.null(names(shval))){
    names(shval) <- paste("X", formatC(1:length(shval), width=nchar(length(shval)), flag="0"), sep="")
  }

  mp <- ggplot() +
    geom_bar(aes(x = names(shval),
                 y = shval),
             stat = "identity") +
    theme_light() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = text_size)
    ) +
    labs(y = "Sh.v.", title = title)+
    geom_hline(yintercept = 0)

  plt <- ggarrange(mp, pp , ncol = 1, heights = c(1,3))
}

