#' Creates correlation heatmap between PCA basis vectors and feature matrix
#'
#' Creates a colour-coded correlation matrix between the PCA basis vectors (determined using cum_var) and the matrix of features (determined manually after using amplitude_plots).  Color-coding is done for only those correlations that achieve statistical significance after adjusting for multiple comparisons.
#' @param corr.scores.pval A matrix containing the p-values for each correlation coefficient after being adjusted for multiple comparisons; columns reflect plasticity features while rows reflect PCA basis vectors.
#' @param corr.scores.rval A matrix containing the correlation coefficients between each feature and PCA basis vector; columns reflect plasticity features while rows reflect PCA basis vectors.
#' @param thresh A threshold which each correlation coefficient’s p-value must be below to meet significance and retain their color-code (valid entries range from 0 – 1).
#' @return A colour-coded correlation matrix displaying the correlations between features (X-axis) and PCA basis vectors (Y-axis). The cells in the matrix are only color-coded and labelled with the coprresponding Pearson's R correlation coefficient when statistically significant. Red tiles reflect significant negative correlations while green tiles reflect significant positive correlations. Grey tiles are not statistically significant.
#' @examples
#' feature_matrix(
#'      corr.scores.pval = corr.scores.bfpval,
#'      corr.scores.rval = corr.scores.rval,
#'      thresh = 0.05
#'       )
#' @export

feature_matrix <- function(corr.scores.pval,
                           corr.scores.rval,
                           thresh) {

  library(tidyverse)

  corr.scores.pval2 <- corr.scores.pval
  corr.scores.pval2[corr.scores.pval2 > (thresh)] <- NA

  while (T) {

    for (i in 1:nrow(corr.scores.pval)){

      for (j in 1:ncol(corr.scores.pval)) {

        if (is.na(corr.scores.pval2[i,j]) == TRUE ) {

          corr.scores.rval[i,j] <- NA

        } else {

          corr.scores.rval[i,j] <- round(corr.scores.rval[i,j],3)

        }

      }

    }

    break

  }

  corr.scores.rval <- as.data.frame(corr.scores.rval)

  corr.scores.rval$dim <- rownames(corr.scores.rval)

  corr.scores.rval2 <- gather(corr.scores.rval, feature, value, -dim)

  names(corr.scores.rval2) <- c("Dimension", "variable", "value")

  index.ord <- colnames(corr.scores.rval)[-ncol(corr.scores.rval)]

  ggplot(corr.scores.rval2,
         aes(x = variable,
             y = Dimension)) +
    geom_tile(data = corr.scores.rval2,
              aes(fill = value,
                  width = 0.95,
                  height = 0.95))+
    scale_fill_continuous(low = "red",
                          high = "darkgreen",
                          limits = c(-1, 1),
                          breaks = c(-1,0,1),
                          labels = c(-1,0,1),
                          na.value = "grey",
                          name = "Correlation")+
    theme(panel.background = element_rect(fill = 'gray95'),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14,face = "bold"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          plot.margin = margin(2, 2, 2, 2, "cm"))+
    scale_x_discrete(expand = c(0,0),
                     name = "\nProtein Indices" ,
                     limits = index.ord
                     )+
    scale_y_discrete(expand = c(0,0),
                     limits = rev(unique(sort(corr.scores.rval2$Dimension))),
                     name = "Dimension\n")+
    coord_fixed(ratio = 1) +
    geom_text(aes(label = value))

}
