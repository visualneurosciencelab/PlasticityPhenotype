#' Creates correlation heatmap between PCA basis vectors and feature matrix
#'
#' Creates a correlation matrix between the identified PCA basis vectors (using cum_var) and the matrix of features (assigned by user), then plots only the significant correlations in a heatmap
#' @param corr.scores.pval The name of the matrix containing bonferroni adjusted p-values
#' @param corr.scores.rval The name of the matrix containing the correlation coefficients
#' @param thresh The significance threshold value, which identifies the significance value of Pearson's R correlations to be inclded in the matrix
#' @return Correlation matrix plotted as a heatmap, with filled cells representing only the significantly correlated features
#' @export

feature_matrix <- function(corr.scores.pval,
                           corr.scores.rval,
                           thresh) {

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
          axis.text.x = element_text(angle = 65,hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_x_discrete(expand = c(0,0),
                     name = "Protein Indices" ,
                     limits = index.ord
                     )+
    scale_y_discrete(expand = c(0,0),
                     limits = rev(unique(sort(corr.scores.rval2$Dimension)))
                     )+
    coord_fixed(ratio = 1) +
    geom_text(aes(label = value
    ))+
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14,
                                    face = "bold"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          plot.margin = margin(0, 0, 0, 0, "cm"))

}
