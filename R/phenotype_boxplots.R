#' Create colour-coded boxplots for each feature across every group
#'
#' Plot individual boxplots for each feature, with color codes corresponding to the output of the plasticity_phenotype function.
#' @param feature_df The dataframe of features used to create phenotypes, containing group identifiers and the values of each feature.
#' @param phenotype_cols Taken from the plasticity_phenotype function; dataframe containing hexadecimal color-code for each feature.
#' @param group_label X-axis title.
#' @return A series of boxplots depicting the expression of each feature across groups, color-code matches phenotypes.
#' @export

phenotype_boxplots <- function(feature_df = feature_df,
                               phenotype_cols = phenotype.cols,
                               group_label = group_label){

  while (T) {

    feat_names <- colnames(feature_df)[2:ncol(feature_df)]

    bp_list <- list()
    plot_list <- list()


    for (i in 2:ncol(feature_df)){

      bp_list[[i-1]] <- feature_df[,c(1,i)]

      current_bp <- bp_list[[i-1]]

      colvect <- NULL

      for (j in 1:ncol(phenotype_cols)) {

        colvect <- c(colvect,phenotype_cols[i-1,j])

      }

      if(missing(group_label)){
        group_label <- "Experimental Conditions"
      } else {
        group_label <- group_label
      }

      colnames(current_bp)[1] <- group_label
      colnames(current_bp)[2] <- "feature"

      x  <-      ggplot(current_bp,
                        aes(x = !!ensym(group_label) ,
                            y = feature))+
        geom_boxplot(fill = colvect)+
        theme_classic()+
        ylab(feat_names[i -  1])+
        xlab(group_label)+
        scale_x_discrete(label=abbreviate)


      plot_list[[i-1]] <- x

    }

    plot_list <<- plot_list

    break

  }}
