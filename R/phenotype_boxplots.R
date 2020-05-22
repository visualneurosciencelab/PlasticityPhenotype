#' Create colour-coded boxplots for each feature across every group
#'
#' Plot individual boxplots for each feature, with color codes corresponding to the output of the plasticity_phenotype function.
#' @param feature_df The dataframe of features used to create phenotypes, containing group identifiers and the values of each feature.
#' @param phenotype_cols Taken from the plasticity_phenotype function; dataframe containing hexadecimal color-code for each feature.
#' @param first_index_column A numeric value that indicates which column number contains the first index feature (column number -1). All subsequent feature columns are assumed to be indices. All preceding feature columns are assumed to be sums.
#' @param group_label X-axis title.
#' @param max_sum A vector that indicates the y-axis maximum value for each feature sum.
#' @param point_size A numeric value that indicates the size of data points (Default value = 0.5 pt)
#' @param point_alpha A numeric value that indicates the transparency (alpha) of data points (Default value = 0.5 pt)
#' @param aspect_ratio A numeric value that indicates the aspect ratio of each plot (Default value= 5/7 or 0.714)
#' @param text_size A numeric value indicating text size
#' @return A series of individual boxplots depicting the expression of each feature across groups, colour-code can optionally be set to match phenotypes.
#' @examples
#' phenotype_boxplots(
#'      feature_df = merged.data2[,c("Condition",ordered_nf_colnames2)],
#'      phenotype_cols = phenotype.cols,
#'      max_sum = c(1.5,
#'           max(merged.data2$'VGLUT\n'),
#'           max(merged.data2$'Receptors Sum\n'),
#'           max(merged.data2$'Scaffolding Sum\n'),
#'           2),
#'      group_label = "\nRearing Conditions",
#'      first_index_column = 6,
#'      point_size = 1.5,
#'      point_alpha = 1,
#'      aspect_ratio = 5/7,
#'      text_size = 8
#'      )
#'
#' @export

phenotype_boxplots <- function(feature_df = feature_df,
                               phenotype_cols = phenotype.cols,
                               first_index_column = first_index_column,
                               group_label = group_label,
                               max_sum = max_sum,
                               point_size = point_size,
                               point_alpha = point_alpha,
                               aspect_ratio = aspect_ratio,
                               text_size = text_size){



  last_sum_column = first_index_column
  first_index_column = first_index_column + 1


  max.sum <- max(sapply(feature_df[,2:last_sum_column],max))
  min.sum <- min(sapply(feature_df[,2:last_sum_column],min))

  max.ind <- max(sapply(feature_df[,first_index_column:ncol(feature_df)],max))
  min.ind <- min(sapply(feature_df[,first_index_column:ncol(feature_df)],min))

  feat_names <- colnames(feature_df)[2:ncol(feature_df)]


  # ##
  #
  # if(missing(group_label)){
  #   max_sum2 <- max(feature_df[,2:last_sum_column])
  # } else {
  #   max_sum2 <-  max_sum[i-1]
  # }
  #
  # ##

  while (T) {

    #feature_df = feats_df3
    #first_index_column = 4
    #phenotype_cols = phenotype.cols.hum




    library(ggpubr)
    library(tidyverse)

    bp_list <- list()
    plot_list <- list()


    for (i in 2:last_sum_column){
      #i = 2
      bp_list[[i-1]] <- feature_df[,c(1,i)]

      current_bp <- bp_list[[i-1]]

      colvect <- NULL

      #feature_df = feats_df4
      #phenotype_cols = phenotype.cols.hum
      #first_index_column = 4

      for (j in 1:ncol(phenotype_cols)) {
        # j = 1

        colvect <- c(colvect,phenotype_cols[i-1,j])

      }

     if(missing(group_label)){
        group_label <- "Experimental Conditions"
     } else {
        group_label <- group_label
     }

      colnames(current_bp)[1] <- group_label
      colnames(current_bp)[2] <- "feature"

      ##

      if(missing(max_sum)){
        max_sum2 <- max(current_bp$feature)
      } else {
        max_sum2 <-  max_sum[i-1]
      }

      ##

      if(missing(aspect_ratio)){
        aspect_ratio = 5/7
      } else{
        aspect_ratio = aspect_ratio
      }

      if(missing(point_size)){
        point_size = 0.5
      } else{
        point_size = point_size
      }

      if(missing(point_alpha)){
        point_alpha = 0.5
      } else{
        point_alpha = point_alpha
      }

      if(missing(text_size)){
        text_size = 8
      } else{
        text_size = text_size
      }

      x  <-      ggplot(current_bp,
                        aes(x = !!ensym(group_label) ,
                            y = feature))+
        geom_boxplot(fill = colvect,
                     outlier.size = point_size,
                     outlier.alpha = point_alpha)+
        theme_classic()+
        ylab(feat_names[i -  1])+
        xlab(group_label)+
       # scale_x_discrete(label=abbreviate)+
        geom_jitter(alpha = point_alpha, size = point_size,
                    position = position_jitter(height = 0, width = 0))+
        coord_cartesian(ylim = c(0, max_sum2))+
        theme(axis.text.y = element_text(size = text_size),
              axis.text.x = element_text(size = text_size, angle = 45, hjust = 1),
          aspect.ratio = aspect_ratio,
          axis.title.x = element_text(size = 9#, face = 'bold'
          ),
          axis.title.y = element_text(size = 9# ,face = 'bold'
          )
          )+
        scale_y_continuous(labels = scales::number_format(accuracy = 0.1))


      plot_list[[i-1]] <- x

    }
    for (i in first_index_column:ncol(feature_df)){

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
        geom_boxplot(fill = colvect,
                     outlier.size = point_size,
                     outlier.alpha = point_alpha)+
        theme_classic()+
        ylab(feat_names[i -  1])+
        xlab(group_label)+
       # scale_x_discrete(label=abbreviate)+
        geom_jitter(alpha = point_alpha, size = point_size,
                    position = position_jitter(height = 0, width = 0))+
        theme(axis.text.y = element_text(size = text_size),
              axis.text.x = element_text(size = text_size ,angle = 45, hjust = 1),
              axis.title.x = element_text(size = 9#, face = 'bold'
                                          ),
              axis.title.y = element_text(size = 9# ,face = 'bold'
                                          ),
          aspect.ratio = aspect_ratio
         # aspect.ratio=4/3
          )+
        #coord_cartesian(ylim = c(min.ind, max.ind))
        scale_y_continuous(labels = scales::number_format(accuracy = 0.1))

      plot_list[[i-1]] <- x

    }

    plot_list <<- plot_list

    break

  }


  }

