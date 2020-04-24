#' Create plasticity phenotypes
#'
#' Plots a series of stacked, colour-coded horizontal bars for each group. Color of each bar represents the group mean for the designated feature. Produces a phenotype for each group. Note: For all arguments matrix column numbers are identified as N-1 (e.g. column 4 is indexed as 3).
#' @param means A data frame with the first column containing group labels and the remaining columns containing the average values of each significant plasticity feature across all groups; this function requires all sums to be grouped together and all indices to be grouped together. Note: Features that describe variable sums must precede features that describe variable indices (e.g. sums = columns 2-4, indices = columns 5 – 10).
#' @param first_mean_indices_column Numerical value indicating the column number that marks the first feature index (plotted as green-yellow-red); all feature columns before this should be feature sums (plotted as grey-black).
#' @param feature_order Numeric vector that describes the order of features to be plotted. The order of values in the vector represent the column numbers for each feature (c(1,3,4,5,2,6,7)). Default order is taken from the data matrix.
#' @param grouping_var X-axis title. Default is "Experimental Conditions".
#' @param group_select Numberic vector that restricts the phenotype to plot only a subset of group phenotypes. Groups phenotypes to be plotted are indicated using their row numbers (c(1,3,2,4). Default order taken from means dataframe row number (C(1,2,3,4)), unless groups are assigned as factor.
#' @param indices_colors Character vector containing the color codes for indices color palette. Order is low, middle, high. Default is c("red", "yellow", "green").
#' @param sums_colors Character vector containing the color codes for sums color palette. Order is low, middle, high. Default is c( "#d3d3d3","grey","black").
#' @return Plots phenotype for each group. Stores object plasticity.cols in your Global Environment. This object will be consulted using the phenotype_boxplots function for each feature across each group.
#' @export


plasticity_phenotype <- function(means = means,
                                 first_mean_indices_column =  first_mean_indices_column,
                                 feature_order = feature_order,
                                 grouping_var = grouping_var,
                                 group_select = group_select,
                                 indices_colors = indices_colors,
                                 sums_colors = sums_colors
                                ){

  maxes <- sapply(means[,2:ncol(means)],max)

  maxes

  mins <- sapply(means[,2:ncol(means)],min)

  mins

  last_mean_sums_column <- first_mean_indices_column

  first_mean_indices_column <- first_mean_indices_column + 1

  grouping_mean_column <- 1

  first_mean_sums_column  <- 2

  first_max_sums_vector <- first_mean_sums_column - 1

  last_max_sums_vector <-  last_mean_sums_column - 1

  ##
  if(missing(indices_colors)) {
    indices_colors <- c("red", "yellow", "green")
  } else {
    indices_colors <- indices_colors
  }
  ##
  ##
  if(missing(sums_colors)) {
    sums_colors <- c( "#d3d3d3","grey","black")
  } else {
    sums_colors <- sums_colors
  }
  ##

  while (T) {

    library(tidyverse)

    nfeats <- length(maxes)

    feat.cols <- list()
    phenotype.cols <- data.frame(matrix(ncol = length(means[,grouping_mean_column]),
                                        nrow = ncol(means[,2:ncol(means)])))

    rownames(phenotype.cols) <- colnames(means[,
                                               2:ncol(means)])

    colnames(phenotype.cols) <- means[,grouping_mean_column]

    for (i in 1:length(colnames(means[,
                                      first_mean_indices_column:ncol(means)]))) {

      feat.min<-mins[last_max_sums_vector+i]
      feat.max<-maxes[last_max_sums_vector+i]
      feat.range<-(feat.max - feat.min)/2
      feat.mid<- feat.max-feat.range
      feat.col<-scale_colour_gradient2(
        low  = indices_colors[1],
        mid  = indices_colors[2],
        high = indices_colors[3],
        midpoint=feat.mid,
        breaks=c(feat.min,
                 feat.mid,
                 feat.max),
        labels=c("Below Normal",
                 "Normal",
                 "Above Normal"),
        limits=c(feat.min,feat.max))

      phenotype.cols[last_max_sums_vector+i,] <- feat.col$map(
        means[,i +
                last_mean_sums_column])

    }

    break

  }

  maxes <- maxes[first_max_sums_vector:last_max_sums_vector]

  mins <- mins[first_max_sums_vector:last_max_sums_vector]

  while (T) {


    for (i in 1:length(colnames(means[,first_mean_sums_column:last_mean_sums_column]))) {

      feat.min<-mins[i]
      feat.max<-maxes[i]
      feat.range<-(feat.max - feat.min)/2
      feat.mid<- feat.max-feat.range
      feat.col<-scale_colour_gradient2(
        low  = sums_colors[1],
        mid  = sums_colors[2],
        high = sums_colors[3],
        midpoint=feat.mid,
        breaks=c(feat.min,
                 feat.mid,
                 feat.max),
        labels=c("Below Normal",
                 "Normal",
                 "Above Normal"),
        limits=c(feat.min,
                 feat.max))


      phenotype.cols[i,] <- feat.col$map(means[,
                                               i+grouping_mean_column])

    }

    print(phenotype.cols)
    phenotype.cols <<- phenotype.cols
    break
    #colnames(phenotype.cols)

  }



  ##

  if(missing(group_select)){
    group_select <- 1:ncol(phenotype.cols)
  } else {
    group_select <- group_select
  }

  rown <-   rownames(phenotype.cols)
  coln <- colnames(phenotype.cols)[group_select]

  phenotype.cols <- as.data.frame(phenotype.cols[,group_select])

  rownames(phenotype.cols) <- rown
  colnames(phenotype.cols) <- coln
  ##
  phenotype.cols$index <- row.names(phenotype.cols)

  # ##
  if(missing(feature_order)) {
    feature_order <- 1:nfeats
  } else {
    feature_order <- feature_order
  }

  phenotype.cols2 <- phenotype.cols

  ##
  ##
  if(missing(grouping_var)) {
    grouping_var2 <- "\nExperimental Conditions"
  } else {
    grouping_var2 <-   grouping_var
  }
  ##

  long.phenotype.cols <- gather(phenotype.cols2,
                                key = agebin,
                                value = color,
                                -index )

  long.phenotype.cols$agebin <- factor(long.phenotype.cols$agebin,levels = means[,grouping_mean_column])


  ggplot(long.phenotype.cols,
         aes(x = agebin ,
             y = index)) +
    geom_tile(data = long.phenotype.cols,
              width=0.95,
              height=0.95,
              aes(fill = color))+
    scale_fill_identity()+
    theme(panel.background = element_rect(fill = 'white'),
          axis.text.x = element_text(angle=45,
                                     vjust=0.5,
                                     size=11),
          axis.text.y = element_text(angle=0,
                                     vjust=0.5,
                                     size=11),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size=24),
          axis.title.y = element_text(size=24),
          plot.margin = margin(2, 2, 2, 2, "cm")
    )+
    scale_x_discrete(expand=c(0,0),
                     name=grouping_var2,
                     position="bottom")+
    scale_y_discrete(expand=c(0,0),
                     name="Plasticity Features\n",
                     limits= rev(colnames(means)[
                       feature_order+1]),
                     position="left")+
    coord_fixed(ratio=0.25)

}

