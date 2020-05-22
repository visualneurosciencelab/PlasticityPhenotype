#' Create plasticity phenotypes
#'
#' Plots a series of stacked, colour-coded horizontal bars for each group. Color of each bar represents the group mean for the designated feature. Produces a phenotype for each group. Note: For all arguments matrix column numbers are identified as N-1 (e.g. column 4 is indexed as 3).
#' @param phenotype_data A data frame with the first column containing group labels and the remaining columns containing the average values of each significant plasticity feature across all groups; this function requires all sums to be grouped together and all indices to be grouped together. Note: Features that describe variable sums must precede features that describe variable indices (e.g. sums = columns 2-4, indices = columns 5 – 10).
#' @param first_index_column Numerical value indicating the column number that marks the first feature index (plotted as green-yellow-red); all feature columns before this should be feature sums (plotted as grey-black).
#' @param feature_order Numeric vector that describes the order of features to be plotted. The order of values in the vector represent the column numbers for each feature (c(1,3,4,5,2,6,7)). Default order is taken from the data matrix.
#' @param group_label X-axis title. Default is "Experimental Conditions".
#' @param group_select Numberic vector that restricts the phenotype to plot only a subset of group phenotypes. Groups phenotypes to be plotted are indicated using their row numbers (c(1,3,2,4). Default order taken from means dataframe row number (C(1,2,3,4)), unless groups are assigned as factor.
#' @param indices_colors Character vector containing the color codes for indices color palette. Order is low, middle, high. Default is c("red", "yellow", "green").
#' @param sums_colors Character vector containing the color codes for sums color palette. Order is low, middle, high. Default is c( "#d3d3d3","grey","black").
#' @param translation Assign min/max data values for colour scales (feature indices only). Four options include: 'local', 'global', 'absolute', 'scaled'.
#' local = local min/max determined within each feature. All features will have an absolute min (default=red) and max (default=green).
#' global= global min/max determined across ALL features. There will only be one occurence of the global min (default=red) AND the global max (default=green) across all features.
#' absolute= global absolute min/max determined across ALL features. There will only be one occurence of the global min (default=red) OR the global max (default=green) across all features. Min/max determined by greatest ABSOLUTE value.
#' scaled= scaled min/max are set to -1 and +1 respectively. This setting ignores min/max in the data.
#' Default= 'local'.
#' @return Plots phenotype for each group. Stores object plasticity.cols in your Global Environment. This object will be consulted using the phenotype_boxplots function for each feature across each group. Stores object plas.phen in your Global Environment. plas.phen contains ggplot visualization of each comparison.
#' @examples
#'plasticity_phenotype(
#'     phenotype_data = meds,
#'     first_mean_index_column =  6,
#'     group_label = "\nRearing Conditions",
#'     translation = 'absolute',
#'     group_select = 1:4,
#'     indices_colors = c('red','yellow','green'),
#'     sums_colors = c('white','grey','black'),
#'     feature_order = 1:9
#'     )
#' @export


plasticity_phenotype <- function(phenotype_data = phenotype_data,
                                 first_index_column =  first_index_column,
                                 feature_order = feature_order,
                                 group_label = group_label,
                                 group_select = group_select,
                                 indices_colors = indices_colors,
                                 sums_colors = sums_colors,
                                 translation = 'local'
                                ){

  if (translation == "local"){

  maxes <- sapply(phenotype_data[,2:ncol(phenotype_data)],max)

  maxes

  mins <- sapply(phenotype_data[,2:ncol(phenotype_data)],min)

  mins

  last_mean_sums_column <- first_index_column

  first_index_column <- first_index_column + 1

  grouping_mean_column <- 1

  first_mean_sums_column  <- 2

  first_max_sums_vector <- first_mean_sums_column - 1

  last_max_sums_vector <-  last_mean_sums_column - 1
  #
  max.sum <- max(sapply(phenotype_data[,first_mean_sums_column:last_mean_sums_column],max))
  min.sum  <- 0
  #
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
    phenotype.cols <- data.frame(matrix(ncol = length(phenotype_data[,grouping_mean_column]),
                                        nrow = ncol(phenotype_data[,2:ncol(phenotype_data)])))

    rownames(phenotype.cols) <- colnames(phenotype_data[,
                                               2:ncol(phenotype_data)])

    colnames(phenotype.cols) <- phenotype_data[,grouping_mean_column]

    for (i in 1:length(colnames(phenotype_data[,
                                      first_index_column:ncol(phenotype_data)]))) {

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
        phenotype_data[,i +
                last_mean_sums_column])

    }

    break

  }

  maxes <- maxes[first_max_sums_vector:last_max_sums_vector]

  mins <- mins[first_max_sums_vector:last_max_sums_vector]

  while (T) {


    for (i in 1:length(colnames(phenotype_data[,first_mean_sums_column:last_mean_sums_column]))) {

      #i = 1

      #feat.min<-mins[i]
      feat.max<-maxes[i]

      # feat.max <- max.sum
      feat.min <- min.sum
      feat.range<-(feat.max - feat.min)/2
      feat.mid<- feat.max-feat.range
      #feat.mid <- max.sum/2
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


      phenotype.cols[i,] <- feat.col$map(phenotype_data[,
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
  if(missing(group_label)) {
    group_label2 <- "\nExperimental Conditions"
  } else {
    group_label2 <-   group_label
  }
  ##

  long.phenotype.cols <- gather(phenotype.cols2,
                                key = agebin,
                                value = color,
                                -index )

  long.phenotype.cols$agebin <- factor(long.phenotype.cols$agebin,levels = phenotype_data[,grouping_mean_column])



  }
  else if (translation == 'global'){
    ########################################



    ###

    maxes <- sapply(phenotype_data[,2:ncol(phenotype_data)],max)

    mins <- sapply(phenotype_data[,2:ncol(phenotype_data)],min)

    ###

    last_mean_sums_column <- first_index_column

    first_index_column <- first_index_column + 1

    grouping_mean_column <- 1

    first_mean_sums_column  <- 2

    first_max_sums_vector <- first_mean_sums_column - 1

    last_max_sums_vector <-  last_mean_sums_column - 1

    ###

    max.ind <- max(sapply(phenotype_data[,first_index_column:ncol(phenotype_data)],max))
    mid.ind = 0
    min.ind <- min(sapply(phenotype_data[,first_index_column:ncol(phenotype_data)],min))

    max.sum <- max(sapply(phenotype_data[,first_mean_sums_column:last_mean_sums_column],max))
    min.sum  <- 0

    ###

    #indices_colors <- c("red", "yellow", "green")
    #sums_colors <- c( "white","grey","black")

    ###

    ##
    if(missing(indices_colors)) {
      indices_colors <- c("red", "yellow", "green")
    } else {
      indices_colors <- indices_colors
    }
    ##
    ##
    if(missing(sums_colors)) {
      sums_colors <- c( "white","grey","black")
    } else {
      sums_colors <- sums_colors
    }
    ##
    ###

    while (T) {

      library(tidyverse)

      nfeats <- length(maxes)

      feat.cols <- list()
      phenotype.cols <- data.frame(matrix(ncol = length(phenotype_data[,grouping_mean_column]),
                                          nrow = ncol(phenotype_data[,2:ncol(phenotype_data)])))

      rownames(phenotype.cols) <- colnames(phenotype_data[,
                                                 2:ncol(phenotype_data)])

      colnames(phenotype.cols) <- phenotype_data[,grouping_mean_column]

      for (i in 1:length(colnames(phenotype_data[,
                                        first_index_column:ncol(phenotype_data)]))) {

        #i = 1

        ##feat.range<-(feat.max - feat.min)/2
        ##feat.mid<- feat.max-feat.range


        #feat.min<-mins[last_max_sums_vector+i]
        #feat.max<-maxes[last_max_sums_vector+i]
        #feat.mid <-  max.ind - ((max.ind - min.ind)/2)

        feat.max <- 1
        feat.min  <- -1
        feat.mid <- mid.ind

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
          #limits=c(feat.min,feat.max)
          limits=c(min.ind,max.ind)
        )

        phenotype.cols[last_max_sums_vector+i,] <- feat.col$map(
          phenotype_data[,i +
                  last_mean_sums_column])

        if(max.ind %in% phenotype_data[,i + last_mean_sums_column]){
          max_rowname <-  colnames(phenotype_data)[i + last_mean_sums_column]
          max_colname <- as.character(phenotype_data[,1][match(max.ind,phenotype_data[,i + last_mean_sums_column])])
          #i + last_mean_sums_column
          phenotype.cols[max_rowname,max_colname] <- "green"

        } else if (min.ind %in% phenotype_data[,i + last_mean_sums_column]){

          min_rowname <-  colnames(phenotype_data)[i + last_mean_sums_column]
          min_colname <- as.character(phenotype_data[,1][match(min.ind,phenotype_data[,i + last_mean_sums_column])])
          #i + last_mean_sums_column
          phenotype.cols[min_rowname,min_colname] <- "red"

        } else {
          next
        }

      }



      break

    }

    ###

    maxes <- maxes[first_max_sums_vector:last_max_sums_vector]

    mins <- mins[first_max_sums_vector:last_max_sums_vector]

    while (T) {


      for (i in 1:length(colnames(phenotype_data[,first_mean_sums_column:last_mean_sums_column]))) {

        #i = 1

        #feat.min<-mins[i]
        feat.max<-maxes[i]

       # feat.max <- max.sum
        feat.min <- min.sum
        feat.range<-(feat.max - feat.min)/2
        feat.mid<- feat.max-feat.range
        #feat.mid <- max.sum/2
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


        phenotype.cols[i,] <- feat.col$map(phenotype_data[,
                                                 i+grouping_mean_column])

      }

      print(phenotype.cols)
      phenotype.cols <<- phenotype.cols
      break
      #colnames(phenotype.cols)

    }


    ###


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
    if(missing(group_label)) {
      group_label2 <- "\nExperimental Conditions"
    } else {
      group_label2 <-   group_label
    }
    ###

    long.phenotype.cols <- gather(phenotype.cols2,
                                  key = agebin,
                                  value = color,
                                  -index )

    long.phenotype.cols$agebin <- factor(long.phenotype.cols$agebin,levels = phenotype_data[,grouping_mean_column])


  }
  else if (translation == 'absolute') {

    ###
    ###

    maxes <- sapply(phenotype_data[,2:ncol(phenotype_data)],max)

    mins <- sapply(phenotype_data[,2:ncol(phenotype_data)],min)


    max.ind <- abs(max(maxes[first_index_column:length(maxes)]))
    min.ind <- abs(min(mins[first_index_column:length(mins)]))
    med.ind <- 0

    while (T) {



    if(max.ind > min.ind){

      feat.max <- max.ind
      feat.min  <- -1*max.ind
      feat.mid <- med.ind

    }

    else if(max.ind < min.ind){

      feat.max <- min.ind
      feat.min  <- -1*min.ind
      feat.mid <- med.ind

    }

      break

    }

    ###

    last_mean_sums_column <- first_index_column

    first_index_column <- first_index_column + 1

    grouping_mean_column <- 1

    first_mean_sums_column  <- 2

    first_max_sums_vector <- first_mean_sums_column - 1

    last_max_sums_vector <-  last_mean_sums_column - 1

    ###


    max.sum <- max(sapply(phenotype_data[,first_mean_sums_column:last_mean_sums_column],max))
    min.sum  <- 0

    ###

    #indices_colors <- c("red", "yellow", "green")
    #sums_colors <- c( "white","grey","black")

    ###

    ##
    if(missing(indices_colors)) {
      indices_colors <- c("red", "yellow", "green")
    } else {
      indices_colors <- indices_colors
    }
    ##
    ##
    if(missing(sums_colors)) {
      sums_colors <- c( "white","grey","black")
    } else {
      sums_colors <- sums_colors
    }
    ##
    ###

    while (T) {

      library(tidyverse)

      nfeats <- length(maxes)

      feat.cols <- list()
      phenotype.cols <- data.frame(matrix(ncol = length(phenotype_data[,grouping_mean_column]),
                                          nrow = ncol(phenotype_data[,2:ncol(phenotype_data)])))

      rownames(phenotype.cols) <- colnames(phenotype_data[,
                                                 2:ncol(phenotype_data)])

      colnames(phenotype.cols) <- phenotype_data[,grouping_mean_column]

      for (i in 1:length(colnames(phenotype_data[,
                                        first_index_column:ncol(phenotype_data)]))) {

        #i = 1

        ##feat.range<-(feat.max - feat.min)/2
        ##feat.mid<- feat.max-feat.range


        #feat.min<-mins[last_max_sums_vector+i]
        #feat.max<-maxes[last_max_sums_vector+i]
        #feat.mid <-  max.ind - ((max.ind - min.ind)/2)

        #feat.max <- 1
        #feat.min  <- -1
        #feat.mid <- mid.ind

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
          limits=c(feat.min,feat.max)
          #limits=c(min.ind,max.ind)
        )

        phenotype.cols[last_max_sums_vector+i,] <- feat.col$map(
          phenotype_data[,i +
                  last_mean_sums_column])


      }



      break

    }

    ###

    maxes <- maxes[first_max_sums_vector:last_max_sums_vector]

    mins <- mins[first_max_sums_vector:last_max_sums_vector]

    while (T) {


      for (i in 1:length(colnames(phenotype_data[,first_mean_sums_column:last_mean_sums_column]))) {

        #i = 1

        #feat.min<-mins[i]
        feat.max<-maxes[i]

        # feat.max <- max.sum
        feat.min <- min.sum
        feat.range<-(feat.max - feat.min)/2
        feat.mid<- feat.max-feat.range
        #feat.mid <- max.sum/2
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


        phenotype.cols[i,] <- feat.col$map(phenotype_data[,
                                                 i+grouping_mean_column])

      }

      print(phenotype.cols)
      phenotype.cols <<- phenotype.cols
      break
      #colnames(phenotype.cols)

    }


    ###


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
    if(missing(group_label)) {
      group_label2 <- "\nExperimental Conditions"
    } else {
      group_label2 <-   group_label
    }
    ###

    long.phenotype.cols <- gather(phenotype.cols2,
                                  key = agebin,
                                  value = color,
                                  -index )

    long.phenotype.cols$agebin <- factor(long.phenotype.cols$agebin,levels = phenotype_data[,grouping_mean_column])


  }
  else if (translation == 'scaled'){

    ###
    ###

    maxes <- sapply(phenotype_data[,2:ncol(phenotype_data)],max)

    mins <- sapply(phenotype_data[,2:ncol(phenotype_data)],min)


    feat.max <- 1
    feat.mid <- 0
    feat.min <- -1

    ###

    last_mean_sums_column <- first_index_column

    first_index_column <- first_index_column + 1

    grouping_mean_column <- 1

    first_mean_sums_column  <- 2

    first_max_sums_vector <- first_mean_sums_column - 1

    last_max_sums_vector <-  last_mean_sums_column - 1

    ###


    max.sum <- max(sapply(phenotype_data[,first_mean_sums_column:last_mean_sums_column],max))
    min.sum  <- 0

    ###

    #indices_colors <- c("red", "yellow", "green")
    #sums_colors <- c( "white","grey","black")

    ###

    ##
    if(missing(indices_colors)) {
      indices_colors <- c("red", "yellow", "green")
    } else {
      indices_colors <- indices_colors
    }
    ##
    ##
    if(missing(sums_colors)) {
      sums_colors <- c( "white","grey","black")
    } else {
      sums_colors <- sums_colors
    }
    ##
    ###

    while (T) {

      library(tidyverse)

      nfeats <- length(maxes)

      feat.cols <- list()
      phenotype.cols <- data.frame(matrix(ncol = length(phenotype_data[,grouping_mean_column]),
                                          nrow = ncol(phenotype_data[,2:ncol(phenotype_data)])))

      rownames(phenotype.cols) <- colnames(phenotype_data[,
                                                 2:ncol(phenotype_data)])

      colnames(phenotype.cols) <- phenotype_data[,grouping_mean_column]

      for (i in 1:length(colnames(phenotype_data[,
                                        first_index_column:ncol(phenotype_data)]))) {



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
          limits=c(feat.min,feat.max)
          #limits=c(min.ind,max.ind)
        )

        phenotype.cols[last_max_sums_vector+i,] <- feat.col$map(
          phenotype_data[,i +
                  last_mean_sums_column])

      }



      break

    }

    ###

    maxes <- maxes[first_max_sums_vector:last_max_sums_vector]

    mins <- mins[first_max_sums_vector:last_max_sums_vector]

    while (T) {


      for (i in 1:length(colnames(phenotype_data[,first_mean_sums_column:last_mean_sums_column]))) {

        #i = 1

        #feat.min<-mins[i]
        feat.max<-maxes[i]

        # feat.max <- max.sum
        feat.min <- min.sum
        feat.range<-(feat.max - feat.min)/2
        feat.mid<- feat.max-feat.range
        #feat.mid <- max.sum/2
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


        phenotype.cols[i,] <- feat.col$map(phenotype_data[,
                                                 i+grouping_mean_column])

      }

      print(phenotype.cols)
      phenotype.cols <<- phenotype.cols
      break
      #colnames(phenotype.cols)

    }


    ###


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
    if(missing(group_label)) {
      group_label2 <- "\nExperimental Conditions"
    } else {
      group_label2 <-   group_label
    }
    ###

    long.phenotype.cols <- gather(phenotype.cols2,
                                  key = agebin,
                                  value = color,
                                  -index )

    long.phenotype.cols$agebin <- factor(long.phenotype.cols$agebin,levels = phenotype_data[,grouping_mean_column])


  }

  plas.phen <<-  ggplot(long.phenotype.cols,
                        aes(x = agebin ,
                            y = index)) +
    geom_tile(data = long.phenotype.cols,
              width=0.95,
              height=0.95,
              aes(fill = color))+
    scale_fill_identity()+
    theme(panel.background = element_rect(fill = 'white'),
          axis.text.x = element_text(#angle=45,
                                     vjust=0.5,
                                     size=11),
          axis.text.y = element_text(angle=0,
                                     vjust=0.5,
                                     size=11),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size=14, face = 'bold'),
          #aspect.ratio = 4/6,
          axis.title.y = element_text(size=14, face = 'bold')#,
          #plot.margin = margin(2, 2, 2, 2, "cm")
    )+
    scale_x_discrete(expand=c(0,0),
                     name=group_label2,
                     # labels=factor(phenotype_data[,grouping_mean_column])[group_select],
                     position="bottom")+
    scale_y_discrete(expand=c(0,0),
                     name="Plasticity Features\n",
                     limits= rev(colnames(phenotype_data)[#2:ncol(phenotype_data)
                       feature_order+1]),
                     position="left")+
    coord_fixed(ratio=0.25)

  return(plas.phen)


}

