#' Create plasticity phenotypes
#'
#' Plots a series of stacked, colour-coded horizontal bars for each group. Color of each bar represents the group mean for the designated feature. Produces a phenotype for each group. Note: For all arguments matrix column numbers are identified as N-1 (e.g. column 4 is indexed as 3).
#' @param df_list A list of data frames with the rows containing group labels and the remaining columns containing the average/median values of each significant plasticity feature across all groups; this function requires all sums to be grouped together and all indices to be grouped together. Note: Features that describe variable sums must precede features that describe variable indices (e.g. sums = columns 1-4, indices = columns 5 – 10).
#' @param first_index_column Numerical value indicating the column number that marks the first feature index (plotted as green-yellow-red); all feature columns before this should be feature sums (plotted as grey-black).
#' @param feature_order Numeric vector that describes the order of features to be plotted. The order of values in the vector represent the column numbers for each feature (c(1,3,4,5,2,6,7)). Default order is taken from the data matrix.
#' @param group_label A list containing the x-axis title for each corresponding data frame in "df_list". Default is "Experimental Conditions".
#' @param indices_colors Character vector containing the color codes for indices color palette. Order is low, middle, high. Default is c("red", "yellow", "green").
#' @param sums_colors Character vector containing the color codes for sums color palette. Order is low, middle, high. Default is c( "#d3d3d3","grey","black").
#' @param translation Assign min/max data values for colour scales (feature indices only). Four options include: 'local', 'absolute'.
#' local = local min/max determined within each feature. All features will have an absolute min (default=red) and max (default=green).
#' absolute= global absolute min/max determined across ALL features. There will only be one occurence of the global min (default=red) OR the global max (default=green) across all features. Min/max determined by greatest ABSOLUTE value.
#' Default= 'local'.
#' @return Plots phenotype for each group. Stores object plasticity.cols in your Global Environment. This object will be consulted using the phenotype_boxplots function for each feature across each group. Stores object plas.phen in your Global Environment. plas.phen contains ggplot visualization of each comparison.
#' @examples
#'plasticity_phenotype(
#'     df_list = df_list,
#'     first_mean_index_column =  6,
#'     group_label = c("\nRearing Conditions"),
#'     translation = 'absolute',
#'     indices_colors = c('red','yellow','green'),
#'     sums_colors = c('white','grey','black'),
#'     feature_order = 1:9
#'     )
#' @export

plasticity_phenotype <- function(df_list = df_list,
                                 first_index_column = first_index_column,      
                                 group_label = group_label,
                                 indices_colors = indices_colors,
                                 sums_colors = sums_colors,
                                 feature_order = feature_order,
                                 translation = 'local'){
  
  # Load tidyverse
  
  library(tidyverse)

  # Initialize 2 empty lists:
  #
  #       1.) plas.phen -- will store the actual plasticity phenotypes
  #       2.) phen.cols -- will store the colour-code of the phenotypes
  
  plas.phen <<- list()
  phen.cols <<- list()
  
  # If the "group_label" parameter isn't explicitly stated,
  #       then its default value is "\nExperimental Conditions".
  
  if(missing(group_label)) {
    group_label2 <- rep("\nExperimental Conditions",length(df_list))
  } else {
    group_label2 <-   group_label
  }
  
  # Store the column names of each index feature from 'df_list[[1]]' in the 
  #       object, 'ind.names'.
  
  ind.names <- colnames(df_list[[1]])[first_index_column:ncol(df_list[[1]])]
  
  # From every phenotype date frame in 'df_list', retain only the columns that
  #       contain data for the index features. Store these subset data frames
  #       in the object, 'new_list'. Note, the order of the list
  #       is the same as it was in 'df_list'. 
  
  new_list <- lapply(df_list, function(x) x%>% select(ind.names))
  
  # Convert every column across each item in a list from a positive or negative to an absolute value.
  #       Store this new absolute value list in 'new_list2'. 
  
  new_list2 <- lapply(new_list, abs)
  
  # Calculate the largest absolute value in each item of the 'new_list2' object. Store each
  #       maximum value in the same order, but a new list called 'new_list3'.
  
  new_list3 <- lapply(new_list2,max)
  
  # Calculate the largest absolute value from the 'new_list3' object and store it
  #       in 'new_list4'. 
  
  new_list4 <- max(unlist(new_list3))
  
  # Assign the value of 'first_index_column' to 'last_mean_sums_column'.
  
  last_mean_sums_column <- first_index_column
  
  # Add '1' to the current value of 'first_index_column' and reassign it to itself.
  
  first_index_column <- first_index_column + 1
  
  # Assign the object, 'grouping_mean_column' a value of '1'.
  
  grouping_mean_column <- 1
  
  # Assign 'first_mean_sums_column' a value of '2'.
  
  first_mean_sums_column  <- 2
  
  # Assign 'first_max_sums_vector ' a value of 'first_mean_sums_column - 1'.
  
  first_max_sums_vector <- first_mean_sums_column - 1
  
  # Assign 'last_max_sums_vector ' a value of 'last_mean_sums_column - 1'.
  
  last_max_sums_vector <-  last_mean_sums_column - 1
  
  
  for (a in 1:length(df_list)){
    
    # Assign the currently-indexed entry in 'df_list[[a]] to the 
    #       object, 'phenotype_data'.
   
    phenotype_data <- df_list[[a]]
    
    # Store the column headers of 'phenotype_data' in the object,
    #       'phen.colheads'.
    
    phen.colheads <- colnames(phenotype_data)
    
    # Store the row names of 'phenotype_data' in a  new column called
    #       'phenotype_data$Group.1'.
    
    phenotype_data$Group.1 <-  rownames(phenotype_data) 
    
    # Reorder the 'phenotype_data' object to have the group names
    #       and the feature headers in that order. Reassign that data 
    #       to the object, 'phenotype_data'. 
    
    phenotype_data <- phenotype_data[,c('Group.1',phen.colheads)]
    
    # if the 'translation' parameter is set to 'local'...
    
    if (translation == "local"){
      
      # Calculate the maximum values of across every feature column 
      #       in 'phenotype_data'. Assign this vector to the object,
      #       'maxes'.
      
      maxes <- sapply(phenotype_data[,2:ncol(phenotype_data)],max)
      
      maxes
      
      # Calculate the minimum values of across every feature column 
      #       in 'phenotype_data'. Assign this vector to the object,
      #       'mins'.
      
      mins <- sapply(phenotype_data[,2:ncol(phenotype_data)],min)
      
      mins
      
      # Calculate the maximum values of every feature sum in 'phenotype_data' and
      #       store these values in the object, 'max.sum'. Assign a value of 0 to 
      #       the object 'min.sum'. 
      
      max.sum <- max(sapply(phenotype_data[,first_mean_sums_column:last_mean_sums_column],max))
      min.sum  <- 0

      # If the 'indices_colors' parameter is empty, then assign a default 
      #       colour-palette. 
      
      if(missing(indices_colors)) {
        indices_colors <- c("red", "yellow", "green")
      } else {
        indices_colors <- indices_colors
      }

      # If the 'sums_colors' parameter is empty, then assign a default 
      #       colour-palette. 
      
      if(missing(sums_colors)) {
        sums_colors <- c( "white","grey","black")
      } else {
        sums_colors <- sums_colors
      }
      
      # While this loop is 'true'...
      
      while (T) {
        
        # Assign the number corresponding to the number of entries in 'maxes', to the object 'nfeats'.
        
        nfeats <- length(maxes)
        
        # Initialize the empty list, 'feat.cols'. 
        
        feat.cols <- list()
        
        # Initialize an empty data frame called 'phenotype.cols'. This data frame has a number of columns equal to the
        #       the number of experimental groups, while the number of rows correspond to the number of plasticity
        #       features. 
        
        phenotype.cols <- data.frame(matrix(ncol = length(phenotype_data[,grouping_mean_column]),
                                            nrow = ncol(phenotype_data[,2:ncol(phenotype_data)])))
        
        # Assign the names of the plasticity features as the row names of the 'phenotype.cols' object.
        
        rownames(phenotype.cols) <- colnames(phenotype_data[,2:ncol(phenotype_data)])
        
        # Assign the names of the experimental groups as the column headers of the 'phenotype.cols' object.
        
        colnames(phenotype.cols) <- phenotype_data[,grouping_mean_column]
        
        # For i in 1 to the number of feature indices (and not sums)...
        
        for (i in 1:length(colnames(phenotype_data[,first_index_column:ncol(phenotype_data)]))) {
          
          # For the current value of i, add it to the value of 'last_max_sums_vector'. This will provide an number
          #       that can be used to select one of the feature indices (and skip over any sums). Select the maximum and
          #       and minimum values of the feature indices from the 'mins' and 'maxes' vectors, and store them in the 
          #       object 'feat.min' and 'feat.max' respectively.
          
          feat.min <- mins[last_max_sums_vector+i]
          
          feat.max <- maxes[last_max_sums_vector+i]
          
          # Calculate the range of the previously-indexed max and min values, and store in 'feat.range'.
          
          feat.range <- (feat.max - feat.min)/2
          
          # Subtract the value of 'feat.range' from 'feat.max' to yield a midpoint value. Store this
          #       midpoint value in the object 'feat.mid'.
          
          feat.mid <- feat.max-feat.range
          
          # Establish the colour-scale for the each feature index on a feature-by-feature basis. Set the colour for the lowest,
          #       median, and highest values to be the first, second and third entries in 'indices_colors'. As the values of 
          #       feat.min, feat.max, and feat.mid change on a feature-by-feature basis, the colour-scale is set such that ecery feature
          #       has its own max and min set to the 'green' and 'red', while the median value is 'yellow'. Store this infomation
          #       in the object, 'feat.col'.
          
          feat.col<- scale_colour_gradient2(low = indices_colors[1],
                                            mid = indices_colors[2],
                                            high = indices_colors[3],
                                            midpoint = feat.mid,
                                            breaks = c(feat.min,
                                                       feat.mid,
                                                       feat.max),
                                            labels = c("Below Normal",
                                                       "Normal",
                                                       "Above Normal"),
                                            limits = c(feat.min,feat.max))
          
          # This line maps the previously calculated colour-scale for te current feature to the actual values of the feature across
          #       every experimental group. Store this vector of colour-codes to the properly indexed row of 'phenotype.cols'.
          
          phenotype.cols[last_max_sums_vector+i,] <- feat.col$map(phenotype_data[,i + last_mean_sums_column])
          
        }
        
        # Break out of the while-loop. 
        
        break
        
      }
      
      # Reassign the maximum values of the feature sums (not indices) to the 'maxes' object.
      
      maxes <- maxes[first_max_sums_vector:last_max_sums_vector]
      
      # Assign the minimum values for the feature sumes to be '0'.
      
      mins <- 0
      
      # While true...
      
      while (T) {
        
        # For i in 1 to the number of feature sums (and not feature indices)...
        
        for (i in 1:length(colnames(phenotype_data[,first_mean_sums_column:last_mean_sums_column]))) {
          
          # Assign the currently-indexed feature sum's max value in the object, 'feat.max'.

          feat.max <- maxes[i]

          # Assign the minimum value for each sum (which should be 0) to the object, 'feat.min'.
          
          feat.min <- mins
          
          # Calculate half the difference of 'feat.max' and 'feat.min' and store in object, 'feat.range'.
          
          feat.range <- (feat.max - feat.min)/2
          
          # Subtract the value of 'feat.range' from 'feat.max', and store in object, 'feat.mid'.
          
          feat.mid <- feat.max - feat.range
          
          # Determine the colour-codes for each feature sum on a feature-by-feature basis. Minimum values are assigned the color 'white', mid values are assigned
          #       a 'grey', and maximum values are assigned a 'black'. The resulting colour-scale for the current feature is stored in the object 'feat.col'.
          
          feat.col<-scale_colour_gradient2(low  = sums_colors[1],
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
          
          # Convert the numeric values for a given sum across all experimental groups to hex-code and store the vector of hex-codes in 
          #     a correspondingly-indexed row of 'phenotyp.cols'.
          
          phenotype.cols[i,] <- feat.col$map(phenotype_data[,i+grouping_mean_column])
          
        }
        
        # Reassign 'phenotype.cols' to itself and break out of this while-loop.
        
        phenotype.cols <<- phenotype.cols
        
        break
        
      }
      
    }
    
    # Else, if the 'translation' parameter is set to 'absolute'...
    
    else if (translation == 'absolute') {
      
      # Store a vector of all maximum values across all features in the object, 'maxes'. 
      
      maxes <- sapply(phenotype_data[,2:ncol(phenotype_data)], max)
      
      # Assign a value of '0' to 'mins'.
      
      mins <- 0
      
      # Assign 'feat.mid' a value of '0'. 
      
      feat.mid <- 0
      
      # Assign the object, 'new_list4' (which contains the largest absolute value across all data frames in 
      #       'df_list'), to the object 'feat.max'.
      
      feat.max <- new_list4
      
      # Assign 'feat.min' the negative of 'new_list4'.
      
      feat.min  <- -1*new_list4
      
      # If the 'indices_colors' parameter is empty, then assign a default 
      #       colour-palette.
      
      if(missing(indices_colors)) {
        indices_colors <- c("red", "yellow", "green")
      } else {
        indices_colors <- indices_colors
      }
      
      # If the 'sums_colors' parameter is empty, then assign a default 
      #       colour-palette.
      
      if(missing(sums_colors)) {
        sums_colors <- c( "white","grey","black")
      } else {
        sums_colors <- sums_colors
      }
      
      # While true... 
      
      while (T) {
        
        # Store the total number of features (both sums and indices) in the object,
        #       'nfeats'.
        
        nfeats <- ncol(df_list[[a]])
        
        # Initialize an empty list called 'feat.cols'. 

        feat.cols <- list()
        
        # Initialize an empty data frame called 'phenotype.cols'. This data frame has a number of columns equal to the
        #       the number of experimental groups, while the number of rows correspond to the number of plasticity
        #       features.
        
        phenotype.cols <- data.frame(matrix(ncol = length(phenotype_data[,grouping_mean_column]),
                                            nrow = ncol(phenotype_data[,2:ncol(phenotype_data)])))
        
        # Assign the names of the plasticity features as the row names of the 'phenotype.cols' object.
        
        rownames(phenotype.cols) <- colnames(phenotype_data[,2:ncol(phenotype_data)])
        
        # Assign the names of the experimental groups as the column headers of the 'phenotype.cols' object.
        
        colnames(phenotype.cols) <- phenotype_data[,grouping_mean_column]
        
        # For i in 1 to the number of feature indices (and not feature sums)... 
        
        for (i in 1:length(colnames(phenotype_data[,first_index_column:ncol(phenotype_data)]))) {
          
          # Set the colour-scale for the feature indices only. the min, mid and max get 'red', 'yellow' and 'green', respectively.
          #       The value of 'feat.max' is the largest absolute values across all data frames in 'df_list', 'feat.min' is the 
          #       the value of 'feat.min' multiplied by -1, and 'feat.mid' is assigned a value of '0'. These limits are set for every single
          #       feature index across every data frame entered in 'df_list'.
          
          feat.col <- scale_colour_gradient2(low  = indices_colors[1],
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
          
          # This line is how the same, constant colour-scale is applied to each feature across each data frame in 'df_list'. For a given data frame
          #       that was assigned to 'phenotype_data', take all values of a feature index across all experimental conditions and convert them
          #       to hex-code values that are based on the previously established colour-scale. Assign these hex-values in their proper
          #       order to their properly indexed row in 'phenotype.cols'. 
          
          phenotype.cols[last_max_sums_vector+i,] <- feat.col$map(phenotype_data[,i +last_mean_sums_column])
          
        }
        
        # Break out of this while-loop. 
        
        break
        
      }
      
      # Reassign the maximum values of the feature sums (and not the feature indices) to the object, 'maxes'.
      
      maxes <- maxes[first_max_sums_vector:last_max_sums_vector]
      
      # Assign 'mins' a value of '0'. 
      
      mins <- 0
      
      # While true...
      
      while (T) {
        
        # For i in 1 to the number of feature sums (and not feature indices)... 
        
        for (i in 1:length(colnames(phenotype_data[,first_mean_sums_column:last_mean_sums_column]))) {
          
          # Assign the currently-indexed entry of 'maxes' to the object, 'feat.max'. 
          
          feat.max <- maxes[i]
          
          # Assign the value of 'mins' to the object, 'feat.min'.
          
          feat.min <- mins 
          
          # Subtract the value of 'feat.min' from 'feat.max', take half of that
          #       and store in the object, 'feat.range'.
          
          feat.range <- (feat.max - feat.min)/2
          
          # Subtract 'feat.range' from 'feat.max' and store in 'feat.mid'.
          
          feat.mid <- feat.max - feat.range
          
          # Determine the colour-codes for each feature sum on a feature-by-feature basis. Minimum values are assigned the color 'white', mid values are assigned
          #       a 'grey', and maximum values are assigned a 'black'. The resulting colour-scale for the current feature is stored in the object 'feat.col'. 
          
          feat.col <- scale_colour_gradient2(low  = sums_colors[1],
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
          
          # Convert the numeric values for a given sum across all experimental groups to hex-code and store the vector of hex-codes in 
          #     a correspondingly-indexed row of 'phenotyp.cols'.
          
          phenotype.cols[i,] <- feat.col$map(phenotype_data[,i+grouping_mean_column])
          
        }
        
        # Reassign 'phenotype.cols' to itself.
        
        phenotype.cols <<- phenotype.cols
        
        # Break out of this while-loop.
        
        break
        
      }
      
      
    }
    
    # Print the current value of 'phenotype.cols'. 
    
    print(phenotype.cols)
    
    # Assign the current instance of 'phenotype.cols' as an indexed entry into the 
    #       'phen.cols' list. 
    
    phen.cols[[a]] <<- phenotype.cols
    
    # Populate a new column in 'phenotype.cols' called 'index' with the row names of 
    #       'phenotype.cols'. 
  
    phenotype.cols$index <- row.names(phenotype.cols)
    
    # If the parameter 'feature_order' is empty, populate it with default order of the features. 
    
    if(missing(feature_order)) {
      feature_order <- 1:nfeats
    } else {
      feature_order <- feature_order
    }
    
    # Assign 'phenotype.cols' to 'phenotype.cols2' to ensure that the original 'phenotype.cols' is
    #       is not lost/overwritten. 
    
    phenotype.cols2 <- phenotype.cols
    
    # Convert 'phenotype.cols2' from wide to long format, and store in the object
    #       'long.phenotype.cols'.
    
    long.phenotype.cols <- gather(phenotype.cols2,
                                  key = agebin,
                                  value = color,
                                  -index )
    
    # Factor the 'agebin' column of 'long.phenotype.cols' order age-bins ascendingly. 
    
    long.phenotype.cols$agebin <- factor(long.phenotype.cols$agebin,
                                         levels = phenotype_data[,grouping_mean_column])
    
    # Create a phenotype and store it as an indexed entry of the list, 'plas.phen'.
    
    plas.phen[[a]] <<-  ggplot(long.phenotype.cols,
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
                       name=group_label2[a],
                       # labels=factor(phenotype_data[,grouping_mean_column])[group_select],
                       position="bottom")+
      scale_y_discrete(expand=c(0,0),
                       name="Plasticity Features\n",
                       limits= rev(colnames(phenotype_data)[#2:ncol(phenotype_data)
                         feature_order+1]),
                       position="left")+
      coord_fixed(ratio=0.25)
    
    
  }
  
  return(plas.phen)
  return(phen.cols)
  
}
