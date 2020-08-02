#' Performs an overrepresentation analysis (ORA) across all groups and within each feature, then visualize results in a phenotype.
#'
#' ORA between groups and visualizes the significantly different features in a plasticity phenotype.
#' @param features_df_row A dataframe that contains all features to be explored using an ORA. Specify row names for each observation must contain the experimental group.
#' @param condition_list A list of factors indicating the group identifiers as they appear in the data frame (row names)
#' @param group_label X-axis title.
#' @param reference_group A character string of the reference group as it appears in the row name. This indicates the group against which an ORA is performed.
#' @param percentiles A list of two values (in ascending order) that act as thresholds for the ORA. These values range between [0,1], with default values of c(0.25,0.75). The first number acts as the lower percentile for the experimental group, while the second acts as the upper percentile. 
#' @return Plots an ORA phenotype for each group, where coloured boxes indicate that a feature is outside of the 90% confidence interval --  default = c(0.05,0.95) -- of the reference group. The color "yellow" reflects when an experimental group is overrepresented (its lower percentile is larger than the reference groups' upper percentile ),  while "purple" is underrepresented  (its upper percentile is less than the reference groups' lower percentile ). Stores the object, 'ora.phen' in your global environment. The 'ora.phen' object contains ggplot visualization of each comparison.
#' @examples
#' ORA_phenotype(features_df_row  = tsne.processed.4[,2:10], # Data frame for bootstrap analysis
#'              condition_list = as.list(c('Normal 1\n C,P,M',
#'                                          'LT BV 1\n C,P,M',
#'                                          'MD 1\n P,M',
#'                                          'LT BV 5\n P,M',
#'                                          'LT BV 4\n P,M',
#'                                          'RO 2\n C,P,M',
#'                                          'ST BV 3\n C,P,M',
#'                                          'MD 3\n C,P',
#'                                          'ST BV 5\n P',
#'                                          'ST BV 1\n C,P,M',
#'                                          'LT BV 6\n P',
#'                                          'BD 3\n C,P,M',
#'                                          'BD 6\n P')), # List of subclusters as they appear in row names of "features_df_row"
#'               reference_group =  'Normal 1\n C,P,M', # Name of reference group as it appears in the row names of "features_df_row"
#'               group_label  = "\nSubclusters",
#'               percentiles = c(0.1,0.9)
#' )
#' @export
#'

ORA_phenotype <- function(features_df_row = features_df_row,
                          condition_list = condition_list,
                          group_label = group_label,
                          reference_group = reference_group,
                          percentiles = percentiles){
  
  
  
  # Load the tidyverse.
  
  library(tidyverse)
  
  # If the 'percentiles' parameter is empty, assign it the default values of '0.25' & '0.75'.
  
  if(missing(percentiles)){
    percentiles = c(0.25,0.75)
  } else {
    percentiles <- percentiles
  }
  
  # Convert the decimals in 'percentiles' into percentage characters, and assign to 
  #       'low.per' and 'high.per' for the low and high values, respectively.
  
  low.per  <- paste0(as.character(percentiles[1]*100),"%")
  high.per <- paste0(as.character(percentiles[2]*100),"%")
  
  # Assign a vector of 1 to the number of features to 'feature_order'. 
  
  feature_order = 1:ncol(features_df_row)
  
  ### This multi-nested while-loop calculates the mean, SD, and number of data points
  ###       for each experimental condition in 'condition_list' across each plasticity feature
  ###       in 'features_df_row'. Store each mean, SD and length that is calculated in 
  ###       3 different lists: mean_list, sd_list, and length_list. 
  
  # While True...
  
  while (T) {
    
    # Initialize 3 empty lists:
    #
    #       1.) 'mean_list' -- A list of the average value of every feature across all experimental groups
    #       2.) 'sd_list' -- A list of the standard deviations of every feature across all experimental groups
    #       3.) 'length_list' -- A list of the number of data points for every feature across all experimental groups
    
    mean_list <- list()
    sd_list <- list()
    length_list <- list()
    
    # Initialize a counter variable to the value of '0'.
    
    counter <- 0
    
    # For i in 1 to the number of experimental conditions...
    
    for (i in 1:length(condition_list)) {
    
      # For j in 1 to the number of features
        
      for (j in 1:ncol(features_df_row)) {
        
        # Increment the counter variable by '1' every time 'j' updates.
        
        counter <- counter + 1
        
        # For the currently indexed feature in 'features_df_row', retain only the 
        #       data points for the currently-indexed experimental condition in 'condition_list'.
        #       Store all such data points for the current experimental group at the current
        #       feature in an object called 'SubsetCondition'. 
        
        SubsetCondition <- as.data.frame(subset(features_df_row[,j],
                                                grepl(condition_list[i],
                                                      rownames(features_df_row),
                                                      fixed = T)))
        
        # Convert 'SubsetCondition' into a data frame and store it in the object 'mean.Condition'.
        #
        # Assign 'mean.Condition' a column name that denotes what calculation this number reflects --
        #       in the case of 'mean.Condition', it is the an average value -- and also denotes what experimental group
        #       and feature this value summarizes. 
        #
        # Convert 'mean.Condition' to a data frame and assign it as an instance of the 'mean_list' 
        #       list. It's indexed according to whatever the current value of 'counter' is. 
        
        mean.Condition <- as.data.frame(mean(SubsetCondition[,1], na.rm=TRUE))
        colnames(mean.Condition) <- paste0("Average -- ",colnames(features_df_row)[j]," (",condition_list[i],")")
        mean_list[[counter]] <- as.data.frame(mean.Condition)
        
        # Convert 'SubsetCondition' into a data frame and store it in the object 'length.Condition'.
        #
        # Assign 'length.Condition' a column name that denotes what calculation this number reflects --
        #       in the case of 'length.Condition', it is the number of data points for an experimental group at a given feature -- 
        #       and also denotes what experimental group and feature this value summarizes. 
        #
        # Convert 'length.Condition' to a data frame and assign it as an instance of the 'length_list' 
        #       list. It's indexed according to whatever the current value of 'counter' is. 
        
        
        length.Condition <-  as.data.frame(length(SubsetCondition[,1]))
        colnames(length.Condition) <- paste0("Length -- ",colnames(features_df_row)[j]," (",condition_list[i],")")
        length_list[[counter]] <- as.data.frame(length.Condition)
        
        # Convert 'SubsetCondition' into a data frame and store it in the object 'sd.Condition'.
        #
        # Assign 'sd.Condition' a column name that denotes what calculation this number reflects --
        #       in the case of 'sd.Condition', it is the standard deviation -- and also denotes what experimental group
        #       and feature this value summarizes. 
        #
        # Convert 'sd.Condition' to a data frame and assign it as an instance of the 'sd_list' 
        #       list. It's indexed according to whatever the current value of 'counter' is. 
        
        sd.Condition <-  as.data.frame(sd(SubsetCondition[,1]))
        colnames(sd.Condition) <- paste0("Std.Dev -- ",colnames(features_df_row)[j]," (",condition_list[i],")")
        sd_list[[counter]] <- as.data.frame(sd.Condition)
        
      }
      
    }
    
    # Break out of this while-loop
    
    break
    
  }
  
  ### Initialize an empty data frame with the same number of rows as there are items in the 
  ###       'condition_list', and a number columns equal to 1 plus the number of columns in
  ###       the 'features_df_row'. Call this empty data frame 'cols_df'.
  
  cols_df <- data.frame(matrix(ncol = ncol(features_df_row)+1,
                               nrow = length(condition_list)))
  
  ### Assign proper column headers to 'cols_df'. The first column will be 'Conditions',
  ###       while the remaining headers will be the same headers as 'features_df_row'.
  
  colnames(cols_df) <- c("Conditions",
                         colnames(features_df_row))
  
  ### This while-loop calculates the quantiles specified in 'percentiles' for the 
  ###       currently indexed experimental group, and the (5%, 95%) quantiles for 
  ###       the reference/normal group. This is done for the same feature across
  ###       the normals and experimental groups. It then assesses whether the experimental
  ###       group is over-represented/under-represented/similar to the 
  ###       quantiles of the reference group. 
  
  # While true...
  
  while (T) {
    
    # Initialize an empty 'counter' variable to '0'. 
    
    counter <- 0
  
    # For i in 1 to the number of experimental groups...
      
    for (i in 1:(length(condition_list))) {
      
      # Assign the currently-indexed row of the 'Conditions' column
      #       of 'cols_df' the name of the current experimental group
      #       being analyzed. 
      
      cols_df[i,1] <- condition_list[i]
      
      # For j in 1 to the number of features...
      
      for (j in 1:ncol(features_df_row)) {
        
        # Increment 'counter' by 1 every time 'j' updates. 
        
        counter <- counter + 1        
        
        # For the currently indexed feature in 'features_df_row', retain only the 
        #       data points for the currently-indexed experimental condition in 'condition_list'.
        #       Store all such data points for the current experimental group at the current
        #       feature in an object called 'SubsetCondition'. 
        
        SubsetCondition <- as.data.frame(subset(features_df_row[,j],
                                                grepl(condition_list[i],
                                                      rownames(features_df_row),
                                                      fixed = T)))
        
        # Simulate a normal distribution with N = 1 million, and ensure that the 
        #       the mean and SD of this distribution are of the same feature in the 
        #       normal condition. Assign these 1 million data points
        #       to the object 'simulated.Norms'. 
        
        simulated.Norms <- rnorm(1000000 #simulate population of MD that has ONE MILLION =N
                         , mean = mean_list[[j]][[1]] #using the mean of MD
                         , sd =  sd_list[[j]][[1]] #and the stdev of MD
        )
        
        
        # Resample with replacement from 'simulated.Norms' a total of 100'000 times. Each time, resample as many
        #       data points as there are exist for the reference/normal group. Store this list 
        #       of 100'000 instances of x many data points in 'resampled.Norms'.
        
        
        resampled.Norms <- lapply(1:100000 #pull from a data set ONE HUNDRED THOUSAND times
                                      ,function(i) sample(simulated.Norms           #the data set to pull from is the simulated MD population
                                                          , length_list[[1]][[1]]       #pull as many times as there are samples in NORMAL
                                                          , replace = T)   #replace the samples after you pull them so each scenario has same P of happening
        )
        
        # Calculate the average of each of the 100'000 groups of resampled data points, and store in 'mean.resampled.Norms'. 
        
        mean.resampled.Norms <- sapply(resampled.Norms, mean)
        
        # Calculate the 5% and 95% quantiles on the data points stored in 'mean.resampled.Norms'. Assign these points to 'quart.Norms'. 
        
        quart.Norms <- quantile(mean.resampled.Norms, c(0.05,0.95))
        
        # Calculate the specified 'percentiles' of the data points stored in 'SubsetCondition'. 
        
        quart.Comps <- quantile(SubsetCondition[,1],percentiles)
        
        # If the 5% quantile of 'quart.Norms' is greater than the highest percentile of 'quart.Comps',
        #       then assign the relevant cell in 'cols_df' the colour "purple" --> 'under-represented'.
        #
        # Else, if the 95% quantile of 'quart.Norms' is less than the lowest percentile of 'quart.Comps',
        #       then assign the relevant cell in 'cols_df' the colour "yellow" --> 'over-represented'.
        #
        # Else, in all other cases,assign the relevant cell in 'cols_df' the colour "grey" 
        
        cols_df[i,1+j] <- if(quart.Norms["5%"] > quart.Comps[high.per]){
                                    #print( "purple")
                                    "purple"
                                  } else if (quart.Norms["95%"] < quart.Comps[low.per]) {
                                    #print("yellow")
                                    "yellow"
                                  } else {
                                    #print('grey')
                                    "#dce0e5"
                                  }
                        
        
        
        
      }
      
    }
    
    # Break out of the while-loop. 
    
    break
  }
  
  # Order the levels of 'cols_df$Conditions' based on the predefined 
  #       order of the 'conditions_list'. This is the order that
  #       conditions will appear on the x-axis of the phenotype.
  
  cols_df$Conditions <-  factor(cols_df$Conditions,
                                levels = unlist(condition_list),
                                ordered = T)
  
  # Convert the 'cols_df' object from wide to long format
  
  long.phenotype.cols <- gather(cols_df, feature, color, - Conditions )
  
  # Create the ORA phenotype, and assign it to an object, 'ora.phen'.
  
  ora.phen <<- ggplot(long.phenotype.cols,
                      aes(x =  Conditions,
                          y = feature)) +
    geom_tile(data = long.phenotype.cols,
              width=0.95,
              height=0.95,
              color = "black",
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
          axis.title.y = element_text(size=14,  face = 'bold')
    )+
    scale_x_discrete(expand=c(0,0),
                     name=group_label,
                     
                     position="bottom")+
    scale_y_discrete(expand=c(0,0),
                     name="Plasticity Features\n",
                     limits= rev(colnames(cols_df)[
                       feature_order+1]),
                     position="left")+
    coord_fixed(ratio=0.25)
  
  print(ora.phen)
  
}
