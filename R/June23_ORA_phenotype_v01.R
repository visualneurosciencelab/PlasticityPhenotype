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

  
  library(tidyverse)
  
  if(missing(percentiles)){
    percentiles = c(0.25,0.75)
  } else {
    percentiles <- percentiles
  }
  
  low.per <- paste0(as.character(percentiles[1]*100),"%")
  high.per <- paste0(as.character(percentiles[2]*100),"%")
  
  feature_order = 1:ncol(features_df_row)
  
  while (T) {
    
    colnames(features_df_row)
    mean_list <- list()
    sd_list <- list()
    length_list <- list()
    counter <- 0
    
    for (i in 1:length(condition_list)) {
      
      for (j in 1:ncol(features_df_row)) {
        
        counter <- counter+1
       
        SubsetCondition <- as.data.frame(subset(features_df_row[,j],
                                                grepl(condition_list[i],
                                                      rownames(features_df_row),
                                                      fixed = T)))
        
        mean.Condition <- as.data.frame(mean(SubsetCondition[,1], na.rm=TRUE))
        colnames(mean.Condition) <- paste0("Average -- ",colnames(features_df_row)[j]," (",condition_list[i],")")
        mean_list[[counter]] <- as.data.frame(mean.Condition)
        
        
        length.Condition <-  as.data.frame(length(SubsetCondition[,1]))
        colnames(length.Condition) <- paste0("Length -- ",colnames(features_df_row)[j]," (",condition_list[i],")")
        length_list[[counter]] <- as.data.frame(length.Condition)
        
        sd.Condition <-  as.data.frame(sd(SubsetCondition[,1]))
        colnames(sd.Condition) <- paste0("Std.Dev -- ",colnames(features_df_row)[j]," (",condition_list[i],")")
        sd_list[[counter]] <- as.data.frame(sd.Condition)
        
      }
      
    }
    
    
    break
    
  }
  
  
  cols_df <- data.frame(matrix(ncol = ncol(features_df_row)+1,
                               nrow = length(condition_list)))
  
  colnames(cols_df) <- c("Conditions",colnames(features_df_row))
  
  while (T) {
    
    counter <- 0

    for (i in 1:(length(condition_list)
    )) {
      
      cols_df[i
              ,1] = condition_list[i
                                   ]
   
      for (j in 1:ncol(features_df_row)) {
        
        counter <- counter + 1
        mean.Normal <- mean_list[[j]][[1]]
        mean.MD1 = mean_list[[counter]][[1]]
        adj.pval = 0.05/length_list[[1]][[1]]
        
         SubsetCondition <- as.data.frame(subset(features_df_row[,j],
                                                 grepl(condition_list[i],
                                                       rownames(features_df_row),
                                                       fixed = T)))
        
        
        sim.MD1 <- rnorm(1000000 #simulate population of MD that has ONE MILLION =N
                         , mean = mean_list[[j]][[1]] #using the mean of MD
                         , sd =  sd_list[[j]][[1]] #and the stdev of MD
        )
        
        
        
        
        
        resamples.NormalMD1 <- lapply(1:100000 #pull from a data set ONE HUNDRED THOUSAND times
                                      ,function(i) sample(sim.MD1           #the data set to pull from is the simulated MD population
                                                          , length_list[[1]][[1]]       #pull as many times as there are samples in NORMAL
                                                          , replace = T)   #replace the samples after you pull them so each scenario has same P of happening
        )
        
        
        mean.resamples.NormalMD1 <- sapply(resamples.NormalMD1, mean)
        
        CI.NormalMD1 <- quantile(mean.resamples.NormalMD1, c(0.05,0.95))
        
        quart.comp <- quantile(SubsetCondition[,1],percentiles)
        
        
        cols_df[i
                ,1+
                  j] <- if(CI.NormalMD1["5%"] > quart.comp[high.per]){
                    #print( "purple")
                    "purple"
                  } else if (CI.NormalMD1["95%"] < quart.comp[low.per]) {
                    #print("yellow")
                    "yellow"
                  } else {
                    #print('grey')
                    "#dce0e5"
                  }
        
        
        
        
      }
      
    }
    break
  }
  
  cols_df$Conditions <-  factor(cols_df$Conditions,
                                levels = unlist(condition_list),
                                ordered = T)
  
  long.phenotype.cols <- gather(cols_df, feature, color, - Conditions )
  

  
  
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
