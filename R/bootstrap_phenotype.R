#' Perform a bootstrap analysis across groups and within each feature, then visualize results in a phenotype.
#'
#' Bootstrap comparison between groups and visualizes the significantly different features in a plasticity phenotype.
#' @param features_df_row A dataframe that contains all features to be analyszed with bootstrap comparison. Specify row names for each observation must contain the experimental group.
#' @param condition_list A list of factors indicating the group identifiers as they appear in the data frame (row names)
#' @param group_label X-axis title.
#' @param reference_group A character string of the reference group as it appears in the row name. This indicates the group against which all significance comparisons are made.
#' @return Plots significance phenotype for each group, where colored groups are significantly different than the reference group. Red is significant and greater than normal, blue is significant and less than the reference group. Stores object sig.phen in your Global Environment. sig.phen contains ggplot visualization of each comparison.
#' @examples
#'bootstrap_phenotype(
#'     features_df_row = NewFeatures[,ordered_nf_colnames],
#'     condition_list = c('  normal  ',
#'                        '  1wk MD  ',
#'                        '  fluoxetine + 1wk MD  ',
#'                        '  fluoxetine  '),
#'     reference_group = '  normal  ',
#'     group_label = "\nRearing Conditions"
#'     )
#' @export
#'
bootstrap_phenotype <- function(features_df_row = features_df_row,
                                condition_list = condition_list,
                                group_label = group_label,
                                reference_group = reference_group){
  
  ### Load the tidyverse package (necessary for this function)
  
  library(tidyverse)
  
  ### Assign 'feature_order' a vector of the numbers 1 to the number of columns in 'features_df_row'
  
  feature_order <- 1:ncol(features_df_row)
  
  ### This multi-nested while-loop calculates the mean, SD, and number of data points
  ###       for each experimental condition in 'condition_list' across each plasticity feature
  ###       in 'features_df_row'. Store each mean, SD and length that is calculated in 
  ###       3 different lists: mean_list, sd_list, and length_list. 
  
  while (T) {
    
    counter <- 0
    
    mean_list <- list()
    sd_list <- list()
    length_list <- list()
 
    
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
  
  ### Initialize an empty data frame with the same number of rows as there are items in the 
  ###       'condition_list', and a number columns equal to 1 plus the number of columns in
  ###       the 'features_df_row'. Call this empty data frame 'cols_df'.
  
  cols_df <- data.frame(matrix(ncol = ncol(features_df_row)+1,
                               nrow = length(condition_list)))
  
  ### Assign proper column headers to 'cols_df'. The first column will be 'Conditions',
  ###       while the remaining headers will be the same headers as 'features_df_row'.
  
  colnames(cols_df) <- c("Conditions",
                         colnames(features_df_row))
  
  ### This multi-nested while-loop runs for each item in 'condition_list' across
  ###       each column in 'features_df_row'. For each combination, it performs a 
  ###       bootstrap + Monte Carlo simulation that compares the 'reference' group's
  ###       average value on a given plasticity feature against a 'comparison' group on
  ###       the same feature. And stores the result as a colour-code in the 'cols_df' object.
  ###       If the comparison group's mean is bigger than the reference
  ###       mean, and this difference is significant, then it is coloured red on the phenotype.
  ###       If the reference mean is larger than the comparison mean (and this difference is significant),
  ###       then it appears as 'blue' on the phenotype. If it is a non-significant difference, it 
  ###       appears as grey on the phenotype. In all other cases, it appears as black on the 
  ###       phenotype, suggesting an error. 
  
  while (T) {
    
    counter <- 0
    
    for (i in 1:(length(condition_list))) {
      
      cols_df[i,1] <- condition_list[i]
      
      for (j in 1:ncol(features_df_row)) {
        
        counter <- counter + 1
        
        sim.comp <- rnorm(1000000, #simulate population of MD that has ONE MILLION =N
                          mean = mean_list[[counter]][[1]], #using the mean of MD
                          sd =  sd_list[[counter]][[1]]) #and the stdev of MD
        
        resamples.Normal.comp <- lapply(1:100000 #pull from a data set ONE HUNDRED THOUSAND times
                                      ,function(i) sample(sim.comp           #the data set to pull from is the simulated MD population
                                                          ,length_list[[1]][[1]]       #pull as many times as there are samples in NORMAL
                                                          ,replace = T)   #replace the samples after you pull them so each scenario has same P of happening
        )
        
        mean.resamples.Normal.comp <- sapply(resamples.Normal.comp,
                                             mean)
        
        Norm.v.comp.pval <- if(mean_list[[counter]][[1]]<mean_list[[j]][[1]]){

          1-(sum(mean.resamples.Normal.comp<mean_list[[j]][[1]])/100000)
          
        } else if (mean_list[[counter]][[1]]>mean_list[[j]][[1]]) {
          
          1-(sum(mean.resamples.Normal.comp>mean_list[[j]][[1]])/100000)
          
        } else {
          
          1-(sum(mean.resamples.Normal.comp == mean_list[[j]][[1]])/100000)
          
        }
        
        cols_df[i,1+j] <- if((Norm.v.comp.pval<(0.05/(ncol(features_df_row)*(length(condition_list)-1)))) && (mean_list[[counter]][[1]]<mean_list[[j]][[1]])) {
                    
                    "blue"
                    
                  } else if ((Norm.v.comp.pval<(0.05/(ncol(features_df_row)*(length(condition_list)-1)))) && (mean_list[[counter]][[1]]>mean_list[[j]][[1]])) {
                    
                    "red"
                    
                  } else if ((Norm.v.comp.pval>(0.05/(ncol(features_df_row)*(length(condition_list)-1))))) {
                    
                    "#dce0e5"
                    
                  } else {
                    
                    "black"
                    
                  }
        
      }
      
    }
    
    break
  }
  
  ### Assign a factor order to the 'conditions' column in 'cols_df'
  
  cols_df$Conditions <-  factor(cols_df$Conditions,
                                levels = unlist(condition_list),
                                ordered = T)
  
  ### Convert 'cols_df' from wide to long format
  
  long.phenotype.cols <- gather(cols_df,
                                feature,
                                color,
                                - Conditions)
  
  ### Visualize the bootstrap phenotype
  
  sig.phen <<- ggplot(long.phenotype.cols,
                      aes(x =  Conditions,
                          y = feature)) +
    geom_tile(data = long.phenotype.cols,
              width=0.95,
              height=0.95,
              color = "black",
              aes(fill = color))+
    scale_fill_identity()+
    theme(panel.background = element_rect(fill = 'white'),
          axis.text.x = element_text(vjust=0.5,size=11),
          axis.text.y = element_text(angle=0,vjust=0.5,size=11),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size=14, face = 'bold'),
          axis.title.y = element_text(size=14,  face = 'bold'))+
    scale_x_discrete(expand=c(0,0),
                     name=group_label,
                     position="bottom")+
    scale_y_discrete(expand=c(0,0),
                     name="Plasticity Features\n",
                     limits= rev(colnames(cols_df)[feature_order+1]),
                     position="left")+
    coord_fixed(ratio=0.25)
  
  ### Print the significance phenotype 
  
  print(sig.phen)
  
}
