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

library(tidyverse)

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


  #print(sd_list)

  break

}


cols_df <- data.frame(matrix(ncol = ncol(features_df_row)+1,
                             nrow = length(condition_list)))

colnames(cols_df) <- c("Conditions",colnames(features_df_row))

while (T) {
  counter <- 0
#i = 1
  for (i in 1:(length(condition_list)#-1
               )) {

    cols_df[i#+1
            ,1] = condition_list[i#+1
                                    ]
    #j = 2
    for (j in 1:ncol(features_df_row)) {

      counter <- counter + 1
      mean.Normal <- mean_list[[j]][[1]]
      mean.MD1 = mean_list[[counter]][[1]]
      adj.pval = 0.05/length_list[[1]][[1]]


      sim.MD1 <- rnorm(1000000 #simulate population of MD that has ONE MILLION =N
                       , mean = mean_list[[counter]][[1]] #using the mean of MD
                       , sd =  sd_list[[counter]][[1]] #and the stdev of MD
      )

      resamples.NormalMD1 <- lapply(1:100000 #pull from a data set ONE HUNDRED THOUSAND times
                                    ,function(i) sample(sim.MD1           #the data set to pull from is the simulated MD population
                                                        , length_list[[1]][[1]]       #pull as many times as there are samples in NORMAL
                                                        , replace = T)   #replace the samples after you pull them so each scenario has same P of happening
      )

      mean.resamples.NormalMD1 <- sapply(resamples.NormalMD1, mean)

      CI.NormalMD1 <- quantile(mean.resamples.NormalMD1, c(0.001, 0.01,0.05,0.95, 0.99, 0.999))

      Norm.v.MD1.pval<- if(mean.MD1<mean.Normal){
        #print(1-(sum(mean.resamples.NormalMD1<mean.Normal)/100000))
        1-(sum(mean.resamples.NormalMD1<mean.Normal)/100000)
      } else if (mean.MD1>mean.Normal) {
        #print(1-(sum(mean.resamples.NormalMD1>mean.Normal)/100000))
        1-(sum(mean.resamples.NormalMD1>mean.Normal)/100000)
      } else {
        1-(sum(mean.resamples.NormalMD1 == mean.Normal)/100000)
      }


      cols_df[i#+1
              ,1+
                j] <- if( (Norm.v.MD1.pval<adj.pval) && (mean.MD1<mean.Normal) ) {
         #print("blue")
         "blue"
      } else if ((Norm.v.MD1.pval<adj.pval) && (mean.MD1>mean.Normal)) {
         #print("red")
         "red"
      } else if ((Norm.v.MD1.pval>adj.pval)) {
        #print("grey")
         "#dce0e5"
      } else {
        print("error")
        "black"
      }


    }

  }
  break
}

  
    cols_df$Conditions <-  factor(cols_df$Conditions,
                                levels = unlist(condition_list),
                                ordered = T)
  
long.phenotype.cols <- gather(cols_df, feature, color, - Conditions )



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
        axis.title.y = element_text(size=14,  face = 'bold')#,
        #plot.margin = margin(2, 2, 2, 2, "cm")
  )+
  scale_x_discrete(expand=c(0,0),
                   name=group_label,
                   # labels=factor(means[,grouping_mean_column])[group_select],
                   position="bottom")+
  scale_y_discrete(expand=c(0,0),
                   name="Plasticity Features\n",
                   limits= rev(colnames(cols_df)[#2:ncol(means)
                     feature_order+1]),
                   position="left")+
  coord_fixed(ratio=0.25)

print(sig.phen)

}
