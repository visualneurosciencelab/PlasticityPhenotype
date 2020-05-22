#' Plot amplitudes of each variable for a specified number of PCA basis vectors
#'
#' Plots the amplitude of each variable for the specified number of PCA basis vectors identified using the cum_var function.  These plots are used to identify novel plasticity features by comparing the amplitudes of each variable about te basis vectors (e.g. to create variable sums/indices).
#' @param cum.var The positive integer value outputted by the cum_var function (can also be a user-assigned positive integer ranging from 0 – the total number of variables within the dataset).
#' @param pca.var.coord Taken from the PCA object created using the FactoMineR package. Used to extract the vector for each variable.
#' @return A series of plots depicting the amplitude about a specified number of PCA basis vectors. Each PCA basis vectors is plotted separately, and each variable is represented along the x-axis. The user should use these plots to identify variables that might interact to create new features (e.g. such as indices  or sums of variables). The features identified after using the amplitude plots will be saved as a dataframe, which will serve as the input for feature_matrix.
#' @examples
#' amplitude_plots(
#'      cum.var = 4,
#'      pca.var.coord = pca.scaled$var$coord
#'      )
#' @export

amplitude_plots <- function(cum.var,pca.var.coord){


  cum.var = cum.var

  pca.var.coord = pca.var.coord

  while (T) {

    for (i in 1:cum.var){

      name <- paste0("Amplitude (Basis Vector ",i,")\n")

      VarCoordDim1<-data.frame(pca.var.coord[,i])

      setDT(VarCoordDim1, keep.rownames = TRUE)[]

      xyz <- ggplot(data=VarCoordDim1,
                    aes(rn,pca.var.coord...i.))+
        geom_col(colour="black")+
        scale_y_continuous(expand =c(0,0),name=name)+
        scale_x_discrete(limits=VarCoordDim1$rn)+
        theme(axis.line.y=element_line(),
              axis.line.x=element_line(),
              panel.grid=element_blank(),
              axis.text.x = element_text(angle=90,hjust =1,size=26),
              axis.text.y = element_text(angle=0,vjust=0.5,size=26),
              axis.title.x = element_text(size=28,face="bold", color = 'black'),
              axis.title.y = element_text(size=28,face="bold", color = 'black'),
              panel.grid.major = element_blank(),
              plot.margin = margin(2, 2, 2, 2, "cm"),
              panel.grid.minor = element_blank())+
        xlab("\nProteins")

      print(xyz)

    }

    break

  }

}
