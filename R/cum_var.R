#' Determine the number of PCA basis vectors that capture a set value of cumulative variance (threshold)
#'
#' Calculates the number of sequential PCA basis vectors required to explain a desired amount of the total cumulative variance.

#' @param pca.eig.3 Consults the pca$eig$3 product contained in the PCA output from the FactoMineR package.
#' @param thresh A desired threshold value of the percent total variance wished to be captured by sequential PCA basis vectors (values range from 0 – 100).
#' @return Only those PCA basis vectors that contain the set threshold value of cumulative variance.
#' @export

cum_var <- function(pca.eig.3, thresh){

  pca.eig.3 = pca.eig.3
  thresh = thresh

  while(T){

    cum.sum.counter <- 0

    for (i in 1:length(pca.eig.3)) {

      if(pca.eig.3[i] >= thresh){

        cum.sum.counter <- i

        break

      } else if (thresh == 100){

        cum.sum.counter <- length(pca$eig[,3])

      } else {

        next

      }

    }

    break

  }

  return(cum.sum.counter)

}
