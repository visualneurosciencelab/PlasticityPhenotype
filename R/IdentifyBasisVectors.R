#' Identify PCA basis vectors that describe a thresholdedd value of cumulative variance
#'
#' Identifies and stores the PCA basis vectors that contain the desired amount of cumulative variance
#' @param pca.eig.3 Consults the pca.eig.3 product contained in the PCA output from the FactoMineR package
#' @param thresh Assign the threshold value of cumulative variance captured by sequential PCA basis vectors
#' @return Only those PCA basis vectors that contain the set threshold value of cumulative variance
#' @export


cum_var <- function(pca.eig.3, thresh){

  pca.eig.3 = pca.eig.3
  thresh = thresh

  while(T){

    cum.sum.counter <- 0

    for (i in 1:length(pca.eig.3)) {

      if(pca.eig.3[i] >= thresh){ ## You can change the value of 80
        ## to whatever minimum variance threshold you desire.

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
