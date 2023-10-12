
#' This function is used for the second step of filtering which removes OTUs based on a threshold.

#'@param zero.per.group a n x 9 matrix. Columns are number of zero in control groups, proportion of zeros in control group,
#' number of subject in control group, number of zero in treated groups, proportion of zeros in treated group,
#' number of subject in treated group, total number of zeros, proportion of zeros in total, number of subject
#'@param Micro.mat OTU matrix (rows are otus, columns are subjects)
#'@param threshold user can choose. For instance, if threshold is 0.7, the function will remove OTUs having at least 70\% of zeros in one of two groups
#'@param week A specific time point. To use when having different time points in the dataset.
#'@return A smaller microbiome matrix.
#' \item{Micro.mat.new}{an smaller OTU matrix with less OTUs}

#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#'
#' @seealso\code{\link[MicrobiomeSurv]{SecondFilter}}
#' @examples
#' # Read dataset
#' data(Week3_otu)
#' Week3_otu = data.frame(Week3_otu)
#' otu_mat_w3 = t(data.matrix(Week3_otu[ , 1:2720]))
#'
#' # Import dataset from the result of zero_per_group
#' data(data_zero_per_group_otu_w3)
#'
#' # Using the function
#' otu_trim_w3 = SecondFilter(zero.per.group = data_zero_per_group_otu_w3,
#'                            Micro.mat = otu_mat_w3, threshold = 0.7, week = 3)

#' @import stats



#' @seealso \code{\link[MicrobiomeSurv]{SecondFilter}}
#' @import utils
#' @import stats
#' @import Biobase

#' @export SecondFilter




SecondFilter = function(zero.per.group, Micro.mat, threshold = 0.7, week = 0) {

  otu.exclude.indices = array(0, nrow(zero.per.group))
  zeros = data.frame(zero.per.group)

  j=1
  for(i in 1:nrow(zero.per.group))
  {
    if(zeros$propzero.ctrl[i] >= threshold || zeros$propzero.Treated[i] >= threshold){

      otu.exclude.indices[j] = i
      j = j+1
    }

  }

  first.zero.temp = match(0, otu.exclude.indices)
  otu.exclude.indices = otu.exclude.indices[1:(first.zero.temp-1)]

  Micro.mat.new = Micro.mat[-otu.exclude.indices,]

  return(Micro.mat.new)

}



