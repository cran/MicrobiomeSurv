
#' This function is used for the first step of filtering which removes OTUs having all zeros (inactive OTUs).
#' The input is an OTU matrix with rows are OTUs and columns are subjects.
#' @param Micro.mat A large or small microbiome matrix.
#' A matrix with microbiome profiles where the number of rows should be equal to the number of taxa and number of columns should be equal to number of patients.
#' @return A smaller microbiome matrix.
#' \item{Micro.mat.trim}{The OTU matrix after removing all inactive OTUs}

#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{FirstFilter}}
#' @examples
#' # Preparing data for analysis at OTU level
#' data(Week3_otu)
#' Week3_otu = data.frame(Week3_otu)
#' otu_mat_w3 = t(data.matrix(Week3_otu[ , 1:2720]))
#' colnames(otu_mat_w3) = Week3_otu$SampleID

#' # Filtering first step
#' otu_w3 = FirstFilter(Micro.mat = otu_mat_w3)


#' @import stats


#' @export FirstFilter



FirstFilter = function(Micro.mat){
  zero.index.temp = array(0,nrow(Micro.mat))

  Micro.mat.new = Micro.mat

  j=1

  for(i in 1:nrow(Micro.mat.new)){

    if(min(Micro.mat.new[i,]) == max(Micro.mat.new[i,])){

      zero.index.temp[j] = i
      j = j+1
    }
  }

  first.zero = match(0, zero.index.temp)
  zero.index = zero.index.temp[1:first.zero-1]

  Micro.mat.trim = Micro.mat.new[-zero.index, ]

  return(Micro.mat.trim)

}

