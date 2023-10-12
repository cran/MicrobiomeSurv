#' This function returns a matrix with rows are Micros and 9 columns containing number and the proportion of zeros per groups of treatments and in total.

#'@param Micro.mat Micro matrix (rows are Micros, columns are subjects)
#'@param groups Treatment groups or groups of any binary variables
#'@param week A specific time point. To use when having different time points in the dataset.
#'@param n.obs Number of patients.
#'@param n.control Number of patients in control group or in the first group.
#'@param n.treated Number of patients in treated group or in the second group.
#'@param n.mi Number of taxa.
#'@param plot A boolean parameter indicating if the plot should be shown. Default is FALSE.

#'@return A matrix with information of number and the proportion of zeros per groups.
#' \item{zero.per.group}{A matrix with rows are Micros and 9 columns containing number and the proportion of zeros per groups of treatments and in total.}
#' \item{plot}{Plot percentage of zeros per group}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{ZerosPerGroup}}
#' @examples
#' \donttest{
#' # Preparing data for analysis at OTU level
#' data(Week3_otu)
#' data(Week3_response)
#' Week3_otu = data.frame(Week3_otu)
#' otu_mat_w3 = t(data.matrix(Week3_otu[ , 1:2720]))
#' n_obs = dim(otu_mat_w3)[2]
#' n_control = table(Week3_response$Treatment_new)[1]
#' n_treated = table(Week3_response$Treatment_new)[2]
#' n_otu = dim(otu_mat_w3)[1]

#' # Calculate zeros per groups
#' zero_per_group_otu_w3 = ZerosPerGroup(Micro.mat = otu_mat_w3,
#'                                      groups = Week3_response$Treatment_new,
#'                                      week = 3,
#'                                      n.obs = n_obs,
#'                                      n.control = n_control,
#'                                      n.treated = n_treated,
#'                                      n.mi = n_otu,
#'                                      plot = TRUE)
#' }
#' @import stats
#' @import ggplot2

#' @export ZerosPerGroup

ZerosPerGroup = function(Micro.mat, groups, week = 0,
                         n.obs = n.obs, n.control = n.control,
                         n.treated = n.treated, n.mi = n.mi,
                         plot = FALSE){

  subjects = colnames(Micro.mat)

  Micro.treatment = data.frame(rbind(Micro.mat, treatment = groups)) ## 0 - Control, 1 - Treated

  Control = Treated = array (0,nrow(Micro.mat))

  for(i in 1:(nrow(Micro.treatment)-1)) {


    for(j in 1:ncol(Micro.treatment)) {

      if(Micro.treatment[i,j] == 0 && Micro.treatment["treatment",j] == 0){

        Control[i] = Control[i] + 1

      } else if(Micro.treatment[i,j] == 0 && Micro.treatment["treatment",j] == 1){

        Treated[i] = Treated[i] + 1

      }
    }
  }

  zero.per.group = matrix(0, nrow(Micro.mat), 9)
  rownames(zero.per.group) = rownames(Micro.mat)
  colnames(zero.per.group) = c("zero.ctrl", "propzero.ctrl", "nCtrl",
                                "zero.Treated", "propzero.Treated", "nTreated",
                                "zero.total", "propzero.total","nTotal")

  zero.per.group[,1] = Control
  zero.per.group[,2] = Control/n.control
  zero.per.group[,3] = n.control
  zero.per.group[,4] = Treated
  zero.per.group[,5] = Treated/n.treated
  zero.per.group[,6] = n.treated
  zero.per.group[,7] = Control + Treated
  zero.per.group[,8] = (Control + Treated)/n.obs
  zero.per.group[,9] = n.obs
  n = rep(NA, 2*n.mi)
  pp = rep(NA, 2*n.mi)
  Treatment = rep(NA, 2*n.mi)
  prop.zero.group = data.frame(Microu = c(row.names(zero.per.group), row.names(zero.per.group)),
                               n = 1: n.mi,
                               pp = c(zero.per.group[,2], zero.per.group[,5]),
                               Treatment = c(rep("3PATCON", n.mi), rep("3PAT", n.mi)))

  # Produce plots if requested
  if (plot == TRUE){
    plot = ggplot2::ggplot(prop.zero.group, ggplot2::aes(x=n, y=pp, color= Treatment))+
      ggplot2::geom_point()+
      ggplot2::facet_grid( ~ Treatment)+
      ggplot2::theme(legend.position="bottom", legend.direction="horizontal",
            legend.title = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::scale_color_manual("Treatment", values = c("3PAT"="red4","3PATCON"="blue"),
                         labels=c("3PAT"="3PAT","3PATCON"="3PATCON"),
                         guide = ggplot2::guide_legend(direction="horizontal")) +
      ggplot2::xlab(" ") + ggplot2::ylab("% of 0s")
    plot
  }

  return(list(zero.per.group, plot))
}
