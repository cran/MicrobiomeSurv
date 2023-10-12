
#' Classification, Survival Estimation and Visualization
#'
#' The function classifies subjects into Low and High risk groups using the risk scores based on the cut-off point which is the mean of the risk score.
#' Also visualize survival fit along with HR estimates.
#'
#' The risk scores obtained using the taxa is then used to generate the risk group by dividing subjects into low and high risk groups.
#' A Cox model is then fitted with the risk group as covariate in the presence or absence of  prognostic factors and or treatment effect.
#' The extent of survival in the risk groups is known
#' @param Risk.Scores A vector of risk scores with size equals to number of subjects obtained from (\code{\link[MicrobiomeSurv]{Lasoelascox}}).
#' @param Data.Survival A dataframe in which the first column is the Survival and the second column is the Censoring indicator for each subject.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect
#' @param Plots A boolean parameter indicating if plots should be shown. Default is FALSE.
#' @param Mean The cut off value for the classifier, default is the mean cutoff
#' @param Quantile If user want to use quantile as cutoff point. They need to specify Mean = FALSE and a quantile that they want to use. The default is the median cutoff
#' @return An object of is returned, which is a list with the results of the cox regression and some informative plot concerning survival of the risk group.
#'   \item{SurvResult}{The cox proportional regression result}
#'   \item{Riskgroup}{The riskgroup based on the riskscore and the cut off value and length is equal to number of subjects}
#'   \item{KMplot}{The Kaplan-Meier survival plot of the riskgroup}
#'   \item{SurvBPlot}{The distribution of the survival in the riskgroup}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}}, \code{\link[MicrobiomeSurv]{Lasoelascox}}
#' @examples
#' # Prepare data
#' data(Week3_response)
#' Week3_response = data.frame(Week3_response)
#' surv_fam_shan_w3 = data.frame(cbind(as.numeric(Week3_response$T1Dweek),
#' as.numeric(Week3_response$T1D)))
#' colnames(surv_fam_shan_w3) = c("Survival", "Censor")
#' prog_fam_shan_w3 = data.frame(factor(Week3_response$Treatment_new))
#' colnames(prog_fam_shan_w3) = c("Treatment")
#' data(fam_shan_trim_w3)
#' names_fam_shan_trim_w3 =
#' c("Unknown", "Lachnospiraceae", "S24.7", "Lactobacillaceae", "Enterobacteriaceae", "Rikenellaceae")
#' fam_shan_trim_w3 = data.matrix(fam_shan_trim_w3[ ,2:82])
#' rownames(fam_shan_trim_w3) = names_fam_shan_trim_w3

#' # Obtaning the risk score and data survival
#' lasso_fam_shan_w3 = Lasoelascox(Survival = surv_fam_shan_w3$Survival,
#'                                 Censor = surv_fam_shan_w3$Censor,
#'                                 Micro.mat = fam_shan_trim_w3,
#'                                 Prognostic = prog_fam_shan_w3,
#'                                 Plots = TRUE,
#'                                 Standardize = TRUE,
#'                                 Alpha = 1,
#'                                 Fold = 4,
#'                                 nlambda = 100,
#'                                 Mean = TRUE)
#'
#' # Using the function
#' est_HR_fam_shan_w3 = EstimateHR(Risk.Scores = lasso_fam_shan_w3$Risk.Scores,
#'                                 Data.Survival = lasso_fam_shan_w3$Data.Survival,
#'                                 Prognostic = prog_fam_shan_w3, Plots = TRUE,
#'                                 Mean = TRUE)

#' @import stats
#' @import ggplot2
#' @import survival
#' @import survminer
#' @import graphics
#' @import base
#' @importFrom abline arrows axis barplot box boxplot legend lines par points

#' @export EstimateHR

EstimateHR = function(Risk.Scores, Data.Survival,
                     Prognostic = NULL, Plots = FALSE,
                     Mean = TRUE, Quantile = 0.50)
{

  if (missing(Data.Survival)) stop('Surival data are missing!')
  if (missing(Risk.Scores)) stop('Risk scores are missing!')

  Survival = Data.Survival[ ,1]
  Censor   = Data.Survival[ ,2]

  n.Patients = length(Risk.Scores)


  # Group by high risk or low risk
  # Based on cut-off value
  if(Mean){
    Cut = mean(Risk.Scores)
  } else{
    Cut = stats::quantile(Risk.Scores, Quantile, type=6)
  }

  Risk.Group = ifelse(Risk.Scores >= Cut,'High risk', 'Low risk')

  # Make a data frame

  Data.Risk.Group = data.frame(Survival, Censor, Risk.Group)

  # Plot the observed Kaplan-Meier curves per group
  # Make a boxplot per group
  # Only if requested

  if (Plots)
  {
    # Kaplan Meier plot
    Data.KM = data.frame(Survival=Survival,Censor=Censor,Risk.Group=Risk.Group)
    KM  = survival::survfit(survival::Surv(Survival,Censor==1)~ Risk.Group,data=Data.KM)
    KMplot = survminer::ggsurvplot(fit=KM, surv.scale ='percent', risk.table=FALSE,
                         pval = TRUE,
                         ggtheme = ggplot2::theme_classic(),
                         palette='Dark2',conf.int=TRUE,legend='right',
                         legend.title='Risk',legend.labs=c('High risk','Low risk'),
                         data=Data.KM)
    KMplot

    # Boxplot
    Data.Boxplot = data.frame(Survival=Survival,Risk.Group=Risk.Group)
    Data.Boxplot$Risk.Group = factor(Data.Boxplot$Risk.Group, levels=c('Low risk','High risk'))
    SurvBPlot = ggplot2::ggplot(Data.Boxplot, ggplot2::aes(y = Survival,x = Risk.Group)) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = Risk.Group),width=0.4) +
      ggplot2::ggtitle('Low risk vs high risk')+
      ggplot2::ylab(expression('Survival time ')) +
      ggplot2::xlab(expression()) +
      ggplot2::theme_bw()+
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::theme(legend.position = 'none', plot.title = ggplot2::element_text(hjust = 0.5))

    SurvBPlot
    }


  # Do a (second) Cox regression including the group allocation
  # Group variable as well as prognostic variables are used

  if (is.null(Prognostic))
  {
    Cox.Fit.Risk.Group = survival::coxph(survival::Surv(Survival, Censor==1) ~ Risk.Group, data = Data.Risk.Group)
  }

  if (is.data.frame(Prognostic))
  {
    Number.Prognostic = ncol(Prognostic)
    Data.Risk.Group = data.frame(Data.Risk.Group,Prognostic)
    Names.Prognostic=colnames(Prognostic)
    eval(parse(text=paste('Cox.Fit.Risk.Group = survival::coxph(survival::Surv(Survival, Censor==1) ~ Risk.Group', paste('+',Names.Prognostic[1:Number.Prognostic],sep='',collapse =''),',data=Data.Risk.Group)' ,sep='')))
  }

  if (Plots){
    return(list(SurvResult=Cox.Fit.Risk.Group, Riskgroup=Risk.Group, KMplot=KMplot, SurvBPlot=SurvBPlot))
  } else{
    return(list(SurvResult=Cox.Fit.Risk.Group, Riskgroup=Risk.Group))
  }
}

