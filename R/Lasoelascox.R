#'   Wapper function for glmnet
#'
#' The function uses the glmnet function to firstly do the variable selection either with Lasso, Elastic net or ridge regressions before the survial analysis.
#' The survival analysis is based on the selected taxa in the presence or absence of prognostic factors.
#'
#' This is a wrapper function for glmnet and it fits models using either Lasso, Elastic net and Ridge regressions.
#' This is done in the presence or absence of prognostic factors.
#' The prognostic factor when available will always be forced to be in the model so no penalty for it.
#' Optimum lambda will be used to select the non-zero shrinkage coefficients, the nonzero selceted taxa will thus be used in the survival analysis and in calculation of the risk scores.

#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Censor A vector of censoring indicator
#' @param Micro.mat A large or small microbiome matrix. A matrix with microbiome profiles where the number of rows is equal to the number of taxa and number of columns is equal to number of patients.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Plots A boolean parameter indicating if plots should be shown. Default is FALSE. If TRUE, the first plot is the partial likelihood deviance against the logarithmn of each lambda while the second is the coefficients versus the lambdas
#' @param Standardize A Logical flag for the standardization of the microbiome matrix, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is standardize=TRUE.
#' @param Alpha The mixing parameter for glmnet (see \code{\link[glmnet]{glmnet}}). The range is 0<= Alpha <= 1. The Default is 1
#' @param Fold number of folds to be used for the cross validation. Its value ranges between 3 and the number of subjects in the dataset
#' @param nlambda The number of lambda values - default is 100 as in glmnet.
#' @param Mean The cut off value for the classifier, default is the mean cutoff
#' @param Quantile If user want to use quantile as cutoff point. They need to specify Mean = FALSE and a quantile that they want to use. The default is the median cutoff
#' @return A object is returned with the following values
#'   \item{Coefficients.NonZero}{The coefficients of the selected taxa}
#'   \item{Selected.Mi}{The selected taxa}
#'   \item{n}{The number of selected taxa}
#'   \item{Risk.scores}{The risk scores of the subjects}
#'   \item{Risk.group}{The risk classification of the subjects based on the specified cutoff point}
#'   \item{SurvFit}{The cox analysis of the riskgroup based on the selected taxa and the prognostic factors}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}}
#' @seealso \code{\link[survival]{coxph}},
#' \code{\link[MicrobiomeSurv]{EstimateHR}},
#' \code{\link[glmnet]{glmnet}},
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

#' # Using the function
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
#' # View the selected taxa
#' lasso_fam_shan_w3$Selected.mi
#'
#' # Number of selected taxa
#' lasso_fam_shan_w3$n
#'
#' # View the classification group of each subject
#' lasso_fam_shan_w3$Risk.Group
#'
#' # View the survival analysis result
#' lasso_fam_shan_w3$SurvFit

#' @import stats
#' @import glmnet
#' @import survival
#' @import graphics
#' @importFrom coef density median p.adjust princomp qnorm quantile
#' @importFrom abline arrows axis barplot box boxplot legend lines par points


#' @export Lasoelascox


Lasoelascox = function (Survival,
                        Censor,
                        Micro.mat,
                        Prognostic,
                        Plots = FALSE,
                        Standardize = TRUE,
                        Alpha = 1,
                        Fold = 4,
                        nlambda = 100,
                        Mean = TRUE,
                        Quantile = 0.5){

  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Micro.mat)) stop("Argument 'Mdata' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")

  n.mi = dim(Micro.mat)[1]
  n.obs = dim(Micro.mat)[2]

  mi.names = rownames(Micro.mat)

  if(is.null(Prognostic)){
    Data = Micro.mat
    Data.Full = t(Micro.mat)
    Penalty = rep(1, row(Data.Full))
  } else{
    Data = t(Micro.mat)
    Data.Full = cbind(Data,Prognostic)
    Penalty = c(rep(1,ncol(Data)), rep(0, ncol(Prognostic)))
  }

  # Survival times must be larger than 0

  Survival[Survival <= 0] = stats::quantile(Survival, probs = 0.01)
  Lasso.Cox.CV = glmnet::cv.glmnet(x = data.matrix(Data.Full),
                                    y = survival::Surv(as.vector(Survival),as.vector(Censor) == 1),
                                    family = 'cox',
                                    alpha = Alpha,
                                    nfolds= Fold,
                                    nlambda = nlambda,
                                    penalty.factor = Penalty,
                                    standardize = Standardize)


  # Results of the cv.glmnet procedure

  Lambda = Lasso.Cox.CV$lambda.min
  Alllamda = Lasso.Cox.CV$lambda
  Coefficients = stats::coef(Lasso.Cox.CV, s =Lambda)
  Coefficients.NonZero = Coefficients[Coefficients[, 1] != 0, ]

  if (!is.null(dim(Coefficients.NonZero)))
  {
    Selected.mi = setdiff(colnames(Data.Full),colnames(Prognostic))
    Coefficients.NonZero = Coefficients[c(Selected.mi,colnames(Prognostic)),]
  }

  if (is.null(Prognostic)){
    Selected.mi = names(Coefficients.NonZero)
  } else{
    Selected.mi = setdiff(names(Coefficients.NonZero), colnames(Prognostic))
  }



  # What to do if no taxa is selected?
  # Decrease lambda until a taxon is selected
  # Going through the list of lambda values and repeat the above procedure

  n = length(Selected.mi)

  Lambda = Lasso.Cox.CV$lambda.min
  Lambda.Sequence = Lasso.Cox.CV$lambda
  Lambda.Index = which(Lambda.Sequence == Lasso.Cox.CV$lambda.min)

  Lambda.Index.Add = 0

  while (n <= 0) {
    Lambda.Index.Add = Lambda.Index.Add + 1
    Coefficients = stats::coef(Lasso.Cox.CV, s = Lambda.Index + Lambda.Index.Add)
    Coefficients.NonZero = Coefficients[Coefficients[, 1] != 0,]

    if (!is.null(dim(Coefficients.NonZero))){
      Selected.mi = setdiff(colnames(Data.Full), colnames(Prognostic))
      Coefficients.NonZero = Coefficients[c(Selected.mi,colnames(Prognostic)),]
    }

    if (is.null(Prognostic)){
      Selected.mi = names(Coefficients.NonZero)
    } else {
      Selected.mi = setdiff(names(Coefficients.NonZero), colnames(Prognostic))
    }
    Lambda = Lambda.Sequence[Lambda.Index + Lambda.Index.Add]
    n = length(Selected.mi)

    if (Lambda.Index.Add == length(Lambda.Sequence) & is.null(Selected.mi) == TRUE) stop("No taxa are selected")
  }

  Risk.Scores = as.vector(Coefficients.NonZero[Selected.mi] %*% t(Data[,Selected.mi]))

  # Estimate the HR by running the appropriate function

  Data.Survival = cbind(Survival,Censor)
  Estimation.HR = EstimateHR(Risk.Scores = Risk.Scores, Data.Survival = Data.Survival,
                              Prognostic=Prognostic, Plots = FALSE, Mean = Mean, Quantile = Quantile)
  Risk.Group = Estimation.HR$Riskgroup
  Cox.Fit.Risk.Group = Estimation.HR$SurvResult

  # Produce plots if requested

  if (Plots)
  {
    #graphics::par(mfrow=c(1,2))
    # Plot 1

    plot(log(Lasso.Cox.CV$lambda), Lasso.Cox.CV$cvm,
         main = expression(paste(alpha,'=', 1 , sep=" ")),
         ylim = c(min(Lasso.Cox.CV$cvlo), max(Lasso.Cox.CV$cvup)),
         xlab = expression(log(lambda)), ylab = 'Partial likelihood deviance',
         pch = 19, col = 'red')

    for (i in 1:length(Lasso.Cox.CV$cvm))
      graphics::lines(log(c(Lasso.Cox.CV$lambda[i], Lasso.Cox.CV$lambda[i])), c(Lasso.Cox.CV$cvlo[i], Lasso.Cox.CV$cvup[i]))

    Lasso.Cox = glmnet::glmnet(x = data.matrix(Data.Full),
                                y = survival::Surv(as.vector(Survival),as.vector(Censor) == 1),
                                family = 'cox',
                                alpha = Alpha,
                                nlambda = 100,
                                penalty.factor = Penalty,
                                standardize = Standardize)

    # Plot 2

    plot(Lasso.Cox, xvar = 'lambda', label = TRUE,xlab=expression(lambda))
    graphics::abline(v = log(Lambda), lwd = 2, lty = 2, col = 'red')

  }
  return(list(Coefficients.NonZero=Coefficients.NonZero,
              Selected.mi = Selected.mi,
              n = n,
              Risk.Scores=Risk.Scores,
              Risk.Group=Risk.Group,
              SurvFit=Cox.Fit.Risk.Group,
              Data.Survival = Data.Survival))
}
