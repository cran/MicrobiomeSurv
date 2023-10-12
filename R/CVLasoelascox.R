#' Cross Validations for Lasso Elastic Net Survival predictive models and Classification
#'
#' The function does cross validation for Lasso, Elastic net and Ridge regressions models before the survial analysis and classification.
#' The survival analysis is based on the selected taxa in the presence or absence of prognostic factors.
#'
#' The function performs the cross validations for Lasso, Elastic net and Ridge regressions models for Cox proportional hazard model.
#' Taxa are selected at each iteration and then use for the classifier.
#' Which implies that predictive taxa is varied from one cross validation to the other depending on selection.
#' The underline idea is to investigate the Hazard Ratio for the train and test data based on the optimal lambda selected for the non-zero shrinkage coefficients, the nonzero selected taxa will thus be used in the survival analysis and in calculation of the risk scores for each sets of data.
#' @param Survival A vector of survival time with length equals to number of subjects.
#' @param Censor A vector of censoring indicator.
#' @param Micro.mat A large or small microbiome profile matrix. A matrix with microbiome profiles where the number of rows is equal to the number of taxa and number of columns is equal to number of patients.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Standardize A Logical flag for the standardization of the microbiome matrix, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is standardize=TRUE.
#' @param Alpha The mixing parameter for glmnet (see \code{\link[glmnet]{glmnet}}). The range is 0<= Alpha <= 1. The Default is 1.
#' @param Fold Number of folds to be used for the cross validation. Its value ranges between 3 and the number of subjects in the dataset.
#' @param Ncv Number of validations to be carried out. The default is 10.
#' @param nlambda The number of lambda values - default is 100 as in glmnet.
#' @param Mean The cut off value for the classifier, default is the mean cutoff.
#' @param Quantile If users want to use quantile as cutoff point. They need to specify Mean = FALSE and a quantile that they wish to use. The default is the median cutoff.

#' @return A object of class \code{\link[MicrobiomeSurv]{cvle}} is returned with the following values
#'   \item{Coef.mat}{A matrix of coefficients with rows equals to number of cross validations and columns equals to number of taxa.}
#'   \item{lambda}{A vector of estimated optimum lambda for each iterations.}
#'   \item{n}{A vector of the number of selected taxa.}
#'   \item{HRTrain}{A matrix of survival information for the training dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.}
#'   \item{HRTest}{A matrix of survival information for the test dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.}
#'   \item{pld}{A vector of partial likelihood deviance at each cross validations.}
#'   \item{Mi.mat}{A matrix with 0 and 1. Number of rows equals to number of iterations and number of columns equals to number of 1 taxon indicates that the particular taxon was selected or had nonzero coefficient and otherwise it is zero.}
#'   \item{Micro.mat}{The Microbiome data matrix that was used for the analysis either same as Mdata or a reduced version.}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}}, \code{\link[MicrobiomeSurv]{EstimateHR}}, \code{\link[glmnet]{glmnet}}, \code{\link[MicrobiomeSurv]{Lasoelascox}}
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
#'
#' # Using the function
#' CV_lasso_fam_shan_w3 = CVLasoelascox(Survival = surv_fam_shan_w3$Survival,
#'                                      Censor = surv_fam_shan_w3$Censor,
#'                                      Micro.mat = fam_shan_trim_w3,
#'                                      Prognostic = prog_fam_shan_w3,
#'                                      Standardize = TRUE,
#'                                      Alpha = 1,
#'                                      Fold = 4,
#'                                      Ncv = 10,
#'                                      nlambda = 100)
#'
#' # Number of selected taxa per CV
#' CV_lasso_fam_shan_w3@n
#'
#' # Get the matrix of coefficients
#' CV_lasso_fam_shan_w3@Coef.mat
#'
#' # Survival information of the train dataset
#' CV_lasso_fam_shan_w3@HRTrain
#'
#' # Survival information of the test dataset
#' CV_lasso_fam_shan_w3@HRTest

#' @export CVLasoelascox
#' @import superpc
#' @import stats
#' @import methods
#' @import grDevices
#' @import graphics
#' @import glmnet
#' @import survival

CVLasoelascox = function(Survival,
                          Censor,
                          Micro.mat,
                          Prognostic,
                          Standardize = TRUE,
                          Alpha = 1,
                          Fold = 4,
                          Ncv = 10,
                          nlambda = 100,
                          Mean = TRUE,
                          Quantile = 0.5){

  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Micro.mat)) stop("Argument 'Micro.mat' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")



  mi.names = rownames(Micro.mat)

  n.mi = nrow(Micro.mat)
  n.patients = ncol(Micro.mat)
  n.train = (n.patients-floor(n.patients/Fold))
  n.test = floor(n.patients/Fold)
  cv.train = matrix(0,Ncv,n.train)
  cv.test  = matrix(0,Ncv,n.test)
  n.mi = rep(0, Ncv)

  #optimum lambda
  lambda = rep(NA, Ncv)
  pld = rep(NA, Ncv) # partial likelihood deviance

  #HR--------
  HRTrain = matrix(NA,nrow=Ncv,ncol=3)
  HRTest = matrix(NA,nrow=Ncv,ncol=3)



  if(is.null(Prognostic)){
    Data = Micro.mat
    Data.Full = t(Micro.mat)
    Penalty = rep(1,row(Data.Full))
  } else{
    Data = t(Micro.mat)
    Data.Full = cbind(Data,Prognostic)
    Penalty = c(rep(1,ncol(Data)),rep(0,ncol(Prognostic)))
  }

  # Survival times must be larger than 0
  #Survival[Survival <= 0] = stats::quantile(Survival, probs = 0.01)
  perPrognostic=NULL
  coef.mat = mi.mat = matrix(0,nrow=Ncv,ncol=nrow(Micro.mat))

  pIndex = c(1:n.patients)

  for (i in 1:Ncv){
    #set.seed(i)
    message('Cross validation loop ',i)

    cv.train[i,] = sort(sample(pIndex,n.train,replace=F) )
    cv.test[i,] = c(1:n.patients)[-c(intersect(cv.train[i,] ,c(1:n.patients)))]

    Stime=Survival[cv.train[i,]]
    sen= Censor[cv.train[i,]]
    Data.Full2 =Data.Full[cv.train[i,],]

    Lasso.Cox.CV = glmnet::cv.glmnet(x = data.matrix(Data.Full2),
                              y = survival::Surv(as.vector(Stime),as.vector(sen) == 1),
                              family = 'cox',
                              alpha = Alpha,
                              nfolds= Fold,
                              nlambda = nlambda,
                              penalty.factor = Penalty,
                              standardize = Standardize)

    # Results of the cv.glmnet procedure

    Lambda = Lasso.Cox.CV$lambda.min
    Alllamda = Lasso.Cox.CV$lambda
    Coefficients = stats::coef(Lasso.Cox.CV, s = Lambda)
    Coefficients.NonZero = Coefficients[Coefficients[, 1] != 0, ]


    pld[i]=Lasso.Cox.CV$cvm[Lasso.Cox.CV$lambda==Lasso.Cox.CV$lambda.min]  # partial likelihood deviance

    if (is.null(Prognostic)){
      Selected.mi = names(Coefficients.NonZero)
    }else{
      Selected.mi = setdiff(names(Coefficients.NonZero), colnames(Prognostic))
    }

    n.mi[i] = length(Selected.mi)

    # What to do if no taxa is selected?
    # Decrease lambda until a taxon is selected
    # Going through the list of lambda values and repeat the above procedure

    Lambda = Lasso.Cox.CV$lambda.min
    Lambda.Sequence = Lasso.Cox.CV$lambda
    Lambda.Index = which(Lambda.Sequence == Lasso.Cox.CV$lambda.min)

    Lambda.Index.Add = 0

    while (n.mi[i] == 0) {
      Lambda.Index.Add = Lambda.Index.Add + 1
      Coefficients = stats::coef(Lasso.Cox.CV, s = Lambda.Sequence[Lambda.Index + Lambda.Index.Add])
      Coefficients.NonZero = Coefficients[Coefficients[, 1] != 0,]

      if (is.null(Prognostic)){
        Selected.mi = names(Coefficients.NonZero)
      } else {
        Selected.mi = setdiff(names(Coefficients.NonZero), colnames(Prognostic))
      }

      Lambda = Lambda.Sequence[Lambda.Index + Lambda.Index.Add]
      n.mi[i] = length(Selected.mi)
      if (Lambda.Index.Add == length(Lambda.Sequence) & is.null(Selected.mi) == TRUE) stop("No taxa are selected")
    }

    lambda[i]=Lambda
    mi.mat[i, is.element(rownames(Micro.mat),Selected.mi)] = 1
    coef.mat[i,is.element(rownames(Micro.mat),Selected.mi)] = Coefficients.NonZero[Selected.mi]

    scores.train = as.vector(Coefficients.NonZero[Selected.mi] %*% t(Data[cv.train[i,], Selected.mi]))
    scores.test = as.vector(Coefficients.NonZero[Selected.mi] %*% t(Data[cv.test[i,],Selected.mi]))


    #######################
    ## train set ###########
    Sdata=data.frame(Survival=Survival[cv.train[i,]],Censor=Censor[cv.train[i,]])

    if (!is.null(Prognostic)) {perPrognostic=as.data.frame(Prognostic[cv.train[i,],])}
    colnames(perPrognostic) = colnames(Prognostic)

    Results1=EstimateHR(Risk.Scores=scores.train, Data.Survival =Sdata, Prognostic = perPrognostic, Plots = FALSE, Mean = TRUE, Quantile = Quantile)

    if (!is.null(Prognostic)) HRTrain[i,]=(summary(Results1$SurvResult)[[8]])[1,][-2]
    if ( is.null(Prognostic)) HRTrain[i,]= summary(Results1$SurvResult)[[8]][-2]

    #######################
    ## test set ###########
    Sdata=data.frame(Survival=Survival[cv.test[i,]],Censor=Censor[cv.test[i,]])

    if (!is.null(Prognostic)) {perPrognostic=as.data.frame(Prognostic[cv.test[i,],])}
    colnames(perPrognostic) = colnames(Prognostic)

    Results2=EstimateHR(Risk.Scores =scores.test, Data.Survival = Sdata, Prognostic = perPrognostic, Plots = FALSE, Mean = TRUE, Quantile = Quantile)

    if (!is.null(Prognostic)) HRTest[i,]=(summary(Results2$SurvResult)[[8]])[1,][-2]
    if ( is.null(Prognostic)) HRTest[i,]= summary(Results2$SurvResult)[[8]][-2]
  }


  return(methods::new("cvle", Coef.mat=coef.mat, lambda=lambda, n=n.mi, mi.mat=mi.mat,
              HRTrain=HRTrain, HRTest=HRTest, pld=pld, Micro.mat=Micro.mat))
}

