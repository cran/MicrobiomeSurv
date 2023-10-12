#'Cross validation for the Taxon specific analysis
#'
#' The function performs cross validation for each taxon depending the number of fold which guides the division into the train and testing dataset.
#' The classifier is then obtained on the training dataset to be validated on the test dataset.
#'
#' This function performs the cross validation for taxon by taxon analysis.
#' The data will firstly be divided into data train dataset and test datset.
#' Furthermore, a taxon-specific model is fitted on train data and a classifier is built.
#' In addition, the classifier is then evaluated on test dataset for each particular taxon.
#' The Process is repeated for all the full or reduced taxa to obtaind the HR statistics of the low risk group.
#' The following steps depends on the number of cross validation specified.
#' @param Fold Number of times in which the dataset is divided. Default is 3 which implies dataset will be divided into three groups and 2/3 of the dataset will be the train datset and 1/3 will be to test the results.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Micro.mat A large or small microbiome profile matrix. A matrix with microbiome profiles where the number of rows should be equal to the number of taxa and number of columns should be equal to number of patients.
#' @param Censor A vector of censoring indicator.
#' @param Reduce A boolean parameter indicating if the microbiome profile matrix should be reduced, default is TRUE and larger microbiome profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of taxa (default is 5) to be selected from supervised PCA. This is valid only if th argument Reduce=TRUE.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Mean The cut off value for the classifier, default is the mean cutoff.
#' @param Quantile If users want to use quantile as cutoff point. They need to specify Mean = FALSE and a quantile that they wish to use. The default is the median cutoff.
#' @param Ncv The Number of cross validation loop. Default is 100.
#' @return A object of class \code{\link[MicrobiomeSurv]{cvmm}} is returned with the following values.
#'   \item{HRTrain}{The Train dataset HR statistics for each taxon by the number of CV.}
#'   \item{HRTest}{The Test dataset HR statistics for each taxon by the number of CV.}
#'   \item{train}{The selected subjects for each CV in the train dataset.}
#'   \item{test}{The selected subjects for each CV in the test dataset.}
#' \item{n.mi}{The number of taxa used in the analysis.}
#'  \item{Ncv}{The number of cross validation performed.}
#'  \item{Rdata}{The Microbiome data matrix that was used for the analysis either same as Micro.mat or a reduced version.}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}},
#' \code{\link[MicrobiomeSurv]{EstimateHR}}, \code{\link[MicrobiomeSurv]{MSpecificCoxPh}},
#' @examples
#' \donttest{
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
#' CVCox_taxon_fam_shan_w3 = CVMSpecificCoxPh(Fold=3,
#'                                            Survival = surv_fam_shan_w3$Survival,
#'                                            Micro.mat = fam_shan_trim_w3,
#'                                            Censor = surv_fam_shan_w3$Censor,
#'                                            Reduce=TRUE,
#'                                            Select=5,
#'                                            Prognostic=prog_fam_shan_w3,
#'                                            Mean = TRUE,
#'                                            Ncv=10)
#'
#' # Get the class of the object
#' class(CVCox_taxon_fam_shan_w3)     # An "cvmm" Class
#'
#' # Method that can be used for the result
#' show(CVCox_taxon_fam_shan_w3)
#' summary(CVCox_taxon_fam_shan_w3)
#' plot(CVCox_taxon_fam_shan_w3)
#' }
#' @import superpc
#' @import stats
#' @import lmtest
#' @import survival
#' @import methods

#' @export CVMSpecificCoxPh

CVMSpecificCoxPh=function(Fold=3,
                             Survival,
                             Micro.mat,
                             Censor,
                             Reduce=TRUE,
                             Select=5,
                             Prognostic=NULL,
                             Mean = TRUE,
                             Quantile = 0.5,
                             Ncv=100)
{
  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Micro.mat)) stop("Argument 'Micro.mat' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")

  n.mi.full = dim(Micro.mat)[1]
  n.obs = dim(Micro.mat)[2]
  coef = exp.coef = p.value.LRT = c(1 : n.mi.full)

  if (Reduce) {
    if (is.null(Prognostic)){
      DataForReduction=list(x=Micro.mat,y=Survival, censoring.status=Censor, mi.names=rownames(Micro.mat))
      TentativeList=names(sort(abs(superpc::superpc.train(DataForReduction, type="survival")$feature.scores),decreasing =TRUE))[1:Select]
      TentativeList

      ReduMicro.mat=Micro.mat[TentativeList, ]
    }

    if (!is.null(Prognostic)){
      if (is.data.frame(Prognostic)) {
        nPrgFac=ncol(Prognostic)
        NameProg=colnames(Prognostic)
      }

      if (dim(Prognostic)[2] == 1){
        cox.prog = eval(parse(text = paste("survival::coxph(survival::Surv(Survival, Censor)", paste(" ~ Prognostic[ ,1])"))))
      } else{
        cox.prog = eval(parse(text = paste("survival::coxph(survival::Surv(Survival, Censor) ~ NameProg[1]", paste("+", NameProg[2:nPrgFac], sep="", collapse =""), paste(")"))))
      }

      for(i in 1 : n.mi.full){
        xi = Micro.mat[i, ]
        datai = data.frame(Survival, Censor, xi, Prognostic)
        modeli = eval(parse(text = paste("survival::coxph(survival::Surv(Survival, Censor) ~ xi", paste("+", NameProg[1:nPrgFac], sep="", collapse =""), ",data=datai)" , sep="" )))
        coef[i] = round(summary(modeli)$coefficients[1,1], 4)
        exp.coef[i] = round(summary(modeli)$coefficients[1,2], 4)
        p.value.LRT[i] = round(lmtest::lrtest(cox.prog, modeli)[2,5], 4)
      }

      p.value = round(stats::p.adjust(p.value.LRT, method = "BH", n = length(p.value.LRT)), 4)
      summary = cbind(coef, exp.coef, p.value.LRT, p.value)
      rownames(summary) = rownames(Micro.mat)
      colnames(summary) = c("coef", "exp.coef", "p.value.LRT", "p.value")
      TentativeList=names(sort(abs(summary[ ,"p.value.LRT"]),decreasing =TRUE))[1:Select]
      TentativeList

      ReduMicro.mat = Micro.mat[TentativeList, ]
      summary.reduced = summary[TentativeList, ]
    }

  } else {
    ReduMicro.mat = Micro.mat
  }


  n.mi=nrow(ReduMicro.mat)
  gNames=rownames(ReduMicro.mat)

  n.train=(n.obs-floor(n.obs/Fold))
  n.test=floor(n.obs/Fold)
  train =matrix(0,Ncv,n.train)
  test  =matrix(0,Ncv,n.test)

  #     THE HAZARD RATIO
  HRTrain=array(NA,dim=c(n.mi,4,Ncv))
  HRTest=array(NA,dim=c(n.mi,4,Ncv))

  #set.seed(123)
  pIndex = c(1:n.obs)

  taxoncox=function(itaxon,Prognostic,Survival,Censor,index){

    if (is.null(Prognostic)) {

      scdata = data.frame(Survival=Survival[index], Censor=Censor[index], itaxon=itaxon[index])
      model1 = survival::coxph(survival::Surv(Survival, Censor==1) ~ itaxon, data=scdata)
    }
    if (!is.null(Prognostic)) {
      if (is.data.frame(Prognostic)) {
        nprog=ncol(Prognostic)
        prognames=colnames(Prognostic)
        scdata = data.frame(Survival[index], Censor[index], itaxon[index], Prognostic[index,])
        colnames(scdata) = c("Survival", "Censor", "itaxon", prognames)
        eval(parse(text=paste( "model1 =survival::coxph(survival::Surv(Survival, Censor==1) ~ itaxon", paste("+",prognames[1:nprog],sep="",collapse =""),",data=scdata)" ,sep="")))
      } else {
        stop(" Argument 'Prognostic' is NOT a data frame ")
      }

    }

    return(list(model1=model1,itaxon=itaxon,scdata=scdata))
  }


  for (j in 1:Ncv){
    message('Cross validation loop ',j)
    train[j,] =sort(sample(pIndex, n.train, replace=F))
    test[j,] =c(1:n.obs)[-c(intersect(train[j, ], c(1:n.obs)))]

    for (i in 1:n.mi){  #---------------------------  STRAT FOR LOOP ------------------

      #training set ---------------------------------------------------------
      TrainTemp=taxoncox(ReduMicro.mat[i,train[j,]], Prognostic, Survival, Censor, train[j,])

      m1 = TrainTemp$model1
      Trtandtaxon=summary(m1)[[7]][c("itaxon"), 1]
      p1.train = Trtandtaxon[1]*ReduMicro.mat[i, train[j,]]
      p1.test = Trtandtaxon[1]*ReduMicro.mat[i, test[j,]]

      Tempitaxon =EstimateHR(p1.train, Data.Survival = data.frame(Survival=Survival[train[j,]], Censor=Censor[train[j,]]),
                              Prognostic = data.frame(Prognostic[train[j,],]),Plots = FALSE, Mean = TRUE, Quantile = Quantile)

      HRTrain[i,,j]=summary(Tempitaxon$SurvResult)[[8]][1,]

      #testing set ---------------------------------------------------------

      TempitaxonTE =EstimateHR(p1.test, Data.Survival=data.frame(Survival=Survival[test[j,]], Censor=Censor[test[j,]]),
                                Prognostic=data.frame(Prognostic[test[j,],]), Plots = FALSE, Mean = TRUE, Quantile = Quantile)
      HRTest[i,,j]=summary(TempitaxonTE$SurvResult)[[8]][1,]

    }#---------------------------  END OF  FOR LOOP for taxon ----------------------

  }#--------------------------END OF loop over CV ---------------------------------------

  return(methods::new("cvmm",HRTrain=HRTrain,HRTest=HRTest,train=train,test=test,n.mi=n.mi,Ncv=Ncv,Rdata = ReduMicro.mat))
}
