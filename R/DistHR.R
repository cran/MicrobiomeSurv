#'Null Distribution of the Estimated HR
#'
#' This function generates the null distribution of the HR by permutation approach either using a large microbiome matrix or a reduced version by supervised pca approach.
#' Several ways of permutation setting can be implemented.
#' That is, the function can be used to generate null distributions for four different validation schemes which are PLS based, PCA based, Majority votes based and Lasso based.
#' Note this function internally calls function  \code{\link[MicrobiomeSurv]{SurvPcaClass}}, \code{\link[MicrobiomeSurv]{SurvPlsClass}}, \code{\link[MicrobiomeSurv]{Majorityvotes}}, and \code{\link[MicrobiomeSurv]{Lasoelascox}}.
#' @param Survival A vector of survival time with length equals to number of subjects.
#' @param Censor A vector of censoring indicator.
#' @param Micro.mat A large or small microbiome profile matrix. A matrix with microbiome profiles where the number of rows should be equal to the number of taxa and number of columns should be equal to number of patients.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Mean The cut off value for the classifier, default is the mean cutoff.
#' @param Quantile If user want to use quantile as cutoff point. They need to specify Mean = FALSE and a quantile that they want to use. The default is the median cutoff.
#' @param Reduce A boolean parameter indicating if the microbiome profile matrix should be reduced, default is TRUE and larger microbiome profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of taxa (default is 5) to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE.
#' @param nperm Number of permutations to be used and default 100.
#' @param Method A multiplicity adjustment Method that user can choose. The default is BH Method.
#' @param case There are seven different ways on how to call this argument:
#' \enumerate{
#' \item{Permute survival only.}
#' \item{Permute survival and rows of data frame of the prognostic factors.}
#' \item{Permute survival, rows of data frame of the prognostic factors, columns of microbiome matrix independently.}
#' \item{Permute microbiome matrix only.}
#' }
#' @param Validation There are four different validation schemes where the null distribution can be estimated. That is c("PLSbased","PCAbased","L1based","MVbased").
#' @return A object of class \code{\link[MicrobiomeSurv]{perm}} is returned with the following values
#'   \item{HRobs}{Estimated HR for low risk group on the original data.}
#'   \item{HRperm}{Estimated HR for low risk group on the permuted data.}
#'   \item{nperm}{Number of permutations carried out.}
#'   \item{Validation}{The validation scheme that was used.}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}}, \code{\link[MicrobiomeSurv]{EstimateHR}}, \code{\link[MicrobiomeSurv]{SurvPcaClass}}, \code{\link[MicrobiomeSurv]{SurvPlsClass}}, \code{\link[MicrobiomeSurv]{Majorityvotes}}, \code{\link[MicrobiomeSurv]{Lasoelascox}}
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
#' DistHR_fam_shan_w3 = DistHR(Survival = surv_fam_shan_w3$Survival,
#'                             Micro.mat = fam_shan_trim_w3,
#'                             Censor = surv_fam_shan_w3$Censor,
#'                             Prognostic=prog_fam_shan_w3,
#'                             Mean = TRUE,
#'                             Quantile=0.5,
#'                             Reduce= FALSE,
#'                             Select = 5,
#'                             nperm=100,
#'                             case=4,
#'                             Method = "BH",
#'                             Validation="PCAbased")
#'
#' # Method that can be used for the result
#' show(DistHR_fam_shan_w3)
#' summary(DistHR_fam_shan_w3)
#' plot(DistHR_fam_shan_w3)
#' }
#' @import stats
#' @import superpc
#' @import lmtest
#' @import base
#' @import methods
#' @import survival
#' @importFrom coef density median p.adjust princomp qnorm quantile
#' @export DistHR

DistHR=function(Survival,
                 Censor,
                 Micro.mat,
                 Prognostic=NULL,
                 Mean = TRUE,
                 Quantile=0.5,
                 Reduce= FALSE,
                 Select = 5,
                 nperm=100,
                 case=2,
                 Method = "BH",
                 Validation=c("PLSbased","PCAbased","L1based","MVbased")

){

  Validation = match.arg(Validation)
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

      p.value = round(stats::p.adjust(p.value.LRT, method = Method, n = length(p.value.LRT)), 4)
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
  HRlowPerm=matrix(NA,nrow=nperm,ncol=3)
  HRlowObs=as.vector(rep(NA,3))


  #set.seed(123)
  ind.s=ind.mi=ind.w=matrix(NA,nrow=nperm,ncol=n.obs)

  for (i in 1:nperm) {
    ind.s[i,]=1:n.obs
    ind.w[i,]=1:n.obs
    ind.mi[i,]=1:n.obs
  }

  switch(case,

         {#case 2: permute survival
           for (i in 1:nperm) {
             ind.s[i,]=sample(c(1:n.obs),replace=F)
           }

         },

         {#case 3: permute survival, prognostic
           for (i in 1:nperm) {
             ind.s[i,]=sample(c(1:n.obs),replace=F)
           }
           ind.w=ind.s

         },


         {#case 4 b: permute survival, prognostic and permute microbiome independently
           for (i in 1:nperm) {
             ind.s[i,]=sample(c(1:n.obs),replace=F)
           }

           ind.w=ind.s
           for (i in 1:nperm) {
             ind.mi[i,]=sample(c(1:n.obs),replace=F)
           }

         },

         {#case 4 : permute  microbiome only
           for (i in 1:nperm) {
             ind.mi[i,]=sample(c(1:n.obs),replace=F)
           }
         }
  )

  perPrognostic=NULL

  if (Validation=="PLSbased") {
    for (i in 1:nperm) {
      message('Permutation loop ',i)
      if (!is.null(Prognostic)) perPrognostic=Prognostic[ind.w[i,],]

      Temp=SurvPlsClass(Survival=Survival[ind.s[i,]],
                         Micro.mat= ReduMicro.mat[,ind.mi[i,]],
                         Censor= Censor[ind.s[i,]],
                         Reduce=Reduce,
                         Select=Select,
                         Prognostic=data.frame(perPrognostic),
                         Plots = FALSE,
                         Mean = TRUE,
                         Quantile = Quantile)

      if (!is.null(Prognostic))  HRlowPerm[i,]=summary(Temp$SurvFit)[[8]][1,][-2]
      if ( is.null(Prognostic))  HRlowPerm[i,]=summary(Temp$SurvFit)[[8]][-2]
    }


    TempObs=SurvPlsClass(Survival,
                          Micro.mat=ReduMicro.mat,
                          Censor,
                          Reduce=Reduce,
                          Select=Select,
                          Prognostic=Prognostic,
                          Plots = FALSE,
                          Quantile = Quantile)

    if (!is.null(Prognostic)) HRlowObs=summary(TempObs$SurvFit)[[8]][1,][-2]
    if ( is.null(Prognostic)) HRlowObs=summary(TempObs$Survfit)[[8]][-2]
  }



  if (Validation=="PCAbased") {
    for (i in 1:nperm) {
      message('Permutation loop ',i)
      if (!is.null(Prognostic)) perPrognostic=Prognostic[ind.w[i,],]

      Temp=SurvPcaClass(Survival=Survival[ind.s[i,]],
                         Micro.mat=ReduMicro.mat[,ind.mi[i,]],
                         Censor=Censor[ind.s[i,]],
                         Reduce=Reduce,
                         Select=Select,
                         Prognostic=data.frame(perPrognostic),
                         Plots = FALSE,
                         Quantile = Quantile)

      if (!is.null(Prognostic))  HRlowPerm[i,]=summary(Temp$SurvFit)[[8]][1,][-2]
      if ( is.null(Prognostic))  HRlowPerm[i,]=summary(Temp$SurvFit)[[8]][-2]

    }
    TempObs=SurvPcaClass(Survival,
                          Micro.mat=ReduMicro.mat,
                          Censor,
                          Reduce=Reduce,
                          Select=Select,
                          Prognostic=Prognostic,
                          Plots = FALSE,
                          Mean = TRUE,
                          Quantile = Quantile)

    if (!is.null(Prognostic)) HRlowObs=summary(Temp$SurvFit)[[8]][1,][-2]
    if ( is.null(Prognostic)) HRlowObs=summary(Temp$SurvFit)[[8]][-2]
  }

  if (Validation=="MVbased") {
    for (i in 1:nperm) {
      message('Permutation loop ',i)
      if (!is.null(Prognostic)) perPrognostic=Prognostic[ind.w[i,],]

      Ana1=MSpecificCoxPh(Survival = Survival[ind.s[i,]],
                           Micro.mat = ReduMicro.mat[,ind.mi[i,]],
                           Censor = Censor[ind.s[i,]],
                           Reduce = Reduce,
                           Select = Select,
                           Prognostic = data.frame(perPrognostic),
                           Mean = TRUE,
                           Quantile = Quantile)

      Temp=Majorityvotes(Ana1,
                          Prognostic = Prognostic[ind.w[i,],],
                          Survival = Survival[ind.s[i,]],
                          Censor = Censor[ind.s[i,]],
                          J = 1)

      if (!is.null(Prognostic)) HRlowPerm[i,]=summary(Temp$Model.result)[[8]][1,][-2]
      if ( is.null(Prognostic)) HRlowPerm[i,]=summary(Temp$Model.result)[[8]][-2]
    }

    Ana2=MSpecificCoxPh( Survival,
                          Micro.mat=ReduMicro.mat,
                          Censor,
                          Reduce = Reduce,
                          Select = Select,
                          Prognostic = Prognostic,
                          Quantile = Quantile)

    TempObs=Majorityvotes(Ana2,
                           Prognostic,
                           Survival,
                           Censor,
                           J=1)

    if (!is.null(Prognostic)) HRlowObs=(summary(TempObs$Model.result)[[8]])[1,][-2]
    if ( is.null(Prognostic)) HRlowObs=(summary(TempObs$Model.result)[[8]])[-2]
  }



  if (Validation=="L1based") {
    for (i in 1:nperm) {
      message('Permutation loop ',i)
      Temp=NA
      if (!is.null(Prognostic)) perPrognostic=Prognostic[ind.w[i,],]

      Temp=Lasoelascox(Survival=Survival[ind.s[i,]],
                           Censor=Censor[ind.s[i,]],
                           Micro.mat[,ind.mi[i,]],
                           Prognostic = data.frame(perPrognostic),
                           Plots = FALSE,
                           Mean = TRUE,
                           Quantile = Quantile,
                           Standardize = TRUE,
                           Alpha=1,
                           Fold = 4,
                           nlambda = 100)

      if ((!is.na(Temp))[1]) {
        if (!is.null(Prognostic)) HRlowPerm[i,]=summary(Temp$SurvFit)[[8]][1,][-2]
        if ( is.null(Prognostic)) HRlowPerm[i,]=summary(Temp$SurvFit)[[8]][-2]
      }
      if ((is.na(Temp))[1])  HRlowPerm[i,]=NA
    }
    TempObs=Lasoelascox(Survival,Censor,
                        Micro.mat,
                        Prognostic=Prognostic,
                        Plots = FALSE,
                        Mean = TRUE,
                        Quantile = Quantile,
                        Standardize = TRUE,
                        Alpha = 1,
                        Fold = 4,
                        nlambda = 100)

    if (!is.null(Prognostic)) HRlowObs=summary(TempObs$SurvFit)[[8]][1,][-2]
    if ( is.null(Prognostic)) HRlowObs=summary(TempObs$SurvFit)[[8]][-2]
  }

  return(methods::new("perm",HRobs=HRlowObs,HRperm=HRlowPerm,nperm=nperm,Validation=Validation))

}
