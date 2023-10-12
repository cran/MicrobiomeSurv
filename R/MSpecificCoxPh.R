#'Taxon by taxon Cox proportional analysis
#'
#' The Function fits cox proportional hazard model and does classification for each taxon separately
#'
#' This function fits  taxon by taxon Cox proportional hazard model and perform the classification based on a microbiome risk score which has been estimated using a single taxon.
#' Function is useful for majority vote classification method and taxon by taxon analysis and also for top K taxa.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Micro.mat A large or small microbiome profile matrix. A matrix with microbiome profiles where the number of rows should be equal to the number of taxa and number of columns should be equal to number of subjects.
#' @param Censor A vector of censoring indicator.
#' @param Reduce A boolean parameter indicating if the microbiome profile matrix should be reduced, default is TRUE and larger microbiome profile matrix is reduced by supervised pca approach.
#' @param Select Number of taxa (default is 5) to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Mean The cut off value for the classifier, default is the mean cutoff.
#' @param Quantile If users want to use quantile as cutoff point. They need to specify Mean = FALSE and a quantile that they wish to use. The default is the median cutoff.
#' @param Method Multiplicity adjustment methods.
#' @return A object of class \code{\link[MicrobiomeSurv]{ms}} is returned with the following values
#'   \item{Result}{The cox proportional regression result for each taxon}
#'   \item{HRRG}{The hazard ratio statistics (Hazard-ratio, Lower confidence interval and upper confidence interval) of the riskgroup based on the riskscore and the cut off value for each taxon}
#'   \item{Group}{The classification of the subjects based on each taxon analysis}
#'   \item{Mi.names}{The names of the taxa for the analysis}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}},  \code{\link[MicrobiomeSurv]{EstimateHR}}
#' @examples
#' # Prepare data
#' data(Week3_response)
#' Week3_response = data.frame(Week3_response)
#' surv_fam_shan_w3 =
#' data.frame(cbind(as.numeric(Week3_response$T1Dweek), as.numeric(Week3_response$T1D)))
#' colnames(surv_fam_shan_w3) = c("Survival", "Censor")
#' prog_fam_shan_w3 = data.frame(factor(Week3_response$Treatment_new))
#' colnames(prog_fam_shan_w3) = c("Treatment")
#' data(fam_shan_trim_w3)
#' names_fam_shan_trim_w3 =
#' c("Unknown", "Lachnospiraceae", "S24.7", "Lactobacillaceae", "Enterobacteriaceae", "Rikenellaceae")
#' fam_shan_trim_w3 = data.matrix(fam_shan_trim_w3[ ,2:82])
#' rownames(fam_shan_trim_w3) = names_fam_shan_trim_w3
#' # Using the function
#' Cox_taxon_fam_shan_w3 = MSpecificCoxPh(Survival = surv_fam_shan_w3$Survival,
#'                                       Micro.mat = fam_shan_trim_w3,
#'                                       Censor = surv_fam_shan_w3$Censor,
#'                                       Reduce=FALSE,
#'                                       Select=5,
#'                                       Prognostic = prog_fam_shan_w3,
#'                                       Mean = TRUE,
#'                                       Method = "BH")
#'
#' # Results
#' show(Cox_taxon_fam_shan_w3)
#' summary(Cox_taxon_fam_shan_w3)

#' @import stats
#' @import superpc
#' @import survival
#' @import lmtest
#' @import methods
#' @export MSpecificCoxPh

MSpecificCoxPh=function(Survival,
                         Micro.mat,
                         Censor,
                         Reduce=FALSE,
                         Select=5,
                         Prognostic=NULL,
                         Mean = TRUE,
                         Quantile = 0.5,
                         Method = "BH"){

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
  mi.name=rownames(ReduMicro.mat)

  HRp = matrix(0, n.mi, 4)
  gr = matrix(0, n.mi, n.obs)
  res =  vector("list", n.mi)

  for (i in 1:n.mi){  #---------------------------  STRAT FOR LOOP ------------------------

    taxoni = ReduMicro.mat[i,]

    if (is.null(Prognostic)) {

      cdata = data.frame(Survival, Censor, taxoni)
      m0 = survival::coxph(survival::Surv(Survival, Censor==1) ~ taxoni, data=cdata)
    }
    if (!is.null(Prognostic)) {
      if (is.data.frame(Prognostic)) {
        nPrgFac=ncol(Prognostic)
        cdata = data.frame(Survival,Censor,taxoni,Prognostic)
        NameProg=colnames(Prognostic)
        eval(parse(text=paste( "m0 =survival::coxph(survival::Surv(Survival, Censor==1) ~ taxoni",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
      } else {

        stop(" Argument 'Prognostic' is NOT a data frame ")
      }

    }
    #risk Score
    beta1 = summary(m0)[[7]][c("taxoni"), 1]
    p1 = beta1*taxoni

    Temptaxoni = EstimateHR(Risk.Scores = p1, Data.Survival =cdata, Prognostic=Prognostic,
                             Plots = FALSE, Mean = TRUE, Quantile = Quantile )
    res[[i]]= Temptaxoni$SurvResult
    HRp[i, c(1, 2, 3)]= summary(Temptaxoni$SurvResult)[[8]][1,c(1,3,4)]
    HRp[i, 4] = summary(Temptaxoni$SurvResult)[[7]][1,5]
    gr[i,]= Temptaxoni$Riskgroup


  }#---------------------------  END OF  FOR LOOP ------------------------

  colnames(HRp) = c("HR","LowerCI","UpperCI", "p-value")
  Mi.names = rownames(as.matrix(ReduMicro.mat))
  rownames(gr) = Mi.names
  colnames(gr) = paste0("Subject", 1: ncol(ReduMicro.mat))

  return(methods::new("ms", Result=res, HRRG=HRp, Group=gr, Mi.names = Mi.names))
}
