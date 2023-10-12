#'Sequential Increase in Taxa for the PCA or PLS classifier
#'
#' The Function fits cox proportional hazard model and does classification by sequentially increasing the taxa
#' using either PCA or PLS based on the topK taxa specified.
#'
#' This function sequentially increase the number of top K taxa to be used in the PCA or PLS methods in order to obtain the risk score.
#' This function internally calls \code{\link[MicrobiomeSurv]{MSpecificCoxPh}} to rank the taxa based on HR for each taxon.
#' Therefore taxa can be ordered based on increasing order of the HR for low risk group.
#' Thereafter, the function takes few top K (5 is the default) to be used in the sequential analysis.
#' @param TopK 	Top K taxa (5 by default) to be used in the sequential analysis.
#' @param Survival A vector of survival time with length equals to number of subjects.
#' @param Micro.mat A large or small microbiome profile matrix. A matrix with microbiome profiles where the number of rows should be equal to the number of taxa and number of columns should be equal to number of patients.
#' @param Censor A vector of censoring indicator.
#' @param Reduce A boolean parameter indicating if the microbiome profile matrix should be reduced, default is TRUE and larger microbiome profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of taxa to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Plot A boolean parameter indicating if Plot should be shown. Default is FALSE.
#' @param DM Dimension reduction method which can either be PLS or PCA.
#' @param ...	 Additinal arguments for plotting and only valid  if Plot=TRUE
#' @return A list containing a data frame with estimated HR along with 95\% CI at each TopK value for the sequential analysis.
#'   \item{Result}{The hazard ratio statistics (HR, Lower confidence interval and upper confidence interval) of the lower riskgroup based for each sequential metabolite analysis}
#'   \item{TopKplot}{A graphical representation of the Result containing the hazard ratio statistics}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}},  \code{\link[MicrobiomeSurv]{EstimateHR}}, \code{\link[MicrobiomeSurv]{MSpecificCoxPh}}, \code{\link[MicrobiomeSurv]{SurvPcaClass}}, \code{\link[MicrobiomeSurv]{SurvPlsClass}}
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
#' SITaxa_fam_shan_w3 = SITaxa(TopK=5,
#'                             Survival = surv_fam_shan_w3$Survival,
#'                             Micro.mat = fam_shan_trim_w3,
#'                             Censor = surv_fam_shan_w3$Censor,
#'                             Reduce=TRUE,
#'                             Select=5,
#'                             Prognostic=prog_fam_shan_w3,
#'                             Plot = TRUE,
#'                             DM="PLS")
#'
#' # For the HR statistics
#' SITaxa_fam_shan_w3$Result
#'
#' # For the graphical output
#' SITaxa_fam_shan_w3$TopKplot

#' @import stats
#' @import superpc
#' @import survival
#' @import lmtest
#' @import methods
#' @import ggplot2
#' @export SITaxa


SITaxa=function(TopK=15,
                 Survival,
                 Micro.mat,
                 Censor,
                 Reduce=TRUE,
                 Select=5,
                 Prognostic=NULL,
                 Plot = FALSE,
                 DM=c("PLS","PCA"),...)
{
  Decrease=FALSE
  DM = match.arg(DM)

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
  object= MSpecificCoxPh(Survival, ReduMicro.mat, Censor, Reduce = FALSE,
                          Select = Select, Prognostic, Mean = TRUE, Quantile = 0.5)

  Names.Ktaxa=object@Mi.names
  index.Top.Ktaxa = order(object@HRRG[ ,1], decreasing = Decrease)
  index.Top.Ktaxa=index.Top.Ktaxa[1:TopK]

  TopSet=Names.Ktaxa[index.Top.Ktaxa]

  Result=matrix(NA,TopK,4)

  for (i in 1:length(TopSet[1:TopK]) ){

    mlist=1:i

    if (DM=="PLS") {

      Temp= SurvPlsClass(Survival, ReduMicro.mat[intersect(rownames(ReduMicro.mat),TopSet[mlist]) , , drop = FALSE],
                          Censor, Reduce = FALSE, Prognostic = Prognostic,
                          Plots = FALSE, Mean = TRUE)
    } else {
      Temp=  SurvPcaClass(Survival, ReduMicro.mat[intersect(rownames(ReduMicro.mat),TopSet[mlist]), , drop = FALSE],
                           Censor, Reduce = FALSE, Prognostic = Prognostic,
                           Plots = FALSE, Mean = TRUE)
    }

    if (is.null(Prognostic)) Result[i,]=c(i,(summary(Temp$SurvFit)[[8]][1,])[-2] )
    if (!is.null(Prognostic)) Result[i,]=c(i,(summary(Temp$SurvFit)[[8]][1,])[-2] )
  }
  HR = LowerCI = UpperCI = NULL
  colnames(Result)=c("Topk","HR","LowerCI","UpperCI")
  Result = data.frame(Result)

  if (Plot) {
    TopKplot = ggplot2::ggplot(data=Result, ggplot2::aes(x=1:nrow(Result)),HR) +
      ggplot2::geom_errorbar(ggplot2::aes(x =1:nrow(Result), ymax = UpperCI, ymin = LowerCI)) +
      ggplot2::ylab(paste("HR Interval Range Based on ",DM,sep="")) +
      ggplot2::xlab("Top K")+ ggplot2::theme_classic() + ggplot2::geom_point(ggplot2::aes(y=Result[,2],x=1:nrow(Result)),colour="red")

    return(list(Result=Result, TopKplot=TopKplot))
  }

  else
    return(Result=Result)

}
