
#'Quantile sensitivity analysis
#'
#' The function performs sensitivity of the cut off quantile for obtaining the risk group obtained
#' under \code{\link[MicrobiomeSurv]{SurvPlsClass}}, \code{\link[MicrobiomeSurv]{SurvPcaClass}}
#' or \code{\link[MicrobiomeSurv]{Lasoelascox}} requires for the survival analysis and classification.
#'
#' This function investigates how each analysis differs from the general median cutoff of 0.5,
#' therefore to see the sensitive nature of the survival result different quantiles ranging from 10th percentile to 90th percentiles were used.
#' The sensitive nature of the quantile is investigated under \code{\link[MicrobiomeSurv]{SurvPlsClass}}, \code{\link[MicrobiomeSurv]{SurvPcaClass}}
#' or \code{\link[MicrobiomeSurv]{Lasoelascox}} while relate to the 3 different Dimension method to select from.
#' @param Survival A vector of survival time with length equals to number of subjects.
#' @param Micro.mat A large or small microbiome profile matrix. A matrix with microbiome profiles where the number of rows should be equal to the number of taxa and number of columns should be equal to number of patients.
#' @param Censor A vector of censoring indicator.
#' @param Reduce A boolean parameter indicating if the microbiome profile matrix should be reduced, default is TRUE and larger microbiome profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of taxa (default is 5) to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Plots A boolean parameter indicating if the graphical represenataion of the analysis should be shown. Default is FALSE and it is only valid for the PCA or PLS dimension method.
#' @param DM The dimension method to be used. PCA implies using the \code{\link[MicrobiomeSurv]{SurvPcaClass}}, PLS uses \code{\link[MicrobiomeSurv]{SurvPcaClass}} while SM uses the \code{\link[MicrobiomeSurv]{Lasoelascox}} which ruses the shrinkage method techniques such as lasso and elastic net.
#' @param Alpha The mixing parameter for glmnet (see \code{\link[glmnet]{glmnet}}). The range is 0<= Alpha <= 1. The Default is 1.
#' @return A Dataframe is returned depending on weather a data reduction method should be used or not. The dataframe contains the HR of the low risk group for each percentile.
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}},\code{\link[MicrobiomeSurv]{EstimateHR}},
#' \code{\link[MicrobiomeSurv]{SurvPcaClass}},
#'  \code{\link[MicrobiomeSurv]{SurvPlsClass}},\code{\link[MicrobiomeSurv]{Lasoelascox}}
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

#' # Using the PCA method
#' QuantileAnalysis_PCA_fam_shan_w3 = QuantileAnalysis(Survival = surv_fam_shan_w3$Survival,
#'                                                     Micro.mat = fam_shan_trim_w3,
#'                                                     Censor = surv_fam_shan_w3$Censor,
#'                                                     Reduce=TRUE,
#'                                                     Select= 5,
#'                                                     Prognostic=prog_fam_shan_w3,
#'                                                     Plots = TRUE,
#'                                                     DM="PCA",
#'                                                     Alpha =1)
#'

#' @import stats
#' @import superpc
#' @import survival
#' @import lmtest
#' @import methods
#' @import gplots
#' @import graphics
#' @export QuantileAnalysis

QuantileAnalysis=function(Survival,
                           Micro.mat,
                           Censor,
                           Reduce=TRUE,
                           Select= 5,
                           Prognostic=NULL,
                           Plots = FALSE,
                           DM=c("PLS","PCA","SM"),
                           Alpha =1){

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
  cutpoint = seq(0.10, 0.9, 0.05)

  Grinding = matrix(NA,ncol=5,nrow=length(cutpoint))
  Grinding2 = list()

  for (i in 1:length(cutpoint)) {
    message('Running analysis for Quantile = ',cutpoint[i])
    if (DM=="PLS") {
      Temp= SurvPlsClass(Survival, Micro.mat, Censor, Reduce = FALSE, Select = Select, Prognostic, Plots = FALSE, Mean = FALSE, Quantile = cutpoint[i])
    } else if (DM=="PCA") {
      Temp=  SurvPcaClass(Survival, Micro.mat, Censor, Reduce = FALSE, Select = Select, Prognostic, Plots = FALSE, Mean = FALSE, Quantile = cutpoint[i])
    } else{
      Temp = Lasoelascox(Survival, Censor, Micro.mat, Prognostic, Mean = FALSE, Quantile = cutpoint[i],
                         Plots = FALSE, Standardize = TRUE, Alpha = Alpha)
    }

      Grinding[i,]=c(summary(Temp$SurvFit)[[8]][1,],cutpoint[i])
      colnames(Grinding)=c("EstimatedHR","IHR","LowerCI","UpperCI","Quantile")
  }

  Grinding[is.infinite(Grinding[,3]), 3] = Grinding[1,1]
  if (Plots) {
      gplots::plotCI(x=Grinding[,1],xaxt="n", ui=Grinding[,4],li=Grinding[,3],  col="black", barcol="red", lwd=2,ylab="HR",xlab="CutOFF (Quantile)",main=paste("Dimension reduction using ",DM,sep=""))
    graphics::axis(side=1,at=1:length(cutpoint),seq(0.10, 0.9, 0.05)*100,cex=0.7)
  }

  data1=Grinding[,-2]
  colnames(data1)=c("EstimatedHR","LowerCI","UpperCI","Quantile")
    return(data1)
}
