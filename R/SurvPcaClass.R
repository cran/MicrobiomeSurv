#'Survival PCA and Classification for microbiome data
#'
#' The function performs principal component analysis (PCA) on microbiome matrix and fit Cox proportional hazard model with covariates using also the first PCA as covariates.
#'
#' This function can handle single and multiple microbiome. For larger microbiome matrix,
#' this function will reduce largermicrobiome matrix to smaller version using supervised pca approach and this is by default done and can be control by using the argument Reduce.
#' Other prognostic factors can be included to the model.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Micro.mat A large or small microbiome profile matrix. A matrix with microbiome profiles where the number of rows should be equal to the number of microbiome and number of columns should be equal to number of patients.
#' @param Censor A vector of censoring indicator
#' @param Reduce A boolean paramier indicating if the microbiome profile matrix should be reduced, default is TRUE and larger microbiome profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of microbiome (default is 15) to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Plots A boolean paramier indicating if the plots should be shown. Default is FALSE
#' @param Mean The cut off value for the classifier, default is the mean cutoff
#' @param Quantile If user want to use quantile as cutoff point. They need to specify Mean = FALSE and a quantile that they want to use. The default is the median cutoff

#' @return A object of class SurvPca is returned with the following values
#'   \item{Survfit}{The cox proportional regression result using the first PCA}
#'   \item{Riskscores}{A vector of risk scores which is equal to the number of patents.}
#'   \item{Riskgroup}{The classification of the subjects based on the PCA into low or high risk group}
#'   \item{pc1}{The First PCA scores based on either the reduced microbiome matrix or the full matrix}
#' \item{KMplot}{The Kaplan-Meier survival plot of the riskgroup}
#'  \item{SurvBPlot}{The distribution of the survival in the riskgroup}
#'  \item{Riskpca}{The plot of Risk scores vs first PCA}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}},
#' \code{\link[MicrobiomeSurv]{EstimateHR}}, \code{\link[stats]{princomp}},
#'  \code{\link[MicrobiomeSurv]{SurvPlsClass}}
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
#' SPCA_fam_shan_w3 = SurvPcaClass(Survival = surv_fam_shan_w3$Survival,
#'                                 Micro.mat = fam_shan_trim_w3,
#'                                 Censor = surv_fam_shan_w3$Censor,
#'                                 Reduce=TRUE,
#'                                 Select=5,
#'                                 Prognostic = prog_fam_shan_w3,
#'                                 Plots = TRUE,
#'                                 Mean = TRUE)
#'
#' # Getting the survival regression output
#' SPCA_fam_shan_w3$SurvFit
#'
#' # Getting the riskscores
#' SPCA_fam_shan_w3$Riskscores
#'
#' # Getting the riskgroup
#' SPCA_fam_shan_w3$Riskgroup
#'
#' # Obtaining the first principal component scores
#' SPCA_fam_shan_w3$pc1

#' @import stats
#' @import superpc
#' @import survival
#' @import lmtest
#' @import ggplot2
#' @export SurvPcaClass



SurvPcaClass = function(Survival,
                 Micro.mat,
                 Censor,
                 Reduce=TRUE,
                 Select=5,
                 Prognostic=NULL,
                 Plots = FALSE,
                 Mean = TRUE,
                 Quantile = 0.5){

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

  f.pca = function (x){
    ca = match.call()
    if (ncol(x) > nrow(x)){
      u = stats::princomp(t(x))
      u$call = ca
      return(u)
    }

    mu = rowMeans(x)
    xb = x - mu
    xb.svd = svd(xb)
    pc = t(xb) %*% xb.svd$u
    dimnames(pc)[[2]] = paste("PC", 1:ncol(pc), sep = "")
    loading = xb.svd$u
    dimnames(loading) = list(paste("V", 1:nrow(loading), sep = ""),
                             paste("Comp.", 1:ncol(loading), sep = ""))
    class(loading) = "loadings"
    sd = xb.svd$d/sqrt(ncol(x))
    names(sd) = paste("Comp.", 1:length(sd), sep = "")
    pc = list(sdev = sd, loadings = loading, center = mu,
              scale = rep(1, length(mu)), n.obs = ncol(x), scores = pc, call = ca)
    class(pc) = "princomp"
    return(pc)
  }

  if (is.matrix(ReduMicro.mat)){
    pc1 = f.pca(as.matrix(ReduMicro.mat))[[6]][,1]
  } else {
    pc1=ReduMicro.mat
  }

  if (is.null(Prognostic)) {

    cdata = data.frame(Survival,Censor,pc1)
    m0 = survival::coxph(survival::Surv(Survival, Censor==1) ~ pc1,data=cdata)
  }
  if (!is.null(Prognostic)) {
    if (is.data.frame(Prognostic)) {
      nPrgFac=ncol(Prognostic)
      cdata = data.frame(Survival,Censor,pc1,Prognostic)
      NameProg=colnames(Prognostic)
      eval(parse(text=paste( "m0 =survival::coxph(survival::Surv(Survival, Censor==1) ~ pc1",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
    } else {

      stop("Argument 'Prognostic' is NOT a data frame")
    }

  }
  Riskscores = Riskgroup = NULL

  #risk Score
  TrtandPC1=summary(m0)[[7]][c("pc1"),1]
  p1 = TrtandPC1*pc1

  TempRes= EstimateHR(Risk.Scores = p1, Data.Survival = cdata, Prognostic = Prognostic,
                       Plots = TRUE, Mean = TRUE, Quantile = 0.5)
  gg = data.frame(Riskscores = p1,Riskgroup = TempRes$Riskgroup,pc1 = pc1)
  ab = ggplot2::ggplot(gg, ggplot2::aes(x=Riskscores, y=pc1, shape=Riskgroup, color=Riskgroup)) + ggplot2::geom_point()

  tempp=list(SurvFit=TempRes$SurvResult,Riskscores = p1,
              Riskgroup=TempRes$Riskgroup, pc1=pc1, ReduMicro.mat=ReduMicro.mat)
  class(tempp)="SPCA"

  temp=list(SurvFit=TempRes$SurvResult,Riskscores = p1, Riskgroup=TempRes$Riskgroup, pc1=pc1,
             KMplot = TempRes$KMplot, SurvBPlot = TempRes$SurvBPlot, Riskpca = ab,
             ReduMicro.mat = ReduMicro.mat)
  class(temp)="SPCA"

  if (Plots){
    return(temp)
  } else{
    return(tempp)
  }
}
