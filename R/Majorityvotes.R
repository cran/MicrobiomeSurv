#'Classifiction for Majority Votes
#'
#' The Function fits cox proportional hazard model and does classification based on the majority votes.
#'
#' The Function fits cox proportional hazard model and does classification based on the majority votes
#' while estimating the Hazard ratio of the low risk group.
#' The function firstly count the number of low risk classification for each subject based on the taxon specific analysis
#' which determines the majority votes.
#' In addition, function visualizes the taxon specific calssification for the subjects. 25 subjects is taken for visualization purpose.
#' @param Result An object obtained from the taxon specific analysis (\code{\link[MicrobiomeSurv]{MSpecificCoxPh}}) which is of class "ms"
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Censor A vector of censoring indicator
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param J The jth set of subjects required for the visualization. The default is J=1 which is the first set of subjects. For visualization, J should be less than the number of subjects divided by 25
#' @return A list is returned with the following values
#'   \item{Model.result}{The cox proportional regression result based on the majority vote classification}
#'   \item{N}{The majority vote for each subject}
#'   \item{Classif}{The majority vote classification for each subjects}
#'   \item{Group}{The classification of the subjects based on each taxon analysis}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{MSpecificCoxPh}}, \code{\link[survival]{coxph}},  \code{\link[MicrobiomeSurv]{EstimateHR}}
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

#' # Running the taxon specific function
#' Cox_taxon_fam_shan_w3 = MSpecificCoxPh(Survival = surv_fam_shan_w3$Survival,
#'                                        Micro.mat = fam_shan_trim_w3,
#'                                        Censor = surv_fam_shan_w3$Censor,
#'                                        Reduce=FALSE,
#'                                        Select=5,
#'                                        Prognostic = prog_fam_shan_w3,
#'                                        Mean = TRUE,
#'                                        Method = "BH")
#'
#' # Using the function
#' Majority_fam_shan_w3 = Majorityvotes(Result = Cox_taxon_fam_shan_w3,
#'                                      Prognostic = prog_fam_shan_w3,
#'                                      Survival = surv_fam_shan_w3$Survival,
#'                                      Censor = surv_fam_shan_w3$Censor,
#'                                      J=1)
#'
#' # The survival analysis for majority vote result
#' Majority_fam_shan_w3$Model.result
#'
#' # The majority vote for each subject
#' Majority_fam_shan_w3$N
#'
#' # The majority vote classification for each subject
#' Majority_fam_shan_w3$Classif
#'
#' # The group for each subject based on the taxon specific analysis
#' Majority_fam_shan_w3$Group

#' @import stats
#' @import survival
#' @import graphics
#' @import base
#' @importFrom coef density median p.adjust princomp qnorm quantile
#' @importFrom abline arrows axis barplot box boxplot legend lines par points

#' @export Majorityvotes


Majorityvotes=function(Result, Prognostic, Survival, Censor, J=1){


  if (inherits(Result, "ms") == FALSE) stop("Invalid class object.")
  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")


  ggr = per.R=per.NR = NULL


  Group = Result@Group
  Group = data.frame(Group)

  for (i in 1:length(Survival)){
    per.R[i] = sum(Group[,i]=="Low risk")
    ggr[i] = ifelse((nrow(Group)-per.R[i]) > per.R[i],"High risk","Low risk")
  }

  ggr=as.factor(ggr)

  if (is.null(Prognostic)) {

    cdata = data.frame(Survival,Censor,ggr)
    m0 = survival::coxph(survival::Surv(Survival, Censor==1) ~ ggr,data=cdata)
  }

  if (!is.null(Prognostic)) {
    if (is.data.frame(Prognostic)) {
      nProg=ncol(Prognostic)
      cdata = data.frame(Survival,Censor,ggr,Prognostic)
      NameProg=colnames(Prognostic)
      eval(parse(text=paste( "m0 =survival::coxph(survival::Surv(Survival, Censor==1) ~ ggr",paste("+",NameProg[1:nProg],sep="",collapse =""),",data=cdata)" ,sep="")))
    } else {

      stop(" Argument 'Prognostic' is NOT a data frame ")
    }

  }

  VoteMat=Group
  ng=nrow(VoteMat)
  np=ncol(VoteMat)
  n0 = np
  Jmax=floor(np/n0)
  if (J<=Jmax) {
    slist=(1+(J-1)*ng):(J*ng)
    #graphics::par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(0,0, xlim=c(0,n0),ylim=c(1,ng),xlab="Subject Index",ylab="Taxa",axes=FALSE,type="n", main="Taxon-Specific Classification of subjects",cex.main=0.9)
    for(i in 1:n0){
      graphics::points(rep(i, ng)[t(VoteMat)[i, slist]=="Low risk"],which(t(VoteMat)[i, slist]=="Low risk"),col="blue", pch=15)
      graphics::points(rep(i, ng)[t(VoteMat)[i, slist]=="High risk"],which(t(VoteMat)[i, slist]=="High risk"),col="yellow", pch=15)
    }
    graphics::axis(1, at = 1:n0);graphics::axis(2,at=1:ng,slist);graphics::box()
    graphics::legend("topright",inset=c(-0.37,0), c("Low Risk","High Risk"),
           pch=c(15,15),col=c("blue", "yellow"), cex=0.6)
  } else {stop("J should be less than or equal to (no.of subjects / number of selected subjects)")
  }
  return(list(Model.result=m0, N= per.R, Classif=ggr,Group=VoteMat))

}
