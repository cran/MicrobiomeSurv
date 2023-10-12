#'Cross validation for majority votes
#'
#' This function does cross validation for the Majority votes based classification which is a cross validated approach to \code{\link[MicrobiomeSurv]{Majorityvotes}}.
#' @param Survival A vector of survival time with length equals to number of subjects.
#' @param Censor A vector of censoring indicator.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Micro.mat A large or small microbiome profile matrix. A matrix with microbiome profiles where the number of rows should be equal to the number of taxa and number of columns should be equal to number of patients.
#' @param Reduce A boolean parameter indicating if the microbiome profile matrix should be reduced, default is TRUE and larger microbiome profile matrix is reduced by supervised pca approach.
#' @param Select Number of taxa (default is 5) to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE.
#' @param Fold Number of times in which the dataset is divided. Default is 3 which implies dataset will be divided into three groups and 2/3 of the dataset will be the train datset and 1/3 will be to train the results.
#' @param Ncv The Number of cross validation loop. Default is 100.
#' @param Mean The cut off value for the classifier, default is the mean cutoff.
#' @param Quantile If users want to use quantile as cutoff point. They need to specify Mean = FALSE and a quantile that they wish to use. The default is the median cutoff.
#' @return A object of class \code{\link[MicrobiomeSurv]{cvmv}} is returned with the following values
#'   \item{HRTrain}{A matrix of survival information for the training dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.}
#'   \item{HRTest}{A matrix of survival information for the test dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.}
#'   \item{Ncv}{The number of cross validation used.}
#'  \item{Micro.mat}{The microbiome data matrix that was used for the analysis either same as Micro.mat or a reduced version.}
#'   \item{Progfact}{The names of prognostic factors used.}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{Majorityvotes}}
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
#' CVMajority_fam_shan_w3 = CVMajorityvotes(Survival = surv_fam_shan_w3$Survival,
#'                                          Micro.mat = fam_shan_trim_w3,
#'                                          Censor = surv_fam_shan_w3$Censor,
#'                                          Reduce=TRUE,
#'                                          Select=5,
#'                                          Mean = TRUE,
#'                                          Prognostic = prog_fam_shan_w3,
#'                                          Fold=3,
#'                                          Ncv=10)
#'
#' # Get the class of the object
#' class(CVMajority_fam_shan_w3)     # An "cvmv" Class
#'
#' # Method that can be used for the result
#' show(CVMajority_fam_shan_w3)
#' summary(CVMajority_fam_shan_w3)
#' plot(CVMajority_fam_shan_w3)

#' @import survival
#' @import superpc
#' @import stats
#' @import lmtest
#' @import methods
#' @importFrom coef density median p.adjust princomp qnorm quantile
#' @export CVMajorityvotes

CVMajorityvotes = function(Survival,
                          Censor,
                          Prognostic=NULL,
                          Micro.mat,
                          Reduce=TRUE,
                          Select=5,
                          Fold=3,
                          Ncv=100,
                          Mean = TRUE,
                          Quantile = 0.5){

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
  n.train=(n.obs-floor(n.obs/Fold))
  n.test=floor(n.obs/Fold)
  ind.train =matrix(0,Ncv,n.train)
  ind.test  =matrix(0,Ncv,n.test)
  res =  vector("list", n.mi)


  HRp.train = matrix(0,Ncv,3)  # Training
  HRp.test =  matrix(0,Ncv,3)     # Testing

  #set.seed(123)
  pIndex = c(1:n.obs)
  res =res1=res2=res3= vector("list", Ncv)

  for (j in 1:Ncv){
    message('Cross validation loop ',j)
    p1=NA
    p2=NA
    gr.train = matrix(0, n.mi,ncol(ind.train))
    gr.test  = matrix(0,n.mi,ncol(ind.test))

    ind.train[j,] =sort(sample(pIndex,n.train,replace=F) )
    ind.test[j,] =c(1:n.obs)[-c(intersect(ind.train[j,] ,c(1:n.obs)))]

    #---------------------------------------  Training  Set -------------------------------------------

    for (i in 1:n.mi){
      taxoni = ReduMicro.mat[i,ind.train[j,]]

      if (is.null(Prognostic)) {

        cdata = data.frame(Survival=Survival[ind.train[j,]], Censor=Censor[ind.train[j,]], taxoni)
        m0 = survival::coxph(survival::Surv(Survival, Censor==1) ~ taxoni, data=cdata)
      }
      if (!is.null(Prognostic)){
        if (is.data.frame(Prognostic)) {
          nPrgFac=ncol(Prognostic)
          NameProg=colnames(Prognostic)
          cdata = data.frame(Survival[ind.train[j,]], Censor[ind.train[j,]], taxoni, Prognostic[ind.train[j,],])
          colnames(cdata) = c("Survival", "Censor", "taxoni", NameProg)
          eval(parse(text=paste( "m0 =survival::coxph(survival::Surv(Survival, Censor==1) ~ taxoni",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
        } else {

          stop(" Argument 'Prognostic' is NOT a data frame ")
        }

      }
      #risk Score
      beta1=summary(m0)[[7]][c("taxoni"),1]
      p1 = beta1*taxoni
      Prognostic.train=as.data.frame(Prognostic[ind.train[j,],])
      colnames(Prognostic.train) = NameProg
      Temptaxoni =EstimateHR(p1,Data.Survival=cdata, Prognostic=Prognostic.train, Plots = FALSE, Mean = TRUE, Quantile = Quantile)
      gr.train[i,]=Temptaxoni$Riskgroup


      #---------------------------------------  Testing  Set -------------------------------------------
      taxonit = ReduMicro.mat[i,ind.test[j,]]
      if (!is.null(Prognostic)) {
        if (is.data.frame(Prognostic)) {
          nPrgFac=ncol(Prognostic)
          NameProg=colnames(Prognostic)
          cdata = data.frame(Survival[ind.test[j,]],Censor[ind.test[j,]],taxonit, Prognostic[ind.test[j,],])
          colnames(cdata) = c("Survival", "Censor", "taxonit", NameProg)
          eval(parse(text=paste( "m0 =survival::coxph(survival::Surv(Survival, Censor==1) ~ taxonit",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
        } else {
    stop(" Argument 'Prognostic' is NOT a data frame ")
        }

      } else{
        cdata = data.frame(Survival=Survival[ind.test[j,]],Censor=Censor[ind.test[j,]],taxonit)
        eval(parse(text=paste("m0 =survival::coxph(survival::Surv(Survival, Censor==1) ~ taxonit",",data=cdata)" ,sep="")))
      }


      beta1=summary(m0)[[7]][c("taxonit"),1]
      p2 = beta1*taxonit
      Prognostic.test=as.data.frame(Prognostic[ind.test[j,],])
      colnames(Prognostic.test) = NameProg
      Temptaxonit =EstimateHR(p2,Data.Survival=cdata, Prognostic=Prognostic.test, Plots = FALSE, Mean = TRUE, Quantile = Quantile)
      gr.test[i,]=Temptaxonit$Riskgroup

    } # END OF LOOP over taxa---------------------------------------------------------------------

    # ------------ count majority votes for jth Cross validation and estimate HR --------------

    ggr.train = per.R = per.NR = NULL
    for (k in 1:ncol(ind.train)){
      per.R[k]=sum(gr.train[,k]=="Low risk")
      ggr.train[k]=ifelse((n.mi-per.R[k])>per.R[k],"High risk","Low risk")
    }



    #-------------------- HR estimation for Training  ----------------------
    GS=as.factor(ggr.train)
    if (is.null(Prognostic)) {

      cdata = data.frame(Survival=Survival[ind.train[j,]],Censor=Censor[ind.train[j,]],GS)
      mTrain = survival::coxph(survival::Surv(Survival, Censor==1) ~ GS,data=cdata)
    }
    if (!is.null(Prognostic)) {
      if (is.data.frame(Prognostic)) {
        nPrgFac=ncol(Prognostic)
        NameProg=colnames(Prognostic)
        cdata = data.frame(Survival[ind.train[j,]],Censor[ind.train[j,]],GS, Prognostic[ind.train[j,],])
        colnames(cdata) = c("Survival", "Censor", "GS", NameProg)
        eval(parse(text=paste( "mTrain =survival::coxph(survival::Surv(Survival, Censor==1) ~ GS",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
      } else {
        stop(" Argument 'Prognostic' is NOT a data frame ")
      }

    }
    HRp.train[j,]=(summary(mTrain)[[8]][1,])[-2]

    #-------------------- HR estimation for Testing  ----------------------
    ggr.test = per.R=per.NR = NULL
    for (k in 1:ncol(ind.test)){
      per.R[k]=sum(gr.test[,k]=="Low risk")
      ggr.test[k]=ifelse((n.mi-per.R[k])>per.R[k],"High risk","Low risk")
    }

    GS=as.factor(ggr.test)
    if (is.null(Prognostic)) {

      cdata = data.frame(Survival=Survival[ind.test[j,]],Censor=Censor[ind.test[j,]],GS)
      mTest = survival::coxph(survival::Surv(Survival, Censor==1) ~ GS,data=cdata)
    }
    if (!is.null(Prognostic)) {
      if (is.data.frame(Prognostic)) {
        nPrgFac=ncol(Prognostic)
        NameProg=colnames(Prognostic)
        cdata = data.frame(Survival[ind.test[j,]], Censor[ind.test[j,]], GS, Prognostic[ind.test[j,],])
        colnames(cdata) = c("Survival", "Censor", "GS", NameProg)
        eval(parse(text=paste( "mTest =survival::coxph(survival::Surv(Survival, Censor==1) ~ GS",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
      } else {
        stop(" Argument 'Prognostic' is NOT a data frame ")
      }

    }
    HRp.test[j,]=(summary(mTest)[[8]][1,])[-2]

  }#---------------------------  END OF  FOR LOOP over Cross Validations ------------------------

  pFactors=NA
  if (!is.null(Prognostic)) pFactors =colnames(Prognostic)

  return(methods::new("cvmv",HRTrain=HRp.train,HRTest=HRp.test,Ncv=Ncv,Micro.mat=ReduMicro.mat, Progfact=pFactors))
  }
