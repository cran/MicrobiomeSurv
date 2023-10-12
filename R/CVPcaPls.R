#'Cross Validations for PCA and PLS based methods
#'
#' This function does cross validation for the analysis performs by \code{\link[MicrobiomeSurv]{SurvPcaClass}}
#' and \code{\link[MicrobiomeSurv]{SurvPlsClass}} functions where the dimension reduction methods can either be PCA and PLS.
#'
#' This function does cross validation for the analysis using two reduction method. The reduction method can be PCA or PLS.
#' If it is PCA then the \code{\link[MicrobiomeSurv]{SurvPcaClass}} is internally used for the cross validation
#' and \code{\link[MicrobiomeSurv]{SurvPlsClass}} otherwise.
#' @param Fold Number of times in which the dataset is divided. Default is 3 which implies dataset will be divided into three groups and 2/3 of the dataset will be the train datset and 1/3 will be to test the results.
#' @param Survival A vector of survival time with length equals to number of subjects.
#' @param Micro.mat A large or small microbiome profile matrix. A matrix with microbiome profiles where the number of rows should be equal to the number of taxa and number of columns should be equal to number of patients.
#' @param Censor A vector of censoring indicator.
#' @param Reduce A boolean parameter indicating if the microbiome profile matrix should be reduced, default is TRUE and larger microbiome profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of taxa (default is 5) to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Ncv The Number of cross validation loop. Default is 100.
#' @param DR The dimension reduction method. It can be either "PCA" for Principle components analysis or "PLS" for Partial least squares.
#' @return A object of class \code{\link[MicrobiomeSurv]{cvpp}} is returned with the following values
#'   \item{Result}{A dataframe containg the estimated Hazard ratio of the test dataset and the training dataset.}
#'   \item{Ncv}{The number of cross validation performed.}
#'   \item{Method}{The dimesion reduction method used.}
#'   \item{CVtrain}{The training dataset indices matrix used for the cross validation.}
#'   \item{CVtest}{The test dataset indices matrix used for the cross validation.}
#'   \item{Select}{The number of taxa used for the dimesion reduction method used.}
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{SurvPlsClass}},
#' \code{\link[MicrobiomeSurv]{SurvPcaClass}}
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
#' CVPls_fam_shan_w3 = CVPcaPls(Fold = 3,
#'                             Survival = surv_fam_shan_w3$Survival,
#'                             Micro.mat = fam_shan_trim_w3,
#'                             Censor = surv_fam_shan_w3$Censor,
#'                             Reduce=TRUE,
#'                             Select=5,
#'                             Prognostic = prog_fam_shan_w3,
#'                             Ncv=10,
#'                             DR = "PLS")
#'
#' # Get the class of the object
#' class(CVPls_fam_shan_w3)     # An "cvpp" Class
#'
#' # Method that can be used for the result
#' show(CVPls_fam_shan_w3)
#' summary(CVPls_fam_shan_w3)
#' plot(CVPls_fam_shan_w3)

#' @import stats
#' @import superpc
#' @import lmtest
#' @import methods
#' @import survival
#' @importFrom coef density median p.adjust princomp qnorm quantile

#' @export CVPcaPls

CVPcaPls=function(Fold = 3,
                   Survival,
                   Micro.mat,
                   Censor,
                   Reduce=TRUE,
                   Select=15,
                   Prognostic=NULL,
                   Ncv=5,
                   DR ="PCA"){

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


  n.mi=nrow(Micro.mat)
  sen=Censor               # Censoring indicator
  surv = Survival

  if (is.data.frame(Prognostic)) {data1=data.frame(sen,surv,Prognostic)
  } else{
    data1=data.frame(sen,surv)
  }


  HRp.train = HRn.train = matrix(0,Ncv,4)  # Training
  HRp.test = HRn.test= matrix(0,Ncv,4)     # Testing

  n.train=(n.obs-floor(n.obs/Fold))
  n.test=floor(n.obs/Fold)
  cv.train =matrix(0,Ncv,n.train)
  cv.test  =matrix(0,Ncv,n.test)

  #----------------------------------- function for PCA ----------------------------------------------
  f.pca = function (x){
    ca <- match.call()
    if (ncol(x) > nrow(x)){
      u = stats::princomp(t(x))
      u$call = ca
      return(u)
    }

    mu = rowMeans(x)
    xb <- x - mu
    xb.svd <- svd(xb)
    pc <- t(xb) %*% xb.svd$u
    dimnames(pc)[[2]] <- paste("PC", 1:ncol(pc), sep = "")
    loading <- xb.svd$u
    dimnames(loading) <- list(paste("V", 1:nrow(loading), sep = ""),
                              paste("Comp.", 1:ncol(loading), sep = ""))
    class(loading) <- "loadings"
    sd = xb.svd$d/sqrt(ncol(x))
    names(sd) <- paste("Comp.", 1:length(sd), sep = "")
    pc <- list(sdev = sd, loadings = loading, center = mu,
               scale = rep(1, length(mu)), n.obs = ncol(x), scores = pc, call = ca)
    class(pc) <- "princomp"
    return(pc)
  }

  #----------------------------------- Intermediate PCA ----------------------------------------------

  IntermediatePCA=function(Micro.mat,
                           Prognostic,
                           Survival,
                           Censor,
                           index){
    if (is.matrix(Micro.mat)) {
      pc1 = f.pca(as.matrix(Micro.mat[ ,index]))[[6]][,1]
    } else {
      pc1=Micro.mat[index]
    }

    if (is.null(Prognostic)) {

      cdata = data.frame(Survival=Survival[index], Censor=Censor[index], pc1)
      m0 = survival::coxph(survival::Surv(Survival, Censor==1) ~ pc1, data=cdata)
    }

    if (!is.null(Prognostic)) {
      if (is.data.frame(Prognostic)) {
        nPrgFac=ncol(Prognostic)
        NameProg=colnames(Prognostic)
        cdata = data.frame(Survival[index], Censor[index], pc1, Prognostic[index,])
        colnames(cdata) = c("Survival", "Censor", "pc1", NameProg)
        eval(parse(text=paste( "m0 =survival::coxph(survival::Surv(Survival, Censor==1) ~ pc1",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
      } else {
        stop(" Argument 'Prognostic' is NOT a data frame ")
      }

    }

    return(list(m0=m0,pc1=pc1,cdata=cdata))
  }


  #---------------------------------------------------------------------------------------------------
  #----------------------------------- Intermediate PLS ----------------------------------------------


  IntermediatePLS=function(Micro.mat,
                           Prognostic,
                           Survival,
                           Censor,
                           index){
    if (is.matrix(Micro.mat)) {

      DataPLS=data.frame(1:length(index))
      DataPLS$g=as.matrix(t(Micro.mat[ ,index]))
      colnames(DataPLS)[1]=c("Survival")
      DataPLS[,1]=Survival[index]
      plsr.1 = pls::plsr(Survival ~ g, method="simpls", ncomp = 2, scale =TRUE,data = DataPLS, validation =  "CV")
      pc1=pls::scores(plsr.1)[,1] # extract the first com
    } else {
      pc1=Micro.mat[index]
    }

    if (is.null(Prognostic)) {
      cdata = data.frame(Survival=Survival[index], Censor=Censor[index], pc1)
      m0 = survival::coxph(survival::Surv(Survival, Censor==1) ~ pc1, data=cdata)
    }
    if (!is.null(Prognostic)) {
      if (is.data.frame(Prognostic)) {
        nPrgFac=ncol(Prognostic)
        NameProg=colnames(Prognostic)
        cdata = data.frame(Survival[index], Censor[index], pc1, Prognostic[index,])
        colnames(cdata) = c("Survival", "Censor", "pc1", NameProg)
        eval(parse(text=paste( "m0 =survival::coxph(survival::Surv(Survival, Censor==1) ~ pc1",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
      } else {
        stop(" Argument 'Prognostic' is NOT a data frame ")
      }
    }
    return(list(m0=m0,pc1=pc1,cdata=cdata))
  }


  #set.seed(123)
  pIndex = c(1:n.obs)
  res =res1=res2=res3= vector("list", Ncv)
  #------------------------  Begin FOR LOOP :1  --------------------------- i=1
  for (i in 1:Ncv){
    message('Cross validation loop ',i)
    p1=NA
    p2=NA

    cv.train[i,] =sort(sample(pIndex,n.train,replace=F) )
    cv.test[i,] =c(1:n.obs)[-c(intersect(cv.train[i, ], c(1:n.obs)))]

    if (DR=="PCA") { #-------------------------------------------------------------------------------------------

  Temp1=IntermediatePCA(Micro.mat,Prognostic,Survival,Censor,as.vector(cv.train[i, ]))
      pc1  = Temp1$pc1
      cdata =Temp1$cdata
      m2 = Temp1$m0

  Temp2=IntermediatePCA(Micro.mat,Prognostic,Survival,Censor,cv.test[i, ])
      pc1test = Temp2$pc1
      ctestdata = Temp2$cdata


      # classification method A1
      TrtandPC1=summary(m2)[[7]][c("pc1"),1]
      p1 =  TrtandPC1[1]*pc1
      p2 =  TrtandPC1[1]*pc1test
    }#-------------------------------------------------------------------------------------------


    if (DR=="PLS") { #-------------------------------------------------------------------------------------------

      Temp3=IntermediatePLS(Micro.mat,Prognostic,Survival,Censor,cv.train[i,])
      pls.comp1  = Temp3$pc1
      cdata =Temp3$cdata
      m3 = Temp3$m0

      Temp4=IntermediatePLS(Micro.mat,Prognostic,Survival,Censor,cv.test[i,])
      pls.comp1.test = Temp4$pc1
      ctestdata = Temp4$cdata

      # classification method A1
      TrtandPLSc1=summary(m3)[[7]][c("pc1"),1]
      p1 = TrtandPLSc1[1]*pls.comp1
      p2 = TrtandPLSc1[1]*pls.comp1.test

    }#-------------------------------------------------------------------------------------------


    #######################
    ## training set ###########

    TrainSet=EstimateHR(p1,Data.Survival=cdata, Prognostic=data.frame(Prognostic[cv.train[i,],]), Plots = FALSE, Mean = TRUE, Quantile = 0.5)
  HRp.train[i,] = summary(TrainSet$SurvResult)[[8]][1,]

    #######################
    ## test set ###########

    mp1 = stats::median(p1)
    TestSet=EstimateHR(p2,Data.Survival=ctestdata, Prognostic=data.frame(Prognostic[cv.test[i,],]), Plots = FALSE, Mean = TRUE, Quantile = 0.5)

    HRp.test[i,] = summary( TestSet$SurvResult)[[8]][1,]


  }
  #------------------------  END of FOR LOOP :1


  Results=data.frame(Training=HRp.train[,1],Test=as.numeric(HRp.test[,1]))

  return(methods::new("cvpp",Results=Results, Ncv=Ncv, Method=DR, CVtrain=cv.train, CVtest=cv.test, Select=n.mi))
}
