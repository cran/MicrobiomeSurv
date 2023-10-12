#' The cvmv Class.
#'
#' Class of object returned by function \code{\link[MicrobiomeSurv]{CVMajorityvotes}}.
#'
#' @name cvmv-class
#' @rdname cvmv-class
#' @exportClass cvmv
#' @param x	 A cvmv class object
#' @param y	 missing
#' @param  object A cvmv class object
#' @param ...	 The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot HRTrain A matrix of survival information for the training dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.
#' @slot HRTest A matrix of survival information for the test dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.
#' @slot Ncv The number of cross validation used.
#' @slot Micro.mat The microbiome data matrix that was used for the analysis either same as Micro.mat or a reduced version.
#' @slot Progfact The names of prognostic factors used.
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{Majorityvotes}}, \code{\link[MicrobiomeSurv]{CVPcaPls}}, \code{\link[MicrobiomeSurv]{SurvPcaClass}}, \code{\link[MicrobiomeSurv]{SurvPlsClass}}
#' @importFrom methods graphics stats setClass setGeneric setMethod setRefClass
#' @importFrom grDevices
#' @importFrom hasArg new
#' @importFrom coef density median na.exclude p.adjust princomp qnorm quantile
#' @importFrom abline arrows axis barplot box boxplot legend lines par points


setClass("cvmv",representation(HRTrain="matrix",HRTest="matrix",Ncv="numeric",Micro.mat="matrix",Progfact="vector"),
         prototype=list(HRTrain=matrix(1,1,1),HRTest=matrix(1,1,1),Ncv=100,Micro.mat=matrix(1,1,1),Progfact=c(NA))
)
#' Method show.
#' @name cvmv
#' @rdname cvmv-class
#' @exportMethod show
#' @rdname cvmv-class
#' @aliases show,cvmv-method
setMethod("show",signature="cvmv"
          , function(object){
            cat("Cross validation for Majority Votes Based Classification Analysis\n")
            cat("Number of cross valdiations used: ", object@Ncv, "\n")
            if(!is.null(object@Progfact)) cat("Prognostic factors used: ",object@Progfact,"\n")
          })





#' Method summary.
#' @name cvmv-class
#' @rdname cvmv-class
#' @exportMethod summary
#' @rdname cvmv-class
#' @aliases summary,cvmv-method
setMethod("summary",signature="cvmv", function(object){
  cat("Summary of majority votes cross validation analysis\n")
  cat("Number of prognostic factor used :",length(object@Progfact),"\n")
  cat("Number of cross validation: ", object@Ncv, "\n")
  cat("Estimated quantiles of the HR in the train dataset \n")
  print(stats::quantile(object@HRTrain[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
  cat("\n")
  cat("Estimated quantiles of the HR in the test dataset \n")
  print(stats::quantile(object@HRTest[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
})


#' Method plot.
#' setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#' @name cvmv-class
#' @rdname cvmv-class
#' @exportMethod plot
#' @rdname cvmv-class
#' @aliases plot,cvmv,ANY-method
#' @aliases cvmv-method
setMethod("plot", signature(x="cvmv"),
          function(x,  y, ...) {
            if (inherits(x, "cvmv") == FALSE) stop("Invalid class object")
            HRTest=x@HRTest
            HRTrain=x@HRTrain
            nCV=x@Ncv
            dotsCall = substitute(list(...))
            ll = eval(dotsCall)
            if(!methods::hasArg("xlab")) ll$xlab = "MCCV index"
            if(!methods::hasArg("ylab")) ll$ylab = "HR estimate"
            ll$main = "Estimated HR on Test Set \n for Low risk group"
            if(!methods::hasArg("cex.lab")) ll$cex.lab = 1.5
            if(!methods::hasArg("cex.main")) ll$cex.main = 1
            if(!methods::hasArg("col")) ll$col = 2

            ll$x=HRTest[,1]
            if(!methods::hasArg("ylim")) ll$ylim = c(0,2) #max(x@HRTrain,x@HRTest)


            #graphics::par(mfrow=c(1,2))
            t1 = which(HRTest[,1]<1)
            do.call(plot,args=ll)
            #plot(HRp.test[,1],ylim=c(0,2),ylab="HR",main="")
            for(i in 1:nCV){
              graphics::lines(c(i,i),HRTest[i,2:3])
            }
            for(i in t1){
              graphics::lines(c(i,i),HRTest[i,2:3],col=2)
            }
            graphics::abline(h=1)


            Results=data.frame(HRTrain=HRTrain[,1],HRTest=as.numeric(HRTest[,1]))
            ll$x=Results
            ll$names=c("Training ","Test ")
            ll$main = "Estimated HR on Training and Test Set \n for Low risk group"
            if(!methods::hasArg("xlab")) ll$xlab = ""
            if(!methods::hasArg("ylab")) ll$ylab = "HR estimate"
            if(!methods::hasArg("cex.lab")) ll$cex.lab = 1.5
            if(!methods::hasArg("cex.main")) ll$cex.main = 1
            if(!methods::hasArg("col")) ll$col = 2:3
            do.call(graphics::boxplot,args=ll)
          })
