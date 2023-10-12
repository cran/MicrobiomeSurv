#' The cvle Class.
#'
#' Class of object returned by function \code{\link[MicrobiomeSurv]{CVLasoelascox}}.
#'
#' @name cvle-class
#' @rdname cvle-class
#' @exportClass cvle
#' @param x	 A cvle class object
#' @param y	 missing
#' @param type Plot type. 1 distribution of the HR under training and test set. 2 HR vs number selected taxa.
#' @param  object A cvle class object
#' @param ...	 The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot Coef.mat A matrix of coefficients with rows equals to number of cross validations and columns equals to number of taxa,
#' @slot lambda A vector of estimated optimum lambda for each iterations.
#' @slot n A vector of the number of selected taxa.
#' @slot mi.mat A matrix with 0 and 1. Number of rows equals to number of iterations and number of columns equals to number of taxa. 1 indicates that the particular taxon was selected or had nonzero coefficient and otherwise it is zero.
#' @slot HRTrain A matrix of survival information for the training dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.
#' @slot HRTest A matrix of survival information for the test dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.
#' @slot pld A vector of partial likelihood deviance at each cross validations.
#' @slot Micro.mat The microbiome matrix that was used for the analysis which can either be the full the full data or a reduced supervised PCA version.
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{EstimateHR}}, \code{\link[glmnet]{glmnet}}, \code{\link[MicrobiomeSurv]{Lasoelascox}}
#' @importFrom methods graphics stats setClass setGeneric setMethod setRefClass
#' @importFrom grDevices
#' @importFrom hasArg new

#' @docType class

setClass("cvle", representation(Coef.mat="matrix", lambda="vector", n="vector", mi.mat="matrix",
                                HRTrain="matrix", HRTest="matrix", pld="vector", Micro.mat="matrix"),
                 prototype=list(Coef.mat=matrix(1,1,1), lambda=c(NA), n=c(NA), mi.mat=matrix(1,1,1),
                                HRTrain=matrix(1,1,1), HRTest=matrix(1,1,1), pld=c(NA), Micro.mat=matrix(1,1,1))
)


#' Method show.
#' @name cvle
#' @rdname cvle-class
#' @exportMethod show
#setGeneric("show", function(object) standardGeneric("show"))

#' @rdname cvle-class
#' @aliases show,cvle-method
setMethod("show", signature ="cvle"
          , function(object){
            cat("Cross Valdiated Results for Lasso and Elastic Net based taxa\n")
            cat("Number of taxa used: ", length(rownames(object@Micro.mat)), "\n")
            cat("Number of CV: ", length(object@lambda), "\n")
          })


#' Method summary.
#' @name cvle-class
#' @rdname cvle-class
#' @exportMethod summary
#setGeneric("summary", function(object,...) standardGeneric("summary"))

#' @rdname cvle-class
#' @aliases summary,cvle-method
setMethod("summary", signature ="cvle"
          , function(object){
            cat("Summary of cross valdiateion for Lasso and Elastic Net based taxa\n")
            cat("Number of taxa used: ", length(rownames(object@Micro.mat)), "\n")
            cat("Estimated  quantiles of HR on test data\n")
            print(stats::quantile(object@HRTest[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
            cat("\n")
            cat("Estimated quantiles of HR on train data\n")
            print(stats::quantile(object@HRTrain[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
            cat("Mostly selected 5 taxa:\n")
            Freq=colSums(object@mi.mat)
            names(Freq)=rownames(object@Micro.mat)
            sFreq=sort(Freq, decreasing = TRUE)
            sFreq=sFreq[sFreq>0]
            maxG=length(sFreq)
            if (maxG>5) maxG=5
print(names(sFreq)[1:maxG])
          })


#' Method plot.
#' setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#' @name cvle-class
#' @rdname cvle-class
#' @exportMethod plot

#' @rdname cvle-class
#' @aliases plot,cvle,missing-method
#' @aliases cvle-method
          setMethod(f ="plot", signature(x="cvle", y="missing"),
                    function(x,  y, type=1, ...) {
                      if (inherits(x, "cvle") == FALSE) stop("Invalid class object")
                      if (type==1) {

                        DistHR=data.frame(HRTrain=x@HRTrain[,1],HRTest=x@HRTest[,1])

                        colnames(DistHR)=c("Training","Test")
                        dotsCall = substitute(list(...))
                        ll = eval(dotsCall)
                        if(!methods::hasArg("xlab")) ll$xlab = ""
                        if(!methods::hasArg("ylab")) ll$ylab = "HR estimate"
                        ll$main = "Distribution of HR on Training and Test Set \n for Low risk group"
                        if(!methods::hasArg("cex.lab")) ll$cex.lab = 1.5
                        if(!methods::hasArg("cex.main")) ll$cex.main = 1
                        if(!methods::hasArg("col")) ll$col = 2:3
                        ll$x=DistHR
                        do.call(graphics::boxplot,args=ll)

                      }


                      if (type==2) {

                        HRTest=x@HRTest[,1]
                        dotsCall = substitute(list(...))
                        ll = eval(dotsCall)
                        if(!methods::hasArg("xlab")) ll$xlab = "Estimated HR on Test Data"
                        if(!methods::hasArg("ylab")) ll$ylab = "Number of non zero coef."
                        ll$main = "HR vs number of taxa"
                        if(!methods::hasArg("cex.lab")) ll$cex.lab = 0.8
                        if(!methods::hasArg("cex.main")) ll$cex.main = 1
                        if(!methods::hasArg("col")) ll$col = 2
                        ll$x=HRTest
                        ll$y=x@n
                        do.call(plot,args=ll)

                      }

                    }
          )
