#' The cvmm Class.
#'
#' Class of object returned by function \code{\link[MicrobiomeSurv]{CVMSpecificCoxPh}}.
#'
#' plot signature(x = "cvmm"): Plots for \code{\link[MicrobiomeSurv]{CVMSpecificCoxPh}} class analysis results.
#'
#' Any parameters of \code{\link[graphics]{plot.default}} may be passed on to this particular plot method.
#'
#' @rdname cvmm-class
#' @exportClass cvmm
#' @param x	 A CVMSpecificCoxPh class object \code{\link[MicrobiomeSurv]{CVMSpecificCoxPh}}
#' @param y	 missing
#' @param object A CVMSpecificCoxPh class object
#' @param which This specify which taxon for which estimated HR information need to be visualized. By default results of the first taxon is used.
#' @param ...	The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot HRTrain A 3-way array, The first dimension is the number of taxa, the second dimension is the HR statistics for the low risk group in the train dataset (HR,1/HR LCI, UCI) while the third dimension is the number of cross validation performed.
#' @slot HRTest A 3-way array, The first dimension is the number of taxa, the second dimension is the HR statistics for the low risk group in the test dataset (HR,1/HR LCI, UCI) while the third dimension is the number of cross validation performed.
#' @slot train The selected subjects for each CV in the train dataset.
#' @slot test The selected subjects for each CV in the test dataset.
#' @slot n.mi The number of taxa used in the analysis.
#' @slot Ncv The number of cross validation performed.
#' @slot Rdata The microbiome data matrix that was used for the analysis either same as Micro.mat or a reduced version
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{CVMSpecificCoxPh}}
#' @importFrom methods graphics stats setClass setGeneric setMethod setRefClass base
#' @importFrom grDevices
#' @importFrom hasArg new
#' @importFrom abline arrows axis barplot box boxplot legend lines par points

setClass("cvmm",slots = representation(HRTrain="array",HRTest="array",train="matrix",test="matrix",n.mi="numeric",Ncv="numeric",Rdata="matrix"),
         prototype=list(HRTrain=array(NA,dim=c(1,1,1)),HRTest=array(NA,dim=c(1,1,1)),train=matrix(0,0,0), test=matrix(0,0,0),n.mi=1,Ncv=3,Rdata=matrix(0,0,0)))

#' Method show.
#' @name cvmm
#' @rdname cvmm-class
#' @exportMethod show
#setGeneric("show", function(object,...) standardGeneric("show"))

#' @rdname cvmm-class
#' @aliases show,cvmm-method
setMethod("show",signature="cvmm"
          , function(object){
            cat("Cross Validation for the Taxon specific analysis\n")
            cat("Number of Taxa used: ", object@n.mi, "\n")
            cat("Number of cross validations  performed: ", object@Ncv, "\n")
          })

#' Method summary.
#' @name cvmm-class
#' @rdname cvmm-class
#' @exportMethod summary
#' @aliases summary,cvmm-method
setMethod("summary",signature="cvmm"
          ,function(object,which=1){
            cat("Summary of Cross Validation for the Taxon specific analysis\n")
            cat("Estimated Median of the cross Validated HR for Taxa: ",which,"\n")
            HRTest  = object@HRTest[which,,][1,]
            HRTrain = object@HRTrain[which,,][1,]

            mean.alpha = stats::median(HRTrain,na.rm=T)
            se.alphal = stats::quantile(HRTrain,na.rm=T,probs = c(0.025))
            se.alphau = stats::quantile(HRTrain,na.rm=T,probs = c(0.975))
            cat("Estimated HR for Train Dataset \n")
            cat(paste(round(mean.alpha,4),"(",round(se.alphal,4)," , ",
                      round(se.alphau,4),")",sep=""))

            cat("\n")
            mean.alpha = stats::median(HRTest,na.rm=T)
            se.alphal = stats::quantile(HRTest,na.rm=T,probs = c(0.025))
            se.alphau = stats::quantile(HRTest,na.rm=T,probs = c(0.975))
            cat("Estimated HR for Test Dataset \n")
            cat(paste(round(mean.alpha,4),"(",round(se.alphal,4)," , ",
                      round(se.alphau,4),")",sep=""))
            cat("\n")
          }
)



#' Method plot.
#' @name cvmm-class
#' @rdname cvmm-class
#' @exportMethod plot

#' @rdname cvmm-class
#' setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#' @aliases plot,cvmm,ANY-method
#' @aliases cvmm-method
setMethod(f="plot", signature = "cvmm",
          definition =  function(x,  y, which=1, ...) {
            HRTest  = x@HRTest[which,,][1,]
            HRTrain = x@HRTrain[which,,][1,]
            Results=data.frame(HRTrain, HRTest)
            dotsCall = substitute(list(...))
            ll = eval(dotsCall)
            if(!methods::hasArg("xlab")) ll$xlab = ""

            if(!methods::hasArg("ylab")) ll$ylab = "HR estimate"
            if(!methods::hasArg("main")) ll$main = paste("Estimated HR of Low risk group for Taxon ", which, "\n Number of CVs = ",x@Ncv,sep="")
            if(!methods::hasArg("cex.lab")) ll$cex.lab = 1.5
            if(!methods::hasArg("cex.main")) ll$cex.main = 1
            if(!methods::hasArg("ylim")) ll$ylim = c(0, 5) #max(max(HRTest), max(HRTrain))
            if(!methods::hasArg("col")) ll$col = 2:3
            if(!methods::hasArg("names"))  ll$names=c("Training","Test")
            ll$x=Results
            do.call(graphics::boxplot,args=ll)
            return(invisible())
          }

)



