#' The cvsit Class.
#'
#' Class of object returned by function \code{\link[MicrobiomeSurv]{cvsit}}.
#'
#' @name cvsit-class
#' @rdname cvsit-class
#' @exportClass cvsit
#' @param x	 A cvsit class object
#' @param y	 missing
#' @param type Plot type. 1 distribution of the HR under test For the Top K taxa using PCA. 2 distribution of the HR under test For the Top K taxa using PLS.
#' @param  object A cvsit class object
#' @param ...	 The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot HRpca A 3-way array in which first, second, and third dimensions correspond to number of taxa, Hazard ratio information (Estimated HR, LowerCI and UpperCI), and number of cross validation respectively. This contains the estimated HR on test data and dimension reduction method is PCA.
#' @slot HRpls A 3-way array in which first, second, and third dimensions correspond to number of taxa, Hazard ratio information (Estimated HR, LowerCI and UpperCI), and number of cross validation respectively. This contains the estimated HR on test data and dimension reduction method is PLS.
#' @slot Ntaxa The number of taxa in the reduced matrix.
#' @slot Ncv The number of cross validation done.
#' @slot Top A sequence of top k taxa considered. Default is Top=seq(5,100,by=5).
#'
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{CVPcaPls}}, \code{\link[MicrobiomeSurv]{SurvPcaClass}}, \code{\link[MicrobiomeSurv]{SurvPlsClass}}

#' @importFrom methods graphics stats setClass setGeneric setMethod setRefClass
#' @importFrom grDevices
#' @importFrom hasArg new
#' @importFrom coef density median p.adjust princomp qnorm quantile
#' @importFrom abline arrows axis barplot box boxplot legend lines par points

setClass("cvsit",representation(HRpca="array",HRpls="array",Ntaxa="numeric",Ncv="numeric",Top="numeric"),
         prototype=list(HRpca=array(NA,dim=c(1,1,1)),HRpls=array(NA,dim=c(1,1,1)),Ntaxa=1,Ncv=3,Top=seq(5,100,by=5))
)

#' Method show.
#' @name cvsit
#' @rdname cvsit-class
#' @exportMethod show
#setGeneric("show", function(object) standardGeneric("show"))

#' @rdname cvsit-class
#' @aliases show,cvsit-method
setMethod("show",signature="cvsit"
          , function(object){
            cat("Cross Validation for sequentially increase taxa Analysis\n")
            cat("Number of Top K taxa used: ", object@Top, "\n")
            cat("Number of cross valdiations used: ", object@Ncv, "\n")
          })

#' Method summary.
#' @name cvsit-class
#' @rdname cvsit-class
#' @exportMethod summary
#setGeneric("summary", function(object,...) standardGeneric("summary"))

#' @rdname cvsit-class
#' @aliases summary,cvsit-method
setMethod("summary",signature="cvsit", function(object){
  cat("Results Based on Test Data\n")
  cat("Summary of Cross Validated Top K Taxa  Analysis\n")
  cat("Estimated Median of the HR for the cross Validated HR for Top K Taxa \n")
  nn=dim(object@HRpca)[3]
  mean.alpha= sapply(1:nn,function(i) stats::median(object@HRpca[,1,i],na.rm=T))
  se.alphal= sapply(1:nn,function(i) stats::quantile(object@HRpca[,1,i],na.rm=T,probs = c(0.025)))
  se.alphau= sapply(1:nn,function(i) stats::quantile(object@HRpca[,1,i],na.rm=T,probs = c(0.975)))
  mx1=paste(round(mean.alpha,3),"(",round(se.alphal,3),"-",round(se.alphau,3),")",sep="")


  mean.alpha= sapply(1:nn,function(i) stats::median(object@HRpls[,1,i],na.rm=T))
  se.alphal= sapply(1:nn,function(i) stats::quantile(object@HRpls[,1,i],na.rm=T,probs = c(0.025)))
  se.alphau= sapply(1:nn,function(i) stats::quantile(object@HRpls[,1,i],na.rm=T,probs = c(0.975)))
  mx3=paste(round(mean.alpha,3),"(",round(se.alphal,3),"-",round(se.alphau,3),")",sep="")


  HR=data.frame(rbind(mx1,mx3))
  colnames(HR)=object@Top
  rownames(HR)=c("(PCA)","(PLS)")
  print(HR)
}
)

#' Method plot.
#' setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#' @name cvsit-class
#' @rdname cvsit-class
#' @exportMethod plot

#' @rdname cvsit-class
#' @aliases plot,cvsit,missing-method
#' @aliases cvsit-method
setMethod("plot", signature(x="cvsit", y="missing"),
          function(x,  y, type=1, ...) {
            if (inherits(x, "cvsit") == FALSE) stop("Invalid class object")
            if (type==1) {
            nn=dim(x@HRpca)[3]
            PC.HRp=x@HRpca[,1,1:nn]
            colnames(PC.HRp)=x@Top

            dotsCall = substitute(list(...))
            ll = eval(dotsCall)
            if(!methods::hasArg("xlab")) ll$xlab = "Top K Taxa"
            if(!methods::hasArg("ylab")) ll$ylab = "Cross Validated HR"
            ll$main = "Estimated HR on Test Data \n for Top K Taxa (PCA)"
            if(!methods::hasArg("cex.lab")) ll$cex.lab = 1.2
            if(!methods::hasArg("cex.main")) ll$cex.main = 1.3
            if(!methods::hasArg("col")) ll$col = 1:nn+1
            if(!methods::hasArg("ylim")) ll$ylim = c(0, 5) #max(PC.HRp)
            ll$x=PC.HRp
            do.call(graphics::boxplot,args=ll)
            }

            if (type==2) {
            nn=dim(x@HRpca)[3]
            PL.HRp=x@HRpls[,1,1:nn]
            colnames(PL.HRp)=x@Top

            dotsCall = substitute(list(...))
            lll = eval(dotsCall)
            if(!methods::hasArg("xlab")) lll$xlab = "Top K Taxa"
            if(!methods::hasArg("ylab")) lll$ylab = "Cross Validated HR"
            lll$main = "Estimated HR Test Data \n for Top K Taxa (PLS)"
            if(!methods::hasArg("cex.lab")) lll$cex.lab = 1.2
            if(!methods::hasArg("cex.main")) lll$cex.main = 1.3
            if(!methods::hasArg("col")) lll$col = 1:nn+1
            if(!methods::hasArg("ylim")) lll$ylim = c(0,5) #max(PL.HRp)
            lll$x=PL.HRp
            do.call(graphics::boxplot,args=lll)
            return(invisible())
          }
}
)
