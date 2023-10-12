#' The perm Class.
#'
#' Class of object returned by function \code{\link[MicrobiomeSurv]{DistHR}}.
#'
#' @name perm-class
#' @rdname perm-class
#' @exportClass perm
#' @param x	 A perm class object
#' @param y	 missing
#' @param  object A perm class object
#' @param ...	 The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot HRobs Estimated HR for low risk group on the original data.
#' @slot HRperm Estimated HR for low risk group on the permuted data.
#' @slot nperm Number of permutations carried out.
#' @slot Validation The validation scheme that was used.
#'
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{DistHR}}, \code{\link[MicrobiomeSurv]{EstimateHR}}, \code{\link[MicrobiomeSurv]{SurvPcaClass}}, \code{\link[MicrobiomeSurv]{SurvPlsClass}}, \code{\link[MicrobiomeSurv]{Majorityvotes}}, \code{\link[MicrobiomeSurv]{Lasoelascox}}
#' @importFrom methods graphics stats setClass setGeneric setMethod setRefClass
#' @importFrom grDevices
#' @importFrom hasArg new
#' @importFrom coef density median na.exclude p.adjust princomp qnorm quantile
#' @importFrom abline arrows axis barplot box boxplot legend lines par points

#' @note The first, third and last vertical line on the plot are the lower,
#' median  and upper CI of the permuted data estimated HR while
#' the red line is the estimated HR of the original data

setClass("perm",representation(HRobs="vector",HRperm="matrix",nperm="numeric",Validation="vector"),
         prototype=list(HRobs=as.vector(rep(NA,3)),HRperm=matrix(1,1,1),nperm=100,Validation=c(NA))
)

#' Method show.
#' @name perm
#' @rdname perm-class
#' @exportMethod show
#setGeneric("show", function(object) standardGeneric("show"))

#' @rdname perm-class
#' @aliases show,perm-method
setMethod("show",signature="perm"
          , function(object){
            cat("Estimated Null Ditribution of the ",object@Validation,"\n")
            cat("Number of Permutations: ", object@nperm, "\n")
          })


#' Method summary.
#' @name perm-class
#' @rdname perm-class
#' @exportMethod summary
#setGeneric("summary", function(object,...) standardGeneric("summary"))

#' @rdname perm-class
#' @aliases summary,perm-method
setMethod("summary",signature="perm", function(object){
  cat("Summary of Permutation Analysis\n")
  cat("validation scheme used :",object@Validation,"\n")
  cat("Number of Permutations: ", object@nperm, "\n")
  cat("\n")
  cat("Estimated  quantiles of the null distribution of HR\n")
  print(stats::quantile(object@HRperm[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
  cat("\n")
  cat("Estimated HR on original data\n")
  ttt=object@HRobs
  names(ttt)=c("Estimate","lower95CI","Upper95CI")
  print(ttt)
  })


#' Method plot.
#' setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#' @name perm-class
#' @rdname perm-class
#' @exportMethod plot

#' @rdname perm-class
#' @aliases plot,perm,ANY-method
#' @aliases perm-method
setMethod("plot", signature("perm"),
          function(x, y, ...) {
            if (inherits(x, "perm") == FALSE) stop("Invalid class object")

            HR=x@HRperm[,1]
            HR=stats::na.exclude(HR)
            n=x@nperm
            vv=x@HRobs[1]
            pvalue=sum(vv>HR)/n

             dotsCall = substitute(list(...))
            ll = eval(dotsCall)
            if(!methods::hasArg("xlab")) ll$xlab = paste("Estimated HR: Emperical p-value = ", pvalue ,sep="")
            if(!methods::hasArg("ylab")) ll$ylab = ""
            ll$main = "Null Distribution of HR on Permuted Data \n for low risk group"
            if(!methods::hasArg("cex.lab")) ll$cex.lab = 0.8
            if(!methods::hasArg("cex.main")) ll$cex.main = 1
            if(!methods::hasArg("col")) ll$col = 1
            if(!methods::hasArg("ylim")) ll$ylim = c(0,4)

            ll$x = stats::density(HR,from=0,to=(max(HR)+0.25))
            do.call(plot,args=ll)
            graphics::abline(v=vv,col=2)
            #CI for permuated cases
            qq=stats::quantile(sort(HR),prob=c(0.05,0.95))
            graphics::abline(v=qq[1],col=3)
            graphics::abline(v=qq[2],col=3)
            graphics::abline(v=stats::median(HR),col=3,lwd=3)
            })

