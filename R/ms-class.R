#' The ms Class.
#'
#' Class of object returned by function \code{\link[MicrobiomeSurv]{MSpecificCoxPh}}.
#' plot signature(x = "ms"): Plots for ms class analysis results
#'
#' Any parameters of \code{\link[graphics]{plot.default}} may be passed on to this particular plot method.
#'
#' show(ms-object)
#' @name ms-class
#' @rdname ms-class
#' @exportClass ms
#' @param x	 A ms class object
#' @param y	 missing
#' @param object	 A ms class object
#' @param ...	The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot Result A list of dataframes of each output object of coxph for the taxa.
#' @slot HRRG A dataframe with estimated taxon-specific HR for low risk group and 95 percent CI.
#' @slot Group A matrix of the classification group a subject belongs to for each of the taxon analysis. The taxa are on the rows and the subjects are the columns
#' @slot Mi.names The names of the taxon for the analysis
#' @author Thi Huyen Nguyen, \email{thihuyen.nguyen@@uhasselt.be}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@@gsk.com}
#' @author Ziv Shkedy
#' @seealso \code{\link[MicrobiomeSurv]{MSpecificCoxPh}}
#' @importFrom methods graphics stats setClass setGeneric setMethod setRefClass
#' @importFrom grDevices rainbow
#' @importFrom hasArg new
#' @importFrom coef density median na.exclude p.adjust princomp qnorm quantile
#' @importFrom abline arrows axis barplot box boxplot legend lines par points
#' @importFrom graphics par

setClass("ms",slots = list(Result="list", HRRG="matrix", Group="matrix", Mi.names="vector"),
         prototype=list(Result=list(1), HRRG=matrix(0,0,0), Group=matrix(0,0,0), Mi.names = vector()))


#' Method show.
#' @name ms
#' @rdname ms-class
#' @exportMethod show
#setGeneric("show", function(object) standardGeneric("show"))

#' @rdname ms-class
#' @aliases show,ms-method
#' @aliases ms,ANY
setMethod("show",signature="ms"
          , function(object){
            cat("Taxon by taxon CoxPh Model\n")
            cat("Number of taxa used: ", length(object@Mi.names), "\n")
          })

#' Method summary.
#' @name ms-class
#' @rdname ms-class
#' @exportMethod summary
# setGeneric("summary", function(object,...) standardGeneric("summary"))

#' @rdname ms-class
#' @aliases summary,ms-method
setMethod("summary",signature="ms",function(object){
  cat("Summary of taxon by taxon CoxPh Models\n")
  cat("Number of taxa used: ", length(object@Mi.names), "\n")
  cat("Top", length(object@Mi.names), " taxa out of ", length(object@Mi.names), "\n")
  cat("Estimated HR for the low risk group\n")
  # top taxa based on upper CI HR GS+
  Names.Ktaxa=object@Mi.names
  index.Top.Ktaxa = order(object@HRRG[,1],decreasing =FALSE)
  index.Top.Ktaxa=index.Top.Ktaxa[1:length(object@Mi.names)]
  Top.Ktaxa.GSplus=data.frame(Mi.names=Names.Ktaxa[index.Top.Ktaxa],object@HRRG[index.Top.Ktaxa,])

  colnames(Top.Ktaxa.GSplus)=c("Mi.names", "HR GS+", "LowerCI", "UpperCI", "p-value")
  # FDR corrected CI for top k taxa
  cilevel = 1-0.05*length(object@Mi.names)/nrow(object@HRRG)

  HRpadj = exp(log(object@HRRG[index.Top.Ktaxa,1]) + log(object@HRRG[index.Top.Ktaxa,c(2,3)])-log(object@HRRG[index.Top.Ktaxa,1])*stats::qnorm(cilevel)/1.96)  #
  res.topKtaxa=data.frame(Top.Ktaxa.GSplus,HRpadj)

  colnames(res.topKtaxa)=c("Mi.names", "HR", "LCI", "UCI", "p-value", "FDRLCI", "FDRUCI")
  print(res.topKtaxa)
})



#' Method plot.
#' setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#' @name ms-class
#' @rdname ms-class
#' @exportMethod plot

#' @rdname ms-class
#' @aliases plot,ms,ANY-method
#' @aliases ms-method
setMethod(f="plot", signature = "ms",
          definition = function(x,y,...){
            object =  x
            Names.Ktaxa=object@Mi.names
            index.Top.Ktaxa = order(object@HRRG[,1],decreasing =FALSE)
            Top.Ktaxa.GSplus=data.frame(Mi.names=Names.Ktaxa, object@HRRG)

            colnames(Top.Ktaxa.GSplus)=c("Mi.names", "HR GS+", "LowerCI", "UpperCI", "p-value")
            # FDR corrected CI for top k taxa
            cilevel = 1-0.05*5/nrow(object@HRRG)

            HRpadj = exp(log(object@HRRG[index.Top.Ktaxa,1]) + ((log(object@HRRG[index.Top.Ktaxa,c(2,3)])-log(object@HRRG[index.Top.Ktaxa,1]))*stats::qnorm(cilevel)/1.96))  #
            res.topKtaxa=data.frame(Top.Ktaxa.GSplus,HRpadj)

            colnames(res.topKtaxa)=c("Mi.names","HR","LCI","UCI","p-value", "FDRLCI","FDRUCI")

            x.axis= 1:length(object@Mi.names)

            old_pars = graphics::par(mfrow=c(3,1))
            on.exit(graphics::par(old_pars))

            Group2 = object@Group
            prop = sapply(1:ncol(Group2),function(k) sum(Group2[,k]=="Low risk")/nrow(Group2))
            graphics::barplot(prop, col=grDevices::rainbow(length(Group2[1,])), xlab = "Subject index", ylab = "Proportion",main = "Proportion of being classified as low risk group",ylim = c(0,max(prop)))


            plot(x=x.axis, y = res.topKtaxa[,2],
                 ylim= c(0,max(res.topKtaxa[,4])),
                 pch=19, xlab="Taxa", ylab="Hazard ratio",
                 main="Hazard ratio plot with unadjusted confidence interval")
            graphics::arrows(x.axis, res.topKtaxa[,3], x.axis, res.topKtaxa[,4], length=0.05, angle=90, code=3)
            graphics::abline(h=1,col="red2",lwd=2.0)

            plot(x=x.axis, y = res.topKtaxa[,2],
                 ylim= c(0,max(res.topKtaxa[,7])),
                 pch=19, xlab="Taxa", ylab="Hazard ratio",
                 main="Hazard ratio plot with adjusted confidence interval")
            graphics::arrows(x.axis, res.topKtaxa[,6], x.axis, res.topKtaxa[,7], length=0.05, angle=90, code=3)
            graphics::abline(h=1,col="red2",lwd=2.0)


            return(invisible())
          })

