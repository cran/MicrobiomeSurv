% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SurvPcaClass.R
\name{SurvPcaClass}
\alias{SurvPcaClass}
\title{Survival PCA and Classification for microbiome data}
\usage{
SurvPcaClass(
  Survival,
  Micro.mat,
  Censor,
  Reduce = TRUE,
  Select = 5,
  Prognostic = NULL,
  Plots = FALSE,
  Mean = TRUE,
  Quantile = 0.5
)
}
\arguments{
\item{Survival}{A vector of survival time with length equals to number of subjects}

\item{Micro.mat}{A large or small microbiome profile matrix. A matrix with microbiome profiles where the number of rows should be equal to the number of microbiome and number of columns should be equal to number of patients.}

\item{Censor}{A vector of censoring indicator}

\item{Reduce}{A boolean paramier indicating if the microbiome profile matrix should be reduced, default is TRUE and larger microbiome profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.}

\item{Select}{Number of microbiome (default is 15) to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE}

\item{Prognostic}{A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.}

\item{Plots}{A boolean paramier indicating if the plots should be shown. Default is FALSE}

\item{Mean}{The cut off value for the classifier, default is the mean cutoff}

\item{Quantile}{If user want to use quantile as cutoff point. They need to specify Mean = FALSE and a quantile that they want to use. The default is the median cutoff}
}
\value{
A object of class SurvPca is returned with the following values
  \item{Survfit}{The cox proportional regression result using the first PCA}
  \item{Riskscores}{A vector of risk scores which is equal to the number of patents.}
  \item{Riskgroup}{The classification of the subjects based on the PCA into low or high risk group}
  \item{pc1}{The First PCA scores based on either the reduced microbiome matrix or the full matrix}
\item{KMplot}{The Kaplan-Meier survival plot of the riskgroup}
 \item{SurvBPlot}{The distribution of the survival in the riskgroup}
 \item{Riskpca}{The plot of Risk scores vs first PCA}
}
\description{
The function performs principal component analysis (PCA) on microbiome matrix and fit Cox proportional hazard model with covariates using also the first PCA as covariates.
}
\details{
This function can handle single and multiple microbiome. For larger microbiome matrix,
this function will reduce largermicrobiome matrix to smaller version using supervised pca approach and this is by default done and can be control by using the argument Reduce.
Other prognostic factors can be included to the model.
}
\examples{
# Prepare data
data(Week3_response)
Week3_response = data.frame(Week3_response)
surv_fam_shan_w3 = data.frame(cbind(as.numeric(Week3_response$T1Dweek),
as.numeric(Week3_response$T1D)))
colnames(surv_fam_shan_w3) = c("Survival", "Censor")
prog_fam_shan_w3 = data.frame(factor(Week3_response$Treatment_new))
colnames(prog_fam_shan_w3) = c("Treatment")
data(fam_shan_trim_w3)
names_fam_shan_trim_w3 =
c("Unknown", "Lachnospiraceae", "S24.7", "Lactobacillaceae", "Enterobacteriaceae", "Rikenellaceae")
fam_shan_trim_w3 = data.matrix(fam_shan_trim_w3[ ,2:82])
rownames(fam_shan_trim_w3) = names_fam_shan_trim_w3
# Using the function
SPCA_fam_shan_w3 = SurvPcaClass(Survival = surv_fam_shan_w3$Survival,
                                Micro.mat = fam_shan_trim_w3,
                                Censor = surv_fam_shan_w3$Censor,
                                Reduce=TRUE,
                                Select=5,
                                Prognostic = prog_fam_shan_w3,
                                Plots = TRUE,
                                Mean = TRUE)

# Getting the survival regression output
SPCA_fam_shan_w3$SurvFit

# Getting the riskscores
SPCA_fam_shan_w3$Riskscores

# Getting the riskgroup
SPCA_fam_shan_w3$Riskgroup

# Obtaining the first principal component scores
SPCA_fam_shan_w3$pc1
}
\seealso{
\code{\link[survival]{coxph}},
\code{\link[MicrobiomeSurv]{EstimateHR}}, \code{\link[stats]{princomp}},
 \code{\link[MicrobiomeSurv]{SurvPlsClass}}
}
\author{
Thi Huyen Nguyen, \email{thihuyen.nguyen@uhasselt.be}

Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@gsk.com}

Ziv Shkedy
}
