% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Top1Uni.R
\name{Top1Uni}
\alias{Top1Uni}
\title{This function finds out the taxon has the smallest p-value, then calculate risk score of patients based on that taxon.
Categorized subjects into high or low risk groups based on the mean of the risk score as a cutoff point
Checking whether the two groups are significant difference in the probability to be survival.}
\usage{
Top1Uni(Result, Micro.mat, Survival, Censor, Plots = FALSE)
}
\arguments{
\item{Result}{A Result statistic of all taxon.}

\item{Micro.mat}{A large or small microbiome matrix. A matrix with microbiome profiles where the number of rows should be equal to the number of taxa and number of columns should be equal to number of patients.}

\item{Survival}{Survival A vector of survival time with length equals to number of subjects}

\item{Censor}{A vector of censoring indicator}

\item{Plots}{A boolean parameter indicating if plots should be shown. Default is FALSE. If TRUE, the first plot is plot of the observed Kaplan-Meier curves per group while the second is boxplot of the two groups.}
}
\value{
A list is returned with the following values
\item{name.top1}{Taxon having the smallest p-value in the univariate coxPH model}
\item{sum.top1}{Result statistic of the taxon containing coefficient, exponential of coefficient, raw p.value using LRT, and p.value after using BH adjustment}
\item{KMplot.top1}{Kaplan-Meier plot}
\item{log.rank.top1}{Log-rank test}
}
\description{
This function finds out the taxon has the smallest p-value, then calculate risk score of patients based on that taxon.
Categorized subjects into high or low risk groups based on the mean of the risk score as a cutoff point
Checking whether the two groups are significant difference in the probability to be survival.
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
# Obtain summary statistics for families
summary_fam_shan_w3 = CoxPHUni(Survival = surv_fam_shan_w3$Survival,
                               Censor = surv_fam_shan_w3$Censor,
                               Prognostic = prog_fam_shan_w3,
                               Micro.mat = fam_shan_trim_w3,
                               Method = "BH")

# Analysis of the taxon having smallest p-value (in the result of using CoxPHUni function)
top1_fam_shan_w3 = Top1Uni(Result = summary_fam_shan_w3,
                           Micro.mat = fam_shan_trim_w3,
                           Survival = surv_fam_shan_w3$Survival,
                           Censor = surv_fam_shan_w3$Censor,
                           Plots = TRUE)
}
\author{
Thi Huyen Nguyen, \email{thihuyen.nguyen@uhasselt.be}

Olajumoke Evangelina Owokotomo, \email{olajumoke.x.owokotomo@gsk.com}

Ziv Shkedy
\code{\link[MicrobiomeSurv]{Top1Uni}}
}
