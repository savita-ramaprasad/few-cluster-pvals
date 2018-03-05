\name{cluster.webb.glm}
\alias{cluster.webb.glm}
\title{Calculate Cluster-Robust p-Values and Confidence Intervals Based on Webb (2013)}
\description{
 Calculate cluster robust p-values and confidence intervals using wild cluster bootstrapped t-statistics based on Webb (2013) which is the prefered method to use when the number of clusters are < 15. Webb, M. D. (2013). <Reworking wild bootstrap based inference for clustered errors (No. 1315). Queen's Economics Department Working Paper>.}
\usage{
cluster.webb.glm(mod, dat, cluster, vars.boot = NULL, 
ci.level = 0.95, impose.null = TRUE, boot.reps = 1000, 
report = TRUE, prog.bar = TRUE, output.replicates = FALSE) 
}
\arguments{
  \item{mod}{A linear model estimated using glm.}
  \item{dat}{the data set used to estimate mod.}
  \item{cluster}{The clustering variable.}
  \item{vars.boot}{The variables to bootstrap over when null is imposed. Default is p-values     for all variables. Displays p-values for all variables with null not imposed.}
  \item{ci.level}{The confidence level of the confidence interval. Reported when impose.null = FALSE}
  \item{impose.null}{Should null be imposed?}
  \item{boot.reps}{The number of bootstrap repititions.}
  \item{report}{Report the result to the console?}
  \item{prog.bar}{Show a progress bar of the bootstrap?}
  \item{output.replicates}{Should the cluster bootstrap replicates be outputted as well?}
}

\value{
	\item{p.values}{A vector of estimated p-values.}
	\item{ci}{A matrix of confidence intervals, reported when null is not imposed.}
}

\author{
	Savita Ramaprasad
}

\note{
	Original code to estimate p-values and ci from GLM wild cluster robust bootstrap-t statistics by Justin Esarey: <https://cran.r-project.org/web/packages/clusterSEs/clusterSEs.pdf> modified in the following ways. 1. Modified code to be based on Webb (2013) (rather than be based on CGM (2004) 2. included an option to obtain bootstrap p-values for specific variables (when null hypothesis is imposed) to reduce run-time. Previously bootstrap p-values were calculated for every variable in the model. 3. Changed the way results are printed (when report = TRUE) to reduce run-time. Results reported look less pretty now. Overall, run time for reduced from 12 mins to 18 secs. 	
}

\references{
	
	Esarey, Justin, and Andrew Menger. 2017. "Practical and Effective Approaches to Dealing with Clustered Data." Political Science Research and Methods forthcoming: 1-35. <URL:http://jee3.web.rice.edu/cluster-paper.pdf>.
	
	Reworking wild bootstrap based inference for clustered errors (No. 1315). Queen's Economics Department Working Paper
	
}
\examples{
}
\keyword{file}
