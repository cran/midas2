% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/platform_midas2.R
\name{platform_midas2}
\alias{platform_midas2}
\title{An Information borrowing Bayesian platform design with subgroup efficacy exploration(midas-2)}
\usage{
platform_midas2(
  seed,
  p,
  p_tox,
  n.burnin = 10000,
  n.iter = 20000,
  n.thin = 2,
  C_T = 0.85,
  C_E1 = 0.15,
  C_E2 = 0.999
)
}
\arguments{
\item{seed}{set a random seed to maintain the repeatability of the simulation results.}

\item{p}{a matrix indicating the efficacy. Row number represents the number of candidate drugs.}

\item{p_tox}{a vector indicating the toxicity.}

\item{n.burnin}{the number of iterations in burn-in period, which means the first \emph{n.burnin} iterations are discarded.}

\item{n.iter}{the number of posterior iterations.}

\item{n.thin}{every \emph{n.thin} iterations after the burn-in period is retained for analysis.}

\item{C_T}{early toxicity stopping threshold, which refers to a predefined threshold used to determine when a clinical trial should be stopped early due to unacceptable levels of toxicity or adverse events in the study participants. This threshold is established to ensure the safety and well-being of the trial participants and to prevent further harm.}

\item{C_E1}{early futility stopping threshold, which refers to a predefined threshold used to determine when a clinical trial should be stopped early due to lack of efficacy or futility. It is established to prevent the continuation of a trial that is unlikely to demonstrate a significant treatment effect, thus saving time, resources, and participant exposure to ineffective treatments.}

\item{C_E2}{early efficacy stopping threshold, which refers to a predefined threshold used to determine when a clinical trial should be stopped early due to the demonstration of significant efficacy or positive treatment effects. This threshold is established to allow for timely decision-making and saves sample size.}
}
\value{
term.tox the indicator of whether early stopping for toxicity

term.fut the indicator of whether early stopping for futility

term.eff the indicator of whether early stopping for efficacy

final.eff a vector of final decision, either efficacy or inefficacy

post.subg subgroup analysis for treatments

post.sign signature analysis for treatments

best selection of best treatment for each subgroup

N sample size, which refers to the number of participants included in a study or experiment.
}
\description{
The rapid screening of effective and optimal therapies from large numbers of candidate combinations, as well as exploring subgroup efficacy, remains challenging, which necessitates innovative, integrated, and efficient trial designs.
MIDAS-2 package enables quick and continuous screening of promising combination strategies and exploration of their subgroup effects within a unified platform design framework. We used a regression model to characterize the efficacy pattern in subgroups. Information borrowing was applied through Bayesian hierarchical model to improve trial efficiency considering the limited sample size in subgroups. MIDAS-2 provides an adaptive drug screening and subgroup exploring framework to accelerate immunotherapy development in an efficient, accurate, and integrated fashion.
}
\details{
The sample size for a particular subgroup is limited. Therefore, it is difficult to accurately estimate treatment-related effects in each experiment arm separately, and the corresponding subgroup power may also be insufficient. Given that all candidate combination therapies in each arm consist of investigational drugs, it is justifiable to leverage information from specific subgroups across arms. Specifically, we adopt the Bayesian hierarchical model and assign half-cauchy prior distributions to the the standard deviation parameters.
}
\examples{


  # Example 1
  p0 <- c(0.1,0.1,0.1,0.1)
  p1 <- c(0.1,0.1,0.1,0.1)

  p <- rbind(p0,p1)
  p_tox <- c(0.1,0.4)

  # consider 1 candidate drugs with 4 subgroups
  result <- platform_midas2(seed=20,p,p_tox,n.burnin=1000,
            n.iter=2000,n.thin=2,C_T=0.85,C_E1=0.15,C_E2=0.999)
  result


  \donttest{
  # Example 2
  p0 <- c(0.05,0.10,0.05,0.10)
  p1 <- c(0.24,0.40,0.12,0.22)
  p2 <- c(0.24,0.40,0.12,0.22)
  p3 <- c(0.12,0.22,0.05,0.10)
  p4 <- c(0.24,0.40,0.12,0.22)
  p5 <- c(0.28,0.45,0.12,0.22)
  p6 <- c(0.24,0.40,0.12,0.22)
  p7 <- c(0.12,0.22,0.05,0.10)

  p <- rbind(p0, p1, p2, p3, p4, p5, p6, p7)
  p_tox <- c(0.10,0.10,0.10,0.10,0.10,0.10,0.15,0.20)

  # consider 7 candidate drugs with 4 subgroups
  result <- platform_midas2(seed=12,p,p_tox,n.burnin=10000,
            n.iter=20000,n.thin=2,C_T=0.85,C_E1=0.15,C_E2=0.999)
  result
  }

}
