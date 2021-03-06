% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parTemper.R
\name{parTemper}
\alias{parTemper}
\title{Parallel tempering for function minimization}
\usage{
parTemper(costFunc, init, ..., control = list())
}
\arguments{
\item{costFunc}{The cost function (often a negative log-likelihood)}

\item{init}{The starting point for sampling or the output from a previous call 
to parTemper}

\item{...}{Further arguments to be passed to costFunc}

\item{control}{A list of control parameters. See Details.}
}
\value{
An object of class \code{parTemper} that is a list containing the samples along 
with summary information
}
\description{
Minimize an input cost function using a parallel tempering algorithm that calls the seminal adaptive Metropolis MCMC sampling algorithm in package samMCMC for within chain sampling
}
\details{
costFunc is a function to be minimized with respect to its first input, x. Additional inputs to costFunc can be passed as named variables. The minimization is initialized using the input variable init in one of two ways

(1) The initial value x = init is directly given
(2) This is a continuation of a previous parallel tempering "run" for which init is the output of that run

parTemper requires three key control variables to run: (1) a vector of temperatures, tempVect, ordered from small to large temperatures; (2) the number of cycles, numCyc, where each cycle consists of a set of within-chain Metropolis samples at each temperature followed by an attempt to swap samples across chains; and (3) the number of within-chain Metropolis samples for each temperature, M.

There are three ways to specify these (and other) control variables: a default, a named variable in the input control list, and (for continued chains) the control variable used in the previous chains, init$control.

Different expectations apply for these inputs for new vs. continuing chains.  For example, for tempVect an input is always preferred over the default for new chains (default < input) and tempVect cannot be reset for continuing chains (an error is raised if this is attempted). These expectations can be summarized thus for the three required control variables:

  tempVect
    new:    default < input
    cont:   must use prev; error if tempVect is input
  numCyc 
    new:    default < input
    cont:   prev < input
  M
    new:    default < input
    cont:   prev < input

In addition to these three, key control variables (tempVect, numCyc, and M), additional control variables can be specified to influence both the behavior of parTemper and the call to samMCMC by saTemper. If the control variable verbose is TRUE [the default is FALSE] diagnostic information is printed as the algorithm runs; the preference ordering for verbose is default < prev < input, where prev is the previous value used if this is a continued run.

In addition to control variables that apply directly to parTemper, the following four control variables for the call to samMCMC an be specified: C_0, t_0, epsilon, and verbose. verbose for the call to samMCMC is set via the control variable metropVerbose since verbose is also a control variable used in saTemper. C_0, t_0, and epsilon can only be set for new chains. If the user attempts to set them for continued chains, an error is raised. If they are not set for a new chain, samMCMC will use default values. metropVerbose can always be set and reset, with the preference ordering default [FALSE] < prev < input.

Finally, we provide a brief summary of the algorithm itself. Further details can  be found in the technical material linked to in the README at https://github.com/santafemha/parTempeR.

[TODO: Ensure that this technical material is linked in the README, or update with arxiv post or peer reviewed article when ready.] 

Let k = 1...K index temperatures in the vector tempVect. Each cycle of the algorithm consists of a set of within-in chain Metropolis samples for each temperature, followed by a single attempt to randomly swap samples across chains for a single, randomly chosen pair of adjacent temperatures (k and k+1). M within chain samples are made by calling the adaptive Metropolis MCMC algorithm implemented by samMCMC in package samMCMC (the main function and package have the same name).

Following these within-chain samples, an attempt is made to swap samples between one pair of adjacent temperatures, also using a Metropolis criterion. The pair for which this is attempted is randomly chosen. Let k be the lower temperature. The weighting accorded the non-swapped configuration is exp(-c_k/T_k)*exp(-c_kp1/T_kp1), where kp1 stands for k+1 (k plus one) and c_k (etc.) is the cost function evaluated for the current value of the chain, X_k. The weighting accorded the swapped configuration is exp(-c_kp1/T_k)*exp(-c_k/T_kp1). The swap is always accepted if the swapped configuration has a higher weighting; otherwise, it is accepted with a probability of the proportion of these weightings. To be precise, the Metropolis acceptance ratio is

a = min(1,exp(-(c_kp1-c_k)*(1/T_k-1/T_kp1)))

The swap is accepted with probability a, and otherwise rejected. The cycle of within-chain sampling followed by a single swap attempt is repeated numCyc times. The minimum value of this procedure can be extracted from the output by calling the function parTemperBest.
}
\author{
Michael Holton Price <MichaelHoltonPrice@gmail.com>
}
