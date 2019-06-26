#' @title Parallel tempering for function minimization
#'
#' @description Minimize an input cost function using a parallel tempering algorithm that calls the seminal adaptive Metropolis MCMC sampling algorithm in package samMCMC for within chain sampling
#'
#' @details costFunc is a function to be minimized with respect to its first input, x. Additional inputs to costFunc can be passed as named variables. The minimization is initialized using the input variable init in one of two ways
#'
#' (1) The initial value x = init is directly given
#' (2) This is a continuation of a previous parallel tempering "run" for which init is the output of that run
#'
#' parTemper requires three key control variables to run: (1) a vector of temperatures, tempVect, ordered from small to large temperatures; (2) the number of cycles, numCyc, where each cycle consists of a set of within-chain Metropolis samples at each temperature followed by an attempt to swap samples across chains; and (3) the number of within-chain Metropolis samples for each temperature, M.
#'
#' There are three ways to specify these (and other) control variables: a default, a named variable in the input control list, and (for continued chains) the control variable used in the previous chains, init$control.
#'
#' Different expectations apply for these inputs for new vs. continuing chains.  For example, for tempVect an input is always preferred over the default for new chains (default < input) and tempVect cannot be reset for continuing chains (an error is raised if this is attempted). These expectations can be summarized thus for the three required control variables:
#'
#'   tempVect
#'     new:    default < input
#'     cont:   must use prev; error if tempVect is input
#'   numCyc 
#'     new:    default < input
#'     cont:   prev < input
#'   M
#'     new:    default < input
#'     cont:   prev < input
#'
#' In addition to these three, key control variables (tempVect, numCyc, and M), additional control variables can be specified to influence both the behavior of parTemper and the call to samMCMC by saTemper. If the control variable verbose is TRUE [the default is FALSE] diagnostic information is printed as the algorithm runs; the preference ordering for verbose is default < prev < input, where prev is the previous value used if this is a continued run.
#'
#' In addition to control variables that apply directly to parTemper, the following four control variables for the call to samMCMC an be specified: C_0, t_0, epsilon, and verbose. verbose for the call to samMCMC is set via the control variable metropVerbose since verbose is also a control variable used in saTemper. C_0, t_0, and epsilon can only be set for new chains. If the user attempts to set them for continued chains, an error is raised. If they are not set for a new chain, samMCMC will use default values. metropVerbose can always be set and reset, with the preference ordering default [FALSE] < prev < input.
#'
#' Finally, we provide a brief summary of the algorithm itself. Further details can  be found in the technical material linked to in the README at https://github.com/santafemha/parTempeR.
#'
#' [TODO: Ensure that this technical material is linked in the README, or update with arxiv post or peer reviewed article when ready.] 
#'
#' Let k = 1...K index temperatures in the vector tempVect. Each cycle of the algorithm consists of a set of within-in chain Metropolis samples for each temperature, followed by a single attempt to randomly swap samples across chains for a single, randomly chosen pair of adjacent temperatures (k and k+1). M within chain samples are made by calling the adaptive Metropolis MCMC algorithm implemented by samMCMC in package samMCMC (the main function and package have the same name).
#'
#' Following these within-chain samples, an attempt is made to swap samples between one pair of adjacent temperatures, also using a Metropolis criterion. The pair for which this is attempted is randomly chosen. Let k be the lower temperature. The weighting accorded the non-swapped configuration is exp(-c_k/T_k)*exp(-c_kp1/T_kp1), where kp1 stands for k+1 (k plus one) and c_k (etc.) is the cost function evaluated for the current value of the chain, X_k. The weighting accorded the swapped configuration is exp(-c_kp1/T_k)*exp(-c_k/T_kp1). The swap is always accepted if the swapped configuration has a higher weighting; otherwise, it is accepted with a probability of the proportion of these weightings. To be precise, the Metropolis acceptance ratio is
#'
#' a = min(1,exp(-(c_kp1-c_k)*(1/T_k-1/T_kp1)))
#'
#' The swap is accepted with probability a, and otherwise rejected. The cycle of within-chain sampling followed by a single swap attempt is repeated numCyc times. The minimum value of this procedure can be extracted from the output by calling the function parTemperBest.
#'
#' @param costFunc The cost function (often a negative log-likelihood)
#' @param init The starting point for sampling or the output from a previous call 
#' to parTemper
#' @param ... Further arguments to be passed to costFunc
#' @param control A list of control parameters. See Details.
#'
#' @return An object of class \code{parTemper} that is a list containing the samples along 
#' with summary information
#'
#' @author Michael Holton Price <MichaelHoltonPrice@@gmail.com>
#' 
#' @export

parTemper <- function(costFunc,init,...,control=list()) {
  # The following control variables for the call to samMCMC can be set
  # initially for new runs, but not modified for continuing runs:
  #
  # C_0
  # t_0
  # epsilon
  #
  # The following control variable for the call to samMCMC can be set each
  # time:
  #
  # verbose [metropVerbose in the parTemper control variable]
  #
  # The following control variables for the call to samMCMC are fixed:
  #
  # thinning = 1
  # numSampBurn = 0

  newRun <- !('parTemper' %in% class(init))
  # K is the number of temperatures
  # M is the samples between swap attempts
  # numCyc is the number of cycles

  # Save the input control as inputControl
  inputControl <- control
  # Create a new control object
  control <- list()
  # Set prevControl (if applicable; if not, set to NA)
  if(!newRun) {
    prevControl <- init$control
  } else {
    prevControl <- NA
  }


  # Set tempVect
  if(newRun) {
    if(!('tempVect' %in% names(inputControl))) {
      control$tempVect <- 200 * .75^(50 - 1:50)
    } else {
      control$tempVect <- inputControl$tempVect
    }
  } else { # !newRun
      if('tempVect' %in% names(inputControl)) {
        stop('tempVect should not be input for continuing chains')
      } else {
        control$tempVect <- prevControl$tempVect
      }
  }

  if(is.unsorted(control$tempVect)) {
    stop('tempVect must be ordered from small to large temperatures')
  }

  # Set M (samples between attempted swap), numCycles, and verbose
  control$M <- samMCMC::chooseControlValue('M',20,inputControl,!newRun,prevControl)
  control$numCycles <- samMCMC::chooseControlValue('numCycles',100,inputControl,!newRun,prevControl)
  control$verbose <- samMCMC::chooseControlValue('verbose',F,inputControl,!newRun,prevControl)

  # Build metropControl, the control variable input to samMCMC
  metropControl <- list()
  # C_0, t_0, and epsilon can only be set for new chains. Throw an error if the
  # user attempts to set them for continuing chains.
  if('C_0' %in% names(inputControl)) {
    if(newRun) {
      metropControl$C_0 <- inputControl$C_0
    } else {
      stop('C_0 should not be input for continuing chains')
    }
  }

  if('t_0' %in% names(inputControl)) {
    if(newRun) {
      metropControl$t_0 <- inputControl$t_0
    } else {
      stop('t_0 should not be input for continuing chains')
    }
  }
  if('epsilon' %in% names(inputControl)) {
    if(newRun) {
      metropControl$epsilon <- inputControl$epsilon
    } else {
      stop('epsilon should not be input for continuing chains')
    }
  }

  # verbose in metropControl can always be set
  if('metropVerbose' %in% names(inputControl)) {
    metropControl$verbose <- inputControl$metropVerbose
  } else {
    if(newRun) {
      metropControl$verbose <- F
    } else {
      metropControl$verbose <- prevControl$metropControl$verbose
    }
  }

  # numSamp in metropControl is M
  metropControl$numSamp <- control$M

  # thinning is 1
  metropControl$thinning <- 1
  # numSampBurn is 0
  metropControl$numSampBurn <- 0
  # direct is FALSE (since a temperature must be used)
  metropControl$direct <- FALSE
  # [temp is reset each time]
  control$metropControl <- metropControl

  if(newRun) {
    chainList <- list()
  } else {
    chainList <- init$chainList
  }
  K <- length(control$tempVect) 

  # Matrix to store swap information
  if(newRun) {
    numSampPrev <- 0
  } else {
    numSampPrev <- ncol(init$chainList[[1]]$X_mat)
  }
  # columns are t [samples] / k / swapped
  swapMat <- matrix(NA,3,control$numCycles)

  # Iterate over number of cycles
  for(cc in 1:control$numCycles) {
    # Extend chains for this cycle
    for(k in 1:K) {
      metropControl_k <- control$metropControl
      metropControl_k$temp <- control$tempVect[k]
      if(newRun && cc == 1) {
        chainList[[k]] <- samMCMC::samMCMC(costFunc,init,...,control=metropControl_k)
      } else {
        chainList[[k]] <- samMCMC::samMCMC(costFunc,chainList[[k]],...,control=metropControl_k)
      }
    }

    if(control$verbose) {
      print('**--**--**')
      print(cc)
      print(min(unlist(lapply(chainList,function(x){min(x$sampFuncVect)}))))
      state0 <- unlist(lapply(chainList,function(x){x$sampFuncVect[length(x$sampFuncVect)]}))
    }

    # Randomly choose an adjacent pair of temperatures to attempt a swap
    k <- sample(1:(K-1),1)
    c_k   <- chainList[[k  ]]$sampFuncVect[length(chainList[[k  ]]$sampFuncVect)]
    c_kp1 <- chainList[[k+1]]$sampFuncVect[length(chainList[[k+1]]$sampFuncVect)]
    T_k   <- control$tempVect[k]
    T_kp1 <- control$tempVect[k+1]
    a_swap <- exp(-(c_kp1-c_k)*(1/T_k-1/T_kp1))
    if(a_swap >= 1) {
      accept <- T
    } else {
      accept <- runif(1) <= a_swap
    }
   
    if(accept) {
      X_k     <- chainList[[k  ]]$X_mat[,ncol(chainList[[k  ]]$X_mat)]
      X_kp1   <- chainList[[k+1]]$X_mat[,ncol(chainList[[k+1]]$X_mat)]
      c_k     <- chainList[[k  ]]$sampFuncVect[length(chainList[[k  ]]$sampFuncVect)]
      c_kp1   <- chainList[[k+1]]$sampFuncVect[length(chainList[[k+1]]$sampFuncVect)]
      acc_k   <- chainList[[k  ]]$acceptVect[length(chainList[[k  ]]$acceptVect)]
      acc_kp1 <- chainList[[k+1]]$acceptVect[length(chainList[[k+1]]$acceptVect)]

      chainList[[k  ]]$X_mat[,ncol(chainList[[k  ]]$X_mat)] <- X_kp1
      chainList[[k+1]]$X_mat[,ncol(chainList[[k+1]]$X_mat)] <- X_k
      chainList[[k  ]]$sampFuncVect[length(chainList[[k  ]]$sampFuncVect)] <- c_kp1
      chainList[[k+1]]$sampFuncVect[length(chainList[[k+1]]$sampFuncVect)] <- c_k
      chainList[[k  ]]$acceptVect[length(chainList[[k  ]]$acceptVect)] <- acc_kp1
      chainList[[k+1]]$acceptVect[length(chainList[[k+1]]$acceptVect)] <- acc_k
    }

    swapMat[1,cc] <- numSampPrev + control$M*cc
    swapMat[2,cc] <- k
    if(accept) {
      swapMat[3,cc] <- 1
    } else {
      swapMat[3,cc] <- 0
    }

    if(control$verbose) {
      state1 <- unlist(lapply(chainList,function(x){x$costVect[length(x$costVect)]}))
      print(rbind(state0,state1))
      print(parTemperBest(chainList,chainsOnly=T))
    }
  }

  out <- list(chainList=chainList,control=control)
  if(newRun) {
    out$swapMat <- swapMat
  } else {
    out$swapMat <- rbind(init$swapMat,swapMat)
  }
  class(out) <- 'parTemper'
  return(out)
}
