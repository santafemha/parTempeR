#' @title Parallel tempering for function minimization
#'
#' @author Michael Holton Price <MichaelHoltonPrice@@gmail.com>
#'
#' @export

saTemper <- function(costFunc,init,...,control=list()) {
  # saTemper needs three control variables to run
  #   tempVect
  #     new:    default < input
  #     cont:   must use original; error if tempVect is input
  #   M
  #     new:    default < input
  #     cont:   original < input
  #   numCycles 
  #     new:    default < input
  #     cont:   original < input
  #
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
  # verbose [metropVerbose in the saTemper control variable]
  #
  # The following control variables for the call to samMCMC are fixed:
  #
  # thinning = 1
  # numSampBurn = 0

  newRun <- !('saTemper' %in% class(init))
  # K is the number of temperatures
  # M is the samples between swap attempts
  # numCycles

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
    if(!('tempVect' %in% inputControl)) {
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
      print(saTemperBest(chainList,chainsOnly=T))
    }
  }

  out <- list(chainList=chainList,control=control)
  if(newRun) {
    out$swapMat <- swapMat
  } else {
    out$swapMat <- rbind(init$swapMat,swapMat)
  }
  class(out) <- 'saTemper'
  return(out)
}
