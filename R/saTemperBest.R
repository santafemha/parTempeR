#' @title Get best value from call to saTemper
#'
#' @description Get the best value (the X with the lowest associated value of the cost function) from the output of saTemper. The optional input chainsOnly [default FALSE] allows a set of chains to be passed to saTemperBest rather than the output of saTemper and is used internally by saTemper.
#'
#' @param input The output of saTemper or a list of chains, chainList
#' @param chainsOnly An optional boolean indicating whether the input is a list of chains only [default FALSE]
#'
#' @return The parameter X with the lowest value of the cost function
#'
#' @author Michael Holton Price <MichaelHoltonPrice@@gmail.com>
#'
#' @export

saTemperBest <- function(input,chainsOnly=F) {
  if(chainsOnly) {
    chainList <- input
  } else {
    chainList <- input$chainList
  }
  bestChain <- which.min(unlist(lapply(chainList,function(x){min(x$sampFuncVect)})))
  theta <- chainList[[bestChain]]$X_mat[,which.min(chainList[[bestChain]]$sampFuncVect)]
  return(theta)
}
