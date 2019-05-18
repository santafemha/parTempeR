#' @title Get best value from call to parTemper
#'
#' @description Get the best value (the X with the lowest associated value of the cost function) from the output of parTemper. The optional input chainsOnly [default FALSE] allows a set of chains to be passed to parTemperBest rather than the output of parTemper and is used internally by parTemper.
#'
#' @param input The output of parTemper or a list of chains, chainList
#' @param chainsOnly An optional boolean indicating whether the input is a list of chains only [default FALSE]
#'
#' @return The parameter X with the lowest value of the cost function
#'
#' @author Michael Holton Price <MichaelHoltonPrice@@gmail.com>
#'
#' @export

parTemperBest <- function(input,chainsOnly=F) {
  if(chainsOnly) {
    chainList <- input
  } else {
    chainList <- input$chainList
  }
  bestChain <- which.min(unlist(lapply(chainList,function(x){min(x$sampFuncVect)})))
  theta <- chainList[[bestChain]]$X_mat[,which.min(chainList[[bestChain]]$sampFuncVect)]
  return(theta)
}
