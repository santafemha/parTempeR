#' @title Get best value from call to saTemper
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
