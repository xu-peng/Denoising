# Heuristic Binary Factorization

bifa <- function(M){
  M[M<=0] <- -1
  
  # Initialization
  krow <- 10
  kcol <- 3
  rowIndex <- sample(1:nrow(M), krow)
  colIndex <- sample(1:ncol(M), kcol)
  
  subM <- M[rowIndex,colIndex]
  
  # Kadane's algorithm
  
  # local search to reach optimum
  
}