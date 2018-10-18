# Anomaly Removal

library(CDM)
library(ltm)

# Utilities
get_all_combination <- function(k){
  z <- rep(0, k)
  result <- do.call(rbind, lapply(0:k, function(i) t(apply(combn(1:k,i), 2, function(k) {z[k]=1;z}))))
}

get_metrics <- function(predicted, ground_truth){
  TP <- sum(predicted[ground_truth==1])
  FP <- sum(predicted[ground_truth==0])
  FN <- sum(predicted[ground_truth==1]==0)
  
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  f1 <- 2 * precision * recall /(precision+recall)
  
  return(c("TP"= TP, "FP"= FP, "FN"=FN, "precision"=precision, "recall"=recall, "f1"=f1))
}

diagnose <- function(R_matrix, Q_matrix, profile_pattern){
  fit <- din(R_matrix, Q_matrix)
  index <- apply(fit$posterior, 1, which.max)
  profile <- profile_pattern[index, ]
  profile
}

# filter rasch
filter_rasch <- function(R_matrix, threshold){
  # See the help file of function rasch from ltm package
  fit<-rasch(R_matrix)
  beta <- fit$coefficients[,"beta"][1]
  X <- fit$coefficients[,"beta.i"]
  temp <- factor.scores(fit)[[1]]
  
  # Since factor.scores only gives all the response patterns, we need to find the index of patterns for every student  
  index <- apply(apply(R_matrix, 1, function(m) colSums(m == t(temp[,1:ncol(R_matrix)]))), 2, which.max)
  Y <- temp[, "z1"][index]
  
  pred_rasch <- t(outer(X, Y, function(x, y) 1/(1+exp(-(beta * (y - x))))))
  before_removal <- round(pred_rasch)
  
  # Anomaly Removal
  index_removal <- (abs(pred_rasch-R_matrix) > threshold)
  pred_rasch[index_removal] <- NA
  after_removal <- round(pred_rasch)
  
  # Return result
  list(before_removal = before_removal, after_removal = after_removal)
}

# filter mirt
filter_mirt <- function(R_matrix, k=1){
  mirt.fit<-mirt(R_matrix, model = k, itemtype = 'Rasch', method="EM")
  coeff<-coef(mirt.fit)
  coeff<-head(unlist(coeff),(3+k)*ncol(R_matrix))
  coeff<-matrix(coeff, nrow=3+k, ncol=ncol(R_matrix))
  coeff<-coeff[1:(k+1),]
  X <- coeff[2,]
  
  score<-fscores(mirt.fit)
  Y <- as.vector(score)
  # score<-cbind(score,rep(1,nrow(R_matrix)))
  #pred_IRT<-1/(1+exp(-(score%*%coeff)))
  #pred_IRT <- round(pred.IRT)
  #pred_IRT
  pred_IRT <- outer(X, Y, function(x, y) 1/(1+exp(-(y+x))))
  t(round(pred_IRT))
}

# filter DINA
filter_DINA <- function(R_matrix, Q_matrix, threshold){
  # R_matrix is the observed result matrix
  # Q_matrix is the item-skill matrix
  
  fit <- din(R_matrix, Q_matrix)
  pred <- fit$posterior %*% t(fit$pjk[,2,])
  before_removal <- round(pred)
  
  # Anomaly Removal
  index_removal <- (abs(pred-R_matrix) > threshold)
  pred[index_removal] <- NA
  after_removal <- round(pred)
  
  # Return result
  list(before_removal = before_removal, after_removal = after_removal)
}

one_sweep <- function(R_matrix, Q_matrix, P_matrix, profile_pattern, threshold){
  
  # Before filtering
  profile_before <- diagnose(R_matrix, Q_matrix, profile_pattern)
  metrics_before <- get_metrics(profile_before, P_matrix)
  
  # Filter rasch
  rasch_result <- filter_rasch(R_matrix, threshold)
    # Before removal
  pred_rasch <- rasch_result$before_removal
  profile_rasch <- diagnose(pred_rasch, Q_matrix, profile_pattern)
  metrics_rasch_before_removal <- get_metrics(profile_rasch, P_matrix)
    # After removal
  pred_rasch <- rasch_result$after_removal
  profile_rasch <- diagnose(pred_rasch, Q_matrix, profile_pattern)
  metrics_rasch_after_removal <- get_metrics(profile_rasch, P_matrix)
  
  # Filter DINA
  dina_result <- filter_DINA(R_matrix, Q_matrix, threshold)
    # Before removal
  pred_dina <- dina_result$before_removal
  profile_dina <- diagnose(pred_dina, Q_matrix, profile_pattern)
  metrics_dina_before_removal <- get_metrics(profile_dina, P_matrix)
    # After removal
  pred_dina <- dina_result$after_removal
  profile_dina <- diagnose(pred_dina, Q_matrix, profile_pattern)
  metrics_dina_after_removal <- get_metrics(profile_dina, P_matrix)
  
  # Assembly the results
  result <- rbind(metrics_before, metrics_rasch_before_removal, metrics_rasch_after_removal, metrics_dina_before_removal, metrics_dina_after_removal)
  result
}


getResult <- function(N=100, z=3, slip=0.1, guess=0.1, l=10, threshold = 0.8){
  # N - number of students
  # z - number of latent skills
  # slip - slip
  # guess - guess
  # l - number of items
  
  # Q_init is the first z rows that keeps the DINA model identifiable
  Q_init <- diag(z)
  
  # Get the pool of q-vectors to choose from, remove the no-skill required item
  Q_vector_candidate <- get_all_combination(z)[-1,]
  
  # Get all profile patterns
  profile_pattern <- get_all_combination(z)
  
  # Initate the result list
  myList <- NULL
  
  # k is the index of list component, i.e the times of runs
  k <- 0
  for (i in 1:10){
    
    # Generate a random index to choose the q-vectors 
    index <- sample(1:(2^z-1), (l-z), replace = TRUE)
    for (j in 1:10){
      k <- k + 1
      
      # Generate Q_matrix, P_matrix, R_matrix
      Q_matrix <- rbind(Q_init, Q_vector_candidate[index,])
      gen <- sim.din(N=N, q.matrix = Q_matrix, guess=rep(guess, nrow(Q_matrix)), slip = rep(slip, nrow(Q_matrix)))
      R_matrix <- gen$dat
      P_matrix <- gen$alpha
      
      # Get result for the k-th run 
      myList[[k]] <- one_sweep(R_matrix, Q_matrix, P_matrix, profile_pattern, threshold)
    }
  }
  myList
}

getDisturbedResult <- function(N=100, z=3, slip=0.1, guess=0.1, l=10, disturbed=1, threshold = 0.1){
  # N - number of students
  # z - number of latent skills
  # slip - slip
  # guess - guess
  # l - number of items
  
  # Q_init is the first z rows that keeps the DINA model identifiable
  Q_init <- diag(z)
  
  # Get the pool of q-vectors to choose from, remove the no-skill required item
  Q_vector_candidate <- get_all_combination(z)[-1,]
  
  # Get all profile patterns
  profile_pattern <- get_all_combination(z)
  
  # Initate the result list
  myList <- NULL
  
  # k is the index of list component, i.e the times of runs
  k <- 0
  for (i in 1:10){
    
    # Generate a random index to choose the q-vectors 
    index <- sample(1:(2^z-1), (l-z), replace = TRUE)
    for (j in 1:10){
      k <- k + 1
      
      # Generate Q_matrix, P_matrix, R_matrix
      Q_matrix <- rbind(Q_init, Q_vector_candidate[index,])
      gen <- sim.din(N=N, q.matrix = Q_matrix, guess=rep(guess, nrow(Q_matrix)), slip = rep(slip, nrow(Q_matrix)))
      R_matrix <- gen$dat
      P_matrix <- gen$alpha
      
      # Get result for the k-th run
      availableIndex <- c(1:length(Q_matrix))[Q_matrix==0 | (rep(rowSums(Q_matrix),3)>1) & (Q_matrix==1)]
      disturbedIndex <- sample(availableIndex,disturbed)
      Q_matrix[disturbedIndex] <- 1 - Q_matrix[disturbedIndex]
      myList[[k]] <- one_sweep(R_matrix, Q_matrix, P_matrix, profile_pattern, threshold)
    }
  }
  myList
}


# 3-skill case

myList <- getResult(N=100, z=3, slip=0.2, guess=0.2, l=10, threshold = 0.3)
myList <- getDisturbedResult(N=100, z=3, slip=0.1, guess=0.1, l=10, disturbed = 1)

# Separate results

metrics_before <- t(sapply(myList, function(m) m["metrics_before", ]))
metrics_rasch_before_removal <- t(sapply(myList, function(m) m["metrics_rasch_before_removal", ]))
metrics_rasch_after_removal <- t(sapply(myList, function(m) m["metrics_rasch_after_removal", ]))
metrics_dina_before_removal <- t(sapply(myList, function(m) m["metrics_dina_before_removal", ]))
metrics_dina_after_removal <- t(sapply(myList, function(m) m["metrics_dina_after_removal", ]))


# hypothesis test
t.test(metrics_before[,"precision"], metrics_rasch_before_removal[,"precision"])
t.test(metrics_before[,"recall"], metrics_rasch_before_removal[,"recall"])
t.test(metrics_before[,"f1"], metrics_rasch_before_removal[,"f1"])

t.test(metrics_before[,"precision"], metrics_rasch_after_removal[,"precision"])
t.test(metrics_before[,"recall"], metrics_rasch_after_removal[,"recall"])
t.test(metrics_before[,"f1"], metrics_rasch_after_removal[,"f1"])

t.test(metrics_before[,"precision"], metrics_dina_before_removal[,"precision"])
t.test(metrics_before[,"recall"], metrics_dina_before_removal[,"recall"])
t.test(metrics_before[,"f1"], metrics_dina_before_removal[,"f1"])

t.test(metrics_before[,"precision"], metrics_dina_after_removal[,"precision"])
t.test(metrics_before[,"recall"], metrics_dina_after_removal[,"recall"])
t.test(metrics_before[,"f1"], metrics_dina_after_removal[,"f1"])

t.test(metrics_rasch_before_removal[,"precision"], metrics_rasch_after_removal[,"precision"])
t.test(metrics_rasch_before_removal[,"recall"], metrics_rasch_after_removal[,"recall"])
t.test(metrics_rasch_before_removal[,"f1"], metrics_rasch_after_removal[,"f1"])

t.test(metrics_dina_before_removal[,"precision"], metrics_dina_after_removal[,"precision"])
t.test(metrics_dina_before_removal[,"recall"], metrics_dina_after_removal[,"recall"])
t.test(metrics_dina_before_removal[,"f1"], metrics_dina_after_removal[,"f1"])
