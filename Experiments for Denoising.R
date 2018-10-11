# Experiments for Denoising

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

# filter_IRP
filter_IRP <- function(R_matrix, Q_matrix, profile_patterns){
  IRP <- t(apply(profile_patterns, 1, function(m) apply(Q_matrix, 1, function(n) prod(as.vector(m)^n))))
  patternIndex <- apply(R_matrix, 1, function(m) which.max(colSums(m == t(IRP))))
  pred_IRP <- IRP[patternIndex, ]
  pred_IRP
}

# filter_knn
filter_knn <- function(R_matrix, k = 3){
  index <- apply(R_matrix, 1, function(m) order(colSums(m == t(R_matrix)), decreasing = TRUE)[2:(k+1)])
  pred <- sapply(1:nrow(R_matrix), function(i) colSums(R_matrix[index[,i],])>floor(k/2))
  t(pred+0)
}

# filter rasch
filter_rasch <- function(R_matrix){
  # See the help file of function rasch from ltm package
  fit<-rasch(R_matrix)
  beta <- fit$coefficients[,"beta"][1]
  X <- fit$coefficients[,"beta.i"]
  temp <- factor.scores(fit)[[1]]
  
  # Since factor.scores only gives all the response patterns, we need to find the index of patterns for every student  
  index <- apply(apply(R_matrix, 1, function(m) colSums(m == t(temp[,1:ncol(R_matrix)]))), 2, which.max)
  Y <- temp[, "z1"][index]

  pred_rasch <- outer(X, Y, function(x, y) 1/(1+exp(-(beta * (y - x)))))
  t(round(pred_rasch))
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
filter_DINA <- function(R_matrix, Q_matrix){
  # R_matrix is the observed result matrix
  # Q_matrix is the item-skill matrix
  
  fit <- din(R_matrix, Q_matrix)
  pred <- round(fit$posterior %*% t(fit$pjk[,2,]))
  pred
}

one_sweep <- function(R_matrix, Q_matrix, P_matrix, profile_pattern){
  
  # Before filtering
  profile_before <- diagnose(R_matrix, Q_matrix, profile_pattern)
  metrics_before <- get_metrics(profile_before, P_matrix)

  # Filter IRP
  pred_IRP <- filter_IRP(R_matrix, Q_matrix, profile_pattern)
  profile_IRP <- diagnose(pred_IRP, Q_matrix, profile_pattern)
  metrics_IRP <- get_metrics(profile_IRP, P_matrix)
  
  # Filter knn
  pred_knn <- filter_knn(R_matrix)
  profile_knn <- diagnose(pred_knn, Q_matrix, profile_pattern)
  metrics_knn <- get_metrics(profile_knn, P_matrix)
  
  # Filter rasch
  pred_rasch <- filter_rasch(R_matrix)
  profile_rasch <- diagnose(pred_rasch, Q_matrix, profile_pattern)
  metrics_rasch <- get_metrics(profile_rasch, P_matrix)
  
  # Filter DINA
  pred_dina <- filter_DINA(R_matrix, Q_matrix)
  profile_dina <- diagnose(pred_dina, Q_matrix, profile_pattern)
  metrics_dina <- get_metrics(profile_dina, P_matrix)
  
  # Assembly the results
  result <- rbind(metrics_before, metrics_IRP, metrics_knn, metrics_rasch, metrics_dina)
  result
}


getResult <- function(N=100, z=3, slip=0.1, guess=0.1, l=10){
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
      myList[[k]] <- one_sweep(R_matrix, Q_matrix, P_matrix, profile_pattern)
    }
  }
  myList
}

getDisturbedResult <- function(N=100, z=3, slip=0.1, guess=0.1, l=10, disturbed=1){
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
      myList[[k]] <- one_sweep(R_matrix, Q_matrix, P_matrix, profile_pattern)
    }
  }
  myList
}


# 3-skill case

myList <- getResult(N=1000, z=4, slip=0.3, guess=0.3, l=10)
myList <- getDisturbedResult(N=100, z=3, slip=0.1, guess=0.1, l=10, disturbed = 1)

# Separate results

metrics_before <- t(sapply(myList, function(m) m["metrics_before", ]))
metrics_IRP <- t(sapply(myList, function(m) m["metrics_IRP", ]))
metrics_knn <- t(sapply(myList, function(m) m["metrics_knn", ]))
metrics_rasch <- t(sapply(myList, function(m) m["metrics_rasch", ]))
metrics_dina <- t(sapply(myList, function(m) m["metrics_dina", ]))


# hypothesis test
t.test(metrics_before[,"precision"], metrics_IRP[,"precision"])
t.test(metrics_before[,"recall"], metrics_IRP[,"recall"])
t.test(metrics_before[,"f1"], metrics_IRP[,"f1"])

t.test(metrics_before[,"precision"], metrics_knn[,"precision"])
t.test(metrics_before[,"recall"], metrics_knn[,"recall"])
t.test(metrics_before[,"f1"], metrics_knn[,"f1"])

t.test(metrics_before[,"precision"], metrics_rasch[,"precision"])
t.test(metrics_before[,"recall"], metrics_rasch[,"recall"])
t.test(metrics_before[,"f1"], metrics_rasch[,"f1"])

t.test(metrics_before[,"precision"], metrics_dina[,"precision"])
t.test(metrics_before[,"recall"], metrics_dina[,"recall"])
t.test(metrics_before[,"f1"], metrics_dina[,"f1"])


