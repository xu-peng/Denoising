# Experiments for Denoising

# Generate Q-matrix and Result matrix

# Notice the identifiability condition
library(CDM)

Q_init <- diag(3)

Q_vector_candidate <- as.matrix(expand.grid(c(0,1),c(0,1),c(0,1)))[-1,]

index <- sample(1:7, 7, replace = TRUE)

Q <- rbind(Q_init, Q_vector_candidate[index,])

profile_pattern <- as.matrix(expand.grid(c(0,1),c(0,1),c(0,1)))

profile_pattern <- profile_pattern[c(1,2,3,5,4,6,7,8),]
R <- sim.din(N=100, q.matrix = Q)

get_metrics <- function(predicted, ground_truth){
  TP <- sum(predicted[ground_truth==1])
  FP <- sum(predicted[ground_truth==0])
  FN <- sum(predicted[ground_truth==1]==0)
  
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  f1 <- 2 * precision * recall /(precision+recall)
  
  return(c("TP"= TP, "FP"= FP, "FN"=FN, "precision"=precision, "recall"=recall, "f1"=f1))
}

diagnose <- function(R_matrix, Q_matrix){
  fit <- din(R_matrix, Q_matrix)
  index <- apply(fit$posterior, 1, which.max)
  profile <- profile_pattern[index, ]
  profile
}


# filter DINA
filter_DINA <- function(R_matrix, Q_matrix){
  # R_matrix is the observed result matrix
  # Q_matrix is the item-skill matrix

  fit <- din(R_matrix, Q_matrix)
  pred <- round(fit$posterior %*% t(fit$pjk[,2,]))
  pred
}

pred_dina <- filter_DINA(R_matrix, Q_matrix)
profile_dina <- diagnose(pred_dina, Q_matrix)
metrics_dina <- get_metrics(profile_dina, P_matrix)
metrics_dina


# filter rasch
filter_rasch <- function(R_matrix){
  # See the help file of function rasch from ltm package
  fit<-rasch(R_matrix)
  beta <- fit$coefficients[,"beta"][1]
  X <- fit$coefficients[,"beta.i"]
  temp <- factor.scores(fit)[[1]]
  
  # Since factor.scores only gives all the response patterns, we need to find the index of patterns for every student  
  index <- apply(apply(R_matrix, 1, function(m) colSums(m == t(temp[,1:10]))), 2, which.max)
  Y <- temp[, "z1"][index]

  pred_rasch <- outer(X, Y, function(x, y) 1/(1+exp(-(beta * (y - x)))))
  t(round(pred_rasch))
}

pred_rasch <- filter_rasch(R_matrix)
profile_rasch <- diagnose(pred_rasch, Q_matrix)
metrics_rasch <- get_metrics(profile_rasch, P_matrix)
metrics_rasch



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

pred_mirt <- filter_mirt(R_matrix)
profile_mirt <- diagnose(pred_mirt, Q_matrix)
metrics_mirt <- get_metrics(profile_mirt, P_matrix)
metrics_mirt


# filter ideal response
all_profile_patterns <- expand.grid(c(0,1), c(0,1), c(0,1))

filter_IRP <- function(R_matrix, Q_matrix){
  IRP <- t(apply(all_profile_patterns, 1, function(m) apply(Q_matrix, 1, function(n) prod(as.vector(m)^n))))
  patternIndex <- apply(R_matrix, 1, function(m) which.max(colSums(m == t(IRP))))
  pred_IRP <- IRP[patternIndex, ]
  pred_IRP
}

profile_IRP <- diagnose(pred_IRP, Q_matrix)
# profile_IRP <- all_profile_patterns[patternIndex, ]
metrics_IRP <- get_metrics(profile_IRP, P_matrix)
metrics_IRP



profile_before <- diagnose(R_matrix, Q_matrix)
metrics_before <- get_metrics(profile_before, P_matrix)
metrics_before

# filter k-NN

filter_knn <- function(R_matrix, k = 3){
  index <- apply(R_matrix, 1, function(m) order(colSums(m == t(R_matrix)), decreasing = TRUE)[2:(k+1)])
  
  pred <- sapply(1:nrow(R_matrix), function(i) colSums(R_matrix[index[,i],])>floor(k/2))
  t(pred+0)
}

pred_knn <- filter_knn(R_matrix)
profile_knn <- diagnose(pred_knn, Q_matrix)
metrics_knn <- get_metrics(profile_knn, P_matrix)
metrics_knn


one_sweep <- function(R_matrix, Q_matrix, P_matrix){
  
  # Before filtering
  profile_before <- diagnose(R_matrix, Q_matrix)
  metrics_before <- get_metrics(profile_before, P_matrix)

  # Filter IRP
  pred_IRP <- filter_IRP(R_matrix, Q_matrix)
  profile_IRP <- diagnose(pred_IRP, Q_matrix)
  metrics_IRP <- get_metrics(profile_IRP, P_matrix)
  
  # Filter knn
  pred_knn <- filter_knn(R_matrix)
  profile_knn <- diagnose(pred_knn, Q_matrix)
  metrics_knn <- get_metrics(profile_knn, P_matrix)
  
  # Filter rasch
  pred_rasch <- filter_rasch(R_matrix)
  profile_rasch <- diagnose(pred_rasch, Q_matrix)
  metrics_rasch <- get_metrics(profile_rasch, P_matrix)
  
  # Filter DINA
  pred_dina <- filter_DINA(R_matrix, Q_matrix)
  profile_dina <- diagnose(pred_dina, Q_matrix)
  metrics_dina <- get_metrics(profile_dina, P_matrix)
  
  # Assembly the results
  result <- rbind(metrics_before, metrics_IRP, metrics_knn, metrics_rasch, metrics_dina)
  result
}

x <- one_sweep(R_matrix, Q_matrix, P_matrix)

# 3-skill, 100 student case
myList <- NULL
k <- 0
for (i in 1:10){
  index <- sample(1:7, 7, replace = TRUE)
  for (j in 1:10){
    k <- k + 1
    Q_matrix <- rbind(Q_init, Q_vector_candidate[index,])
    gen <- sim.din(N=100, q.matrix = Q_matrix)
    R_matrix <- gen$dat
    P_matrix <- gen$alpha
    myList[[k]] <- one_sweep(R_matrix, Q_matrix, P_matrix)
  }
}

# 3-skill, 1000 student case
myList_3_1000 <- NULL
k <- 0
for (i in 1:10){
  index <- sample(1:7, 7, replace = TRUE)
  for (j in 1:10){
    k <- k + 1
    Q_matrix <- rbind(Q_init, Q_vector_candidate[index,])
    gen <- sim.din(N=100, q.matrix = Q_matrix)
    R_matrix <- gen$dat
    P_matrix <- gen$alpha
    myList_3_1000[[k]] <- one_sweep(R_matrix, Q_matrix, P_matrix)
  }
}


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

# 1000 students
metrics_before <- t(sapply(myList_3_1000, function(m) m["metrics_before", ]))
metrics_IRP <- t(sapply(myList_3_1000, function(m) m["metrics_IRP", ]))
metrics_knn <- t(sapply(myList_3_1000, function(m) m["metrics_knn", ]))
metrics_rasch <- t(sapply(myList_3_1000, function(m) m["metrics_rasch", ]))
metrics_dina <- t(sapply(myList_3_1000, function(m) m["metrics_dina", ]))
