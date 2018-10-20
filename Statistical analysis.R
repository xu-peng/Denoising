# statistical analysis of multivariate discrete data

library(CDM)
library(MASS)

myData <- sim.dina
myData[] <- lapply(myData, factor)
fit <- mca(myData)

library(FactoMineR)
library(factoextra)
data(poison)
poison.active <- poison[1:55, 5:15]
for (i in 1:4) {
  plot(poison.active[,i], main=colnames(poison.active)[i],
       ylab="Count", col="steelblue", las=2)
}
res.mca <- MCA(poison.active, graph=TRUE)
print(res.mca)

# eig
eig.val <- get_eigenvalue(res.mca)
fviz_screeplot(res.mca, addlabels = TRUE, ylim=c(0,45))

# biplot
fviz_mca_biplot(res.mca, repel=TRUE, ggthme=theme_minimal())
