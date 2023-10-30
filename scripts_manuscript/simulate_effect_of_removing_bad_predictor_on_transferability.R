
## Simulate the effect of removing discordant predictors on prediction accuracy
## (Reviewer 2 suggestion)

## Load necessary packages
library(MASS)
library(glmnet)
library(ggplot2)

## Set seed for reproducibility
set.seed(6495)

## Number of observations and variables
n<-1000
p<-10

## Initialize covariance mantrix filled with zeros
Sigma<-matrix(0, nrow=p, ncol=p)

## Fill "nearby" entries of correlation matrix
for(i in 1:p) {
  for(j in 1:p) {
    ## Check if nearby
    if (abs(i-j) <= 2) {
      ## Fill entry
      Sigma[i,j]<-0.5
    }
  }
}

## Set diagonals to 1
diag(Sigma)<-1

## Generate multivariate normal data
true_data<-mvrnorm(n, mu=rep(0,p), Sigma=Sigma)

## Create a dataset with a discordant variable
discordant_data<-true_data

## Get residuals of the last column after regressing out other columns
residuals<-lm(discordant_data[,p] ~ discordant_data[,-p])$residuals

## Create the discordant column in discordant_data as a linear combination of columns 1:3 and the residuals
discordant_data[,p]<-0.9 * rowMeans(discordant_data[,1:3]) + residuals

## Create response variable with some noise
y<-true_data %*% rnorm(p) + rnorm(n)

## Fit model
fit1<-cv.glmnet(true_data, y, alpha=0.1)

## Calcualte predictions
pred1<-predict(fit1, s="lambda.min", newx=true_data)
pred2<-predict(fit1, s="lambda.min", newx=discordant_data)

## Calculate root mean squared error
rmse1<-sqrt(mean((y - pred1)^2))
rmse2<-sqrt(mean((y - pred2)^2))

## Print RMSE values
print(paste("RMSE without discordant variable: ", round(rmse1, 3)))
print(paste("RMSE with discordant variable: ", round(rmse2, 3)))

## Print correlation difference
print(paste("Correlation difference: ", round(cor(y, pred1)-cor(y, pred2), 3)))

## Create data.frame to hold predictions
df<-data.frame(True = c(y, y), 
               Predicted = c(pred1, pred2), 
               Model = rep(c("Without discordant", "With discordant"), each=n))

## Plot true vs predicted values
ggplot(df, aes(x=True, y=Predicted, color=Model)) + 
  geom_point() + 
  geom_smooth(method="lm", se=FALSE) +
  ggtitle("True vs Predicted values") + 
  theme_minimal()
ggsave("figures/adj_transferability_3_extra_simulation.pdf", width=8, height=6)
