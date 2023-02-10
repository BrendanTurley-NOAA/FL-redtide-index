# https://daviddalpiaz.github.io/r4sl/ensemble-methods.html#tuning-1

gc()

library(caret)
library(fields)
library(lubridate)
library(mgcv)
library(party)
library(partykit)
library(pROC)
library(randomForest)
library(scales)

setwd('~/Documents/nasa/data/lowres_4km')
# write.csv(habs_covar_agg,'habs_covariates_agg.csv',row.names = F)
habs_covar_agg <- read.csv('habs_covariates_agg.csv')
habs_covar_agg$date <- ymd(habs_covar_agg$date)

### random forest
# https://www.r-bloggers.com/2021/04/random-forest-in-r/
# https://stats.stackexchange.com/questions/41443/how-to-actually-plot-a-sample-tree-from-randomforestgettree
# https://cran.rstudio.com/web/packages/randomForestExplainer/vignettes/randomForestExplainer.html
set.seed(222)
ind <- sample(2, nrow(habs_covar_agg), replace = TRUE, prob = c(.55, .45))
# ind_rm <- c(1:3,6:8)
ind_rm <- c(1:8,22:25) # alternative removal of 'climatological' variables
names(habs_covar_agg)[ind_rm]
habs_covar_agg$pa100k <- as.factor(habs_covar_agg$pa100k)
train <- habs_covar_agg[ind==1,-ind_rm] # remove all superfluous variables
test <- habs_covar_agg[ind==2,-ind_rm]
rf <- randomForest(pa100k~., data=train, proximity=T, importance=T) # norm.votes = F for combining trees

print(rf)

hist(treesize(rf),main = "No. of Nodes for the Trees",col = "green")

# https://topepo.github.io/caret/measuring-performance.html
p1 <- predict(rf, train)
confusionMatrix(p1, train$pa100k,positive='1')
p2 <- predict(rf, test)
tabs <- addmargins(table(p2,test$pa100k))
tabs
error_mat <- confusionMatrix(p2, test$pa100k, positive='1', mode='everything')
error_mat
error_mat$byClass # F1 out of 1; https://en.wikipedia.org/wiki/F-score; https://en.wikipedia.org/wiki/Sensitivity_and_specificity
### 2022/12/07 - there is a high specificity (few false positives) and low sensitivity (many false negatives); the opposite of what is desired
tabs[1,2]/tabs[3,2] # FNR or 1 - sensitivity
tabs[2,1]/tabs[3,1] # FPR or 1 - specificity


######## logistic ######## 
log_mod <- glm(pa100k~., data=train, family='binomial')

log_preds <- predict(log_mod,train,type='response')

### ROC analysis
p3 <- predict(log_mod, train, type='response')
temproc <- roc(train$pa100k , p3, plot=TRUE, grid=TRUE)
# CALCULATE AREA UNDER THE CURVE
temproc$auc  

roctable <- cbind(temproc$sensitivities, temproc$specificities, temproc$thresholds, 
                  temproc$sensitivities+temproc$specificities)
Threshold <- roctable[roctable[,4] == max(roctable[,4]),][3]
Threshold 

p2.1 <- ifelse(log_preds>Threshold,1,0)
p2.1 <- as.factor(p2.1)
error_mat2 <- confusionMatrix(p2.1, train$pa100k, positive='1', mode='everything')
error_mat2
######## logistic ######## 


error_analysis <- data.frame(true=test$pa100k,prediction=p2)
error_analysis$diff <- ifelse(error_analysis$true==error_analysis$prediction,1,0)
errors <- cbind(error_analysis,test)
errors <- errors[,-c(1:2,17)]
errors <- cbind(error_analysis,test)
error_tree <- partykit::ctree(diff~., data=errors)

er_cl <- rep(NA,nrow(error_analysis))
er_cl[which(error_analysis$true==0 & error_analysis$prediction==0)] <- 'TN'
er_cl[which(error_analysis$true==0 & error_analysis$prediction==1)] <- 'FP'
er_cl[which(error_analysis$true==1 & error_analysis$prediction==0)] <- 'FN'
er_cl[which(error_analysis$true==1 & error_analysis$prediction==1)] <- 'TP'
errors <- cbind(error_analysis,er_cl,test)
errors <- errors[,-c(1:3,18)]
errors$er_cl <- as.factor(errors$er_cl)
error_tree <- partykit::ctree(er_cl~., data=errors)

setwd('~/Downloads')
pdf('errors_box.pdf',width=7,height=11,pointsize = 12,useDingbats = T)
par(mfrow=c(2,2))
for(i in 2:15){
  boxplot(errors[,i]~errors$er_cl,ylab=paste(names(errors)[i]))
}
dev.off()

error_tree
setwd('~/Downloads')
png('error_tree_3.png',width=80,height=15,units='in',res=300)
plot(error_tree,type='simple')
dev.off()

plot(error_tree,type='simple')
node_boxplot(error_tree)

library(rpart)
library(rpart.plot)

tree <- rpart(survived~., data=TitanicData, cp=.02)
error_tree2 <- rpart(diff~., data=errors, cp=.001)
error_tree2 <- rpart(er_cl~., data=errors, cp=.002)

plot(error_tree2)
summary(error_tree2)

setwd('~/Downloads')
png('error_tree_4.png',width=20,height=15,units='in',res=300)
rpart.plot(error_tree2)
dev.off()


fn <- errors[which(errors$er_cl=='FN'),]

plot(rf,log='y')
legend('topright',c('OOB','Neg','Pos'),col=c(1,2,3),lty=1)

### ROC analysis
p3 <- predict(rf, test, type='prob')
temproc <- roc(test$pa100k , p3[,2], plot=TRUE, grid=TRUE)
# CALCULATE AREA UNDER THE CURVE
temproc$auc  

yr <- 2005
subset <- test[which(test$year==yr ),]
phat1 <- p3[,2][which(test$year==yr )]

par(mar=c(5,5,1,6))
plot(subset$LONGITUDE,subset$LATITUDE,asp=1)
quilt.plot(subset$LONGITUDE,subset$LATITUDE,phat1,col=plasma(60),asp=1,add=T)


var_imp <- importance(rf,scale=F)
var_imp
par(mar=c(4,7,1,1))
barplot(sort(var_imp[,1]),las=1,horiz=T,main='Importance (absence)')
barplot(sort(var_imp[,2]),las=1,horiz=T,main='Importance (presence)')
barplot(t(var_imp[,1:2]),las=1,horiz=T,col=c('gray20','gray80'),beside=F)
legend('bottomright',c('absence','presence'),fill=c('gray20','gray80'),bty='n')
plot(var_imp[,3],var_imp[,4],typ='n',xlab='MeanDecreaseAccuracy',ylab='MeanDecreaseGini')
text(var_imp[,3],var_imp[,4],row.names(var_imp),cex=.8)

png('rti_rf_varimp.png',width=9,height=7,pointsize=12,unit='in',res=300)
varImpPlot(rf,
           sort = T, scale = F,
           main = "Variable Importance")
dev.off()

var_imp2 <- varImp(rf,scale=F)
barplot(sort(var_imp2[,2]),names.arg = rownames(var_imp2)[order(var_imp2[,2])],las=2,horiz=T)

par(mfrow=c(2,2))
plot(var_imp[,3],var_imp[,4],typ='n',xlab='MeanDecreaseAccuracy',ylab='MeanDecreaseGini')
text(var_imp[,3],var_imp[,4],row.names(var_imp),cex=.8)

plot(var_imp[,3],var_imp2[,2],typ='n',xlab='MeanDecreaseAccuracy',ylab='VarImp (caret)')
text(var_imp[,3],var_imp2[,2],row.names(var_imp),cex=.8)

plot(var_imp2[,2],var_imp[,4],typ='n',xlab='VarImp (caret)',ylab='MeanDecreaseGini')
text(var_imp2[,2],var_imp[,4],row.names(var_imp),cex=.8)

ind_var <- rownames(var_imp)[order(var_imp[,4],decreasing=T)]
# pdf('rti_rf_partial.pdf',width=9,height=7,pointsize=12)
par(mfrow=c(3,3),mar=c(4,4,1,1))
for (i in seq_along(ind_var)) {
  partialPlot(rf, train, ind_var[i], xlab=ind_var[i],
              main=paste("Partial Dependence on", ind_var[i]),rug=T)
}
# dev.off()




set.seed(222)
ind <- sample(2, nrow(habs_covar_agg), replace = TRUE, prob = c(.55, .45))
names(habs_covar_agg)[c(1:3,6:8)]
habs_covar_agg$pa100k <- as.factor(habs_covar_agg$pa100k)
train <- habs_covar_agg[ind==1,-c(1:3,6:8)] # remove all superfluous variables
test <- habs_covar_agg[ind==2,-c(1:3,6:8)]
# rf2 <- cforest(pa100k~., data=train) # apparently this take awhile
rf2 <- cforest(pa100k~ LATITUDE + LONGITUDE + chlor_a + chl_anom + nflh + nflh_anom + rrs_667 + ssnlw488 + carder_bbp + morel_bbp + cm_bbp + abi + rbd + kbbi + sst + year + month + yday + week + pa100k + depth_m,
               data=train)

plot(rf2)

# https://topepo.github.io/caret/measuring-performance.html
# https://stats.stackexchange.com/questions/342650/cforest-runs-out-of-ram-when-running-predict-function
p1 <- predict(rf2)
confusionMatrix(p1, train$pa100k, positive='1')
p2 <- predict(rf2, test, OOB = TRUE)
tabs <- addmargins(table(p2,test$pa100k))
tabs
error_mat <- confusionMatrix(p2, test$pa100k, positive='1', mode='everything')
error_mat
error_mat$byClass # F1 out of 1; https://en.wikipedia.org/wiki/F-score; https://en.wikipedia.org/wiki/Sensitivity_and_specificity
### 2022/12/07 - there is a high specificity (few false positives) and low sensitivity (many false negatives); the opposite of what is desired
tabs[1,2]/tabs[3,2] # FNR or 1 - sensitivity
tabs[2,1]/tabs[3,1] # FPR or 1 - specificity


### ROC analysis
p3 <- predict(rf2, test, type='prob')
temproc <- roc(test$pa100k , p3[,2], plot=TRUE, grid=TRUE)
# CALCULATE AREA UNDER THE CURVE
temproc$auc  

yr <- 2005
subset <- test[which(test$year==yr),]
phat1 <- p3[,2][which(test$year==yr)]

par(mar=c(5,5,1,6))
plot(subset$LONGITUDE,subset$LATITUDE,asp=1)
quilt.plot(subset$LONGITUDE,subset$LATITUDE,phat1,col=plasma(60),asp=1,add=T)

boxplot(phat1~subset$month,varwidth=T)

boxplot(p3[,2]~test$month)
boxplot(p3[,2]~test$year,las=2)
boxplot(p3[,2]~test$week,las=2)


varimp(rf2)
# LATITUDE  LONGITUDE    chlor_a   chl_anom       nflh  nflh_anom    rrs_667 
# 0.05664329 0.02495758 0.01955025 0.02565890 0.04490371 0.01369197 0.01595241 
# ssnlw488 carder_bbp  morel_bbp     cm_bbp        abi        rbd       kbbi 
# 0.06584009 0.01429602 0.03173409 0.01395783 0.04463399 0.05242650 0.01081208 
# sst       year      month       yday       week    depth_m 
# 0.04396409 0.03376888 0.01981681 0.02185371 0.02037945 0.05045956

imp <- as.numeric(unlist(strsplit(paste('0.05664329 0.02495758 0.01955025 0.02565890 0.04490371 0.01369197 0.01595241 0.06584009 0.01429602 0.03173409 0.01395783 0.04463399 0.05242650 0.01081208 0.04396409 0.03376888 0.01981681 0.02185371 0.02037945 0.05045956',sep=','),' ')))
names(imp) <- names(train[-20])

par(mar=c(4,6,1,1))
barplot(imp[order(imp)],horiz = T,names.arg = names(imp)[order(imp)],las=2)


