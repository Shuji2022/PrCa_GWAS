library(fmsb)
library(rcompanion)
library(pROC)
library(ggplot2)


files <- dir(pattern = ".pheno$")
files <- files[order(nchar(files), files)]
data <- list()

#susceptibility
results <- data.frame(matrix(data=NA, nrow=length(files), ncol=5))
names(results) <- c("PRS.OR","PRS.CI","PRS.p","PRS.auc","PRS.auc.CI")
data <- list()


for(i in 1:length(files)){
  data[[i]] <- read.delim(files[i], header=TRUE, sep=" ")  
  names(data[[i]])[6] <- "SCORE"
  #data[[i]] <- na.omit(data[[i]])
  data[[i]]$newpheno <- as.factor(as.character(data[[i]]$newpheno))
  #calculate the PRS along
  mylogit <- glm(newpheno ~ SCORE, data = data[[i]], family = "binomial")
  results[i,"PRS.p"] <- formatC(summary(mylogit)$coefficients["SCORE","Pr(>|z|)"],format = "e",digits = 2)
  prob=predict(mylogit,type=c("response"))
  data[[i]]$prob.prs=prob
  
  ## Confidence intervals (CIs) using profiled log-likelihood
  results[i,"PRS.OR"] <- formatC(exp(summary(mylogit)$coefficients["SCORE","Estimate"]),format = "f",digits = 3)
  lower <- formatC(exp(confint(mylogit))["SCORE","2.5 %"],format = "f",digits = 3)
  upper <- formatC(exp(confint(mylogit))["SCORE","97.5 %"],format = "f",digits = 3)
  results[i,"PRS.CI"] <- paste0(lower, "-", upper)
  ngR2 <- nagelkerke(mylogit,restrictNobs = FALSE)
  results[i,"PRSonlyPseudoR2"] <- ngR2$Pseudo.R.squared.for.model.vs.null["Nagelkerke (Cragg and Uhler)",]  
}

auc.prs <- list()
auc.cov <- list()
auc.prs.cov <- list()
for(i in 1:length(files)){
  ##calculate the PRS along
  auc.prs[[i]] <- roc(newpheno ~ prob.prs, data = data[[i]])
  results[i,"PRS.auc"] <- round(auc(auc.prs[[i]]),digits = 3)
  lower <- formatC(ci(auc.prs[[i]])[1],format = "f",digits = 3)
  upper <- formatC(ci(auc.prs[[i]])[3],format = "f",digits = 3)
  results[i,"PRS.auc.CI"] <- paste0(lower, "-", upper)
}


#generate the OR plots
tiff("AUC.OR.tiff", width = 28, height = 20, units = "in",res=300)
par(mfcol=c(4,5), cex=1.5) 
for(i in 1:length(files)){
  title <- paste("p <",pval[i,1])
  plot(auc.prs[[i]], main = title )
  #plot(auc.prs.cov[[i]], add= TRUE, col="red")
}
dev.off()

results$name <- files
results
write.csv(results, "ORandAUC.csv")


