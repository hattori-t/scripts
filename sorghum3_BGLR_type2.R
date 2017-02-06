setwd("/Users/tomo/Dropbox/sorghum3")

## Type:2 Cross-validation with BGLR ##
# parameters
data <- commandArgs(trailingOnly = T)[1]
type <- commandArgs(trailingOnly = T)[2]
repeatNo <- commandArgs(trailingOnly = T)[3]

## data
geno <- read.csv("data/amat_GATK_all.csv", row.names = 1)
geno_hetero <- read.csv("data/dmat_GATK_all.csv", row.names = 1)
pheno <- read.csv(paste("data/",data,"_all.csv",sep=""), row.names=1)
original <- read.csv(paste("data/",data,"_",type,".csv",sep=""), row.names=1)

pheno <- pheno[,1:8]
original <- original[,1:8]

rownames(pheno) <- gsub("B2/","B2.",rownames(pheno))
rownames(pheno) <- gsub("B31/","B31.",rownames(pheno))
rownames(original) <- gsub("B2/","B2.",rownames(original))
rownames(original) <- gsub("B31/","B31.",rownames(original))


pheno_trim <- na.omit(pheno)
line <- intersect(rownames(pheno_trim),colnames(geno))
Pheno <- pheno_trim[line,]
geno_trim <- geno[line,line]
Geno <- t(geno_trim)
geno_trim_hetero <- geno_hetero[line,line]
Geno_hetero <- t(geno_trim_hetero)
phenolist <- colnames(Pheno)

# partition
CreateRandomPartition<-function(N, Nfold, Nrepeat){
  #N: number of lines
  #Nfold: number of folds of CV
  #Nrepeat: number of repeats of CV

  for(r in 1:Nrepeat){
    Partition <- sample(1:N, N, replace = F)
    Output <- paste(Nfold, "fold.N", N, ".repeat", r, ".txt", sep = "")
    print(write(c(Nfold, ceiling(N/Nfold)), Output, ncol = 2))
    Partition <- c(Partition, rep(-9, Nfold*ceiling(N/Nfold)-N))
    write(matrix(Partition, ncol = ceiling(N/Nfold), nrow = Nfold), Output, ncol = Nfold, append = TRUE)
  }
}

#CreateRandomPartition(nrow(Pheno),10,10) #make partition only one time!

Partition <- as.matrix(read.table(paste("data/partition/10fold.N", nrow(Pheno), ".repeat", repeatNo, ".txt", sep = ""), skip = 1))

#############
## BGLR
## Additive with Dominance (AD) ##
setwd(paste("/Users/tomo/Dropbox/sorghum3/type2/", data, "_", type, "/AD/fold", repeatNo, sep="" ))
Prediction.BGLR_AD <- function(Za, Zd, Pheno, Partition){

  Nl <- nrow(Pheno)
  stopifnot(Nl == nrow(Za))
  stopifnot(Nl == nrow(Zd))
  Ntrait <- ncol(Pheno)
  require(BGLR)

  Partition[Partition == -9] <- 0
  Nfold <- ncol(Partition)
  Predictions <- matrix(0, nc = Ntrait, nr = Nl)
  for(trait in 1:Ntrait){
    for (fold in 1:Nfold){
      cat("trait",trait,"fold",fold,"\n")
      Test <- Partition[,fold]
      train <- Pheno
      train[Test, trait] <- NA
      ETA <- list(Additive = list(K=Za, model = "RKHS"), Dominance = list(K=Zd, model = "RKHS"))
      Result <- BGLR(y = train[,trait], ETA = ETA, verbose = F)
      Predictions[Test,trait] <- as.vector(Result$yHat[Test])
    }
  }
  dimnames(Predictions) <- dimnames(Pheno)
  return(Predictions)
}

Predictedvalues.BGLR_AD <- Prediction.BGLR_AD(Geno, Geno_hetero, Pheno, Partition)

#line selecting
select <- intersect(rownames(original),rownames(Predictedvalues.BGLR_AD))
Predictedvalues.BGLR_AD <- Predictedvalues.BGLR_AD[select,]
Original <- original[select,]

#plot
cor_BGLR_AD <- NULL
rmse_BGLR_AD <- NULL
Ntrait <- ncol(Original)

for(trait in 1:Ntrait){
  pdf(paste("res_",data,"_",type,"_",repeatNo,"_", phenolist[trait], "_BGLR_AD.pdf", sep = ""))
  plot(Original[,trait], Predictedvalues.BGLR_AD[,trait], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_BGLR_AD",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Original[,trait], Predictedvalues.BGLR_AD[,trait], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Original[,trait] - Predictedvalues.BGLR_AD[,trait])^2) / length(Original[,trait]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_BGLR_AD <- rbind(cor_BGLR_AD, Core)
  rmse_BGLR_AD <- rbind(rmse_BGLR_AD,rmse)
  dev.off()
}

dimnames(Predictedvalues.BGLR_AD) <- dimnames(Original)
write.csv(Predictedvalues.BGLR_AD,paste("res_",data,"_",type,"_",repeatNo,"_Predictedvalues_BGLR_AD.csv", sep = ""))
rownames(cor_BGLR_AD) <- colnames(Original)
write.csv(cor_BGLR_AD,paste("res_",data,"_",type,"_",repeatNo,"_cor_BGLR_AD.csv", sep = ""))
rownames(rmse_BGLR_AD) <- colnames(Original)
write.csv(rmse_BGLR_AD,paste("res_",data,"_",type,"_",repeatNo,"_rmse_BGLR_AD.csv", sep = ""))


#######################
## Additive only (A) ##
setwd(paste("/Users/tomo/Dropbox/sorghum3/type2/", data, "_", type, "/A/fold", repeatNo, sep="" ))
Prediction.BGLR_A <- function(Za, Pheno, Partition){

  Nl <- nrow(Pheno)
  stopifnot(Nl == nrow(Za))
  Ntrait <- ncol(Pheno)
  require(BGLR)

  Partition[Partition == -9] <- 0
  Nfold <- ncol(Partition)
  Predictions <- matrix(0, nc = Ntrait, nr = Nl)
  for(trait in 1:Ntrait){
    for (fold in 1:Nfold){
      cat("trait",trait,"fold",fold,"\n")
      Test <- Partition[,fold]
      train <- Pheno
      train[Test, trait] <- NA
      ETA <- list(Additive = list(K=Za, model = "RKHS"))
      Result <- BGLR(y = train[,trait], ETA = ETA, verbose = F)
      Predictions[Test,trait] <- as.vector(Result$yHat[Test])
    }
  }
  dimnames(Predictions) <- dimnames(Pheno)
  return(Predictions)
}

Predictedvalues.BGLR_A <- Prediction.BGLR_A(Geno, Pheno, Partition)

#line selecting
select <- intersect(rownames(original),rownames(Predictedvalues.BGLR_A))
Predictedvalues.BGLR_A <- Predictedvalues.BGLR_A[select,]
Original <- original[select,]

#plot
cor_BGLR_A <- NULL
rmse_BGLR_A <- NULL
Ntrait <- ncol(Original)

for(trait in 1:Ntrait){
  pdf(paste("res_",data,"_",type,"_",repeatNo,"_", phenolist[trait], "_BGLR_A.pdf", sep = ""))
  plot(Original[,trait], Predictedvalues.BGLR_A[,trait], xlab = "Observed Value", ylab = "Predicted Value", main = paste(phenolist[trait],"_BGLR_A",sep = ""))
  abline(0, 1, lty = "dotted")
  Cor <- cor(Original[,trait], Predictedvalues.BGLR_A[,trait], use="pair")
  Core <- sprintf("%.2f", Cor)
  mse <- round(sum((Original[,trait] - Predictedvalues.BGLR_A[,trait])^2) / length(Original[,trait]), 2)
  rmse <- round(sqrt(mse), 2)
  legend("bottomright", legend = paste("r=", Core, " rmse=", rmse, sep = ""), bty="n")
  cor_BGLR_A <- rbind(cor_BGLR_A, Core)
  rmse_BGLR_A <- rbind(rmse_BGLR_A,rmse)
  dev.off()
}

dimnames(Predictedvalues.BGLR_A) <- dimnames(Original)
write.csv(Predictedvalues.BGLR_A,paste("res_",data,"_",type,"_",repeatNo,"_Predictedvalues_BGLR_A.csv", sep = ""))
rownames(cor_BGLR_A) <- colnames(Original)
write.csv(cor_BGLR_A,paste("res_",data,"_",type,"_",repeatNo,"_cor_BGLR_A.csv", sep = ""))
rownames(rmse_BGLR_A) <- colnames(Original)
write.csv(rmse_BGLR_A,paste("res_",data,"_",type,"_",repeatNo,"_rmse_BGLR_A.csv", sep = ""))
