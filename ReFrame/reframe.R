# Set working directory
path <- getwd()

# load packages
library(openxlsx)
library(pracma)

# fit function of main peak ratio depending on micro satellite length
data.sigmoid <- read.xlsx(paste(path, "MPR_fit.xlsx", sep="/"),1)
optim.logistic.lss <- function(data, par) with(data, sum((par[1] + ((1-par[1])/(1 + exp(-par[2] * (x - par[3])))))-y)^2)
m.s <- optim(par=c(0, -1, 9), optim.logistic.lss, data=data.sigmoid)
logistic.fun <- function(x, par) par[1] + ((1-par[1])/(1 + exp(-par[2] * (x - par[3]))))
invlogistic.fun <- function(y, par) (-1/par[2])*log(((1-par[1])/(y-par[1])-1)*(1/exp(par[2]*par[3])))

# plot(data.sigmoid$x, data.sigmoid$y, xlim=c(0, max(data.sigmoid$x)), ylim=c(-0.5, 2))
# x <- seq(0, max(data.sigmoid$x), length.out = 500)
# lines(x, logistic.fun(x, m.s$par), col="red")


#### GENERAL INPUT
# Define formula to adjust MS length            
formel <- function(x) logistic.fun(x, m.s$par)

# Inverse of formula
formel.invers <- function(y) invlogistic.fun(y, m.s$par)

#Read in list of references
ref.all <- read.xlsx(paste(path, "reflist.xlsx", sep="/"),1)

#List all files in folder "in"
files<-list.files(path = paste(path, "in", sep="/"), pattern = ".xlsx")


#### Run through files in "in"-folder
for (ii in 1:length(files)){
  act.file <- files[ii]
  print(act.file)
  gene <- gsub("_.*$","",act.file, perl = TRUE)
  act.file_results <- paste0(paste(path, "out/", sep="/"),gene,"_results.xlsx")
  
  # Get references for current gene
  ref <- ref.all[which(ref.all$Marker==gene),3:11]
  ref[is.na(ref)] <- 0
 
  # Convert reference into percentage values
  refP <- ref/sum(ref[!is.na(ref)])

  # Read in input file for gene
  data <- read.xlsx(paste(path, "in", act.file, sep="/"),"Heights", colNames=TRUE)

  # compute weights
  results <- data
  results[,3:11] <- t(apply(data, 1, get.weights, refP=refP))
  colnames(results) <- c("TU_ID","Lauf_ID","l4","l3","l2","l1","wt","r1","r2","r3","r4")


  #### OUTPUT
  results$gene <- gene
  results <- results[,c(1,12,2:11)]
  
  # Write original data to file
  write.xlsx(data, file=act.file_results, sheetName="Heights", row.names=FALSE, col.names=TRUE, showNA=FALSE)
  # Write results to same file
  write.xlsx(results, file=act.file_results, sheetName="Result_of_algorithm", row.names=FALSE, col.names=TRUE, append=TRUE, showNA=FALSE)

}#end of loop

###########################################################################################################################################
# function to set up matrix for linear equation system
setup.matrix <- function(refP){
  #Compute fictive length of microsatellite
  fictiv.length <- formel.invers(refP$MP)
  
  #Compute values for main diagonal in equation system
  #vec<-seq(-4,4,by=1)
  hd<-formel(fictiv.length+seq(-4,4,by=1))
  
  #Compute Cofactors
  S<-(1-hd)/(1-refP$MP)
  
  #Define matrix for equation system
  A <- matrix(c( 
    hd[1], refP$L1*S[2], refP$L2*S[3], refP$L3*S[4], refP$L4*S[5],            0,            0,            0,            0,
    refP$R1*S[1],        hd[2], refP$L1*S[3], refP$L2*S[4], refP$L3*S[5], refP$L4*S[6],            0,            0,            0,
    refP$R2*S[1], refP$R1*S[2],        hd[3], refP$L1*S[4], refP$L2*S[5], refP$L3*S[6], refP$L4*S[7],            0,            0,
    refP$R3*S[1], refP$R2*S[2], refP$R1*S[3],        hd[4], refP$L1*S[5], refP$L2*S[6], refP$L3*S[7], refP$L4*S[8],            0,
    refP$R4*S[1], refP$R3*S[2], refP$R2*S[3], refP$R1*S[4],        hd[5], refP$L1*S[6], refP$L2*S[7], refP$L3*S[8], refP$L4*S[9],
    0, refP$R4*S[2], refP$R3*S[3], refP$R2*S[4], refP$R1*S[5],        hd[6], refP$L1*S[7], refP$L2*S[8], refP$L3*S[9], 
    0,            0, refP$R4*S[3], refP$R3*S[4], refP$R2*S[5], refP$R1*S[6],        hd[7], refP$L1*S[8], refP$L2*S[9], 
    0,            0,            0, refP$R4*S[4], refP$R3*S[5], refP$R2*S[6], refP$R1*S[7],        hd[8], refP$L1*S[9], 
    0,            0,            0,            0, refP$R4*S[5], refP$R3*S[6], refP$R2*S[7], refP$R1*S[8],        hd[9]),
    nrow=9, ncol=9,byrow=T)
  
  return(A)
}


# function to solve linear equation system with contraints (weights in [0,1] and sum(weights)=1)
get.weights <- function(x, refP){
  heights<-as.numeric(x[3:11])                          # Get peak heights
  heights[is.na(heights)] <- 0                # Replace NA's by 0
  heights<-heights/sum(heights)               # Convert heights into procentage values
  
  A <- setup.matrix(refP)                     # Set up matrix for linear equation system
  result <- lsqlincon(C=A, d=as.numeric(heights), 
                      Aeq=matrix(rep(1, 9), nrow=1), beq=1, 
                      lb=rep(0, 9), ub=rep(1, 9))
  result
}

