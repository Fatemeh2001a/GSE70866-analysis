############################################################################################################################################################
############################################################################################################################################################
## Example 1
##############################################################################
rm(list = ls())

library(bnlearn)

set.seed(123)

##############################################################################
## Parameters
##############################################################################

p <- 100
groups <- 2
n.group <- 50
n <- groups*n.group
upsilon <- 1
shape <- 1
rate <- 0.5
NREP <- 10

FP.BGe <- numeric(NREP)
FP.BGeCM <- numeric(NREP)
FP.Res <- numeric(NREP)

##############################################################################
## Repeat simulation
##############################################################################
for(rep in 1:NREP){
  
  ##############################
  ## Generate data
  ##############################
  
  psi <- 1/rgamma(p,shape=shape,rate=rate)
  
  B <- matrix(0,p,groups)
  
  for(i in 1:p){
    B[i,] <- rnorm(groups,0,sqrt(psi[i]))
  }
  
  X <- matrix(0,n,p)
  
  for(i in 1:p){
    X[1:n.group,i] <- B[i,1] +rnorm(n.group,0,sqrt(psi[i]))
    X[(n.group+1):100,i] <- B[i,2] +rnorm(n.group,0,sqrt(psi[i]))
  }
  
  ##############################################################################
  ## BGe
  ##############################################################################
  
  data.bge <- as.data.frame(X)
  
  names(data.bge) <- paste0("G",1:p)
  
  dag <- hc(data.bge,score="bge")
  
  FP.BGe[rep] <- narcs(dag)
  
  ##############################################################################
  ## BGeCM
  ##############################################################################
  
  Q <- model.matrix(~factor(rep(1:2,each=50))-1)
  
  J <- diag(n) -Q %*%solve(upsilon*diag(ncol(Q))+t(Q)%*%Q) %*%t(Q)
  
  L <- t(chol(J))
  
  X.star <- t(L)%*% X
  
  data.bgecm <- as.data.frame(X.star)
  
  names(data.bgecm) <- paste0("G",1:p)
  
  dag.bgecm <- hc(data.bgecm,score="bge")
  
  FP.BGeCM[rep] <- narcs(dag.bgecm)
  
  ##############################################################################
  ## Residual approach
  ##############################################################################
  
  QR <- qr(Q)
  
  P <- qr.Q(QR,complete=TRUE)[,(ncol(Q)+1):n]
  
  Y <- t(P)%*%X
  
  data.res <- as.data.frame(Y)
  
  names(data.res) <- paste0("G",1:p)
  
  dag.res <- hc(data.res,score="bge")
  
  FP.Res[rep] <- narcs(dag.res)
  
  cat(rep,"\n")
  
}

##############################################################################
## Results
##############################################################################

table=data.frame(Method=c("BGe","BGeCM","Residual"),
                 mean.FP=c(mean(FP.BGe),mean(FP.BGeCM),mean(FP.Res)),
                 sd.FP=c(sd(FP.BGe),sd(FP.BGeCM),sd(FP.Res)))
table
############################################################################################################################################################
############################################################################################################################################################
## Example 2
##############################################################################

rm(list = ls())
set.seed(123)

##############################################################################
## Parameters
##############################################################################
n <- 10          # observations
p <- 20          # variables
m <- 3           # columns of design matrix
K <- 100          # number of simulated datasets


## Variance parameters
psi <- numeric(p)
for(i in 1:18){
  psi[i] <- 1/rgamma(1,shape=1,scale=0.5)}
psi[19] <- 1/rgamma(1,shape=2,scale=0.5)
psi[20] <- 1/rgamma(1,shape=1.5,scale=0.5)


## Regression coefficients
gamma19 <- rnorm(2,mean=0,sd=sqrt(psi[19]))
gamma20 <- rnorm(1,mean=0,sd=sqrt(psi[20]))


## Random effects
B <- matrix(0,p,m)
for(i in 1:p){
  B[i,] <- rnorm(m,mean=0,sd=sqrt(psi[i]))}



## Design matrix Q
Q <- matrix(rnorm(n*m), nrow=n, ncol=m)
QB <- Q %*% t(B)


## Simulate datasets
DataList <- vector("list",K)
for(k in 1:K){
  X <- matrix(0,n,p)
  
  ## Variables 1-18
  for(i in 1:18){
    X[,i] <- QB[,i] +rnorm(n,mean=0,sd=sqrt(psi[i]))}
  
  ## Variable 19
  X[,19] <- QB[,19] +gamma19[1]*X[,1] +gamma19[2]*X[,2] +rnorm(n,mean=0,sd=sqrt(psi[19]))
  
  ## Variable 20
  X[,20] <- QB[,20] +gamma20*X[,19] +rnorm(n,mean=0,sd=sqrt(psi[20]))
  
  colnames(X) <- paste0("X",1:p)
  DataList[[k]] <- as.data.frame(X)
}


## Check
#length(DataList)
#dim(DataList[[1]])
#dim(DataList[[1]])
#(DataList[[1]])
##############################################################################
## True edge
##############################################################################

library(bnlearn)

data=DataList

dag=model2network("[X1][X2][X3][X4][X5][X6][X7][X8][X9][X10][X11][X12][X13][X14][X15][X16][X17][X18][X19|X1:X2][X20|X19]")

plot(dag)

##############################################################################
## BGe
##############################################################################
library(bnlearn)

edge.bge <- numeric(K)

correct.bge  <- numeric(K)

spurious.bge <- numeric(K)

for(k in 1:K){
  
  data.bge <- as.data.frame(DataList[[k]])
  
  dag.bge <- hc(data.bge, score = "bge", restart = 50, perturb = 5)
  
  edge.bge[k] <- narcs(dag.bge)
  
  ## Learned edges
  learn.arc <- arcs(dag.bge)
  
  ## True edges
  true.arc <- arcs(dag)
  
  ## Convert to character strings
  learn.edge <- paste(learn.arc[,1], learn.arc[,2], sep = "->")
  true.edge  <- paste(true.arc[,1],  true.arc[,2],  sep = "->")
  
  ## Number of correct edges
  correct.bge[k] <- sum(learn.edge %in% true.edge)
  
  ## Number of spurious edges
  spurious.bge[k] <- sum(!(learn.edge %in% true.edge))
}

##############################################################################
## BGeCM
##############################################################################

nu <- 1

H <- Q %*%solve(nu*diag(ncol(Q))+t(Q)%*%Q) %*%t(Q)

J <- diag(n)-H

L <- t(chol(J))

#بررسی فرضیات ماتریس J

#ابعاد
#dim(J)

#بررسی تقارن (مقدار خروجی درست نشان دهنده تقارن ماتریس است.)
#all.equal(J,t(J))

#بررسی معین مثبت بودن(در صورتی که تمام مقادیرویژه مثبت باشند ماتریس معین مثبت است.)
#eigen(J)$values

#نداشتن مقادیر گمشده(در صورتی که مقدار خروجی غلط باشد کماتریس مقادیر گمشده ندارد.)
#any(is.na(J))

#بررسی تجزیه
#all.equal(L%*%t(L),J)

edge.bgecm <- numeric(K)

correct.bgecm  <- numeric(K)

spurious.bgecm <- numeric(K)

for(k in 1:K){
  
  Xstar <- t(L) %*% as.matrix(DataList[[k]])
  
  data.bgecm <- as.data.frame(Xstar)
  
  dag.bgecm <- hc(data.bgecm, score = "bge", restart = 50, perturb = 5)
  
  edge.bgecm[k] <- narcs(dag.bgecm)
  
  ## Learned edges
  learn.arc <- arcs(dag.bgecm)
  
  ## True edges
  true.arc <- arcs(dag)
  
  ## Convert to character strings
  learn.edge <- paste(learn.arc[,1], learn.arc[,2], sep = "->")
  true.edge  <- paste(true.arc[,1],  true.arc[,2],  sep = "->")
  
  ## Number of correct edges
  correct.bgecm[k] <- sum(learn.edge %in% true.edge)
  
  ## Number of spurious edges
  spurious.bgecm[k] <- sum(!(learn.edge %in% true.edge))
}


##############################################################################
## Residual approach
##############################################################################


QR <- qr(Q)

P <- qr.Q(QR,complete=TRUE)[,(ncol(Q)+1):n]

#بررسی ویژگی اول: مساوی با ماتریس همانی باشد
#a=t(P)%*%P
#isTRUE(all.equal(a, diag(nrow(a))))

#بررسی ویژگی دوم: مساوی با صفر باشد
#b=t(P)%*%Q
#all.equal(unname(b), matrix(0, nrow(b), ncol(b)))

#بررسی ویژگی سوم 
#cc=P%*%t(P)
#d=diag(nrow(Q))-Q%*%solve(t(Q)%*%Q)%*%t(Q)
#all.equal(unname(cc), unname(d))

edge.res <- numeric(K)

correct.res  <- numeric(K)

spurious.res <- numeric(K)

for(k in 1:K){
  
  ## Residual data
  y <- t(P) %*% as.matrix(DataList[[k]])
  data.res <- as.data.frame(y)
  
  ## Learn Bayesian network
  dag.res <- hc(data.res, score = "bge", restart = 50, perturb = 5)
  
  edge.res[k] <- narcs(dag.res)
  
  ## Learned edges
  learn.arc <- arcs(dag.res)
  
  ## True edges
  true.arc <- arcs(dag)
  
  ## Convert to character strings
  learn.edge <- paste(learn.arc[,1], learn.arc[,2], sep = "->")
  true.edge  <- paste(true.arc[,1],  true.arc[,2],  sep = "->")
  
  ## Number of correct edges
  correct.res[k] <- sum(learn.edge %in% true.edge)
  
  ## Number of spurious edges
  spurious.res[k] <- sum(!(learn.edge %in% true.edge))
}

##############################################################################
## Results
##############################################################################
cat("BGe:\n",
    "Correct Edges= ",mean(correct.bge), " (", sd(correct.bge), ")\n",
    "Spurious Edges= ",mean(spurious.bge), " (", sd(spurious.bge), ")\n",
    "BGeCM:\n",
    "Correct Edges= ",mean(correct.bgecm), " (", sd(correct.bgecm), ")\n",
    "Spurious Edges= ",mean(spurious.bgecm), " (", sd(spurious.bgecm), ")\n",
    "Residual Approach:\n",
    "Correct Edges= ",mean(correct.res), " (", sd(correct.res), ")\n",
    "Spurious Edges= ",mean(spurious.res), " (", sd(spurious.res), ")\n")

table=data.frame(Method=c("BGe","BGeCM","Residual"),
                 edge.mean=c(mean(edge.bge),mean(edge.bgecm),mean(edge.res)),
                 Correct.Edges.mean=c(mean(correct.bge),mean(correct.bgecm),mean(correct.res)),
                 Correct.Edges.sd=c(sd(correct.bge),sd(correct.bgecm),sd(correct.res)),
                 Spurious.Edges.mean=c(mean(edge.bge),mean(edge.bgecm),mean(edge.res)),
                 Spurious.Edges.sd=c(sd(edge.bge),sd(edge.bgecm),sd(edge.res)))
table
##############################################################################
##رسم شبکه ها 
##############################################################################
par(mfrow = c(1,3))

plot(dag.bge,main = "BGe")

plot(dag.bgecm,main = "BGeCM")

plot(dag.res,main = "Residual approach")

# بازگشت به حالت پیش فرض
par(mfrow = c(1,1))
############################################################################################################################################################
############################################################################################################################################################
##Real Data
##############################################################################
#نصب و فراخوانی کتابخانه ها
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")

if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma")

if (!requireNamespace("Rgraphviz", quietly = TRUE))
  BiocManager::install("Rgraphviz")

if (!requireNamespace("bnlearn", quietly = TRUE))
library(GEOquery)
library(limma)
library(bnlearn)
library(Rgraphviz)

##############################################################################

##data
##############################################################################

rm(list = ls())

#فراخوانی داده ها
gse1=getGEO(filename="C:/Users/USER/Desktop/GSE70866-GPL14550_series_matrix.txt.gz", GSEMatrix = TRUE, getGPL = FALSE)
gse2=getGEO(filename="C:/Users/USER/Desktop/GSE70866-GPL17077_series_matrix.txt.gz", GSEMatrix = TRUE, getGPL = FALSE)

# تعداد پلتفرم ها
length(gse1)
length(gse2)

#نام ژن ها
featureData(gse1)
featureData(gse2)

#استخراج ماتریس بیان ژن
expr1=exprs(gse1)
expr2=exprs(gse2)   #این فایل داده های عددی نرمال شده را شامل می شود

#مشاهده ابعاد
dim(expr1)
dim(expr2)

#مشاهده قسمتی از داده ها
head(expr1)
head(expr2)

#خلاصه ای از داده ها
summary(expr1)
summary(expr2)

#استخراج اطلاعات نمونه ها
pheno1 <- pData(gse1)
pheno2 <- pData(gse2)

#مشاهده ستون ها
colnames(pheno1)
colnames(pheno2)


#انتخاب ژن های مشترک
identical(rownames(expr1), rownames(expr2))   #اگر خروجی TRUE باشد یعنی دو ماتریس دقیقا ژن های یکسان و با ترتیب یکسان دارند و آماده ادغام هستند
commonGenes=intersect(rownames(expr1), rownames(expr2))
length(commonGenes)   #هر دو ماتریس شامل 20330 ژن مشترک هستند
expr1=expr1[commonGenes, ]
expr2=expr2[commonGenes, ]


#تجمیع دو گروه
expr.all=cbind(expr1, expr2)
pheno=rbind(pheno1, pheno2)

#دیدن levelهای هر متغیر
sapply(pheno, function(x) length(unique(x)))
lapply(pheno,unique)

#بدیل به dataframe
Data=as.data.frame(t(expr.all))
dim(Data)


#حذف مقادیر گمشده
Data=na.omit(Data)
dim(Data)
dim(expr1)
dim(expr2)
n=nrow(Data)
p=ncol(Data)
##############################################################################
# نمونه گیری از کل ژن ها
##############################################################################

#انتخاب ژن های با بیشترین واریانس
gene.var=apply(Data, 2, var)
top30=names(sort(gene.var, decreasing = TRUE))[1:30]
Data.sample=Data[, top30]


##############################################################################
##BGe
##############################################################################

dag.bge=hc(Data.sample,score="bge", restart = 50, perturb = 5)

##############################################################################
##BGeCM
##############################################################################

colnames(pheno)[colnames(pheno) == "age:ch1"] <- "age"
colnames(pheno)[colnames(pheno) == "cohort:ch1"] <- "cohort"
colnames(pheno)[colnames(pheno) == "diagnosis:ch1"] <- "diagnosis"
colnames(pheno)[colnames(pheno) == "cell type:ch1"] <- "cell_type"
colnames(pheno)[colnames(pheno) == "sex (0=female, 1=male):ch1"] <- "sex"

sapply(pheno[, c("platform_id", "sex", "diagnosis", "cohort", "age", "cell_type")], 
       function(x) length(unique(na.omit(x))))


Q=model.matrix(~ age + sex, data = pheno)
Q=model.matrix(~ platform_id + cohort , data = pheno)
Q=model.matrix(~ platform_id + age + sex  + cohort , data = pheno)
head(Q)
dim(Q)

data=as.matrix(Data.sample)

n=nrow(data)

nu=1000

H=Q %*%solve(nu * diag(ncol(Q)) + t(Q) %*% Q) %*%t(Q)

J=diag(n) - H

L=t(chol(J))

Xstar=L %*% data

data.bgecm=as.data.frame(Xstar)

dag.bgecm=hc(data.bgecm,score = "bge", restart = 50, perturb = 5)


##############################################################################
## Residual approach
##############################################################################

QR=qr(Q)

P=qr.Q(QR,complete=TRUE)[,(ncol(Q)+1):n]

y=t(P) %*% data

data.res=as.data.frame(y)

dag.res=hc(data.res, score = "bge", restart = 50, perturb = 5)

##############################################################################
#result
##############################################################################

cat("Number of edges in BGe    :", narcs(dag.bge), "\n",
    "Number of edges in BGeCM  :", narcs(dag.bgecm), "\n",
    "Number of edges in Residual approach  :", narcs(dag.res), "\n")


table=data.frame(Method=c("BGe","BGeCM","Residual"),
                 edge=c(narcs(dag.bge),narcs(dag.bgecm),narcs(dag.res)))
table


##############################################################################
##رسم شبکه ها 
##############################################################################
par(mfrow = c(1,3))

plot(dag.bge,main = "BGe")

plot(dag.bgecm,main = "BGeCM")

plot(dag.res,main = "Residual approach")

# بازگشت به حالت پیش فرض
par(mfrow = c(1,1))
##############################################################################
##رسم شبکه ها
##############################################################################
library(Rgraphviz)

par(mfrow = c(1,3))

graphviz.plot(dag.bge, main = "BGe")

graphviz.plot(dag.bgecm,main = "BGeCM")

graphviz.plot(dag.res,main = "Residual approach")
par(mfrow = c(1,1))
############################################################################################################################################################
############################################################################################################################################################

