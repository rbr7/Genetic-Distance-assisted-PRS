library(data.table)
library(dplyr)
library(stringr)
library(rio)
library(optparse)

## define functions ##

Euclidean <- function(vec1, vec2, n) {
sqr.sum = 0
for (i in seq(n)) {sqr.sum <- sqr.sum + (vec1[i]-vec2[i])^2}
dis <- sqrt(sqr.sum)
return(dis)}

Euclidean.matrix <- function(mat, vec, n) {
N=nrow(mat)
sqr.sum = matrix(rep(0, N), nrow = N, byrow=TRUE)
for (i in seq(n)) {
PC <- as.double(vec[i])
sqr.sum <- sqr.sum + (mat[,i]-PC)^2}
dis <- sqrt(sqr.sum)
return(dis)}

shrinkage.vector <- function(med, n) {
N <- nrow(med)
dis.list <- vector()
for (i in seq(N*N)){
row = (i-1) %/% N +1
col = (i-1) %% N +1 
dis.list[i] <- Euclidean(med[row,], med[col,], n)
}
dis.list <- unlist(dis.list)
A <- matrix(dis.list, nrow=N, byrow=TRUE)
A.inv <- solve(A)
a <- A.inv %*% rep(1,N)
return(a)}
## finish defining function ##

## read options from command line ##
option_list = list(
  make_option(c("-m", "--med.file"), type="character", default=NULL,
              help="the file that contains the median of PC of the fine-tuning cohort. Have header. 1st row is cohort name",
              metavar="character"),
  make_option(c("-p", "--pca.file"), type="character", default=NULL,
              help="the file that contains PC of the testing cohort. Have header. 1st row is individual ID. Other rows are top PCs in ascending order.",
              metavar="character"),
  make_option(c("--prs.list"), type="character", default=NULL,
              help="the list of prs files. In the same order as in the med.file",
              metavar="character"),
  make_option(c("-A", "--a.list"), type="character", default=NULL,
              help="the list of weight. Separate items with ,", metavar="character"),
  make_option(c("-s", "--select.col"), type="character", default=NULL,
              help="the IID and PRS column in the prs files. Separate items with , ", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="the output name prefix", metavar="character"),
  make_option(c("--regress.PCA"), type="logical", default=TRUE,
              help="choose weather regress out PCA from the PRS before combining them", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

med.file <- opt$med.file
pca.file <- opt$pca.file
prs.list <- opt$prs.list
a.list <- opt$a.list
select.col <- opt$select.col
out <- opt$out
regress.PCA <- as.logical(opt$regress.PCA)

cat("\ncheck the input:\n")
Good.Input=TRUE

cat("\ncheck the number of PCA:\n")
Med <- data.frame(fread(med.file, header = T))
npca1 = ncol(Med) - 1
target.pca <- data.frame(fread(pca.file, header = T))
npca2 = ncol(target.pca) - 1

if(min(npca1) < 5) {
cat("\nFewer than 5 top PCs in med.file. May lead to inefficient genetic distance calculation! \n")}
if(min(npca2) < 5) {
cat("\nFewer than 5 top PCs in prs.file. May lead to inefficient genetic distance calculation! \n")}
npca = min(npca1, npca2)
print(paste("Calculating PRS based on", npca, "PCs"))
if ( npca < 5) {
cat("\nToo few top PCs. Please check the input \n")
Good.Input=F}

cat("\ncheck the number of fine-tuning cohorts\n")
med <- Med[, seq(2,npca+1)]
base.list <- as.matrix(Med)[, 1]
N <- length(base.list)
cat(paste("\nThe combined PRS based on PRS fine-tuned in", N, "cohorts.\n"))
print(base.list)
if (N<=1) {
cat("\nPlease provide 2 or more PRS to combine\n")
Good.Input=F}

target.pca <- target.pca[, seq(npca+1)]

if (is.null(prs.list)) {
cat("\nPlease provide list of PRS to combine\n")
Good.Input=F } else {
prs.list <- unlist(strsplit(prs.list, ","))}

cat("\ncheck the sele-defined shrinkage parameters\n")
if (is.null(a.list)) { A <- rep(1, N)} else {
A <- as.numeric(unlist(strsplit(a.list, ",")))}
if (any(is.na(A)) | any(A< 0) ) {
cat("Error: Unexpected value in the a.list. Please check input \n")
Good.Input=F}

if (length(prs.list)!=N | length(A)!=N) {
cat("Error: Inconsistent numbers of fine tuning cohorts. Please check input\n")
Good.Input=F}

cat("\ncheck the header to select from the PRS files\n")
select.col <- unlist(strsplit(select.col, ","))
if (length(select.col) == 2*N) {select.col <- select.col
} else if (length(select.col) == 2) {select.col <- rep(select.col, N)
} else {cat("\nPlease check the --selectcol. Incorrect input\n")
Good.Input=F}

if (is.na(regress.PCA)) {
cat("\n--regress.PCA has invalid value. The correct input should be TRUE, T or FALSE, F. Please check.\n" )
Good.Input=F}

cat("\n##### Finish checking the Input #####\n")

if (Good.Input==F) {
cat("The input is incorrect and the program is stopped. Please re-check your input!")
} else {
cat("\nStart to calculate...\n")

cat("\ngenerate the shrinkage parameter depends on the genetic distance of fine-tuning cohorts\n")
d <- shrinkage.vector(med, npca)

cat("\nCalculate genetic distance between fine-tuning cohorts and individuals in testing cohort\n")
colnames(target.pca) <- c("IID", paste0("PC", seq(npca1)))
mat <- target.pca[, 2:ncol(target.pca)]

dis.matrix <- data.frame(matrix(,nrow=nrow(target.pca),ncol=0))
for (i in seq(N)) {
dis <- Euclidean.matrix(mat, med[i,], npca)
dis.matrix <- cbind(dis.matrix, dis)}
colnames(dis.matrix) <- c(base.list)

cat("\nCalculate weight to combine PRS for individuals in testing cohort\n") 
w.matrix <- data.frame(matrix(,nrow=nrow(target.pca),ncol=0))
for (i in seq(N)) {
w <- d[i]*A[i]/dis.matrix[,i]
w.matrix <- cbind(w.matrix, w)}
w.sum <- rowSums(w.matrix)

IID <- target.pca$IID
a.matrix <- data.frame(matrix(,nrow=nrow(target.pca),ncol=0))
for (i in seq(N)) {
a <- as.numeric(w.matrix[,i]/w.sum)
a.matrix <- cbind(a.matrix, a)}
a.matrix <- cbind(IID, a.matrix)
colnames(a.matrix) <- c("IID", paste0(base.list, ".a"))
a.matrix$IID <- as.character(a.matrix$IID)

cat("\nGather PRS to combine\n")
prs.matrix <- fread(prs.list[1], select=select.col[c(1,2)])
colnames(prs.matrix) <- c("IID", base.list[1])
for (i in seq(2, N)) {
score <- fread(prs.list[i], select=select.col[c(2*i-1, 2*i)])
colnames(score) <- c("IID", base.list[i])
prs.matrix <- merge(prs.matrix, score)}
Nprs <- nrow(prs.matrix)
cat(paste("\nRead PRS for", Nprs, "individuals\n"))
Npca <- nrow(target.pca)
cat(paste("\nRead PCA for", Npca, "individuals\n"))

dat <- merge(prs.matrix, target.pca, by="IID")
Ndat <- nrow(dat)
cat(paste("\n", Ndat, "individuals have PRS and PCA information\n"))
  
if (Ndat < 0.4*min(Nprs, Npca)){
cat("\nWarning: the overlap between PRS and PCA information is too small. Please check input")}
if (Ndat==0) {
stop("No overlap between PRS and PCA information")} else {
cat("\n Calculate combined PRS\n")

if (regress.PCA){
IID <- dat$IID
prs.res.matrix <- data.frame(matrix(,nrow=nrow(dat),ncol=0))
for (i in seq(N)) {
prs.res <- scale(resid(glm(lm(formula = as.formula(paste(base.list[i], " ~", paste(paste0("PC", seq(10), collapse = "+")))),data = dat))))
prs.res.matrix <- cbind(prs.res.matrix, prs.res)}
prs.res.matrix <- cbind(IID, prs.res.matrix)
colnames(prs.res.matrix) <- c("IID", paste0(base.list, ".res"))
prs.res.matrix$IID <- as.character(prs.res.matrix$IID)

dat <- merge(prs.res.matrix, a.matrix, by= "IID") 
} else {
dat <- merge(prs.matrix, a.matrix, by= "IID")} 
print(summary(dat))

prs <- rep(0, nrow(dat))
for (i in seq(N)) {
print(paste(colnames(dat)[1+i], "x", colnames(dat)[1+N+i]))
prs <- prs + dat[, 1+i] * dat[, 1+N+i]
}
out.dat <- cbind(IID, prs)
colnames(out.dat) <- c("IID", "PRS")
export(out.dat, paste0(out, ".tsv.gz"), quote = F)
cat("\nPRS combined! \n  Y(^ W ^)Y  \n")}
cat(paste0("Result is written to ", out, ".tsv.gz\n"))
}
