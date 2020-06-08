############################Working directory and library loading####################
setwd("C:/xampp/htdocs/asrpro/rcode") # Set the working directory of created folder
library(e1071)
library(Biostrings)
library(R2HTML)
if (file.exists("asrpro.txt")) file.remove("asrpro.txt")
x <- readAAStringSet("test.fasta")
nx <- as.character(names(x))
xx <- toupper(as.character(as.character(x)))
#Checking of standard residues#
std <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
zx <- xx
z <-  sapply(zx, function(s) strsplit(s, split=""))
zz <- lapply(z, table)
zz1 <- unique(as.character(unlist(lapply(zz, names))))
if(length(union(zz1, std))>20){
pp <- data.frame("Dataset contain non-standard residues,so kindly submit sequences having standard residues only")
write.table(pp,"asrpro.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
}else{
###########################################ACC features########################
tst <- readAAStringSet("test.fasta")
setwd("C:/xampp/htdocs/asrpro/rcode/pscdnc") # Set the working directory of pscdnc
writeXStringSet(tst,"file name.fasta") # User has to supply the test sequences in FASTA format with name 

shell("acc.py test.fasta test_acc_3.txt Protein ACC -lag 3 -f tab", intern=TRUE)
shell("acc.py test.fasta test_acc_4.txt Protein ACC -lag 4 -f tab", intern=TRUE)
shell("acc.py test.fasta test_acc_10.txt Protein ACC -lag 10 -f tab", intern=TRUE)

ac3 <-read.table("test_acc_3.txt") 
ac4 <- read.table("test_acc_4.txt") 
ac10 <-read.table("test_acc_10.txt") 
 
setwd("C:/xampp/htdocs/asrpro/rcode") # Set the working directory of created folder
write.table(names(tst),"name_tst.txt",row.names=F, col.names=F)
write.table(ac3,"test_acc_3.txt", sep="\t", col.names=FALSE, row.names=FALSE)
write.table(ac4,"test_acc_4.txt", sep="\t", col.names=FALSE, row.names=FALSE)
write.table(ac10,"test_acc_10.txt", sep="\t", col.names=FALSE, row.names=FALSE)
#################################################PTM Features###############################################

if (file.exists("output.csv")) file.remove("output.csv")
if (file.exists("ptm_test.txt")) file.remove("ptm_test.txt")


system("ModPred_win64.exe -i test.fasta -o output.csv")

x <- read.csv("output.csv")
kk <- as.character(x[,1])
zz <- as.character(sapply(kk, function(s)substr(s,1,1)))
mm <- which(zz==">")
str <- c(1,mm+1)
edn <- c(mm-1,length(kk))
pp <- c("Acetylation","ADP-ribosylation","Amidation","Carboxylation","C-linked_glycosylation","Disulfide_linkage","GPI_anchor_amidation","Hydroxylation","Methylation","N-linked_glycosylation","N-terminal_acetylation","O-linked_glycosylation","Palmitoylation","Phosphorylation","Proteolytic_cleavage","PUPylation","Pyrrolidone_carboxylic_acid","Sulfation","SUMOylation","Ubiquitination","Farnesylation","Geranylgeranylation","Myristoylation")
my_ptm <- matrix(0,nrow=length(str), ncol=23)

for(i in 1:length(str)){
s1 <- str[i]
s2 <- edn[i]
ss <- kk[s1:s2]
my_mat <- function(pm){ unlist(strsplit(as.character(pm), split="\t"))[c(2,5)]}
zm <- t(sapply(ss,my_mat))
zmd <- as.matrix(zm)
zmm <- table(as.character(zmd[ which(as.numeric(zmd[,2])>=0.5),][,1]))
zn <- names(zmm)
my_ptm[i,which(pp%in%zn)]<-as.numeric(zmm)
}
write.table(my_ptm,"ptm_test.txt", row.names=FALSE, col.names=FALSE, sep="\t")


####################################################################end###########################################
len <- width(tst)

tst_ptm <- read.table("ptm_test.txt")[,1:20]*(100/len)
tst_acc_3 <- read.table("test_acc_3.txt")
tst_acc_4 <- read.table("test_acc_4.txt")
tst_acc_10 <- read.table("test_acc_10.txt")


###################reading and combining features##########################

zz_n_1 <- cbind(tst_acc_3, tst_ptm)
zz_n_2 <- cbind(tst_acc_4, tst_ptm)
zz_n_3 <- cbind(tst_acc_10, tst_ptm)

########################End############################################

res <- matrix(0,nrow=length(tst), ncol=8)
colnames(res)<- c("COLD","DROUGHT","HEAT","LIGHT","OXIDATIVE","PHYTO","SALT","WOUND")
ksm <- c("cold","drought","heat","light","oxidative","phyto","salt","wound")
atr_pp <- c("Cold","Drought","Heat","Light","Oxidative","Phyto","Salt","Wound")

tsm <- c("zz_n_2", "zz_n_1", "zz_n_3","zz_n_3","zz_n_1","zz_n_3","zz_n_1","zz_n_3")

for(j in 1:8){

ann_trn <- readRDS(paste(ksm[j],".rds", sep=""))
pred_p <- predict (ann_trn, newdata=eval(parse(text=tsm[j])), probability=TRUE)
kk <- round(attr(pred_p,"prob")[,"Y"], 3)
res[,j]<- paste(atr_pp[j],"(",kk,")", sep="")
}
###############################result processing#########################
kk <- res
zst <- as.character(read.delim("name_tst.txt", header=FALSE)[,1])
sink("asrpro.txt")
cat("===================================================================================\n")
cat("Probability of sequence being predicted in different stress classes. Value inside brackets denotes probability of prediction.\n")
cat("===================================================================================\n")
for(i in 1:length(zst)){
cat("-----------------------------------------------------------------------------------\n")
#cat(paste(rep("-",nchar(zst[i])),collapse=""), sep="\n")
cat("[",i,"]",">",zst[i], sep="")
cat("\n")
#cat(paste(rep("-",nchar(zst[i])),collapse=""), sep="\n")
cat("-----------------------------------------------------------------------------------\n")
for(j in 1:8){
writeLines(as.character(kk[i,j]))
}

}
sink()

###################################END#######################################
}

