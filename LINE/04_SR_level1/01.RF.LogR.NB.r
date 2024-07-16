#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
set.seed(1)
datpath<-args[1]
ver<-args[2]
strand<-args[3]
sub<-args[4]
masterpath<-args[5]
set.seed(1)

library('randomForest')
library('e1071')
library('glmnet')
### read in the RF, LogR and NB models ###
filename <- paste(as.character(masterpath),"/LINE/04_SR_level1/RFII_1.rds", sep="")
rf2 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/04_SR_level1/RFII_2.rds", sep="")
rf3 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/04_SR_level1/LogRII.rds", sep="")
logr <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/04_SR_level1/NBII.rds", sep="")
nb <- readRDS(filename)

    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_",as.character(strand),"/",as.character(sub),".sr.LINE.matrix", sep="")
    clone1 <- read.table(filename, sep="\t", header=T)
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_0/",as.character(sub),".sr.LINE.matrix", sep="")
    clone2 <- read.table(filename, sep="\t", header=T)
    clone  <- rbind(clone1, clone2)
    subclone <- subset(clone, (clone[,7] == 1 | clone[,8] == 1) & ref == 0)
    subclone$gap <- subclone$gap/(abs(subclone$refpos1 - subclone$refpos2))

    test2 <- subset(subclone, oldSR==1)
    nx2 <- test2[,c("seg", "map", "depth","map_size","map_ratio", "short","end3","end5","direction", "refpos", "dist", "A_pair", "A_insert", "A_mm", "A_MapQ", "A_AS", "A_XS")]
    pred.RF2 <- predict(rf2,newdata=nx2,type='prob')[,2]
    pred2 <- data.frame(cbind(as.character(test2$chr), test2$cord1, test2$cord2, as.character(test2$read), pred.RF2))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_",as.character(strand),"/LINE/",as.character(sub),".sr.pred.G1.txt", sep="")
    write(t(pred2), file=filename, ncol=5, sep="\t")

    test3 <- subset(subclone, oldSR==0)
    nx3 <- test3[,c("seg", "map", "depth","map_size","map_ratio", "short","end3","end5","direction", "refpos", "dist", "A_pair", "A_insert", "A_mm", "A_MapQ", "A_AS", "A_XS")]
    pred.RF3 <- predict(rf3,newdata=nx3,type='prob')[,2]
    pred3 <- data.frame(cbind(as.character(test3$chr), test3$cord1, test3$cord2, as.character(test3$read), pred.RF3))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_",as.character(strand),"/LINE/",as.character(sub),".sr.pred.G2.txt", sep="")
    write(t(pred3), file=filename, ncol=5, sep="\t")

    test4 <-subclone
    nx4 <- test4[,c("seg", "map", "depth","map_size","map_ratio", "short","end3","end5","direction", "refpos", "dist", "A_pair", "A_insert", "A_mm", "A_MapQ", "A_AS", "A_XS", "oldSR")]
    pred.RF4 <- predict(nb,newdata=nx4, type="raw")[,2]
    pred4 <- data.frame(cbind(as.character(test4$chr), test4$cord1, test4$cord2, as.character(test4$read), pred.RF4))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_",as.character(strand),"/LINE/",as.character(sub),".sr.pred.NB.txt", sep="")
    write(t(pred4), file=filename, ncol=5, sep="\t")

    pred.RF5 <- predict(logr,newdata=nx4, type="response")
    pred5 <- data.frame(cbind(as.character(test4$chr), test4$cord1, test4$cord2, as.character(test4$read), pred.RF5))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_",as.character(strand),"/LINE/",as.character(sub),".sr.pred.logR.txt", sep="")
    write(t(pred5), file=filename, ncol=5, sep="\t")
