#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
set.seed(1)
datpath<-args[1]
ver<-args[2]
strand<-args[3]
sub<-args[4]
masterpath<-args[5]

library('randomForest')
library('e1071')
library('glmnet')

### read in the RF, LogR and NB models ###
filename <- paste(as.character(masterpath),"/ALU/04_SR_level1/RFIV.rds", sep="")
rf1 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/ALU/04_SR_level1/LogRIV.rds", sep="")
logr <- readRDS(filename)
filename <- paste(as.character(masterpath),"/ALU/04_SR_level1/NBIV.rds", sep="")
nb <- readRDS(filename)

    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_0/",as.character(sub),".sr.ALU.matrix", sep="")
    clone1 <- read.table(filename, sep="\t", header=T)
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_",as.character(strand),"/",as.character(sub),".sr.ALU.matrix", sep="")
    clone2 <- read.table(filename, sep="\t", header=T)
    clone  <- rbind(clone1, clone2)
    subclone <- subset(clone, (clone[,7] == 1 | clone[,8] == 1) & ref == 0)
    subclone$refpos <- (subclone$refpos1 + subclone$refpos2)/2

    test1 <- subclone
    nx1 <- test1[,c("seg", "map", "depth","map_size", "map_ratio", "end3", "end5","direction", "refpos", "dist", "A_pair", "A_insert", "A_mm", "A_MapQ", "A_AS", "A_XS")]
    pred.RF1 <- predict(rf1,newdata=nx1,type='prob')[,2]
    pred1 <- data.frame(cbind(as.character(test1$chr), test1$cord1, test1$cord2, as.character(test1$read), pred.RF1))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_",as.character(strand),"/ALU/",as.character(sub),".sr.pred.T1.txt", sep="")
    write(t(pred1), file=filename, ncol=5, sep="\t")

    test2 <-subclone
    nx2 <- test2[,c("pos", "seg", "map", "depth","map_size", "map_ratio", "end3", "end5","direction", "refpos", "dist", "A_pair", "A_insert", "A_mm", "A_MapQ", "A_AS", "A_XS", "oldSR")]
    pred.RF2 <- predict(nb,newdata=nx2,type="raw")[,2]
    pred2 <- data.frame(cbind(as.character(test2$chr), test2$cord1, test2$cord2, as.character(test2$read), pred.RF2))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_",as.character(strand),"/ALU/",as.character(sub),".sr.pred.NB.txt", sep="")
    write(t(pred2), file=filename, ncol=5, sep="\t")

    pred.RF3 <- predict(logr,newdata=nx2,type="response")
    pred3 <- data.frame(cbind(as.character(test2$chr), test2$cord1, test2$cord2, as.character(test2$read), pred.RF3))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_",as.character(strand),"/ALU/",as.character(sub),".sr.pred.logR.txt", sep="")
    write(t(pred3), file=filename, ncol=5, sep="\t")
