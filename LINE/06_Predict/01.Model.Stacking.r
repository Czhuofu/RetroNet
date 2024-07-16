#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
set.seed(1)
datpath<-args[1]
ver<-args[2]
sub<-args[3]
masterpath<-args[4]
library('randomForest')

filename <- paste(as.character(masterpath),"/LINE/06_Predict/RFXIV_1.rds", sep="")
rf1 <- readRDS(filename)
filename <- paste(as.character(masterpath),"/LINE/06_Predict/RFXIV_2.rds", sep="")
rf2 <- readRDS(filename)

    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_0/LINE/",as.character(sub),".pairs", sep="")
    clone1 <- read.table(filename, sep="\t", header=T)
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_1/LINE/",as.character(sub),".pairs", sep="")
    clone2 <- read.table(filename, sep="\t", header=T)
    clone  <- rbind(clone1, clone2)

    test1 <- subset(clone, !is.na(id))
    nx1 <- test1[,c( "insertion", "strand", "read1", "PE1","up1","read1_RF","read1_LR","read1_NB", "read2","PE2","up2","read2_RF","read2_LR","read2_NB","dist1","dist2","dist3", "dist4", "id", "pA1", "pA2")]
    pred.RF1 <- predict(rf1,newdata=nx1,type='prob')[,2]
    pred1 <- data.frame(cbind(as.character(test1$insertion), test1$strand, pred.RF1, as.character(test1$read1), test1$PE1, test1$up1, test1$read1_RF, as.character(test1$read2), test1$PE2, test1$up2, test1$read2_RF, test1$dist5, test1$id))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_1/LINE/",as.character(sub),".pred.P1.txt", sep="")
    write(t(pred1), file=filename, ncol=13, sep="\t")

    test2 <- subset(clone, is.na(id))
    nx2 <- test2[,c("insertion", "strand", "read1", "PE1","up1","read1_RF","read1_LR","read1_NB", "read2","PE2","up2","read2_RF","read2_LR","read2_NB","dist1","dist2","dist3", "dist4", "pA1", "pA2")]
    pred.RF2 <- predict(rf2,newdata=nx2,type='prob')[,2]
    pred2 <- data.frame(cbind(as.character(test2$insertion), test2$strand, pred.RF2, as.character(test2$read1), test2$PE1, test2$up1, test2$read1_RF, as.character(test2$read2), test2$PE2, test2$up2, test2$read2_RF, test2$dist5, test2$id))
    filename <- paste(as.character(datpath),"/",as.character(sub),"/retro_v",as.character(ver),"_1/LINE/",as.character(sub),".pred.P2.txt", sep="")
    write(t(pred2), file=filename, ncol=13, sep="\t")

