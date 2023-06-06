###############################################################################
# stats and plots after SUBSAMPLING analysis:
###############################################################################

# Inputs::
args = commandArgs(trailingOnly=TRUE)
X = read.table(args[1],sep="\t",header=TRUE)
RL= as.vector(unique(X$RL))

dirO=paste0(dirname(args[1]),"/")

options(scipen = 999)

#############################################################
# 1) Plotabsolute number of reads
#############################################################
png(paste0(dirO,"1.",RL,"_AbsoluteNumberReads_MappedFiltPromBackground.png"),res=150, width = 2400,height = 1200)
par(mfrow=c(1,2))
plot(X$TotalReads, X$Mapped, type="b",pch=19,bty="n",main=RL, 
	xlim=c(min(X$TotalReads),max(X$TotalReads)),col="darkred",
	ylim=c(0,max(c(X$Mapped))) ,
	xlab = "Total Reads Subsampled", ylab = "Reads ... ")
points(X$TotalReads, X$Filt,type="b",pch=19,col="red")
points(X$TotalReads, X$ReadsInProm,type="b",pch=19,col="black")
points(X$TotalReads, X$ReadsInBackground,type="b",pch=19,col="purple")
  
plot.new()
# legend(1, 95, 
legend("left", 
       legend=c("Mapped","Filtered","ReadsInBackground","ReadsInProm"),
       col=c("darkred", "red","purple", "black"), lty=1, cex=1.5,lwd = 4,box.lty = 0)
dev.off()
#############################################################


#############################################################
# Plotabsolute number of reads
#############################################################
png(paste0(dirO,"2.",RL,"_AbsoluteNumberReads_in_peaksANDpromoters.png"),res=150, width = 2400,height = 1200)
par(mfrow=c(1,2))
plot(X$TotalReads, X$ReadsInProm,type="b",pch=19,bty="n",main=RL, 
	xlim=c(min(X$TotalReads),max(X$TotalReads)), col="black",
	ylim=c(0,max(X$ReadsInBAMPEpeaks)) ,
	xlab = "Total Reads Subsampled", ylab = "Reads in ...")
points(X$TotalReads, X$ReadsInNFRpeaks,type="b",pch=19,col="dodgerblue2")
points(X$TotalReads, X$ReadsInNFRpeaksOverProm,type="b",pch=19,col="lightskyblue")
points(X$TotalReads, X$ReadsInBAMPEpeaks,type="b",pch=19,col="coral3")
points(X$TotalReads, X$ReadsInBAMPEpeaksOverProm,type="b",pch=19,col="darksalmon")
# points(X$TotalReads, X$ReadsInBackground,type="b",pch=19,col="purple")
plot.new()
# legend(1, 95, 
legend("left", 
       legend=c("ReadsInProm", "ReadsNFRpeaks","ReadsInNFRpeaksOverProm", "ReadsInBAMPEpeaks","ReadsInBAMPEpeaksOverProm"),
       col=c("black", "dodgerblue2","lightskyblue", "coral3","darksalmon"), lty=1, cex=1.5,lwd = 4,box.lty = 0)
dev.off()
#############################################################

######################################################################
# Plot realtive number of reads per every subsampled fraction
######################################################################
png(paste0(dirO,"3.",RL,"_ReadsRelativeToFilteredReads_peaksANDpromoters.png"),res=150, width = 2400,height = 1200)
par(mfrow=c(1,2))
plot(X$TotalReads, X$ReadsInProm/X$Filt*100,type="b",pch=19,bty="n",main=RL,
	xlim=c(min(X$TotalReads),max(X$TotalReads)),col="black",
	ylim=c(0,35) ,
	xlab = "Total Reads Subsampled", ylab = "% Reads in ...")
points(X$TotalReads, X$ReadsInNFRpeaks/X$Filt*100,type="b",pch=19,col="dodgerblue2")
points(X$TotalReads, X$ReadsInNFRpeaksOverProm/X$Filt*100,type="b",pch=19,col="lightskyblue")
points(X$TotalReads, X$ReadsInBAMPEpeaks/X$Filt*100,type="b",pch=19,col="coral3")
points(X$TotalReads, X$ReadsInBAMPEpeaksOverProm/X$Filt*100,type="b",pch=19,col="darksalmon")
plot.new()
# legend(1, 95, 
legend("left", title ="% over filtered reads -every subsampled-",
       legend=c("%ReadsInProm", "%ReadsNFRpeaks","%ReadsInNFRpeaksOverProm", "%ReadsInBAMPEpeaks","%ReadsInBAMPEpeaksOverProm"),
       col=c("black", "dodgerblue2","lightskyblue", "coral3","darksalmon"), lty=1, cex=1.3,lwd = 4,box.lty = 0)
dev.off()
#############################################################


######################################################################
# Plot realtive number of reads per every maximun
######################################################################
png(paste0(dirO,"4.",RL,"_ReadsRelativeToMaxSubsampledFilteredReads_peaksANDpromoters.png"),res=150, width = 2400,height = 1200)
par(mfrow=c(1,2))
maxSubsampled= X$Filt[X$proportion == max(X$proportion)]
plot(X$TotalReads, X$ReadsInProm/maxSubsampled*100,type="b",pch=19,bty="n",main=RL,
	xlim=c(min(X$TotalReads),max(X$TotalReads)),col="black",
	ylim=c(0,40) ,
	xlab = "Total Reads Subsampled", ylab = "% Reads in ...")
points(X$TotalReads, X$ReadsInNFRpeaks/maxSubsampled*100,type="b",pch=19,col="dodgerblue2")
points(X$TotalReads, X$ReadsInNFRpeaksOverProm/maxSubsampled*100,type="b",pch=19,col="lightskyblue")
points(X$TotalReads, X$ReadsInBAMPEpeaks/maxSubsampled*100,type="b",pch=19,col="coral3")
points(X$TotalReads, X$ReadsInBAMPEpeaksOverProm/maxSubsampled*100,type="b",pch=19,col="darksalmon")
plot.new()
legend("bottom", title ="% over filtered reads\n -max sumbsampled-",
       legend=c("%ReadsInProm", "%ReadsNFRpeaks","%ReadsInNFRpeaksOverProm", "%ReadsInBAMPEpeaks","%ReadsInBAMPEpeaksOverProm"),
       col=c("black", "dodgerblue2","lightskyblue", "coral3","darksalmon"), lty=1, cex=1.4,lwd = 4,box.lty = 0)
dev.off()
#############################################################



######################################################################
# Plot realtive number of peaks in Promoters
######################################################################
png(paste0(dirO,"5.",RL,"_PeaksPromRel2TotalPeaks.png"),res=150, width = 2400,height = 1200)
par(mfrow=c(1,2))
plot(X$TotalReads,X$NFRPeaksInPromoters/X$NFRpeaks*100,type="b",pch=19,bty="n",main=RL,
	xlim=c(min(X$TotalReads),max(X$TotalReads)),col="dodgerblue3",
	ylim=c(0,100) ,
	xlab = "Total Reads Subsampled", ylab = "% Peaks")
points(X$TotalReads, X$BAMPEPeaksInPromoters/X$BAMPEPeaks*100,type="b",pch=19,col="coral2")
plot.new()
# legend(1, 95, 
legend("center", title ="",
       legend=c("% NFRpeaksProm / NFRpeaks", "% BAMPEpeaksProm / BAMPEpeaks"),
       col=c("dodgerblue1","coral2"), lty=1, cex=1.4,lwd = 4,box.lty = 0)
dev.off()
#############################################################

#########################################################################################
# Plot realtive number of peaks in Promoters over maximun number of peaks in promoter 
#########################################################################################
png(paste0(dirO,"6.",RL,"_PeaksPromRel2MaxPeaks.png"),res=150, width = 2400,height = 1200)
par(mfrow=c(1,2))
maxNFRinProm = max(X$NFRPeaksInPromoters)
maxBAMPEPinProm = max(X$BAMPEPeaksInPromoters)
plot(X$TotalReads,X$NFRPeaksInPromoters/maxNFRinProm*100,type="b",pch=19,bty="n",main=RL,
	xlim=c(min(X$TotalReads),max(X$TotalReads)), col="dodgerblue3",
	ylim=c(0,100) ,
	xlab = "Total Reads Subsampled", ylab = "% Peaks")
points(X$TotalReads, X$BAMPEPeaksInPromoters/maxBAMPEPinProm*100,type="b",pch=19,col="coral2")
plot.new()
legend("center", title ="relative to max Peaks in Prom\nper every sample",
	legend=c("% NFRpeaksProm / maxNFRpeaksInProm", "% BAMPEpeaksProm / maxBAMPEpeaksInProm"),
	col=c("dodgerblue1","coral2"), lty=1, cex=1.2,lwd = 4,box.lty = 0)
dev.off()
#############################################################


#########################################################################################
# Plot realtive number of peaks in Promoters over maximun number of peaks in promoter 
#########################################################################################
png(paste0(dirO,"7.",RL,"RLsNOVASEQ_PeaksTotalRel2MaxTotalPeaks.png"),res=150, width = 2400,height = 1200)
par(mfrow=c(1,2))
maxNFR = max(X$NFRpeaks)
maxBAMPE = max(X$BAMPEPeaks)
plot(X$TotalReads,X$NFRpeaks/maxNFR*100,type="b",pch=19,bty="n",main=RL,
	xlim=c(min(X$TotalReads),max(X$TotalReads)),col="dodgerblue3",
	ylim=c(0,100) ,
	xlab = "Total Reads Subsampled", ylab = "% Peaks")
points(X$TotalReads, X$BAMPEPeaks/maxBAMPE*100,type="b",pch=19,col="coral2")
plot.new()
legend("center", title ="relative to max Total Peaks \nper every sample",
       legend=c("% NFRpeaks / maxNFRpeaks", "% BAMPEpeaks / maxBAMPEpeaks"),
       col=c("dodgerblue1","coral2"), lty=1, cex=1.2,lwd = 4,box.lty = 0)
dev.off()
#############################################################


#############################################################
# Relative to all samples RL together track 
#############################################################
NovaFilteredReads=1059250812
NovaRLs_ReadsInProm=117906845
(NovaRLs_ReadsInProm / NovaFilteredReads)



######################################################################
# Plot realtive number of reads per every subsampled fraction
######################################################################
png(paste0(dirO,"8.",RL,"_Readsin_NFRpeaksFromAllRLsNovaseqtogether.png"),res=150, width = 2400,height = 1200)
par(mfrow=c(1,2))
maxFilt = max(X$Filt)
plot(X$TotalReads, X$readsin_NFRPeaksAllRLs/X$Filt*100,type="b",pch=19,bty="n",main=RL,
	xlim=c(min(X$TotalReads),max(X$TotalReads)),col="black",
	ylim=c(0,60) ,
	xlab = "Total Reads Subsampled", ylab = "% Reads in ...")
points(X$TotalReads, X$readsin_NFRPeaksAllRLs/maxFilt*100,type="b",pch=19,col="dodgerblue2")
plot.new()
legend("left", title ="% over filtered reads -every subsampled-",
       legend=c("%ReadsIn_allRLsNFRspeaks/FiltReads", "%ReadsIn_allRLsNFRspeaks/maxFiltReads"),
       col=c("black", "dodgerblue2"), lty=1, cex=1,lwd = 4,box.lty = 0)
dev.off()
########################################################################################
