# define your own Rpath if needed

args = commandArgs(trailingOnly=TRUE)
df = read.table(args[1], header=FALSE)	
name=args[2]

png(name,res=100, width = 1200,height = 800)
hist(as.numeric(df$V1),breaks=500, main="BAMPE insert sizes")
abline(v=75,col="red",lty=2)
abline(v=150,col="blue",lty=2)
abline(v=225,col="red",lty=2)
abline(v=300,col="blue",lty=2)
abline(v=375,col="red",lty=2)
abline(v=450,col="blue",lty=2)
dev.off()


