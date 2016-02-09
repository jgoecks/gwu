library(exomeCopy)

args <- commandArgs(TRUE)
input <- args[1]
target.file <- args[2]
output <- args[3]

bam.file <- c(input)

sample.names <- c("sample")
target.df <- read.delim(target.file, header = FALSE)
target <- GRanges(seqname = target.df[, 1], IRanges(start = target.df[,2] + 1, end = target.df[, 3]))
counts <- RangedData(space = seqnames(target), ranges = ranges(target))

for (i in 1:length(bam.file)) {
	counts[[sample.names[i]]] <- countBamInGRanges(bam.file[i],target)
	}

exon_size<-width(counts)
sample_counts<-counts$sample

breaks1=c(0,50,100,150,200,250,300,350,400,450,500,550,600,1000000)
breaks2=c(0,10,20,30,40,50,60,70,80,90,100,1000000)

coverage=(sample_counts/exon_size)*100
plot_coverage=subset(coverage, coverage < 1000)
mean_cov=mean(plot_coverage)
median_cov=median(round((plot_coverage), digits=2))
low_coverage=100*((length(subset(coverage, coverage < 20)))/(length(coverage)))
low_coverage2=round(low_coverage, digits=2)


pdf(file=output)
cols <- c("red", "red")[(plot_coverage < 20)]
hist(plot_coverage,xlim=c(0,500),breaks=50,xlab="Average Coverage of an Exon",ylab="Frequency",main = sample.names[1],col=cols)
abline(v=median_cov,col="purple")
legend("topright",legend=paste("Median Coverage = ", median_cov ,"\n","Less than 20x coverage = ",low_coverage2," %",sep=" "))

dev.off()






