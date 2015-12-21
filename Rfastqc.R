rm(list=ls())
getProgramName<-function(){
	args <- commandArgs(trailingOnly = FALSE)
	sub("--file=", "", args[grep("--file=", args)])
}
program <- getProgramName()
args <- commandArgs(trailingOnly = TRUE)
if (length(args)<2) {
	stop(sprintf("Please input the information of arguments:
	args[1]		outfile
	args[2]		read1(/share/work1/staff/xuxiong/test/13C37198_L7_I012.R1.clean.fastq.gz)
	args[3]		read2
Example:
	Rscript %s out /share/work1/staff/xuxiong/test/13C37198_L7_I012.R1.clean.fastq.gz 
Usage:
	Rscript %s outfile read1 [read2]",program,program)
	)
}
ptm <- proc.time()
options(scipen=999)
outfile<-args[1]
fastq1<-args[2]
fastq2<-ifelse(length(args)>2,args[3],"")
.libPaths("/home/xuxiong/bin/R_installed_package/")

barplot_Data_range<-function(outfile,Data){
	colors <- c('#4682B4','#A0522D','#FF8C00','#87CEEB','#6B8E23','#6A5ACD','#778899','#DAA520','#B22222','#FF6699')
	png(paste(outfile,"_read_freq_count",".png",sep=''),pointsize=18,width=900,height=600)
	Data_range<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,50,100,500)
	Data<-as.numeric(Data)
	Data_range_count<-sapply(Data_range,function(X){
		ifelse(X==max(Data_range),length(Data[Data>=X]),length(Data[Data>=X & Data<Data_range[which(Data_range %in% X)+1]]))
	})
#	cat(Data_range_count,"\n",file=stderr())
	mp <- barplot(Data_range_count,
		beside=TRUE,
		width = 0.5,
		axisnames = F,
		cex.names = 0.8,
		xlab="Hits range",
		ylab="frequency counts",
		col=colors[1],
		ylim=c(0,max(Data_range_count)*1.2),
		axes=TRUE,
		plot=TRUE,
		xpd=FALSE,
		pch=15,
	#	args.legend = list(x = "topleft"),
		main=list("Fastq Hits distribution"),
	)
#	cat(mp,"\n",axTicks(1),"\n",axTicks(2),"\n",file=stderr())
#	cat(sprintf("%.3f%%\n%.3f%%\n",sum(Data_range_count)/sum(Data)*100,(1-Data_range_count[1]/sum(Data))*100),file=stderr())
#	cat(paste(round(Data_range_count/sum(Data)*100,3),"%",sep=""),"\n",file=stderr())
	text(mp,Data_range_count,
#		labels = sapply(Data_range_count,function(x){sprintf("%d\n%.2f%%",x,x/sum(Data)*100)}),
		labels = sapply(Data_range_count,function(x){sprintf("%d",x)}),
		adj=c(0.5,-0.5), cex=0.6,xpd = TRUE
	)
	text(mp, par("usr")[3],
		labels = sapply(Data_range,function(X){
			ifelse(X==max(Data_range),
				paste(">=",max(Data_range),sep=""),
				ifelse(Data_range[which(Data_range %in% X)+1]-X>1,paste('[',X,',',Data_range[which(Data_range %in% X)+1],')',sep=''),X)
			)
		}),
		srt = 45, adj = c(1,1),cex=0.8,xpd = TRUE
	)
	#arrows(x0=mp,y0=as.numeric(Data_range_count)*0.95,x1=mp,y1=as.numeric(Data_range_count)*1.05,angle=90,code=3,length=0.04,lwd=1)
	box()
	pic=dev.off()
	cat("Done barplot Data_range\n",file=stderr())
}

plot_Data_range<-function(outfile,Data){
	colors <- c('#4682B4','#A0522D','#FF8C00','#87CEEB','#6B8E23','#6A5ACD','#778899','#DAA520','#B22222','#FF6699')
	png(paste(outfile,"_dup_level",".png",sep=''),pointsize=18,width=900,height=600)
#	Data_range<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,50,100,500)
	Data_range<-c(1,2,3,4,5,6,7,8,9,10)
	Data_range_count<-sapply(Data_range,function(X){
		ifelse(X==max(Data_range),
			sum(Data[Data>=X]),
			sum(Data[Data>=X & Data<Data_range[which(Data_range %in% X)+1]]))
	})
	Data_range_count_norm<-Data_range_count/Data_range_count[1]*100
#	cat(Data_range_count_norm,"\n",file=stderr())
	plot(seq_along(Data_range),Data_range_count_norm,
		type="l",
		xlab="Duplicate level",
		ylab="Percentage(%)",
		col=colors[1],
		axes=TRUE,
	#	xaxs='r',
		xaxt='n',
		xpd=FALSE,
		ylim=c(0,max(Data_range_count_norm)*1.2), 
		pch=1,
		font.lab=1.2,font.main=2,font.axis=1,lwd=2,
		cex.lab=1.5,cex.main=1.5,cex.axis=1,cex.sub=1,cex=0.5,
		main=list("Sequence duplication level")
	)
	axis(1,at=Data_range,
		label=sapply(Data_range,function(X){
			ifelse(X==max(Data_range),
				paste(">=",max(Data_range),sep=""),
				ifelse(Data_range[which(Data_range %in% X)+1]-X>1,paste('[',X,',',Data_range[which(Data_range %in% X)+1],')',sep=''),X)
			)
		}),
		srt = 45
	)
	legend("topright",
		legend=sprintf("Dup%%: %.3f%%",(1-Data_range_count[1]/sum(Data))*100),
		cex=0.8,
		inset = 0.01
	)
	box()
	pic=dev.off()
	cat("Done plot Data_range\n",file=stderr())
}

plot_GC_density0<-function(outfile,Data,maxLen){
	colors <- c('#4682B4','#A0522D','#FF8C00','#87CEEB','#6B8E23','#6A5ACD','#778899','#DAA520','#B22222','#FF6699')
	png(paste(outfile,"_GC_density",".png",sep=''),pointsize=18,width=900,height=600)
#	cat("Min: ",min(Data),"\tMax: ",max(Data),sprintf("\tMean GC%%: %.2f%%",mean(Data)*100),"\n",file=stderr())
	D_Data<-density(Data*100,n=maxLen)
#	cat(D_Data$x,"\n",D_Data$y,"\n",file=stderr())
	plot(D_Data,
		type="o",
		xlab="GC(%)",
		ylab="density",
		col=colors[1],
		axes=TRUE,
	#	xaxs='r',
	#	xaxt='n',
		xpd=FALSE,
#		ylim=c(0,max(Data)*1.2),
		pch=1,
		font.lab=1.2,font.main=2,font.axis=1,lwd=2,
		cex.lab=1.5,cex.main=1.5,cex.axis=1,cex.sub=1,cex=0.5,
		main=list("GC density distribution")
	)
	segments(D_Data$x[D_Data$y==max(D_Data$y)],max(D_Data$y),D_Data$x[D_Data$y==max(D_Data$y)],0,col="black",lwd=1,lty=2)
	legend("topright",
		legend=c(sprintf("Mean GC%%: %.2f%%",mean(Data)*100),sprintf("Max density GC%%: %.2f%%",D_Data$x[D_Data$y==max(D_Data$y)])),
		cex=0.8,
		inset = 0.01
	)
	box()
	pic=dev.off()
	cat("Done plot_GC_density\n",file=stderr())
}

plot_GC_density<-function(outfile,Data,maxLen){
	colors <- c('#4682B4','#A0522D','#FF8C00','#87CEEB','#6B8E23','#6A5ACD','#778899','#DAA520','#B22222','#FF6699')
	png(paste(outfile,"_GC_density",".png",sep=''),pointsize=18,width=900,height=600)
#	cat("Min: ",min(Data),"\tMax: ",max(Data),sprintf("\tMean GC%%: %.2f%%",mean(Data)*100),"\n",file=stderr())
	Dens<-density(Data*100,n=maxLen)
	x<-Dens$x
	y<-Dens$y
	tab<-data.frame(x=x,y=y)
	nlmod <- nls(y ~ k/(sqrt(2*pi)*sigma)*exp(-1/2*(x-mu)^2/sigma^2), start=c(mu=50,sigma=100,k=0.1) , data = tab)
	v <- summary(nlmod)$parameters[,"Estimate"]
	plot(tab,
		type="o",
		xlab="GC(%)",
		ylab="Count",
		col=colors[1],
		axes=TRUE,
	#	xaxs='r',
	#	xaxt='n',
		xpd=FALSE,
#		ylim=c(0,max(Data)*1.2),
		pch=1,
		font.lab=1.2,font.main=2,font.axis=1,lwd=2,cex.lab=1.5,cex.main=1.5,cex.axis=1,cex.sub=1,cex=0.5,
		main=list("GC density distribution")
	)
	plot(function(x) v[3]/(sqrt(2*pi)*v[2])*exp(-1/2*(x-v[1])^2/v[2]^2),type='o',pch=15,
		font.lab=1.2,font.main=2,font.axis=1,lwd=2,cex.lab=1.5,cex.main=1.5,cex.axis=1,cex.sub=1,cex=0.5,
		col=colors[2],add=T,xlim=range(tab$x) 
	)
	#segments(D_Data$x[D_Data$y==max(D_Data$y)],max(D_Data$y),D_Data$x[D_Data$y==max(D_Data$y)],0,col="black",lwd=1,lty=2)
	legend("topright",
		legend=c("GC count per read","Theoretical Distribution"),
		cex=0.8,
		#fill=TRUE,
		col=colors,
		lty=c(1,1),
		pch=c(1,12),
		inset = 0.01
	)
	box()
	pic=dev.off()
	cat("Done plot_GC_density\n",file=stderr())
}

matrix.axes <- function(data) {
	x <- (1:dim(data)[1] - 1) / (dim(data)[1] - 1);
	axis(side=1, at=x, labels=1:dim(data)[1], las=2,cex.axis=0.6);
	x <- (1:dim(data)[2] - 1) / (dim(data)[2] - 1);
	axis(side=2, at=x, labels=1:dim(data)[2], las=2,cex.axis=0.6);
	grid(nx=(dim(data)[1]-1), ny=(dim(data)[2]-1), col="black", lty=1,lwd=0.6);
}

plot_quality<-function(outfile,Data){
	library(gplots)
	colors <- c('#4682B4','#A0522D','#FF8C00','#87CEEB','#6B8E23','#6A5ACD','#778899','#DAA520','#B22222','#FF6699')
	png(paste(outfile,"_quality",".png",sep=''),pointsize=18,width=900,height=600)
#	Index<-seq_along(1:nrow(Data))[apply(Data,1,function(X){!all(X==0)})]
	M<-t(Data[seq(34,75),])
	filled.contour(M, plot.axes=matrix.axes(M), 
		col=colorpanel(30, "white", "red"), 
		nlevels=20,
		main=paste(outfile,"_quality",sep='')
	)
	pic=dev.off()
	cat("Done plot_quality\n",file=stderr())
}

plot_quality2<-function(outfile,Data){
	library(lattice)
	colors <- c('#4682B4','#A0522D','#FF8C00','#87CEEB','#6B8E23','#6A5ACD','#778899','#DAA520','#B22222','#FF6699')
	trellis.device(device="png", filename=paste(outfile,"_quality2",".png",sep=''),pointsize=18,width=900,height=600)
	M<-t(Data[seq(34,75),])
	rgb.palette <- colorRampPalette(c("white", "blue"), space = "rgb")
	my.plot <-levelplot(M, main=paste(outfile,"_quality",sep=''),
		xlab="cycle",
		ylab="quality score",
		col.regions=rgb.palette(1200), cuts=1000, at=seq(0,max(M),l=1000),
		pretty=TRUE
	)
	print(my.plot)
	dev.off()
	cat("Done plot_quality2\n",file=stderr())
}

plot_boxplot<-function(outfile,Data){
	colors <- c('#4682B4','#A0522D','#FF8C00','#87CEEB','#6B8E23','#6A5ACD','#778899','#DAA520','#B22222','#FF6699')
	png(paste(outfile,"_boxplotquality",".png",sep=''),pointsize=18,width=900,height=600)
	Index<-seq_along(1:nrow(Data))[apply(Data,1,function(X){!all(X==0)})]
	M<-Data[seq(34,75),]
	sumQ30=sum(as.numeric(Data[Index[Index>63],]))
	sumQ20=sum(as.numeric(Data[Index[Index>53],]))
	sumQ=sum(as.numeric(M))
#	cat(sprintf("Q30: %.3f%%(%.0f/%.0f)\nQ20: %.3f%%\n",100*sumQ30/sumQ,sumQ30,sumQ,100*sumQ20/sumQ),file=stderr())
	cat(sprintf("Q30: %.3f%%\nQ20: %.3f%%\n",100*sumQ30/sumQ,100*sumQ20/sumQ),file=stderr())
	rownames(M)<-seq(0,41)
#	print(Index)
	boxplot(apply(M,2,function(X){rep(seq(34,75),X%/%100)}),
		main=paste(outfile,"_quality",sep=''),
		xlab="cycle",
		ylab="quality score",
#		names=colnames(RPKM),
		show.names=T,
		varWidth=T,
		notch=F,
		outline=F,
		col=colors[1],
		plot=TRUE,
		las=0,
		boxwex=0.75,
		ylim=c(34,75),
#		ylim=range(Index),
		cex.lab=1,cex.main=1.5,cex.axis=0.8,cex.sub=1,cex=0.5,
		font.lab=1,font.main=2,font.axis=0.8,lwd=1,
		pch=1
	)
	pic=dev.off()
	cat("Done boxplot_quality\n",file=stderr())
}

draw_length_distribution<-function(outfile_prefix,Len){
	png(paste(outfile_prefix,"len.png",sep='_'),pointsize=18,width=900,height=600)
	colors <- c('#4682B4','#87CEEB','#6B8E23','#A0522D','#FF8C00','#6A5ACD','#778899','#DAA520','#B22222','#FF6699')
	totalBase<-sprintf("Total base: %.0f",sum(as.numeric(as.numeric(names(Len))*Len)))
	totalReads<-sprintf("Total reads: %.0f",sum(as.numeric(Len)))
	meanLen<-sprintf("Mean length: %.1f",weighted.mean(as.numeric(names(Len)), Len))
	cat(totalBase,totalReads,meanLen,"\n",file=stderr())
#	print(Len)
	mp <- barplot(Len,
	#	beside=TRUE,
		width = 1,
		names.arg=names(Len),
		axisnames = T,
		cex.names = 0.8,
		cex.axis = 0.8,
		xlab="length(bp)",
		ylab="Counts",
		col=colors[1],
		ylim=c(0,max(Len)*1.4),
		axes=TRUE,
		plot=TRUE,
		xpd=FALSE,
		pch=15,
		args.legend = names(Len),
		main=list("Length distribution"),
	)
#	text(mp,Len,
#		labels = sprintf("%d\n(%.2f%%)",Len,Len/sum_count*100),
#		adj=c(0.5,-0.5),
#		cex=0.75,
#		xpd = TRUE
#	)
#	arrows(x0=mp,y0=Len*0.95,x1=mp,y1=Len*1.05,angle=90,code=3,length=0.04,lwd=0.4)
#	text(mp, par("usr")[3],labels = names(Len),srt = 35, cex=0.6,adj = c(1,1),xpd = TRUE)
	legend("topleft",
		legend=c(totalBase,meanLen,totalReads),
		cex=0.8,
		col=colors[1],
#		lty=1,
#		lwd=3,
#		fill=colors[1],
	#	pch=c(1,20),
		inset = 0.01
	)
	box()
	pic=dev.off()
	cat("Done Length distribution\n",file=stderr())
}

draw_nucleotide_distribution<-function(outfile_prefix,Data){
	colors <- c('#4682B4','#87CEEB','#6B8E23','#A0522D','#FF8C00','#6A5ACD','#778899','#DAA520','#B22222','#FF6699')
	png(paste(outfile_prefix,'_nucleotide.png',sep=''),pointsize=18,width=900,height=600)
	par(mar=c(5.1,4.1,4.1,2.1))
	plot(1:ncol(Data),Data[1,],
		type="n",
		xlab="Cycle",
		ylab="Counts",
		col=colors,
		axes=TRUE,
#		xaxs='r',
#		xaxt='n',
		xpd=T,
		ylim=c(0,max(Data)*1.5), 
#		xlim=c(0,100),
		pch=1,
		font.lab=1.2,font.main=2,font.axis=1,lwd=2,
		cex.lab=1.5,cex.main=1.5,cex.axis=0.8,cex.sub=1,cex=0.5,
		main=list("Nucleotide Content Distribution")
	)
	sapply(1:nrow(Data),function(X){
		lines(1:ncol(Data),Data[X,],lwd=2,col=colors[X],type = "o",pch=20)
	})
	legend("topright",
		legend=c("T","C","A","G","N"),
		cex=0.8,
		col=colors,
		lty=1,
		lwd=3,
#		fill=colors[c(1,2)],
		pch=20,
		inset = 0.01
	)
	pic=dev.off()
	cat("Done nucleotide distribution\n",file=stderr())
}

dyn.load(paste(dirname(program),ifelse(grepl('Linux|Darwin',Sys.info()["sysname"],perl=T),"Rgzfastq_uniq_3.so","Rgzfastq_uniq_3.dll"),sep="/"))
List<-.Call("qsort_hash_count",fastq1,fastq2)
barplot_Data_range(outfile,List[[1]])
plot_Data_range(outfile,List[[1]])
plot_GC_density(paste(outfile,"R1",sep=""),List[[2]],max(as.numeric(names(List[[5]][List[[5]]>0]))))

List[[3]]<-List[[3]][,apply(List[[3]],2,function(Y){!all(Y==0)})]

#plot_quality(paste(outfile,"R1",sep=""),List[[3]])
plot_quality2(paste(outfile,"R1",sep=""),List[[3]])
plot_boxplot(paste(outfile,"R1",sep=""),List[[3]])

List[[4]]<-List[[4]][,apply(List[[4]],2,function(Y){!all(Y==0)})]
draw_nucleotide_distribution(paste(outfile,"R1",sep=""),List[[4]])

draw_length_distribution(paste(outfile,"R1",sep=""),List[[5]][List[[5]]>0])
if (length(args)>2) {
	plot_GC_density(paste(outfile,"R2",sep=""),List[[6]],max(as.numeric(names(List[[9]][List[[9]]>0]))));

	List[[7]]<-List[[7]][,apply(List[[7]],2,function(Y){!all(Y==0)})]
#	plot_quality(paste(outfile,"R2",sep=""),List[[7]])
	plot_quality2(paste(outfile,"R2",sep=""),List[[7]])
	plot_boxplot(paste(outfile,"R2",sep=""),List[[7]])

	List[[8]]<-List[[8]][,apply(List[[8]],2,function(Y){!all(Y==0)})]
	draw_nucleotide_distribution(paste(outfile,"R2",sep=""),List[[8]])

	draw_length_distribution(paste(outfile,"R2",sep=""),List[[9]][List[[9]]>0])
}
