plotHeatmap <- function(paramfiles, concatfiles, dist) {
	#library(Heatplus)
	suppressMessages(library(gplots))
	for (i in 1:length(paramfiles)) {
	print(i)
	paramfile <- paramfiles[i]
	param <- read.table(paramfile, header=T, sep="\t")
	num.cluster <- length(na.omit(param$prop))
	dim = ncol(param)-5
	concatfile <- concatfiles[i]
	concat <- read.table(concatfile,header=T)
	names <- names(concat)[1:dim]
	locations<- matrix(ncol = dim+1, nrow = num.cluster, dimnames = list(c(1:num.cluster),c(names,"weight")))
	locations[,1:dim] <- matrix(t(matrix(param$mod,dim,num.cluster)),ncol = dim, nrow = num.cluster, dimnames = list(c(1:num.cluster),c(names)))
	weights <- matrix(na.omit(param$prop))*100
	locations[,(dim+1)] = weights
	
	dummies <- matrix(ncol = dim+1, nrow = 3, dimnames = list(c(1:3),c(names,"weight")))
	for (d in 1:dim) {
		dummies[1,d] = quantile(concat[,d],0.1)
		dummies[2,d] = quantile(concat[,d],0.5)
		dummies[3,d] = quantile(concat[,d],0.9)
	}
	dummies[(1:3),(dim+1)] = 0 
	heatmap<- rbind(locations,dummies)
	fname <- strsplit(paramfile,"\\.parameters")[[1]][1]
	write.table(locations,file = paste(fname,"locations.txt",sep='.'),sep="\t",quote=F,row.names=F)	
	heatmap <- heatmap[,1:dim]
	row.names(heatmap) <- c(paste(1:num.cluster, round(weights,1), sep = ":"),paste("ref: ",c("10","50","90"),"%",sep=''))
	colnames(heatmap) <- names

    if (.Platform$OS.type == "windows")
    {
	    png(filename = paste(fname, "heatmap.png",sep='.'),height = 960*1.5, width = 960*1.5)
	}
	else
	{
        library(Cairo, lib.loc=Sys.getenv("R_LIBS"))
	    CairoPNG(filename = paste(fname, "heatmap.png",sep='.'),height = 960*1.5, width = 960*1.5)
    }

    if(.Platform$OS.type == "unix")
	{
	    heatmap.colors <- redblue(75)
	}
	else
	{
	    heatmap.colors <- bluered(75)
	}

	heatmap.2(heatmap, Colv=F, Rowv=F, col = bluered(75),scale = "column", dendrogram="none", key=T,
	keysize = F,symkey=F,density.info = "none",cellnote=round(heatmap,1),notecex=2,notecol="black",trace="none",rowsep=nrow(heatmap)-3,sepcolor="white", main = paste(fname,"cluster log intensities",sep=''))
	dev.off()
	}
}