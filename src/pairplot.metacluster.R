pairplotMetaClusters2D <- function(concatfiles, dim, maxclus){

for (i in 1:length(concatfiles)) {
	print(i)
	filename <- strsplit(concatfiles[i], "\\.member.txt")
	datafile <- data.frame(read.table(concatfiles[i],header=T,sep="\t"))
	colors <- rainbow(maxclus)
	if (.Platform$OS.type == "windows")
    {
	    png(filename = paste(filename, "pairplots.png", sep = "."), width = 960, height = 960)
	}
	else
	{
        library(Cairo, lib.loc=Sys.getenv("R_LIBS"))
	    CairoPNG(filename = paste(filename, "pairplots.png", sep = "."), width = 960, height = 960)
	}
	pairs(datafile[1:dim], main = paste(filename, sep = " "), pch = ".",
		labels = names(datafile)[1:dim], col = colors[datafile$cluster])
	dev.off()
}
}

