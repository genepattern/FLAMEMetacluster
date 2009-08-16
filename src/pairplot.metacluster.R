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
        library(Cairo)
	    CairoPNG(filename = paste(filename, "pairplots.png", sep = "."), width = 960, height = 960)
	}
	pairs(datafile[1:dim], main = paste(filename, sep = " "), pch = ".",
		labels = names(datafile)[1:dim], col = colors[datafile$cluster])
	dev.off()
}

create.pairplot.legend(maxclus)

}


create.pairplot.legend <- function(num.colors)
{
    ###############
    ## COLOR KEY ##
    ###############

    if (.Platform$OS.type == "windows")
    {
	    png(filename = "pairplots_legend.png", width = 640, height = 500)
	}
	else
	{
        CairoPNG("pairplots_legend.png", width=640, height=500)
    }

    colors <- rainbow(num.colors)
    vec <- c(1:num.colors)
    x = vec%%25
    x <- replace(x, which(x==0), 25)

    n.row = trunc((num.colors+25)/25)
    y = rep(n.row:1, each = 25)[1:num.colors]

    plot(x, y, asp=1, col=colors[1:num.colors], pch=15, cex=3.9, xlab='',ylab='', main=paste("Cluster Color Legend"), axes=F)
    text(cex=0.98, x[seq(1, num.colors, by=3)], y[seq(1, num.colors, by=3)], label=seq(1, num.colors, by=3))

    dev.off()
}


