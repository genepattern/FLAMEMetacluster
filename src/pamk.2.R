pamk.2 <- function (classname, data, krange = 2:10, scaling = FALSE, diss = inherits(data, 
    "dist"), ...) 
{
    require(cluster)
    if (scaling) 
        sdata <- scale(data, scale = scaling)
    else sdata <- data
    asw <- numeric(max(krange))
    pams <- list()
    for (k in krange) {
        pams[[k]] <- pam(sdata, k, ...)
        asw[k] <- pams[[k]]$silinfo$avg.width
    }
    print(asw)
    write.table(asw,paste(classname,"asw.txt",sep='.'),sep='\t')

    if (.Platform$OS.type == "windows")
    {
        png(paste(classname,"asw.png",sep='.'),width=960,height=960)
    }
    else{
        library(Cairo, lib.loc=Sys.getenv("R_LIBS"))
        CairoPNG(paste(classname,"asw.png",sep='.'),width=960,height=960)    
    }
    plot(krange,asw[krange],type='b',xlab='k',ylab='silhouette width',main=paste(classname,"asw",sep='.'))
    dev.off()
    k.best <- which.max(asw)
    out <- list(pamobject = pams[[k.best]], nc = k.best)
    out
}
