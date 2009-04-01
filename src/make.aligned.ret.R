makeAlignedRet <- function(metaparamfiles) {
	
for (i in 1:length(metaparamfiles)) {
	paramfile <- metaparamfiles[i]
	filename <- strsplit(paramfile,".aligned.parameters.txt")[[1]][1]
	param <- read.table(paramfile,header=T)
	dim <- ncol(param)-5
	param <- subset(param, param$mus != 0.0000)
	num.component <- length(na.omit(param$props))

#make ret
	pro <- as.vector(na.omit(param$props))
	mu <- matrix(param$mus,ncol = num.component, nrow = dim)
	sigma <- array(0, c(dim, dim, num.component))
	for (c in 1:num.component) {
		sigma[,,c] = as.matrix(param[((c-1)*dim+1):(c*dim),4:7])
	}
	dof <- as.vector(na.omit(param$df))
	delta <- matrix(param$alpha, ncol = num.component, nrow = dim)
	result <- c()
	result$pro = pro
	result$mu = mu
	result$sigma = sigma
	result$dof = dof
	result$delta = delta
	dput(result, paste(filename,"aligned.ret",sep='.'))
}

}