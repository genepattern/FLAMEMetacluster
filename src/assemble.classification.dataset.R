makeDataset<-function(stage, classes, all.k, dim, num.samples, fileset.name) {
allsamples.param <- c()
files.names <- c()
input.classes <- c()
for (c in 1:length(classes)) {
	class.param<-c()
	pfiles <- dir("./", pattern = classes[c])
	pfiles <- pfiles[grep("aligned.parameters.txt",pfiles)]
	for (f in 1:length(pfiles)) {
		print(paste("class",c,"file",f,sep = '.'))
		file.name <- paste(classes[c],f,sep = ".")
		files.names <- c(files.names, file.name)
		input.class <- classes[c]
		input.classes <- c(input.classes, input.class)
		parfile <- read.table(pfiles[f], header = T, sep = "\t")
		c.props = matrix(na.omit(parfile$props))
		#c.props = matrix(parfile$props[dim*c((1:all.k)-1)+1])
		#make sure props add up to 1
		non0prop = length(c.props[c.props!=0])
		if (sum(c.props) != 1) {
			residue = (1-sum(c.props))/non0prop
			c.props[c.props!=0] = c.props[c.props!=0]  + residue
		}
		c.mus = matrix(parfile$mod)
		library(compositions)
		vars<-matrix(nrow=all.k,ncol=dim*(dim+1)/2)	
		diags<-c()
		for (ii in 1:all.k) { #assemble variance matrix where each row is each component and each  column is one of 10 variance values
			tmp.vars<-as.matrix(parfile[(dim*ii-(dim-1)):(dim*ii),4:(dim+3)]) 
			vars[ii,]<-tmp.vars[lower.tri(tmp.vars,diag=T)]
			diags<-c(diags,ceiling(geometricmean(diag(tmp.vars))))
		}
		c.vars <- c()
		for (ii in 1:all.k) {
			c.vars <- c(c.vars, vars[ii,])
		} #collect vars
		c.alpha <- parfile[,(dim(parfile)[2]-1)]
		##find non-empty mclusters##
		e <- c()
		for (z in 1:all.k) {
			prop <- parfile[(z-1)*dim+1,1]
			if (prop != 0) {e <- c(e, z)}
		}

		library(mclust)
		decomp.array<-array(dim=c(dim,dim,all.k))
		sample.sigma.matrix<-as.matrix(parfile[1:(dim*all.k),4:(dim+3)])
		for (k in 0:(all.k-1)) {
			decomp.array[,,(k+1)]<-sample.sigma.matrix[(k*dim+1):(k*dim+dim),1:dim]
		}
		decomp.array<-decomp.array[,,e]
		decomp<-sigma2decomp(decomp.array)
		c.scale<-c()
		c.shape<-c()
		c.orientation<-c()
		G = 0
		for (ii in 1:all.k) {
			if (any(e == ii)) {
				G <- G+1
				personal.scale<-decomp[[2]]$scale[G] #fill in scale info
				personal.shape<-c(decomp[[2]]$shape[,G])
				personal.orientation<-c(decomp[[2]]$orientation[,,G])	
			} else {
			personal.scale <- as.vector(0)
			personal.shape <- as.vector(rep(0, dim))
			personal.orientation <- as.vector(rep(0, dim*dim))
			}
		c.scale<-c(c.scale, personal.scale)
		c.shape<-c(c.shape, personal.shape)
		c.orientation<-c(c.orientation, personal.orientation)
		}
		cluster.param <- matrix(ncol = 1, nrow = all.k*(dim*(dim+1)/2+2*dim+1+1+dim+dim^2)+1)
		cluster.param <- c(c.props, c.mus, c.vars, c.alpha, c.scale, c.shape, c.orientation)
		class.param <- cbind(class.param, cluster.param)
	}#collect one class
	allsamples.param <- cbind(allsamples.param, class.param)
}#collect all classes


###make full dataset###
naming.props<-c()  #initialize naming variables (row names)
naming.means<-c()
naming.vars<-c()
naming.alpha<-c()
naming.scale<-c()
naming.shape<-paste("shape",1:(dim*all.k))
naming.orientation<-paste("orientation",1:(dim*dim*all.k))
cov<-c()

for (i in 1:all.k) {
	naming.props<-c(naming.props,paste("prop",i,sep=""))
	naming.means<-c(naming.means,paste("mus",1:dim,".",i,sep=""))
	cov<-c()
	for (k in 1:dim) {cov<-c(cov,paste(k:dim,k,sep=""))}
	naming.vars<-c(naming.vars,paste("vars",cov,".",i,sep=""))
	naming.alpha<-c(naming.alpha,paste("alpha",1:dim,".",i,sep=""))
	naming.scale<-c(naming.scale,paste("scale",i,sep=""))
}
naming<-as.vector(c(naming.props,naming.means,naming.vars,naming.alpha,naming.scale,naming.shape,naming.orientation))

colnames <- c("Feature Index", "Description", files.names)
dataset <- matrix(ncol = num.samples + 2, 
nrow=all.k*(dim*(dim+1)/2+2*dim+1+1+dim+dim^2),
dimnames = list (c(1:(all.k*(dim*(dim+1)/2+2*dim+1+1+dim+dim^2))), colnames))
dataset[,1]<-naming
dataset[,2]<-c(1:length(naming))
dataset[,3:(num.samples+2)] <- allsamples.param

#make xls file#
write.table(dataset, file = paste(fileset.name, "allparameters.xls", sep = "."),
sep = "\t", row.names = F, quote =F)

#make gct file#
dataset <- data.frame(dataset)
cat("#1.2","\n",file=paste(fileset.name,"mixture.output", stage,"gct",sep="."))
cat(dim(allsamples.param),"\n",sep="\t",file=paste(fileset.name,"mixture.output", stage,"gct",sep="."),append=T)
write.table(dataset,sep="\t",file=paste(fileset.name,"mixture.output", stage,"gct",sep="."),append=T,row.names=F)

#make cls file#
input.names<-vector(length=num.samples)
number.people<-num.samples
number.classes<-length(classes)
cat(number.people,number.classes,1,sep=" ",file=paste(fileset.name,"mixture.output", stage,"cls",sep="."))
cat("\n",file=paste(fileset.name, "mixture.output", stage,"cls",sep="."),append=T)
cat("#"," ", classes, file=paste(fileset.name, "mixture.output", stage,"cls",sep="."),append=T)
cat("\n",file=paste(fileset.name,"mixture.output", stage,"cls",sep="."),append=T)
cat(input.classes,sep=" ",file=paste(fileset.name,"mixture.output", stage,"cls",sep="."),append=T)

###END MAKING DATASET FILE FOR CLASSIFICATION###
}