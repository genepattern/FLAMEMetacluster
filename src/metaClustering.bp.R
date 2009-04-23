#LoadPackages <- function(package.names) {
#source("http://bioconductor.org/biocLite.R")
#for (i in 1:length(package.names)) {
#	package.name = package.names[i]
#	installed <- installed.packages()[,1]
#	if (!any(installed == package.name)) {
#		#install.packages(package.name, repos = "http://cran.r-project.org")
#		biocLite(package.name)
#	}
#}
#}

parseCmdLine <- function(...)
{
    options(warn=-1)
    suppressMessages(metaCluster(...))
}

cleanup <- function()
{
    file.remove(temp.files)

    files <- list.files()
    for (i in 1:length(files))
    {
        if(regexpr(paste(".zip","$",sep=""), tolower(files[[i]]))[[1]] == -1
            && tolower(files[[i]]) != "stderr.txt" && tolower(files[[i]]) != "cmd.out"
            && tolower(files[[i]]) != "stdout.txt")
        {
            file.remove(files[[i]])
        }
    }
}

metaCluster <- function(
libdir, #full pathname to where FLAME program files are kept
OptimalG, #full pathname to OptimalG,zip
#class.difference, #{T,F} whether samples belong to different classes
step=0.1,
output.intermediate.results, #default is no
output.prefix, #<studyname_dataname>
sample.class=NA #optional; full pathname to supporting file describing class membership of samples
){

zip.ext <- regexpr(paste(".zip","$",sep=""), tolower(OptimalG))
if(zip.ext[[1]] == -1)
{
    stop("Input file must be of type zip ")
}

source(paste(libdir,"common.R",sep='/'))
source(paste(libdir,"metacluster.bipartite.R",sep='/'))
source(paste(libdir,"pairplot.metacluster.R",sep='/'))
source(paste(libdir,"crossclass.bipartite.R",sep='/')) 
source(paste(libdir,"rearrange.split.merge.R",sep='/'))
source(paste(libdir,"make.aligned.ret.R",sep='/'))
source(paste(libdir,"pamk.2.R",sep='/'))
source(paste(libdir,"EmSkew.R",sep='/'))
source(paste(libdir,"unzip.R",sep='/'))
source(paste(libdir,"assemble.classification.dataset.R",sep='/'))
source(paste(libdir,"plot.heatmap.R",sep='/'))
source(paste(libdir,"zip.R",sep='/'))

if(is.na(sample.class) || sample.class == '')
{
    class.difference <- "F"
}
else
{    
    class.difference <- "T"
}

on.exit(cleanup())
wkdir <- getwd()
setwd(libdir)

isWindows <- Sys.info()[["sysname"]]=="Windows"
if(isWindows)
{
    file.copy("emskew.dll", to = wkdir)
}
else
{
    file.copy("emskew.so", to = wkdir)
}
setwd(wkdir)

if(libdir!='')
{
    setLibPath(libdir)
    install.required.packages(libdir)
}

#metacluster.packages <- function() {
#	packages <- c("lpSolve","cluster","calibrate","sn","compositions","mclust","gplots")
#	LoadPackages(packages)
#}
#source("http://bioconductor.org/biocLite.R")
#biocLite("lpSolve")
#biocLite("cluster")

#unzip preprocessed data
#optimalGfiles <- unzip.file(OptimalG, getwd())@extracted

temp.dir <- paste(wkdir, "temp", sep="/")
dir.create(temp.dir)
unzip.file(OptimalG, temp.dir)
temp.files <<- list.files(temp.dir)
unlink(temp.dir, recursive = TRUE)

unzip.file(OptimalG, getwd())

specs <- dget(dir("./",pattern="OptimalGSpecs.ret"))
bestG.range <- specs$bestG.range
dist = specs$dist
dim = specs$dim

#######################
### FIND CLASS INFO ###
#######################
#made global for removal in cleanup function
allparamfiles <<- dir("./", pattern = "parameters.txt")
allconcatfiles <<- dir("./", pattern = "membership.txt")
alllocationsfiles <<- dir("./", pattern = "locations.txt")

num.samples <- length(allparamfiles)
if (class.difference == "T") {
	class.info <- as.matrix(read.table(sample.class, header=T))
      classes <- as.vector(unique(class.info[,2]))
      num.classes = length(classes)
#	gather files for each class
	all.concat <- array(NA,dim=c(num.samples,1,num.classes))
	all.param <- array(NA,dim=c(num.samples,1,num.classes))
	all.locations <- array(NA,dim=c(num.samples,1,num.classes))
	for (c in 1:length(classes)){
		class <- classes[c]
		#print(class)
		samples <- subset(class.info, class.info[,2] == class)[,1]
		class.concat<- c()
		class.param <- c()
		class.locations <- c()
		for (i in 1:length(samples)){
			this.concat <- allconcatfiles[grep(samples[i],allconcatfiles)]
			class.concat <- c(class.concat, this.concat)
			this.param <- allparamfiles[grep(samples[i],allparamfiles)]
			class.param <- c(class.param, this.param)
			this.locations <- alllocationsfiles[grep(samples[i],alllocationsfiles)]
			class.locations <- c(class.locations, this.locations)
		}
		all.concat[,,c] <- as.array(rep(class.concat,len=num.samples))
		all.param[,,c] <- as.array(rep(class.param,len=num.samples))
		all.locations[,,c] <- as.array(rep(class.locations,len=num.samples))
	}
}

if (class.difference == "F") {
	classes = "matched"
	num.classes = 1
	all.concat <- array(allconcatfiles,dim=c(num.samples,1,num.classes))
	all.param <- array(allparamfiles,dim=c(num.samples,1,num.classes))
	all.locations <- array(alllocationsfiles,dim=c(num.samples,1,num.classes))
}

##########################
###### META-CLUSTER ######
##########################
kmin = as.numeric(min(bestG.range))
kmax = as.numeric(max(bestG.range))
if (class.difference == "F") {
       n.metaclusters <- metacluster.bp(classname = classes[1],
		concatfiles = as.matrix(allconcatfiles),
		density = dist,
		dim = dim,
		kmin = kmin,
		kmax = kmax,
		classified = F)$pam.k
} else {#class.difference == T
       #metacluster within each class
       pam.k<-c()
       matrices <- c()
	 numbersamples <- c()
       for (c in 1:num.classes) {
		metaresult <- metacluster.bp(classname = classes[c],
			concatfiles = as.matrix(unique(all.concat[,,c])),
			density = dist,
			dim = dim,
			kmin = kmin,
			kmax = kmax,
			classified=F,
			classnum=c)
		pam.k <- c(pam.k,metaresult$pam.k)
		matrices <- c(matrices, metaresult$matrix)
		numbersamples <- c(numbersamples, metaresult$numberfiles)
       }
       #metacluster across classes
	all.class.concat <- dir("./", pattern = "metacluster.membership")
	n.metaclusters <- crossclass.metacluster(classname = "final",
		concatfiles = matrix(all.class.concat),
		classes=classes,
		density = dist,
		dim = dim,
		pam.k = pam.k,
		matrices = matrices,
		numbersamples = numbersamples)
}

#######################
### REARRANGE PARAM ###
#######################
all.memberfiles <- dir(pattern="metacluster.membership")
all.memberfiles <- all.memberfiles [grep("matched.",all.memberfiles)]
all.paramfiles <- dir(pattern="parameters.txt")
all.medoids <- dir(pattern="medoids.txt")

rearrange.parameters(all.memberfiles,all.paramfiles,all.medoids,n.metaclusters,dist,class.difference,step,OS) 

#################################
##### PLOT ALIGNED CLUSTERs #####
#################################
concatfiles <- dir("./", pattern = "aligned.member.txt")
concatfiles <- concatfiles[grep("matched",concatfiles)]
###user input 2 dimensions to plot###
pairplotMetaClusters2D(concatfiles, dim,maxclus = n.metaclusters)

########################
### MAKE ALIGNED RET ###
########################
metaparamfiles <- dir("./",pattern = "aligned.parameters.txt")
makeAlignedRet (metaparamfiles = metaparamfiles)

#########################
## MOVE & REMOVE FILES ##
#########################
#file.remove(optimalGfiles)

#############################
### PLOT ALIGNED HEATMAPS ###
#############################
metaparamfiles <- dir("./",pattern = "aligned.parameters.txt")
metaconcatfiles <- dir("./",pattern = "aligned.member.txt")
metaconcatfiles <- metaconcatfiles[grep("member",metaconcatfiles)]
plotHeatmap(paramfiles = metaparamfiles,
concatfiles = metaconcatfiles,
dist = dist)

#########################################
### MAKE DATSET FOR FEATURE SELECTION ###
#########################################
num.samples=length(dir("./", pattern= "aligned.parameters.txt"))
makeDataset(stage="train",
classes = classes,
all.k = n.metaclusters,
dim = dim,
num.samples=num.samples,
fileset.name=output.prefix)

zip.file(libdir = libdir,files =  "*legend.png",outfile = paste(output.prefix,"FinalAlignedClusters.zip",sep='.'))
zip.file(libdir = libdir,files =  "*.pairplots.png",outfile = paste(output.prefix,"FinalAlignedClusters.zip",sep='.'))
zip.file(libdir = libdir,files =  "*.heatmap.png",outfile = paste(output.prefix,"FinalAlignedClusters.zip",sep='.'))
zip.file(libdir = libdir,files =  "*.locations.txt",outfile = paste(output.prefix,"FinalAlignedClusters.zip",sep='.'))
zip.file(libdir = libdir,files =  "*.parameters.txt",outfile = paste(output.prefix,"FinalAlignedClusters.zip",sep='.'))
zip.file(libdir = libdir,files =  "*.member.txt",outfile = paste(output.prefix,"FinalAlignedClusters.zip",sep='.'))
zip.file(libdir = libdir,files =  "*.aligned.ret",outfile = paste(output.prefix,"FinalAlignedClusters.zip",sep='.'))
zip.file(libdir = libdir,files =  "*.allparameters.xls",outfile = paste(output.prefix,"FinalAlignedClusters.zip",sep='.'))
zip.file(libdir = libdir,files =  "*.output.train.gct",outfile = paste(output.prefix,"FinalAlignedClusters.zip",sep='.'))
zip.file(libdir = libdir,files =  "*.output.train.cls",outfile = paste(output.prefix,"FinalAlignedClusters.zip",sep='.'))

if(output.intermediate.results == "T")
{
    zip.file(libdir = libdir, files = "*.asw.txt", outfile = paste(output.prefix,"MetaclusterIntermediateResults.zip", sep='.'))
    zip.file(libdir = libdir, files = "*.asw.png", outfile = paste(output.prefix,"MetaclusterIntermediateResults.zip", sep='.'))

    zip.file(libdir = libdir, files = "*.metacluster.membership", outfile = paste(output.prefix,"MetaclusterIntermediateResults.zip", sep='.'))
    zip.file(libdir = libdir, files = "*.metacluster.assignments.txt", outfile = paste(output.prefix,"MetaclusterIntermediateResults.zip", sep='.'))
    zip.file(libdir = libdir, files = "*.metacluster.pca.png", outfile = paste(output.prefix,"MetaclusterIntermediateResults.zip", sep='.'))

    zip.file(libdir = libdir, files = "*.pooledlocations.txt", outfile = paste(output.prefix,"MetaclusterIntermediateResults.zip", sep='.'))

    zip.file(libdir = libdir, files = "*.silhouette.png", outfile = paste(output.prefix,"MetaclusterIntermediateResults.zip", sep='.'))
    zip.file(libdir = libdir, files = "*.medoids.txt", outfile = paste(output.prefix,"MetaclusterIntermediateResults.zip", sep='.'))

    if(length(list.files(pattern=".match.txt")) != 0)
    {
        zip.file(libdir = libdir, files = "*.match.txt", outfile = paste(output.prefix,"MetaclusterIntermediateResults.zip", sep='.'))
    }
}
}

install.required.packages <- function(libdir)
{
    if(!is.package.installed(libdir, "mclust"))
    {
		install.package(libdir, "mclust_3.1-10.zip", "mclust_3.1-10.tgz", "mclust_3.1-10.tar.gz")
	}
	if(!is.package.installed(libdir, "lpSolve"))
	{
		install.package(libdir, "lpSolve_5.6.4.zip", "lpSolve_5.6.4.tgz", "lpSolve_5.6.4.tar.gz")
	}
	if(!is.package.installed(libdir, "gtools"))
	{
		install.package(libdir, "gtools_2.5.0.zip", "gtools_2.5.0.tgz", "gtools_2.5.0.tar.gz")
	}
	if(!is.package.installed(libdir, "gdata"))
	{
		install.package(libdir, "gdata_2.4.2.zip", "gdata_2.4.2.tgz", "gdata_2.4.2.tar.gz")
	}	
	if(!is.package.installed(libdir, "gplots"))
	{
		install.package(libdir, "gplots_2.6.0.zip", "gplots_2.6.0.tgz", "gplots_2.6.0.tar.gz")
	}
	if(!is.package.installed(libdir, "scatterplot3d"))
    {
		install.package(libdir, "scatterplot3d_0.3-27.zip", "scatterplot3d_0.3-27.tgz", "scatterplot3d_0.3-27.tar.gz")
	}
	if(!is.package.installed(libdir, "rgl"))
    {
		install.package(libdir,"rgl_0.82.zip", "rgl_0.82.tgz", "rgl_0.82.tar.gz")
	}
	if(!is.package.installed(libdir, "tensorA"))
    {
		install.package(libdir, "tensorA_0.31.zip", "tensorA_0.31.tgz", "tensorA_0.31.tar.gz")
	}
	if(!is.package.installed(libdir, "compositions"))
    {
		install.package(libdir, "compositions_0.91-6.zip", "compositions_0.91-6.tgz","compositions_0.91-6.tar.gz")
	}
}
