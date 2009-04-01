rearrange.parameters <- function(all.memberfiles,all.paramfiles,all.medoids,n.metaclusters,dist,class.difference,step,OS) {
	for (i in 1:length(all.memberfiles)) {
		paramfile <- all.paramfiles[i]
		filename <- strsplit(paramfile,paste('.',dist,sep=''))[[1]][1]
		cat(filename,'\n')
		memberfile <- all.memberfiles[grep(filename,all.memberfiles)]
		if (class.difference=="T"){
			prefix <- strsplit(memberfile,paste('.',filename,sep=''))[[1]][1]
			classname = strsplit(prefix,"matched.")[[1]][2]
		} else { classname="matched" }
		medoidfile <- all.medoids[grep(classname,all.medoids)]
		rearrange(memberfile,paramfile,medoidfile,n.metaclusters,dist,class.difference,classname,filename,OS=OS,step=step)
	}
}

rearrange <- function(memberfile,paramfile,medoidfile,n.metaclusters,
dist = "mst",class.difference,classname,filename,OS="Windows",step=0.1) {

suppressMessages(library(cluster))
member <- read.table(memberfile,header=T)
parameters <- read.table(paramfile,header=T)
medoids <- read.table(medoidfile,header=T)
dim = ncol(parameters)- 5
filename <- strsplit(paramfile, paste(".",dist,sep=''))[[1]][1]
ndist=switch(dist, "mvn"=1,"mvt"=2, "msn"=3, "mst"=4)

#########################################
### ASSIGN MEMBERSHIP AFTER SPLITTING ###
####################################################
# SPLITTING PERFORMED BY EMSKEW USING PAMK MEDOIDS #
####################################################
cluster.assignments <- unique(member[,(dim+1):ncol(member)])
o.clusters <- cluster.assignments[,1]
 
if (class.difference == "T") {
c.clusters <- cluster.assignments[,2:(1+n.metaclusters)]
f.clusters <- cluster.assignments[,(2+n.metaclusters):ncol(cluster.assignments)]

after.split.member <- c()
for (o in sort(o.clusters)){
	print(o)
	index <- which(o.clusters == o)
	member.o <- as.matrix(subset(member,member$o.cluster == o))
	dat <- as.matrix(member.o[,1:dim])
	c <- c.clusters[index,]
	all.c <- c[c!=0]
	if (length(all.c)==1) new.c.cluster=rep(all.c,len=nrow(dat))
	if (length(all.c)>1) {
		k = length(all.c)
		med <- as.matrix(medoids[all.c,])[,1:dim]
#		indices <- c()
#		for (r in 1:nrow(med)) {
#			med.r = med[r,]
#			index = 1:nrow(member)
#			for (d in 1:dim) {
#				index <- index[grep(index,which(abs(dat[,d]-med.r[d]) < 1))]
#			}
#			indices <- c(indices,index[1])
#		}
		clust = pam(dat,k=k)$clustering#,medoids = indices)$clustering
		new.c.cluster <- emskewfit(dat=dat,ng=k,clust=clust,dist=ndist,ncov=3,itmax=1000,epsilon=0.0001,method=1)$clust
		new.c.cluster = all.c[clust]
	}
	after.split.clusters <- c(length=length(new.c.cluster))
	for (m in 1:length(new.c.cluster)) {
		after.split.clusters[m] = f.clusters[which(c.clusters[,1]==new.c.cluster[m])[1],1]
	}
	dat.after.split = cbind(member.o,after.split.clusters)
	after.split.member <- rbind(after.split.member,dat.after.split)
}
}

if (class.difference == "F") {
c.clusters <- cluster.assignments[,2:(1+n.metaclusters)]
after.split.member <- c()
for (o in sort(o.clusters)){
	print(o)
	index <- which(o.clusters == o)
	member.o <- as.matrix(subset(member,member$o.cluster == o))
	dat <- member.o[,1:dim]
	c <- c.clusters[index,]
	all.c <- c[c!=0]
	if (length(all.c)==1) new.c.cluster=rep(all.c,len=nrow(dat))
	if (length(all.c)>1) {
		k = length(all.c)
		med <- as.matrix(medoids[all.c,])[,1:dim]
		indices <- c()
		for (r in 1:nrow(med)) {
			med.r = med[r,]
			index = 1:nrow(member)
			for (d in 1:dim) {
				index <- index[grep(index,which(abs(dat[,d]-med.r[d]) < 1))]
			}
			indices <- c(indices,index[1])
		}
		clust = pam(dat,k=k,medoids = indices)$clustering
		new.c.cluster <- emskewfit(dat=dat,ng=k,clust=clust,dist=dist,ncov=3,itmax=1000,epsilon=0.0001,method=1)$clust
		new.c.cluster = all.c[clust]
	}
	after.split.clusters <- new.c.cluster
	dat.after.split = cbind(member.o,after.split.clusters)
	after.split.member <- rbind(after.split.member,dat.after.split)
}
}

after.split.member <- data.frame(after.split.member)
final.member <- subset(after.split.member,select=c(1:dim,ncol(after.split.member)))
member.names <- c(colnames(member)[1:dim],"cluster")

############################################
### MERGE CLUSTERS, REARRANGE PARAMETERS ###
############################################
num.o.clusters = nrow(parameters)/dim
original.final.assignment <- unique(after.split.member[,c((dim+1),ncol(after.split.member))])
o.param <- array(dim = c(dim,dim+5,num.o.clusters))
for (n in 1:num.o.clusters) {
	o.param[,,n] = as.matrix(parameters[((n-1)*dim+1):(n*dim),])
}
final.parameters <- array(dim = c(dim,dim+5,n.metaclusters))
for (c in 1:n.metaclusters) {
	#refer to no.original clusters making up the final cluster
	num.o.clus <- length(which(original.final.assignment[,2]==c))
	if (num.o.clus == 0) {
		c.param <- matrix(0,ncol=(dim+5),nrow=dim)
		c.param[2:dim,c(1,3)] = NA
		final.parameters[,,c]=c.param
	}
	# if 1, take parameters and relocate accordingly
	if (num.o.clus == 1) {
		o.cluster = original.final.assignment[which(original.final.assignment[,2]==c),1]
		final.parameters[,,c] = o.param[,,o.cluster]
	}
	# if > 1, merge using emskew ng=1
	if (num.o.clus > 1) {
		dat <- as.matrix(subset(after.split.member,after.split.member$after.split.clusters==c)[,1:dim])
		c.param <- matrix(0,ncol=(dim+5),nrow=dim)
		c.param[1,1] = dim(dat)[1]/dim(member)[1]
		obj <- EmSkew(dat=dat,ng=1,dist=ndist,ncov=3,seed=123456,OS=OS)
		if (ndist == 3 | ndist == 4) {
		obb <-  EmSkewMOD(ndist,obj$mu,obj$sigma,step=step,obj$delta,obj$dof)
		obj$mod <- obb$modpts
		}
		if (ndist == 1 | ndist == 2) {
		obj$mod = obj$mu
		}
		c.param[,2] <- obj$mu
		c.param[1,3] <- obj$dof
		c.param[,4:(3+dim)] <- obj $sigma
		c.param[,(4+dim)] <- obj$delta
		c.param[,(5+dim)] <- obj$mod
		c.param[2:dim,c(1,3)]=NA
		final.parameters[,,c] = c.param
	}
}

newparameters <- c()
for (x in 1:n.metaclusters){
	newparameters <- rbind(newparameters,final.parameters[,,x])
}
param.names <- c("props", "mus", "df", paste("Var", 1:dim, sep = ""), "alpha","mod")
if (class.difference == "T") {
	newparamfilename <- paste("matched", classname,filename, dist, n.metaclusters,"aligned.parameters.txt",sep = ".")
	newmemberfilename<- paste("matched", classname,filename, dist, n.metaclusters,"aligned.member.txt",sep = ".")
}
if (class.difference == "F") {
	newparamfilename <- paste("matched", filename, dist, n.metaclusters,"aligned.parameters.txt",sep = ".")
	newmemberfilename<- paste("matched",filename, dist, n.metaclusters,"aligned.member.txt",sep = ".")
}
write.table(newparameters, sep = "\t", quote = F, row.names = F, col.names=param.names,file = newparamfilename)
write.table(final.member,sep="\t",quote=F,row.names=F,col.names = member.names,file=newmemberfilename)
}
