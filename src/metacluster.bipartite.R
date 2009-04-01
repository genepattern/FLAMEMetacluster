
######## Euclidean distance ###########################

d=function(x,y) { sqrt(sum((x-y)^2)) }

######## Integer Programming ##########################

IP = function(dist.ST, cap.S, cap.T, delta, constraint) {

nS = dim(dist.ST)[1]; nT = dim(dist.ST)[2]; n=nS+nT

f.obj = as.vector(t(dist.ST))

f.con.1 = f.con.3 = array(1:nS*nS*nT, dim=c(nS, nS*nT)) * 0
f.con.2 = f.con.4 = array(1:nT*nS*nT, dim=c(nT, nS*nT)) * 0

s=seq(1,nT)
for (i in 1:nS) { f.con.1[i,s]=cap.T; f.con.3[i,s]=rep(1,nT); s=s+nT }

t=seq(1,nS*nT,nT)
for (i in 1:nT) { f.con.2[i,t]=cap.S; f.con.4[i,t]=rep(1,nS); t=t+1 }

if (constraint=="11") {
f.con = rbind(f.con.1,f.con.2,f.con.3,f.con.4)
f.rhs = c(cap.S+delta, cap.T+delta, rep(1, n))
f.dir = c(rep("<=", n), rep(">=", n))
}  

if (constraint=="10") {
f.con = rbind(f.con.1,f.con.2,f.con.3)
f.rhs = c(cap.S+delta, cap.T+delta, rep(1, nS))
f.dir = c(rep("<=", n), rep(">=", nS))
}  

if (constraint=="01") {
f.con = rbind(f.con.3, f.con.4)
f.rhs = c(rep(1, nS),rep(1, nT))
f.dir = rep(">=", n)
}

if (constraint=="00") {
f.con = f.con.3
f.rhs = rep(1, nS)
f.dir = rep(">=", nS)
}


IP=lp ("min", f.obj, f.con, f.dir, f.rhs, all.bin=T)

IP.solution=IP$solution

matrix(IP$solution, nrow=nS, byrow=T)

}

metacluster.bp <- function(classname, concatfiles, density, dim, max.reiter = 2, kmin, kmax, classnum = 0,classified) {
	
	suppressMessages(library(lpSolve))
	suppressMessages(library(cluster))

	##make gfile & list of filenames##
	cat("Class:",classname,'\n')
	cat("Meta-clustering",'  ........\n')
	filenames <- c()
	classfilenames <- c()
	for (i in 1:nrow(concatfiles)) {
		filename <- strsplit(concatfiles[i,1],paste('.',density,sep=''))[[1]][1]
		filenames <- c(filenames, filename)
		classfilename <- paste(classname, i,sep = '.')
		if (classified == T) {
			classfilename <- strsplit(concatfiles[i,1],"\\.")[[1]][1]
		}
		classfilenames <- c(classfilenames, classfilename)
	}
	
	#generate meansmatrix#
	gfile <- c()
	for (i in 1:nrow(concatfiles)) {
		file <- dir("./", pattern = as.character(concatfiles[i,1]))
		datafile <- read.table(file, header = T, sep = "\t")
		g <- length(unique(datafile$cluster))
		gfile <- c(gfile, g)	
	}
	meansmatrix <- matrix(ncol = dim+1, nrow = sum(gfile), dimnames = list(rep(filenames, gfile), c(paste("v",1:dim,sep=""),"props")))
	col.names = colnames(datafile)[1:dim]
	for (i in 1:nrow(concatfiles)) {
		#print(i)
		file <- dir("./", pattern = as.character(concatfiles[i,1]))
		datafile <- read.table(file, header = T, sep = "\t")
		total.size <- dim(datafile)[1]
		g <- unique(datafile$cluster)
		for (m in 1:length(g)) {
			n <- g[m]
			clus<- subset(datafile, datafile$cluster == n)
			size <- dim(clus)[1]
			prop <- size/total.size
			for (d in 1:dim) {
				meansmatrix[(m+sum(gfile[0:(i-1)])), d] <- mean(clus[,d])	
			}
			meansmatrix[(m+sum(gfile[0:(i-1)])), d+1] <- prop	

		}
	}
	allfiles <- rep(as.character(filenames),gfile)
	pooledlocation<-paste(classname,"withprops.pooledlocations.txt", sep = ".")
	write.table(meansmatrix, sep = "\t", quote = F, row.names = allfiles, file = pooledlocation)
	meanfile <- meansmatrix[,1:dim]
	c <- pamk.2(meanfile, k = kmin:kmax,classname=classname)
	pam.k <- c$nc
	pamclus <- pam(meanfile, k = pam.k)

optimal.k = pam.k
Number.samples = length(concatfiles)
k.NN = max(1,trunc(Number.samples * 0.2))
delta = 0.1
modes.pch=20
medoids.pch=22

# Input file
data=meansmatrix
modes = data[,1:dim]
prop = data[,(dim+1)]
sample = allfiles
unique.sample <- unique(allfiles)
pam.out = pam(modes, k=optimal.k)
indices = pam.out$clustering
medoids = pam.out$medoids

median.prop = rep(0, optimal.k)

for (i in 1:optimal.k)
{
pick = which(indices == i)
prop.cluster = prop[pick]
sample.cluster = sample[pick]
modes.cluster = as.matrix(modes[pick,])

dist = rep(Inf, length(pick))
for (j in 1:length(pick)) dist[j] = d(modes.cluster[j,], medoids[i,])

dist.NN = sort(dist)[1:k.NN]
indices.NN = which(dist %in% dist.NN)

median.prop[i] = median(prop.cluster[indices.NN])
}


###### to compare sample S's modes with the medoid

# start with medoid info
ncol=dim+optimal.k+1
matched.points=c(medoids[1,],1,rep(0,optimal.k-1),medoids.pch)
for (i in 2:optimal.k) matched.points=rbind(matched.points,t(c(medoids[i,],i,rep(0,optimal.k-1),medoids.pch)))

for (S in 1:Number.samples) {

delta = 0.05

sample.id = which(sample==unique.sample[S])
cap.S = prop[sample.id]
modes.sample = modes[sample.id,]

# template proportions

cap.T = median.prop

# distances between modes and medoids

nS = length(sample.id); nT = optimal.k

dist.ST = array(1:nS*nT, dim=c(nS,nT))
for (i in 1:nS) { for (j in 1:nT) {
dist.ST[i,j] = d(modes.sample[i,], medoids[j,])
}}

# call to IP

constraint.choices=c("11", "10", "01", "00")
min.obj = Inf

for (constraint in constraint.choices) {
try.match = IP(dist.ST, cap.S, cap.T, delta, constraint)
check = sum(try.match)
if (check) {
try.obj = sum(dist.ST * try.match)
if (try.obj < min.obj) { match = try.match }
}}

# rownames(match) = paste( as.character(S), ".", seq(1:nS) )
# colnames(match) = paste( "T.", seq(1:nT) )

for (k in 1:nS) {

colors = rep(0, nT)

labels=which(match[k,]==1)
l.l=length(labels)

if (l.l>1) {

dist = rep(Inf, l.l)
for (j in 1:l.l) dist[j]=d(modes.sample[k,], medoids[labels[j],])
dist.order=order(dist)
labels=labels[dist.order]
}

colors[1:l.l]=labels

p = t(c(modes.sample[k,],colors,modes.pch))
vec=1:ncol; for (m in 1:ncol) vec[m]=p[[m]]

matched.points=rbind(matched.points,vec)

}
}

write.table(matched.points[1:optimal.k,1:(dim+1)], file=paste(classname,"medoids.txt",sep='.'),quote = FALSE, sep = "\t", row.names = FALSE, col.names = c(col.names,"metacluster#"))


out <- matched.points
out <- subset(out,out[,(dim+optimal.k+1)]!=22)[,1:(dim+optimal.k)]
metafile <- out


mcluster <- (metafile[,(dim+1):(dim+optimal.k)])
metafile <- cbind(meanfile[,1:dim], mcluster)
columns <- paste("member",1:optimal.k,sep='')
write.table (metafile, file = paste(classname,"medoids.metacluster.assignments.txt", sep = "."),quote = F, col.names = c(col.names,columns),row.names =allfiles,sep="\t")
	
	#plot results of metacluster
	#silhouette plot
	if (.Platform$OS.type == "windows")
    {
	    png(filename = paste(classname,"metacluster.silhouette.png",sep = "."),width = 960,height = 960)	
	}
	else
	{
        library(Cairo, lib.loc=Sys.getenv("R_LIBS"))
        CairoPNG(filename = paste(classname,"metacluster.silhouette.png",sep = "."),width = 960,height = 960)	
	}
	plot(pamclus, ask = FALSE, which.plots = 2)
	dev.off()
	
	#clusplot (PCA)
	mtitle <- paste(classname,"metacluster.pca", sep = ".")

	if (.Platform$OS.type == "windows")
    {
	    png(filename = paste(mtitle, "png", sep = "."), width = 960, height = 960)
	}
	else
	{
        library(Cairo, lib.loc=Sys.getenv("R_LIBS"))
	    CairoPNG(filename = paste(mtitle, "png", sep = "."), width = 960, height = 960)
	}
	clusplot(pamclus, main = mtitle, lines = 0, labels = 4)
	dev.off()

	#make new concat file, collect patterns#
	for (i in 1:nrow(concatfiles)) {
		file <- dir("./", pattern = as.character(concatfiles[i,1]))
		tfile <- read.table(file, header = T, sep = "\t")
		g <- unique(tfile$cluster)
		row.start <- sum(as.numeric(gfile[0:(i-1)]))+1
		row.end <- row.start+length(g)-1
		assign <- mcluster[row.start:row.end,1:optimal.k]  #m-clusters assignment for each cluster
		count <- cbind(g, assign)
		newfile <- c()
		for (tclus in g[1:length(g)]) {
			clus<-subset(tfile, tfile$cluster == tclus)
			mclus <- count[,2:(optimal.k+1)][count[,1] == tclus]
			mclus <- matrix(rep(mclus,each=nrow(clus)),ncol=optimal.k)
			clus <- cbind(clus, mclus)
			newfile<-rbind(newfile, clus)
		}
		
		#if (classname == "final" & classified == T) {
		#	colnames = c(names(newfile)[1:dim],"o.cluster", paste("c.cluster",1:optimal.k,sep=''), paste("cluster",1:optimal.k,sep=''))
		#}  else {
			colnames = c(names(newfile)[1:dim], "o.cluster", paste("cluster",1:optimal.k,sep=''))
		#}
		filename <- paste(classname, filenames[i],density, pam.k,"metacluster.membership", sep = ".")
		write.table(newfile, file = filename,
			row.names = F, col.names = colnames, quote = F, sep = "\t")
	}
	cat('\n')
	results <- c()
	results$pam.k <- pam.k
	results$matrix <- pooledlocation
	results$numberfiles <- Number.samples
	return(results)
}