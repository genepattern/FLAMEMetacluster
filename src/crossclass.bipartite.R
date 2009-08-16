
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

#################
## metacluster ##
#################
crossclass.metacluster <- function(classname = "matched",classes,
concatfiles,density,dim,pam.k,classfied,matrices,numbersamples){

suppressMessages(library(lpSolve))
suppressMessages(library(cluster))
#library(scatterplot3d)

optimal.k.1 = pam.k[1]#
optimal.k.2 = pam.k[2]#

Number.samples.1 = numbersamples[1]#
Number.samples.2 = numbersamples[2]#

k.NN.1 = max(1,round(0.2*Number.samples.1))
k.NN.2 = max(1,round(0.2*Number.samples.2))

delta = 0.1

modes.pch=20
medoids.pch=22

###################################

# Input file
data=read.table(matrices[1],skip=1,sep='\t')

modes = as.matrix(data[,2:(1+dim)])#
prop = data[,(2+dim)]#
sample = data[,1]#

pam.out = pam(modes, k=optimal.k.1)
indices = pam.out$clustering
medoids.1 = pam.out$medoids


median.prop.1 = rep(0, optimal.k.1)

for (i in 1:optimal.k.1)
{
pick = which(indices == i)
prop.cluster = prop[pick]
sample.cluster = sample[pick]
modes.cluster = as.matrix(modes[pick,])

dist = rep(Inf, length(pick))
for (j in 1:length(pick)) dist[j] = d(modes.cluster[j,], medoids.1[i,])

dist.NN = sort(dist)[1:k.NN.1]
indices.NN = which(dist %in% dist.NN)

median.prop.1[i] = median(prop.cluster[indices.NN])
}


###############################


# Input file
data=read.table(matrices[2],skip=1,sep='\t')

modes = as.matrix(data[,2:(1+dim)])#
prop = data[,(2+dim)]#
sample = data[,1]#

pam.out = pam(modes, k=optimal.k.2)
indices = pam.out$clustering
medoids.2 = pam.out$medoids


median.prop.2 = rep(0, optimal.k.2)

for (i in 1:optimal.k.2)
{
pick = which(indices == i)
prop.cluster = prop[pick]
sample.cluster = sample[pick]
modes.cluster = as.matrix(modes[pick,])

dist = rep(Inf, length(pick))
for (j in 1:length(pick)) dist[j] = d(modes.cluster[j,], medoids.2[i,])

dist.NN = sort(dist)[1:k.NN.2]
indices.NN = which(dist %in% dist.NN)

median.prop.2[i] = median(prop.cluster[indices.NN])
}


###################################

###### to compare sample S's modes with the medoid

# start with medoid info
ncol=dim+optimal.k.2+1
matched.points=c(medoids.2[1,],1,rep(0,optimal.k.2-1),medoids.pch)
for (i in 2:optimal.k.2) matched.points=rbind(matched.points,t(c(medoids.2[i,],i,rep(0,optimal.k.2-1),medoids.pch)))

delta = 0.05

# template proportions

cap.T.1 = median.prop.1
cap.T.2 = median.prop.2

# distances between medoids of two templates

nT.1 = optimal.k.1 ; nT.2 = optimal.k.2

dist.ST = array(1:nT.1*nT.2, dim=c(nT.1,nT.2))
for (i in 1:nT.1) { for (j in 1:nT.2) {
dist.ST[i,j] = d(medoids.1[i,], medoids.2[j,])
}}

# call to IP

constraint.choices=c("11", "10", "01", "00")
min.obj = Inf

for (constraint in constraint.choices) {
try.match = IP(dist.ST, cap.T.1, cap.T.2, delta, constraint)
check = sum(try.match)
if (check) {
try.obj = sum(dist.ST * try.match)
if (try.obj < min.obj) { match = try.match }
}}

for (k in 1:nT.1) {

colors = rep(0, nT.2)

labels=which(match[k,]==1)
l.l=length(labels)

if (l.l>1) {

dist = rep(Inf, l.l)
for (j in 1:l.l) dist[j]=d(medoids.1[k,], medoids.2[labels[j],])
dist.order=order(dist)
labels=labels[dist.order]
}

colors[1:l.l]=labels

p = t(c(medoids.1[k,],colors,modes.pch))
vec=1:ncol; for (m in 1:ncol) vec[m]=p[[m]]

matched.points=rbind(matched.points,vec)

}

######## plot, view, save output##########
write.table(matched.points[,1:(dim+1)], file="crossclass.bipartite.match.txt", 
quote = F, sep = "\t", row.names = rep(classes[1:2],c(optimal.k.1,optimal.k.2)), 
col.names = FALSE)

if (.Platform$OS.type == "windows")
{
    png("crossclass.bipartite.matching.medoids.png",height=960,width=960)
}
else
{
    library(Cairo)
    CairoPNG("crossclass.bipartite.matching.medoids.png",height=960,width=960)
}

colors.1=rainbow(optimal.k.1)
colors.2=rainbow(optimal.k.2)
pairs(matched.points[,1:dim],col=c(colors.1,colors.2),pch=matched.points[,(dim+max(optimal.k.1,optimal.k.2)+1)],cex=3)
dev.off()
#########################################
#### Assign final cluster membership ####
#########################################
matching.1 = matrix(cbind(c(1:optimal.k.1),matched.points[1:optimal.k.1,(dim+1)]),ncol = 2, nrow = optimal.k.1) #
matching.2 = matrix(cbind(c(1:optimal.k.2),matched.points[(optimal.k.1+1):(optimal.k.1+optimal.k.2),(dim+1)]),ncol = 2, nrow = optimal.k.2) #
matching.1 = rbind(c(0,0),matching.1)
matching.2 = rbind(c(0,0),matching.2)
all.class.memberships <- dir('./',pattern = "metacluster.membership")
class1.memberships <- all.class.memberships[grep(paste(classes[1],'.',sep=''),all.class.memberships)]
class2.memberships <- all.class.memberships[grep(paste(classes[2],'.',sep=''),all.class.memberships)]

for (i in 1:length(class1.memberships)) {
	gc()
	file <- read.table(class1.memberships[i],header=T)
	n.class.metaclusters = ncol(file)-1-dim
	final.cluster <- c()
	for (c in (dim+2):(dim+1+n.class.metaclusters)) {
		#print(c)
		current.member <- file[,c]
		member <- vector(length = length(current.member))
		for (m in 1:length(current.member)) {
		#	print(m)
			member[m] <- matching.1[,2][matching.1[,1]==current.member[m]]
		}
		final.cluster <- cbind(final.cluster,member)
	}
	newfile <- cbind(file,final.cluster)
	names(newfile) = c(names(newfile)[1:(dim+1)],paste("c.cluster",1:n.class.metaclusters,sep=''),paste("cluster",1:n.class.metaclusters,sep=''))
	write.table(newfile,paste("matched",class1.memberships[i],sep='.'),sep='\t',row.names=F)
}

for (i in 1:length(class2.memberships)) {
	gc()
	file <- read.table(class2.memberships[i],header=T)
	n.class.metaclusters = ncol(file)-1-dim
	final.cluster <- c()
	for (c in (dim+2):(dim+1+n.class.metaclusters)) {
		#print(c)
		current.member <- file[,c]
		member <- vector(length = length(current.member))

		for (m in 1:length(current.member)) {
		#	print(m)
		    
			member[m] <- matching.2[,2][matching.2[,1]==current.member[m]]
		}
		final.cluster <- cbind(final.cluster,member)
	}
	newfile <- cbind(file,final.cluster)
	names(newfile) = c(names(newfile)[1:(dim+1)],paste("c.cluster",1:n.class.metaclusters,sep=''),paste("cluster",1:n.class.metaclusters,sep=''))
	write.table(newfile,paste("matched",class2.memberships[i],sep='.'),sep='\t',row.names=F)
}
n.metaclusters <- max(optimal.k.1,optimal.k.2)
return(n.metaclusters)
}


