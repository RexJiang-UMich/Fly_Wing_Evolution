library(phytools)
library(clusterGeneration)
library(mvtnorm)
library(geiger)
library(MBESS)

#read the empirical tree
dist<-read.table("dist.txt",sep="\t")
dm=dist[,2:111]
dm=data.matrix(dm)
colnames(dm)=dist[,1];rownames(dm)=dist[,1]
tr<-vcv2phylo(dm)

#read empirical M matrix
m_empirical<-read.table("hom.txt",sep="\t");m_empirical=data.matrix(m_empirical)
r_empirical<-read.table("r_mat_shape.txt",sep="\t");r_empirical=data.matrix(r_empirical)

#simulation
ntrait=24
Ntest=50
slope_all=matrix(0,nrow=Ntest,ncol=4)
cor_all=matrix(0,nrow=Ntest,ncol=4)

for(test in 1:Ntest){

fvar=rgamma(ntrait,shape=0.025,rate=1/mean(diag(m_empirical))) #shape parameter manually set
f1=cor2cov(rcorrmatrix(ntrait),sd=sqrt(fvar))

m_true=m_empirical
f=f1*(sum(diag(m_true))/sum(diag(f1)))
fweight=0.05 #manually set
r_true=(1-fweight)*m_true+fweight*f

nsample=150
m_dat=rmvnorm(nsample,sigma=m_true)
m_obs=cov(m_dat)

#simulate evolution on the tree (BM model)
#amount of evolution on each branch
evo=list()
for(i in 1:length(tr$edge.length)){
	evo[[i]]=rmvnorm(1,sigma=tr$edge.length[i]*r_true)
}
#convert to phenotype at each node
phe=list()
rdat=rep(0,ntrait)
for(i in 1:nrow(tr$edge)){
	phe[[i]]=matrix(0,nrow=ntrait,ncol=2)
	ances=which(tr$edge[,2]==tr$edge[i,1])
	if(length(ances)>0){
		phe[[i]][,1]=phe[[ances]][,2]
	}else{
		phe[[i]][,1]=rep(0,ntrait)
	}
	phe[[i]][,2]=phe[[i]][,1]+t(evo[[i]][1,])
	descend=which(tr$edge[,1]==tr$edge[i,2])
	if(length(descend)==0){
		rdat=rbind(rdat,t(phe[[i]][,2]))
	}
}
rdat=rdat[2:nrow(rdat),]
rownames(rdat)=tr$tip

r=ratematrix(tr,rdat)

nrank=18
#calculate slopes
e_true=eigen(m_true)$vectors
v1=diag((t(e_true)%*%m_true)%*%e_true)
v2=diag((t(e_true)%*%r_true)%*%e_true)
v1=v1[1:nrank];v2=v2[1:nrank]
slope_all[test,1]=lm(log10(v2)~log10(v1))[[1]][2];cor_all[test,1]=cor(log10(v1),log10(v2))
v2=diag((t(e_true)%*%r)%*%e_true)
v2=v2[1:nrank]
slope_all[test,2]=lm(log10(v2)~log10(v1))[[1]][2];cor_all[test,2]=cor(log10(v1),log10(v2))

em=eigen(m_obs)$vectors
v1=diag((t(em)%*%m_obs)%*%em)
v2=diag((t(em)%*%r)%*%em)
v1=v1[1:nrank];v2=v2[1:nrank]
slope_all[test,3]=lm(log10(v2)~log10(v1))[[1]][2];cor_all[test,3]=cor(log10(v1),log10(v2))

rhat=r*(sum(diag(m_obs))/sum(diag(r)))
h=(m_obs+rhat)/2
eh=eigen(h)$vectors
v1=diag((t(eh)%*%m_obs)%*%eh)
v2=diag((t(eh)%*%r)%*%eh)
v1=v1[1:nrank];v2=v2[1:nrank]
slope_all[test,4]=lm(log10(v2)~log10(v1))[[1]][2];cor_all[test,4]=cor(log10(v1),log10(v2))

}

mean(slope_all[,1]);mean(slope_all[,2]);mean(slope_all[,3]);mean(slope_all[,4])
mean(cor_all[,1]);mean(cor_all[,2]);mean(cor_all[,3]);mean(cor_all[,4])

