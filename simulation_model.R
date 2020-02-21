library(clusterGeneration)
library(MBESS)
library(ggplot2)
library(resample)

Ne=10^6;fmin=1-(1/(2*Ne));sdsize=0.01

ntest=20;T=2000;nrep=50
vm_all=rep(0,ntest);pl=rep(0,ntest);physig=rep(0,ntest)
div=rep(0,ntest)
for(test in 1:ntest){
	
	vm=rgamma(1,shape=0.5,rate=1/400);vm_all[test]=vm #mutational variance of focal trait; to get phenotype level variance, divide it by variance of mutation size
	ntrait=16+as.integer(0.5*log2(vm/10)) #total number of traits
	if(ntrait<=0){
		ntrait=1
	}
	pl[test]=ntrait-1
	cor=rep(0.5,ntrait);cor[1]=1
	alpha=rep(1,ntrait);alpha[1]=0

	value=matrix(0,nrow=nrep,ncol=(T+1))
	for(r in 1:nrep){
	    v=matrix(0,nrow=ntrait,ncol=(T+1)) #matrix of phenotype through time
		for(t in 2:(T+1)){
			v[,t]=v[,(t-1)] #ancestral state
	    	nmut=rpois(n=1,lambda=vm) #number of mutations aquired
	    	if(ntrait>1){
	    	    f0=1-(sum((alpha*v[,(t-1)])^2))^(e/2)
	    	}else{
	    	    f0=1
	    	}
	    	if(nmut>0){
	        	for(x in 1:nmut){
		        	effect=rbinom(n=(ntrait),size=1,prob=cor)*rnorm((ntrait),mean=0,sd=sdsize) #effect on each trait
		        	phe=v[,(t-1)]+effect
		        	if(ntrait>1){
		        	    fmean=1-(sum((alpha*phe)^2))^(e/2) #based on euclidean distance
		        	}else{
		        		fmean=1
		        	}
		        	if(fmean>f0){ #mutations increasing overall fitness would all fix
		        		fixprob=1
		        	}else{
		        		if(f0==0){
		        			fixprob=0
		        		}else{
		        			if((fmean/f0)>fmin){
		        				fixprob=1
		        			}else{
		        				fixprob=0
		        			}
		        		}
		        	}
		        	if(fixprob==1){
		        		v[,t]=v[,t]+effect
		        	}
		    	}
	    	}
    	}
    	value[r,]=v[1,]
    }
    physig[test]=cor(0:(T),colVars(value))
    div[test]=var(value[,(T+1)])
}
