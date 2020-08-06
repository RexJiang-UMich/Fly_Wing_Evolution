library(resample)
library(mvtnorm)
library(MBESS)

#parameters & functions involved in all models
Ne=10^6

#fitness function (input can be number or vector)
fitness<-function(v){
	#Euclidean distance
	d=sqrt(sum(v^2))
	#Gaussian fitness function
	w=dnorm(d,mean=0,sd=1)/dnorm(0,mean=0,sd=1)
	if(w<0){
		w=0
	}
	return(w)
}

#selective coefficient
sc<-function(v1,v2){
	wm=fitness(v1);wa=fitness(v2)
	if(wa>0){
	    s=(wm/wa)-1
	}else{
		if(wm==0){
			s=0
		}else{
			s=999
		}
	}
	return(s)
}

#fixation probability
fp<-function(v1,v2){
	s=sc(v1,v2)
	if(s==0){
		p=1/(2*Ne)
	}else{
		p=(1-exp(-2*s))/(1-exp(-4*Ne*s))
	}
	return(p)
}

#constant effect size distribution and variable mutation rate
sdsize=0.01
ntest=20;T=1000;nrep=50
vm_all=rgamma(ntest,scale=0.05,shape=0.5)
lm_all=vm_all/(sdsize^2)
ntrait_all=as.integer(1+20+1.7*log10(vm_all));pl=ntrait_all-1
physig=rep(0,ntest)
div=rep(0,ntest)
for(test in 1:ntest){
	
	lm=lm_all[test]
	ntrait=ntrait_all[test] #number of correlated traits (including focal)
	if(ntrait<=0){
		ntrait=1
	}
	cor=rep(0.5,ntrait);cor[1]=1
	alpha=rep(1,ntrait);alpha[1]=1e-4

	value=matrix(0,nrow=nrep,ncol=(T+1))

	for(r in 1:nrep){
	    v=matrix(0,nrow=ntrait,ncol=(T+1)) #matrix of phenotype through time
		for(t in 2:(T+1)){
			v[,t]=v[,(t-1)] #ancestral state
	    	nmut=rpois(n=1,lambda=lm) #number of mutations aquired
	    	if(nmut>0){
	        	for(x in 1:nmut){
		        	effect=rbinom(n=(ntrait),size=1,prob=cor)*rnorm((ntrait),mean=0,sd=sdsize) #effect on each trait
		        	phe=v[,(t-1)]+effect
		        	if(ntrait>1){
		        	    fixprob=2*Ne*fp(alpha*phe,alpha*v[,(t-1)])
		        	    if(fixprob>1){
		        	    	fixprob=1
		        	    }
		        	}else{
		        		fixprob=1
		        	}
		        	fix=rbinom(n=1,size=1,prob=fixprob) #determine if fixation happens
		        	if(fix==1){
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

pl;vm_all;div;physig

#constant mutation rate and variable effect size distributions
ntest=20;T=1000;nrep=50
vm_all=rgamma(ntest,scale=0.05,shape=0.5)
lm=500
sdsize_all=sqrt(vm_all/lm)
ntrait_all=as.integer(1+20+1.7*log10(vm_all));pl=ntrait_all-1
physig=rep(0,ntest)
div=rep(0,ntest)
for(test in 1:ntest){
	
	vm=vm_all[test]
	sdsize=sdsize_all[test]
	ntrait=ntrait_all[test] #number of correlated traits (including focal)
	if(ntrait<=0){
		ntrait=1
	}
	cor=rep(0.5,ntrait);cor[1]=1
	alpha=rep(1,ntrait);alpha[1]=0

	value=matrix(0,nrow=nrep,ncol=(T+1))

	for(r in 1:nrep){
	    v=matrix(0,nrow=ntrait,ncol=(T+1)) #matrix of phenotype through time
		for(t in 2:(T+1)){
			v[,t]=v[,(t-1)] #ancestral state
	    	nmut=rpois(n=1,lambda=lm) #number of mutations aquired
	    	if(nmut>0){
	        	for(x in 1:nmut){
	        		effect=c(rnorm(1,mean=0,sd=sdsize),rbinom(n=(ntrait-1),size=1,prob=cor)*rnorm((ntrait-1),mean=0,sd=0.01)) #effect on each trait
		        	phe=v[,(t-1)]+effect
		        	if(ntrait>1){
		        	    fixprob=2*Ne*fp(alpha*phe,alpha*v[,(t-1)])
		        	    if(fixprob>1){
		        	    	fixprob=1
		        	    }
		        	}else{
		        		fixprob=1
		        	}
		        	fix=rbinom(n=1,size=1,prob=fixprob) #determine if fixation happens
		        	if(fix==1){
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

pl;vm_all;div;physig

#slope analysis
v1=log10(vm_all);v2=log10(div/T)
fit=lm(v2~v1)
slope=fit$coef[[2]];cor=cor(v1,v2)
sfit<-summary(fit)
tstats<-(1-slope)/sfit$coefficients[2,2]
pval<-2*pt(abs(tstats),df=df.residual(fit),lower.tail=FALSE)
pval



